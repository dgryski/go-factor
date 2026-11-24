package qs

import (
	"math/bits"
	"math/rand"
)

// This is Montgomery's block Lanczos over GF(2), from "A Block Lanczos
// Algorithm for Finding Dependencies over GF(2)", EUROCRYPT 1995.

// mat64 stores a 64x64 GF(2) matrix row-major: bit j of m[i] is entry (i, j).
type mat64 [64]uint64

func mulMat64(a, b mat64) mat64 {
	var out mat64
	for i := 0; i < 64; i++ {
		row := a[i]
		var acc uint64
		for row != 0 {
			k := bits.TrailingZeros64(row)
			acc ^= b[k]
			row &= row - 1
		}
		out[i] = acc
	}
	return out
}

func addMat64(a, b mat64) mat64 {
	var out mat64
	for i := 0; i < 64; i++ {
		out[i] = a[i] ^ b[i]
	}
	return out
}

// mulVecMat64 computes dst += v * m over GF(2), where v is R x 64 and m is 64 x 64.
// Uses an 8-bit lookup table (Method of Four Russians) for 8x fewer ops per row.
func mulVecMat64(dst, v []uint64, m mat64) {
	var lut [8][256]uint64
	for b := 0; b < 8; b++ {
		base := b * 8
		for val := 1; val < 256; val++ {
			lsb := bits.TrailingZeros(uint(val))
			lut[b][val] = lut[b][val&(val-1)] ^ m[base+lsb]
		}
	}

	for r, x := range v {
		if x == 0 {
			continue
		}
		dst[r] ^= lut[0][byte(x)] ^
			lut[1][byte(x>>8)] ^
			lut[2][byte(x>>16)] ^
			lut[3][byte(x>>24)] ^
			lut[4][byte(x>>32)] ^
			lut[5][byte(x>>40)] ^
			lut[6][byte(x>>48)] ^
			lut[7][byte(x>>56)]
	}
}

// innerProd64 computes a^T * b (64 x 64) for R x 64 matrices a and b.
func innerProd64(a, b []uint64) mat64 {
	var res mat64
	for r := range a {
		ar := a[r]
		br := b[r]
		if ar == 0 || br == 0 {
			continue
		}
		for ar != 0 {
			i := bits.TrailingZeros64(ar)
			res[i] ^= br
			ar &= ar - 1
		}
	}
	return res
}

// findWinv picks a maximal invertible principal submatrix of the symmetric t,
// preferring indices not in sPrev, and returns its inverse and index mask.
func findWinv(t mat64, sPrev uint64) (wInv mat64, sMask uint64) {
	m := t
	var w mat64
	for i := 0; i < 64; i++ {
		w[i] = uint64(1) << uint(i)
	}

	// Process indices outside sPrev first, then those in sPrev.
	var order [64]int
	pos := 0
	for i := 0; i < 64; i++ {
		if (sPrev>>uint(i))&1 == 0 {
			order[pos] = i
			pos++
		}
	}
	for i := 0; i < 64; i++ {
		if (sPrev>>uint(i))&1 != 0 {
			order[pos] = i
			pos++
		}
	}

	var used [64]bool
	for oi := 0; oi < 64; oi++ {
		c := order[oi]
		if used[c] {
			continue
		}

		if (m[c]>>uint(c))&1 != 0 {
			// 1x1 pivot on c.
			used[c] = true
			sMask |= uint64(1) << uint(c)
			mc := m[c]
			wc := w[c]
			maskC := uint64(1) << uint(c)
			for r := 0; r < 64; r++ {
				if r != c && (m[r]&maskC) != 0 {
					m[r] ^= mc
					w[r] ^= wc
				}
			}
			continue
		}

		// Search remaining unused index k with m[c][k] == 1 for a 2x2 pivot.
		kFound := -1
		for oj := oi + 1; oj < 64; oj++ {
			k := order[oj]
			if !used[k] && (m[c]>>uint(k))&1 != 0 {
				kFound = k
				break
			}
		}
		if kFound == -1 {
			continue
		}

		k := kFound
		used[c] = true
		used[k] = true
		sMask |= (uint64(1) << uint(c)) | (uint64(1) << uint(k))

		// Normalize 2x2 block [[0, 1], [1, m_kk]] -> [[1, 0], [0, 1]].
		if (m[k]>>uint(k))&1 != 0 {
			m[k] ^= m[c]
			w[k] ^= w[c]
		}
		m[c], m[k] = m[k], m[c]
		w[c], w[k] = w[k], w[c]

		mc, mk := m[c], m[k]
		wc, wk := w[c], w[k]
		maskC := uint64(1) << uint(c)
		maskK := uint64(1) << uint(k)
		for r := 0; r < 64; r++ {
			if r == c || r == k {
				continue
			}
			if (m[r] & maskC) != 0 {
				m[r] ^= mc
				w[r] ^= wc
			}
			if (m[r] & maskK) != 0 {
				m[r] ^= mk
				w[r] ^= wk
			}
		}
	}

	for i := 0; i < 64; i++ {
		if (sMask>>uint(i))&1 != 0 {
			wInv[i] = w[i] & sMask
		}
	}
	return wInv, sMask
}

// blockLanczos runs Montgomery's block Lanczos on A = B·Bᵀ, where B has the
// given rows, and returns subsets of rows that sum to zero mod 2. It finds up
// to 64, less one for each unit of rank(B) − rank(B·Bᵀ), so B should be the
// output of the structured elimination, where that gap was 0 in testing.
func blockLanczos(rows [][]int32, numCols int, seed int64) [][]int {
	R := len(rows)
	C := numCols
	if R == 0 || C == 0 {
		return nil
	}

	// mulBt sets dst = Bᵀ·v and mulB sets dst = B·t.
	mulBt := func(v, dst []uint64) {
		for i := range dst {
			dst[i] = 0
		}
		for r, cols := range rows {
			vr := v[r]
			if vr == 0 {
				continue
			}
			for _, c := range cols {
				dst[c] ^= vr
			}
		}
	}

	mulB := func(t, dst []uint64) {
		for r, cols := range rows {
			var acc uint64
			for _, c := range cols {
				acc ^= t[c]
			}
			dst[r] = acc
		}
	}

	tmpC := make([]uint64, C)
	mulA := func(v, dst []uint64) {
		mulBt(v, tmpC)
		mulB(tmpC, dst)
	}

	rng := rand.New(rand.NewSource(seed))
	y := make([]uint64, R)
	for r := range y {
		y[r] = rng.Uint64()
	}

	v0 := make([]uint64, R)
	mulA(y, v0)

	x := make([]uint64, R)
	vi := append([]uint64(nil), v0...)
	vPrev := make([]uint64, R)
	vPrev2 := make([]uint64, R)
	vNext := make([]uint64, R)
	pi := make([]uint64, R)

	var wPrev, wPrev2 mat64
	var sPrev uint64

	maxIter := R/32 + 64
	for iter := 0; iter < maxIter; iter++ {
		mulA(vi, pi) // P_i = A * V_i
		ti := innerProd64(vi, pi)
		mi := innerProd64(pi, pi)

		wInv, sMask := findWinv(ti, sPrev)
		if sMask == 0 {
			break
		}

		// Accumulate solution: X += V_i * W_i^inv * (V_i^T * V_0)
		vtV0 := innerProd64(vi, v0)
		mulVecMat64(x, vi, mulMat64(wInv, vtV0))

		// If any index dropped in the previous step failed to enter S_i,
		// stop cleanly after the X update before 3-term orthogonality breaks.
		if iter > 0 && ((^sPrev)&^sMask) != 0 {
			break
		}

		// D_i = I + W_i^inv * (M_i * S_i + T_i)
		var ms mat64
		for r := 0; r < 64; r++ {
			ms[r] = mi[r] & sMask
		}
		di := mulMat64(wInv, addMat64(ms, ti))
		for r := 0; r < 64; r++ {
			di[r] ^= uint64(1) << uint(r)
		}

		// V_{i+1} = A*V_i*S_i + V_i*D_i + V_{i-1}*E_i + V_{i-2}*F_i
		for r := 0; r < R; r++ {
			vNext[r] = pi[r] & sMask
		}
		mulVecMat64(vNext, vi, di)

		if iter >= 1 {
			// E_i = W_{i-1}^inv * T_i * S_i
			var tS mat64
			for r := 0; r < 64; r++ {
				tS[r] = ti[r] & sMask
			}
			ei := mulMat64(wPrev, tS)
			mulVecMat64(vNext, vPrev, ei)

			if iter >= 2 {
				// F_i = W_{i-2}^inv * (I - S_{i-1}) * T_i * S_i
				var notPrevTS mat64
				notPrev := ^sPrev
				for r := 0; r < 64; r++ {
					if (notPrev>>uint(r))&1 != 0 {
						notPrevTS[r] = tS[r]
					}
				}
				fi := mulMat64(wPrev2, notPrevTS)
				mulVecMat64(vNext, vPrev2, fi)
			}
		}

		allZero := true
		for r := 0; r < R; r++ {
			if vNext[r] != 0 {
				allZero = false
				break
			}
		}
		if allZero {
			break
		}

		vPrev2, vPrev, vi, vNext = vPrev, vi, vNext, vPrev2
		wPrev2, wPrev = wPrev, wInv
		sPrev = sMask
	}

	// Nullspace extraction from U = [X ^ Y | V_i] (R x 128).
	xy := make([]uint64, R)
	for r := 0; r < R; r++ {
		xy[r] = x[r] ^ y[r]
	}

	btXY := make([]uint64, C)
	btVi := make([]uint64, C)
	mulBt(xy, btXY)
	mulBt(vi, btVi)

	// Transpose B^T * U into 128 rows of C bits each, and eliminate to find ker(B^T * U).
	cWords := (C + 63) / 64
	elimMat := make([][]uint64, 128)
	comb := make([][2]uint64, 128)
	for k := 0; k < 128; k++ {
		elimMat[k] = make([]uint64, cWords)
		if k < 64 {
			comb[k][0] = uint64(1) << uint(k)
		} else {
			comb[k][1] = uint64(1) << uint(k-64)
		}
	}
	for c := 0; c < C; c++ {
		w0 := btXY[c]
		for w0 != 0 {
			k := bits.TrailingZeros64(w0)
			elimMat[k][c>>6] |= 1 << uint(c&63)
			w0 &= w0 - 1
		}
		w1 := btVi[c]
		for w1 != 0 {
			k := bits.TrailingZeros64(w1)
			elimMat[64+k][c>>6] |= 1 << uint(c&63)
			w1 &= w1 - 1
		}
	}

	piv := 0
	for c := 0; c < C && piv < 128; c++ {
		w := c >> 6
		mask := uint64(1) << uint(c&63)
		p := -1
		for r := piv; r < 128; r++ {
			if elimMat[r][w]&mask != 0 {
				p = r
				break
			}
		}
		if p == -1 {
			continue
		}
		elimMat[piv], elimMat[p] = elimMat[p], elimMat[piv]
		comb[piv], comb[p] = comb[p], comb[piv]
		for r := piv + 1; r < 128; r++ {
			if elimMat[r][w]&mask != 0 {
				for j := w; j < cWords; j++ {
					elimMat[r][j] ^= elimMat[piv][j]
				}
				comb[r][0] ^= comb[piv][0]
				comb[r][1] ^= comb[piv][1]
			}
		}
		piv++
	}

	var deps [][]int

	for k := piv; k < 128; k++ {
		d0, d1 := comb[k][0], comb[k][1]
		var dep []int
		for r := 0; r < R; r++ {
			bit := (bits.OnesCount64(xy[r]&d0) + bits.OnesCount64(vi[r]&d1)) & 1
			if bit != 0 {
				dep = append(dep, r)
			}
		}
		if len(dep) > 0 && isDependency(rows, numCols, dep) {
			deps = append(deps, dep)
		}
	}

	return deps
}

// isDependency reports whether the rows in dep sum to zero mod 2.
func isDependency(rows [][]int32, numCols int, dep []int) bool {
	parity := make([]uint64, (numCols+63)/64)
	for _, r := range dep {
		for _, c := range rows[r] {
			parity[c/64] ^= 1 << (uint(c) % 64)
		}
	}
	for _, w := range parity {
		if w != 0 {
			return false
		}
	}
	return true
}

package qs

import (
	"math/bits"
	"math/rand"
)

// This is Montgomery's block Lanczos over GF(2), from "A Block Lanczos
// Algorithm for Finding Dependencies over GF(2)", EUROCRYPT 1995.

// mat64 is a 64×64 GF(2) matrix: bit j of m[i] is entry (i, j).
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

// mulVecMat64 adds v·m to dst, a byte of v at a time.
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

// innerProd64 returns aᵀ·b.
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

// symPrecond is the invertible symmetric P = Lᵀ·L, with L = U1·Perm·L1 for random
// sparse unit triangular L1 and U1 and a random permutation.
type symPrecond struct {
	c    int
	lIdx [][3]int32
	uIdx [][3]int32
	perm []int32
	buf  []uint64
}

func newSymPrecond(c int, rng *rand.Rand) *symPrecond {
	lIdx := make([][3]int32, c)
	uIdx := make([][3]int32, c)
	perm := make([]int32, c)

	for r := 0; r < c; r++ {
		perm[r] = int32(r)
		lIdx[r] = [3]int32{-1, -1, -1}
		uIdx[r] = [3]int32{-1, -1, -1}

		// Pick up to 3 distinct indices in [0, r).
		for k := 0; k < 3 && k < r; k++ {
			for {
				j := int32(rng.Intn(r))
				if j != lIdx[r][0] && j != lIdx[r][1] {
					lIdx[r][k] = j
					break
				}
			}
		}

		// Pick up to 3 distinct indices in [r+1, c).
		rem := c - 1 - r
		for k := 0; k < 3 && k < rem; k++ {
			for {
				j := int32(r + 1 + rng.Intn(rem))
				if j != uIdx[r][0] && j != uIdx[r][1] {
					uIdx[r][k] = j
					break
				}
			}
		}
	}

	rng.Shuffle(c, func(i, j int) {
		perm[i], perm[j] = perm[j], perm[i]
	})

	return &symPrecond{
		c:    c,
		lIdx: lIdx,
		uIdx: uIdx,
		perm: perm,
		buf:  make([]uint64, c),
	}
}

// apply computes v <- P * v = L1^T * Perm^T * U1^T * U1 * Perm * L1 * v.
func (p *symPrecond) apply(v []uint64) {
	c := p.c
	buf := p.buf

	// 1. v <- L1 * v (iterate r from c-1 down to 0)
	for r := c - 1; r >= 0; r-- {
		e := p.lIdx[r]
		if e[0] >= 0 {
			v[r] ^= v[e[0]]
			if e[1] >= 0 {
				v[r] ^= v[e[1]]
				if e[2] >= 0 {
					v[r] ^= v[e[2]]
				}
			}
		}
	}

	// 2. buf <- Perm * v
	for r := 0; r < c; r++ {
		buf[r] = v[p.perm[r]]
	}

	// 3. buf <- U1 * buf (iterate r from 0 up to c-1)
	for r := 0; r < c; r++ {
		e := p.uIdx[r]
		if e[0] >= 0 {
			buf[r] ^= buf[e[0]]
			if e[1] >= 0 {
				buf[r] ^= buf[e[1]]
				if e[2] >= 0 {
					buf[r] ^= buf[e[2]]
				}
			}
		}
	}

	// 4. buf <- U1^T * buf (iterate r from c-1 down to 0)
	for r := c - 1; r >= 0; r-- {
		br := buf[r]
		if br == 0 {
			continue
		}
		e := p.uIdx[r]
		if e[0] >= 0 {
			buf[e[0]] ^= br
			if e[1] >= 0 {
				buf[e[1]] ^= br
				if e[2] >= 0 {
					buf[e[2]] ^= br
				}
			}
		}
	}

	// 5. v <- Perm^T * buf
	for r := 0; r < c; r++ {
		v[p.perm[r]] = buf[r]
	}

	// 6. v <- L1^T * v (iterate r from 0 up to c-1)
	for r := 0; r < c; r++ {
		vr := v[r]
		if vr == 0 {
			continue
		}
		e := p.lIdx[r]
		if e[0] >= 0 {
			v[e[0]] ^= vr
			if e[1] >= 0 {
				v[e[1]] ^= vr
				if e[2] >= 0 {
					v[e[2]] ^= vr
				}
			}
		}
	}
}

// blockLanczos returns subsets of rows that sum to zero mod 2, up to 64. It
// uses A = B·P·Bᵀ, since with A = B·Bᵀ each unit of rank(B) − rank(B·Bᵀ) costs one.
func blockLanczos(rows [][]int32, numCols int, seed int64) [][]int {
	return blockLanczosCheck(rows, numCols, seed, nil)
}

// blockLanczosCheck is blockLanczos, calling onF with each F matrix and
// Montgomery's form of it, for testing.

func blockLanczosCheck(rows [][]int32, numCols int, seed int64, onF func(iter int, f, fMont mat64)) [][]int {
	R := len(rows)
	C := numCols
	if R == 0 || C == 0 {
		return nil
	}

	rng := rand.New(rand.NewSource(seed))
	pre := newSymPrecond(C, rng)

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
		pre.apply(tmpC)
		mulB(tmpC, dst)
	}

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

	var wPrev, wPrev2, tPrev, mPrev mat64
	var sPrev uint64

	maxIter := R/32 + 64
	for iter := 0; iter < maxIter; iter++ {
		mulA(vi, pi)
		ti := innerProd64(vi, pi)
		mi := innerProd64(pi, pi)

		wInv, sMask := findWinv(ti, sPrev)
		if sMask == 0 {
			break
		}

		vtV0 := innerProd64(vi, v0)
		mulVecMat64(x, vi, mulMat64(wInv, vtV0))

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

		for r := 0; r < R; r++ {
			vNext[r] = pi[r] & sMask
		}
		mulVecMat64(vNext, vi, di)

		if iter >= 1 {
			var tS mat64
			for r := 0; r < 64; r++ {
				tS[r] = ti[r] & sMask
			}
			ei := mulMat64(wPrev, tS)
			mulVecMat64(vNext, vPrev, ei)

			if iter >= 2 {
				var notPrevTS mat64
				notPrev := ^sPrev
				for r := 0; r < 64; r++ {
					if (notPrev>>uint(r))&1 != 0 {
						notPrevTS[r] = tS[r]
					}
				}
				fi := mulMat64(wPrev2, notPrevTS)

				if onF != nil {
					// Montgomery's form: W_{i-2}^inv * (I + T_{i-1}*W_{i-1}^inv) * (M_{i-1}*S_{i-1} + T_{i-1}) * S_i
					var mPrevS mat64
					for r := 0; r < 64; r++ {
						mPrevS[r] = mPrev[r] & sPrev
					}
					rhs := addMat64(mPrevS, tPrev)
					for r := 0; r < 64; r++ {
						rhs[r] &= sMask
					}
					tw := mulMat64(tPrev, wPrev)
					for r := 0; r < 64; r++ {
						tw[r] ^= uint64(1) << uint(r)
					}
					fMont := mulMat64(wPrev2, mulMat64(tw, rhs))
					onF(iter, fi, fMont)
				}

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
		tPrev, mPrev = ti, mi
		sPrev = sMask
	}

	// Final extraction against UNPRECONDITIONED B^T:
	// compute ker(B^T * [X ^ Y | V_i]).
	xy := make([]uint64, R)
	for r := 0; r < R; r++ {
		xy[r] = x[r] ^ y[r]
	}

	btXY := make([]uint64, C)
	btVi := make([]uint64, C)
	mulBt(xy, btXY)
	mulBt(vi, btVi)

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

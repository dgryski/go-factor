package linalg

import (
	"math/bits"
	"math/rand"
	"runtime"
	"sync"
)

// This is Montgomery's block Lanczos over GF(2), from "A Block Lanczos
// Algorithm for Finding Dependencies over GF(2)", EUROCRYPT 1995.

// Mat64 is a 64×64 GF(2) matrix: bit j of m[i] is entry (i, j).
type Mat64 [64]uint64

// MulMat64 multiplies two 64×64 GF(2) matrices.
func MulMat64(a, b Mat64) Mat64 {
	var out Mat64
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

// AddMat64 adds two 64×64 GF(2) matrices.
func AddMat64(a, b Mat64) Mat64 {
	var out Mat64
	for i := 0; i < 64; i++ {
		out[i] = a[i] ^ b[i]
	}
	return out
}

func buildLUT64(m Mat64) [8][256]uint64 {
	var lut [8][256]uint64
	for b := 0; b < 8; b++ {
		base := b * 8
		for val := 1; val < 256; val++ {
			lsb := bits.TrailingZeros(uint(val))
			lut[b][val] = lut[b][val&(val-1)] ^ m[base+lsb]
		}
	}
	return lut
}

// mulVecMat64Range adds v[start:end]·m to dst[start:end].
func mulVecMat64Range(dst, v []uint64, lut *[8][256]uint64, start, end int) {
	for r := start; r < end; r++ {
		x := v[r]
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

// MulVecMat64 adds v·m to dst, a byte of v at a time.
func MulVecMat64(dst, v []uint64, m Mat64) {
	lut := buildLUT64(m)
	mulVecMat64Range(dst, v, &lut, 0, len(v))
}

func innerProd64Range(a, b []uint64, start, end int) Mat64 {
	var res Mat64
	for r := start; r < end; r++ {
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

// InnerProd64 returns aᵀ·b.
func InnerProd64(a, b []uint64) Mat64 {
	return innerProd64Range(a, b, 0, len(a))
}

// FindWinv picks a maximal invertible principal submatrix of the symmetric t,
// preferring indices not in sPrev, and returns its inverse and index mask.
func FindWinv(t Mat64, sPrev uint64) (wInv Mat64, sMask uint64) {
	m := t
	var w Mat64
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

// SymPrecond is the invertible symmetric P = Lᵀ·L, with L = U1·Perm·L1 for random
// sparse unit triangular L1 and U1 and a random permutation.
type SymPrecond struct {
	c    int
	lIdx [][3]int32
	uIdx [][3]int32
	perm []int32
	buf  []uint64
}

// NewSymPrecond constructs a random symmetric preconditioner for c columns.
func NewSymPrecond(c int, rng *rand.Rand) *SymPrecond {
	lIdx := make([][3]int32, c)
	uIdx := make([][3]int32, c)
	perm := make([]int32, c)

	for r := 0; r < c; r++ {
		perm[r] = int32(r)
		lIdx[r] = [3]int32{-1, -1, -1}
		uIdx[r] = [3]int32{-1, -1, -1}

		for k := 0; k < 3 && k < r; k++ {
			for {
				j := int32(rng.Intn(r))
				if j != lIdx[r][0] && j != lIdx[r][1] {
					lIdx[r][k] = j
					break
				}
			}
		}

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

	return &SymPrecond{
		c:    c,
		lIdx: lIdx,
		uIdx: uIdx,
		perm: perm,
		buf:  make([]uint64, c),
	}
}

// Apply computes v <- P * v = L1^T * Perm^T * U1^T * U1 * Perm * L1 * v.
func (p *SymPrecond) Apply(v []uint64) {
	c := p.c
	buf := p.buf

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

	for r := 0; r < c; r++ {
		buf[r] = v[p.perm[r]]
	}

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

	for r := 0; r < c; r++ {
		v[p.perm[r]] = buf[r]
	}

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

type parTask struct {
	n  int
	fn func(wid, start, end int)
}

type lanczosPool struct {
	nw  int
	chs []chan parTask
	wg  sync.WaitGroup
}

func newLanczosPool(nw int) *lanczosPool {
	p := &lanczosPool{
		nw:  nw,
		chs: make([]chan parTask, nw),
	}
	for w := 0; w < nw; w++ {
		ch := make(chan parTask, 1)
		p.chs[w] = ch
		wid := w
		go func() {
			for t := range ch {
				start := wid * t.n / nw
				end := (wid + 1) * t.n / nw
				t.fn(wid, start, end)
				p.wg.Done()
			}
		}()
	}
	return p
}

func (p *lanczosPool) run(n int, fn func(wid, start, end int)) {
	p.wg.Add(p.nw)
	t := parTask{n: n, fn: fn}
	for _, ch := range p.chs {
		ch <- t
	}
	p.wg.Wait()
}

func (p *lanczosPool) close() {
	for _, ch := range p.chs {
		close(ch)
	}
}

// BlockLanczos returns subsets of rows that sum to zero mod 2, up to 64. It
// uses A = B·P·Bᵀ, since with A = B·Bᵀ each unit of rank(B) − rank(B·Bᵀ) costs one.
func BlockLanczos(rows [][]int32, numCols int, seed int64) [][]int {
	return BlockLanczosCheck(rows, numCols, seed, nil)
}

// BlockLanczosCheck is BlockLanczos, calling onF with each F matrix and
// Montgomery's form of it, for testing.
func BlockLanczosCheck(rows [][]int32, numCols int, seed int64, onF func(iter int, f, fMont Mat64)) [][]int {
	R := len(rows)
	C := numCols
	if R == 0 || C == 0 {
		return nil
	}

	rng := rand.New(rand.NewSource(seed))
	pre := NewSymPrecond(C, rng)

	nw := runtime.GOMAXPROCS(0)
	if R < 4096 || nw <= 1 {
		nw = 1
	} else if nw > 16 {
		nw = 16
	}

	var pool *lanczosPool
	var localC [][]uint64
	var localM []Mat64
	if nw > 1 {
		pool = newLanczosPool(nw)
		defer pool.close()
		localC = make([][]uint64, nw)
		buf := make([]uint64, nw*C)
		for w := 0; w < nw; w++ {
			localC[w] = buf[w*C : (w+1)*C]
		}
		localM = make([]Mat64, nw)
	}

	mulBt := func(v, dst []uint64) {
		if nw == 1 {
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
			return
		}
		pool.run(R, func(wid, start, end int) {
			lc := localC[wid]
			for i := range lc {
				lc[i] = 0
			}
			for r := start; r < end; r++ {
				vr := v[r]
				if vr == 0 {
					continue
				}
				for _, c := range rows[r] {
					lc[c] ^= vr
				}
			}
		})
		pool.run(C, func(_ int, start, end int) {
			copy(dst[start:end], localC[0][start:end])
			for w := 1; w < nw; w++ {
				lc := localC[w]
				for c := start; c < end; c++ {
					dst[c] ^= lc[c]
				}
			}
		})
	}

	mulB := func(t, dst []uint64) {
		if nw == 1 {
			for r, cols := range rows {
				var acc uint64
				for _, c := range cols {
					acc ^= t[c]
				}
				dst[r] = acc
			}
			return
		}
		pool.run(R, func(_ int, start, end int) {
			for r := start; r < end; r++ {
				var acc uint64
				for _, c := range rows[r] {
					acc ^= t[c]
				}
				dst[r] = acc
			}
		})
	}

	innerProd := func(a, b []uint64) Mat64 {
		if nw == 1 {
			return innerProd64Range(a, b, 0, R)
		}
		pool.run(R, func(wid, start, end int) {
			localM[wid] = innerProd64Range(a, b, start, end)
		})
		res := localM[0]
		for w := 1; w < nw; w++ {
			for i := 0; i < 64; i++ {
				res[i] ^= localM[w][i]
			}
		}
		return res
	}

	mulVecMat := func(dst, v []uint64, m Mat64) {
		lut := buildLUT64(m)
		if nw == 1 {
			mulVecMat64Range(dst, v, &lut, 0, R)
			return
		}
		pool.run(R, func(_ int, start, end int) {
			mulVecMat64Range(dst, v, &lut, start, end)
		})
	}

	tmpC := make([]uint64, C)
	mulA := func(v, dst []uint64) {
		mulBt(v, tmpC)
		pre.Apply(tmpC)
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

	var wPrev, wPrev2, tPrev, mPrev Mat64
	var sPrev uint64

	maxIter := R/32 + 64
	for iter := 0; iter < maxIter; iter++ {
		mulA(vi, pi)
		ti := innerProd(vi, pi)
		mi := innerProd(pi, pi)

		wInv, sMask := FindWinv(ti, sPrev)
		if sMask == 0 {
			break
		}

		vtV0 := innerProd(vi, v0)
		mulVecMat(x, vi, MulMat64(wInv, vtV0))

		if iter > 0 && ((^sPrev)&^sMask) != 0 {
			break
		}

		var ms Mat64
		for r := 0; r < 64; r++ {
			ms[r] = mi[r] & sMask
		}
		di := MulMat64(wInv, AddMat64(ms, ti))
		for r := 0; r < 64; r++ {
			di[r] ^= uint64(1) << uint(r)
		}

		for r := 0; r < R; r++ {
			vNext[r] = pi[r] & sMask
		}
		mulVecMat(vNext, vi, di)

		if iter >= 1 {
			var tS Mat64
			for r := 0; r < 64; r++ {
				tS[r] = ti[r] & sMask
			}
			ei := MulMat64(wPrev, tS)
			mulVecMat(vNext, vPrev, ei)

			if iter >= 2 {
				var notPrevTS Mat64
				notPrev := ^sPrev
				for r := 0; r < 64; r++ {
					if (notPrev>>uint(r))&1 != 0 {
						notPrevTS[r] = tS[r]
					}
				}
				fi := MulMat64(wPrev2, notPrevTS)

				if onF != nil {
					var mPrevS Mat64
					for r := 0; r < 64; r++ {
						mPrevS[r] = mPrev[r] & sPrev
					}
					rhs := AddMat64(mPrevS, tPrev)
					for r := 0; r < 64; r++ {
						rhs[r] &= sMask
					}
					tw := MulMat64(tPrev, wPrev)
					for r := 0; r < 64; r++ {
						tw[r] ^= uint64(1) << uint(r)
					}
					fMont := MulMat64(wPrev2, MulMat64(tw, rhs))
					onF(iter, fi, fMont)
				}

				mulVecMat(vNext, vPrev2, fi)
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
	seen := make(map[uint64]bool)
	for k := piv; k < 128; k++ {
		d0, d1 := comb[k][0], comb[k][1]
		var dep []int
		h := uint64(14695981039346656037)
		for r := 0; r < R; r++ {
			bit := (bits.OnesCount64(xy[r]&d0) + bits.OnesCount64(vi[r]&d1)) & 1
			if bit != 0 {
				dep = append(dep, r)
				for v, i := uint32(r), 0; i < 4; v, i = v>>8, i+1 {
					h = (h ^ uint64(v&0xff)) * 1099511628211
				}
			}
		}
		if len(dep) == 0 || seen[h] {
			continue
		}
		seen[h] = true
		if IsDependency(rows, numCols, dep) {
			deps = append(deps, dep)
		}
	}

	return deps
}

// IsDependency reports whether the rows in dep sum to zero mod 2.
func IsDependency(rows [][]int32, numCols int, dep []int) bool {
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

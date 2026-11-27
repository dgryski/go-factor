package gnfs

import (
	"math/bits"
	"sync"

	"github.com/dgryski/go-factor/internal/arith"
)

type u128 struct {
	hi, lo uint64
}

func (x u128) isZero() bool { return x.hi == 0 && x.lo == 0 }

func (x u128) bitLen() int {
	if x.hi != 0 {
		return 64 + bits.Len64(x.hi)
	}
	return bits.Len64(x.lo)
}

func cmp128(a, b u128) int {
	if a.hi != b.hi {
		if a.hi < b.hi {
			return -1
		}
		return 1
	}
	if a.lo != b.lo {
		if a.lo < b.lo {
			return -1
		}
		return 1
	}
	return 0
}

func sub128(a, b u128) u128 {
	lo, br := bits.Sub64(a.lo, b.lo, 0)
	hi, _ := bits.Sub64(a.hi, b.hi, br)
	return u128{hi: hi, lo: lo}
}

func absDiff128(a, b u128) u128 {
	if cmp128(a, b) >= 0 {
		return sub128(a, b)
	}
	return sub128(b, a)
}

func ctz128(x u128) int {
	if x.lo != 0 {
		return bits.TrailingZeros64(x.lo)
	}
	return 64 + bits.TrailingZeros64(x.hi)
}

func rsh128(x u128, s uint) u128 {
	if s == 0 {
		return x
	}
	if s >= 64 {
		return u128{hi: 0, lo: x.hi >> (s - 64)}
	}
	return u128{hi: x.hi >> s, lo: (x.lo >> s) | (x.hi << (64 - s))}
}

func gcd128(a, b u128) u128 {
	if a.isZero() {
		return b
	}
	if b.isZero() {
		return a
	}
	shift := uint(min(ctz128(a), ctz128(b)))
	a = rsh128(a, uint(ctz128(a)))
	for !b.isZero() {
		if a.hi == 0 && b.hi == 0 {
			g := arith.Gcd64(a.lo, b.lo)
			if shift == 0 {
				return u128{lo: g}
			}
			if shift >= 64 {
				return u128{hi: g << (shift - 64)}
			}
			return u128{hi: g >> (64 - shift), lo: g << shift}
		}
		b = rsh128(b, uint(ctz128(b)))
		if cmp128(a, b) > 0 {
			a, b = b, a
		}
		b = sub128(b, a)
	}
	if shift == 0 {
		return a
	}
	if shift >= 64 {
		return u128{hi: a.lo << (shift - 64)}
	}
	return u128{hi: (a.hi << shift) | (a.lo >> (64 - shift)), lo: a.lo << shift}
}

func exactDiv128(n, d u128) u128 {
	if d.hi == 0 {
		qh, ql, _ := divRem128(n.hi, n.lo, d.lo)
		return u128{hi: qh, lo: ql}
	}
	s := uint(ctz128(d))
	if s > 0 {
		n = rsh128(n, s)
		d = rsh128(d, s)
		if d.hi == 0 {
			qh, ql, _ := divRem128(n.hi, n.lo, d.lo)
			return u128{hi: qh, lo: ql}
		}
	}
	inv := d.lo
	for i := 0; i < 5; i++ {
		inv *= 2 - d.lo*inv
	}
	return u128{hi: 0, lo: n.lo * inv}
}

// mont128 performs Montgomery arithmetic mod an odd 65 to 128 bit integer n.
type mont128 struct {
	n    u128
	nInv uint64
	one  u128
	r2   u128
}

func newMont128(n u128) mont128 {
	inv := n.lo
	for i := 0; i < 5; i++ {
		inv *= 2 - n.lo*inv
	}
	var m mont128
	m.n = n
	m.nInv = -inv

	bl := n.bitLen()
	if bl == 128 {
		lo, br := bits.Sub64(0, n.lo, 0)
		hi, _ := bits.Sub64(0, n.hi, br)
		m.one = u128{hi: hi, lo: lo}
	} else {
		x := u128{hi: uint64(1) << uint(bl-65), lo: 0}
		for k := bl - 1; k < 128; k++ {
			x = m.add(x, x)
		}
		m.one = x
	}

	z := m.add(m.one, m.one)
	for i := 0; i < 7; i++ {
		z = m.mul(z, z)
	}
	m.r2 = z
	return m
}

func (m *mont128) add(a, b u128) u128 {
	lo, c0 := bits.Add64(a.lo, b.lo, 0)
	hi, c1 := bits.Add64(a.hi, b.hi, c0)
	if c1 != 0 || hi > m.n.hi || (hi == m.n.hi && lo >= m.n.lo) {
		var b0 uint64
		lo, b0 = bits.Sub64(lo, m.n.lo, 0)
		hi, _ = bits.Sub64(hi, m.n.hi, b0)
	}
	return u128{hi: hi, lo: lo}
}

func (m *mont128) sub(a, b u128) u128 {
	lo, b0 := bits.Sub64(a.lo, b.lo, 0)
	hi, b1 := bits.Sub64(a.hi, b.hi, b0)
	if b1 != 0 {
		var c0 uint64
		lo, c0 = bits.Add64(lo, m.n.lo, 0)
		hi, _ = bits.Add64(hi, m.n.hi, c0)
	}
	return u128{hi: hi, lo: lo}
}

// mul computes a * b * 2^-128 mod n via two limb CIOS reduction.
func (m *mont128) mul(a, b u128) u128 {
	h0, t0 := bits.Mul64(a.lo, b.lo)
	h1, l1 := bits.Mul64(a.lo, b.hi)
	t1, c := bits.Add64(l1, h0, 0)
	t2 := h1 + c

	q0 := t0 * m.nInv
	mh0, ml0 := bits.Mul64(q0, m.n.lo)
	_, c0 := bits.Add64(t0, ml0, 0)
	mh1, ml1 := bits.Mul64(q0, m.n.hi)
	t0p, c1 := bits.Add64(t1, mh0, c0)
	t1p, c2 := bits.Add64(t2, mh1, c1)
	u0, d0 := bits.Add64(t0p, ml1, 0)
	u1, d1 := bits.Add64(t1p, 0, d0)
	u2 := c2 + d1

	ah0, al0 := bits.Mul64(a.hi, b.lo)
	ah1, al1 := bits.Mul64(a.hi, b.hi)
	v0, vc0 := bits.Add64(u0, al0, 0)
	s1, sc1 := bits.Add64(al1, ah0, 0)
	s2 := ah1 + sc1
	v1, vc2 := bits.Add64(u1, s1, vc0)
	v2, v3 := bits.Add64(u2, s2, vc2)

	q1 := v0 * m.nInv
	nh0, nl0 := bits.Mul64(q1, m.n.lo)
	_, nd0 := bits.Add64(v0, nl0, 0)
	nh1, nl1 := bits.Mul64(q1, m.n.hi)
	w0, nd1 := bits.Add64(v1, nh0, nd0)
	w1, nd2 := bits.Add64(v2, nh1, nd1)
	w2 := v3 + nd2
	r0, e0 := bits.Add64(w0, nl1, 0)
	r1, e1 := bits.Add64(w1, 0, e0)
	r2 := w2 + e1

	if r2 != 0 || r1 > m.n.hi || (r1 == m.n.hi && r0 >= m.n.lo) {
		var b0 uint64
		r0, b0 = bits.Sub64(r0, m.n.lo, 0)
		r1, _ = bits.Sub64(r1, m.n.hi, b0)
	}
	return u128{hi: r1, lo: r0}
}

func (m *mont128) to64(a uint64) u128 {
	return m.mul(u128{hi: 0, lo: a}, m.r2)
}

func (m *mont128) to(a u128) u128 {
	if cmp128(a, m.n) >= 0 {
		a = sub128(a, m.n)
	}
	return m.mul(a, m.r2)
}

func (m *mont128) pow64(base u128, exp uint64) u128 {
	res := m.one
	for ; exp > 0; exp >>= 1 {
		if exp&1 != 0 {
			res = m.mul(res, base)
		}
		base = m.mul(base, base)
	}
	return res
}

func (m *mont128) pow128(base, exp u128) u128 {
	res := m.one
	for e := exp.lo; e > 0; e >>= 1 {
		if e&1 != 0 {
			res = m.mul(res, base)
		}
		base = m.mul(base, base)
	}
	if exp.hi != 0 {
		// Advance base across any remaining zero bits of exp.lo.
		used := bits.Len64(exp.lo)
		for k := used; k < 64; k++ {
			base = m.mul(base, base)
		}
		for e := exp.hi; e > 0; e >>= 1 {
			if e&1 != 0 {
				res = m.mul(res, base)
			}
			base = m.mul(base, base)
		}
	}
	return res
}

func isPrime128(n u128) bool {
	if n.hi == 0 {
		return arith.IsPrime64(n.lo)
	}
	if n.lo&1 == 0 {
		return false
	}
	for _, p := range [...]uint64{3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37} {
		_, _, r := divRem128(n.hi, n.lo, p)
		if r == 0 {
			return false
		}
	}
	m := newMont128(n)
	minusOne := sub128(n, m.one)
	d := u128{hi: n.hi, lo: n.lo - 1}
	s := ctz128(d)
	d = rsh128(d, uint(s))
	for _, a := range [...]uint64{2, 3, 5, 7, 11} {
		x := m.pow128(m.to64(a), d)
		if x == m.one || x == minusOne {
			continue
		}
		composite := true
		for r := 1; r < s; r++ {
			x = m.mul(x, x)
			if x == minusOne {
				composite = false
				break
			}
		}
		if composite {
			return false
		}
	}
	return true
}

var (
	cof128Once sync.Once
	p1Chunks   []uint64
	ecmPP250   []uint64
	ecmPP800   []uint64
	ecmPP2000  []uint64
)

func initCof128Tables() {
	buildPP := func(b1 uint64) []uint64 {
		var out []uint64
		for _, p64 := range arith.PrimesUpTo(int64(b1)) {
			p := uint64(p64)
			pk := p
			for pk <= b1/p {
				pk *= p
			}
			out = append(out, pk)
		}
		return out
	}
	ecmPP250 = buildPP(250)
	ecmPP800 = buildPP(800)
	ecmPP2000 = buildPP(2000)

	cur := uint64(1)
	for _, pk := range buildPP(300) {
		if cur > (^uint64(0))/pk {
			p1Chunks = append(p1Chunks, cur)
			cur = pk
		} else {
			cur *= pk
		}
	}
	if cur > 1 {
		p1Chunks = append(p1Chunks, cur)
	}
}

func pMinusOne128(m *mont128) u128 {
	cof128Once.Do(initCof128Tables)
	x := m.add(m.one, m.one)
	for i, e := range p1Chunks {
		x = m.pow64(x, e)
		if (i+1)&3 == 0 || i+1 == len(p1Chunks) {
			g := gcd128(m.sub(x, m.one), m.n)
			if cmp128(g, u128{lo: 1}) > 0 {
				if g == m.n {
					return u128{}
				}
				return g
			}
		}
	}
	return u128{}
}

type xzPoint128 struct {
	x, z u128
}

func xDBL128(m *mont128, p xzPoint128, a24, c24 u128) xzPoint128 {
	s := m.add(p.x, p.z)
	s = m.mul(s, s)
	d := m.sub(p.x, p.z)
	d = m.mul(d, d)
	e := m.sub(s, d)
	t1 := m.mul(c24, d)
	t2 := m.mul(a24, e)
	return xzPoint128{
		x: m.mul(c24, m.mul(s, d)),
		z: m.mul(e, m.add(t1, t2)),
	}
}

func xADD128(m *mont128, p, q, diff xzPoint128) xzPoint128 {
	u := m.mul(m.sub(p.x, p.z), m.add(q.x, q.z))
	v := m.mul(m.add(p.x, p.z), m.sub(q.x, q.z))
	upv := m.add(u, v)
	umv := m.sub(u, v)
	return xzPoint128{
		x: m.mul(diff.z, m.mul(upv, upv)),
		z: m.mul(diff.x, m.mul(umv, umv)),
	}
}

func xLadder128(m *mont128, p xzPoint128, k uint64, a24, c24 u128) xzPoint128 {
	if k <= 1 {
		return p
	}
	r0 := p
	r1 := xDBL128(m, p, a24, c24)
	top := 63 - bits.LeadingZeros64(k)
	for bit := top - 1; bit >= 0; bit-- {
		if (k>>uint(bit))&1 != 0 {
			r0 = xADD128(m, r0, r1, p)
			r1 = xDBL128(m, r1, a24, c24)
		} else {
			r1 = xADD128(m, r0, r1, p)
			r0 = xDBL128(m, r0, a24, c24)
		}
	}
	return r0
}

// ecmCurve128 runs Stage 1 on Suyama's Montgomery curve for sigma >= 6
// (Montgomery, Math. Comp. 48, 1987, Section 10.3.1).
func ecmCurve128(m *mont128, sigma uint64, pps []uint64) u128 {
	sig := m.to64(sigma)
	five := m.to64(5)
	u := m.sub(m.mul(sig, sig), five)
	v := m.add(sig, sig)
	v = m.add(v, v)

	u2 := m.mul(u, u)
	u3 := m.mul(u2, u)
	v2 := m.mul(v, v)
	v3 := m.mul(v2, v)

	vmu := m.sub(v, u)
	vmu3 := m.mul(m.mul(vmu, vmu), vmu)
	u3pv := m.add(m.add(u, m.add(u, u)), v)
	a24 := m.mul(vmu3, u3pv)

	u3v := m.mul(u3, v)
	u3v4 := m.add(u3v, u3v)
	u3v4 = m.add(u3v4, u3v4)
	c24 := m.add(u3v4, u3v4)
	c24 = m.add(c24, c24)

	if g := gcd128(c24, m.n); cmp128(g, u128{lo: 1}) > 0 {
		if g != m.n {
			return g
		}
		return u128{}
	}

	pt := xzPoint128{x: u3, z: v3}
	for i, k := range pps {
		pt = xLadder128(m, pt, k, a24, c24)
		if (i+1)&31 == 0 || i+1 == len(pps) {
			g := gcd128(pt.z, m.n)
			if cmp128(g, u128{lo: 1}) > 0 {
				if g != m.n {
					return g
				}
				return u128{}
			}
		}
	}
	return u128{}
}

func ecm128(m *mont128) u128 {
	cof128Once.Do(initCof128Tables)
	for sig := uint64(6); sig < 10; sig++ {
		if g := ecmCurve128(m, sig, ecmPP250); !g.isZero() {
			return g
		}
	}
	for sig := uint64(10); sig < 16; sig++ {
		if g := ecmCurve128(m, sig, ecmPP800); !g.isZero() {
			return g
		}
	}
	for sig := uint64(16); sig < 24; sig++ {
		if g := ecmCurve128(m, sig, ecmPP2000); !g.isZero() {
			return g
		}
	}
	return u128{}
}

func rho128(m *mont128) u128 {
	one := u128{lo: 1}
	for c64 := uint64(1); c64 <= 8; c64++ {
		c := m.to64(c64)
		step := func(v u128) u128 {
			return m.add(m.mul(v, v), c)
		}
		y := m.to64(2)
		q := m.one
		g := one
		r := uint64(1)
		var x, ys u128
		for g == one && r <= 1<<15 {
			x = y
			for i := uint64(0); i < r; i++ {
				y = step(y)
			}
			for k := uint64(0); k < r && g == one; k += 64 {
				ys = y
				for i := uint64(0); i < 64 && k+i < r; i++ {
					y = step(y)
					q = m.mul(q, absDiff128(x, y))
				}
				g = gcd128(q, m.n)
			}
			r <<= 1
		}
		if g == m.n {
			for {
				ys = step(ys)
				g = gcd128(absDiff128(x, ys), m.n)
				if cmp128(g, one) > 0 {
					break
				}
			}
		}
		if cmp128(g, one) > 0 && g != m.n {
			return g
		}
	}
	return u128{}
}

func findFactor128(n u128) u128 {
	if n.hi == 0 {
		return u128{lo: arith.Rho64(n.lo)}
	}
	if n.lo&1 == 0 {
		return u128{lo: 2}
	}
	m := newMont128(n)
	if g := pMinusOne128(&m); !g.isZero() {
		return g
	}
	if g := ecm128(&m); !g.isZero() {
		return g
	}
	return rho128(&m)
}

func sortLPs3(out *[3]uint64, n int) {
	if n >= 2 && out[0] > out[1] {
		out[0], out[1] = out[1], out[0]
	}
	if n >= 3 {
		if out[1] > out[2] {
			out[1], out[2] = out[2], out[1]
		}
		if out[0] > out[1] {
			out[0], out[1] = out[1], out[0]
		}
	}
}

func collectPrimes64(v, pMaxSq, lpBound uint64, out *[3]uint64, count *int, maxLP int) bool {
	if v == 1 {
		return true
	}
	if *count >= maxLP {
		return false
	}
	if v < pMaxSq || arith.IsPrime64(v) {
		if v > lpBound {
			return false
		}
		out[*count] = v
		*count++
		return true
	}
	if *count+2 > maxLP {
		return false
	}
	if *count+2 == maxLP {
		lp1, lp2, ok := arith.SplitCofactor(v, 0, lpBound)
		if !ok || !arith.IsPrime64(lp1) || !arith.IsPrime64(lp2) {
			return false
		}
		out[*count] = lp1
		out[*count+1] = lp2
		*count += 2
		return true
	}
	d := arith.Rho64(v)
	if d == 0 {
		return false
	}
	return collectPrimes64(d, pMaxSq, lpBound, out, count, maxLP) &&
		collectPrimes64(v/d, pMaxSq, lpBound, out, count, maxLP)
}

func collectPrimes128(v u128, pMaxSq, lpBound uint64, out *[3]uint64, count *int, maxLP int) bool {
	if v.hi == 0 {
		return collectPrimes64(v.lo, pMaxSq, lpBound, out, count, maxLP)
	}
	if *count+2 > maxLP {
		return false
	}
	if isPrime128(v) {
		return false
	}
	d := findFactor128(v)
	if d.isZero() {
		return false
	}
	q := exactDiv128(v, d)
	if cmp128(d, q) > 0 {
		d, q = q, d
	}
	return collectPrimes128(d, pMaxSq, lpBound, out, count, maxLP) &&
		collectPrimes128(q, pMaxSq, lpBound, out, count, maxLP)
}

// splitNormCofactor splits a trial divided normState into up to maxLP (1..3) primes <= lpBound.
func splitNormCofactor(ns *normState, pMax, lpBound uint64, maxLP, mfb int) ([3]uint64, bool) {
	if !ns.use128 {
		return [3]uint64{}, false
	}
	if ns.hi == 0 {
		c := ns.lo
		if c == 1 {
			return [3]uint64{}, true
		}
		if mfb > 0 && bits.Len64(c) > mfb {
			return [3]uint64{}, false
		}
		var pMaxSq uint64
		if pMax <= 1<<32-1 {
			pMaxSq = pMax * pMax
		} else {
			pMaxSq = ^uint64(0)
		}
		if c <= lpBound && (c < pMaxSq || arith.IsPrime64(c)) {
			return [3]uint64{c, 0, 0}, true
		}
		if maxLP < 2 {
			return [3]uint64{}, false
		}
		if maxLP == 2 {
			lp1, lp2, ok := arith.SplitCofactor(c, pMax, lpBound)
			if !ok {
				return [3]uint64{}, false
			}
			return [3]uint64{lp1, lp2, 0}, true
		}
		var out [3]uint64
		var count int
		if !collectPrimes64(c, pMaxSq, lpBound, &out, &count, maxLP) {
			return [3]uint64{}, false
		}
		sortLPs3(&out, count)
		return out, true
	}

	if maxLP < 2 {
		return [3]uint64{}, false
	}
	v := u128{hi: ns.hi, lo: ns.lo}
	if mfb > 0 && v.bitLen() > mfb {
		return [3]uint64{}, false
	}
	if maxLP == 2 {
		sqHi, sqLo := bits.Mul64(lpBound, lpBound)
		if cmp128(v, u128{hi: sqHi, lo: sqLo}) > 0 {
			return [3]uint64{}, false
		}
	}
	var pMaxSq uint64
	if pMax <= 1<<32-1 {
		pMaxSq = pMax * pMax
	} else {
		pMaxSq = ^uint64(0)
	}
	var out [3]uint64
	var count int
	if !collectPrimes128(v, pMaxSq, lpBound, &out, &count, maxLP) {
		return [3]uint64{}, false
	}
	sortLPs3(&out, count)
	return out, true
}

func factorU128Into(v u128, dst *[]uint64) bool {
	if v.hi == 0 {
		return factorUint64Into(v.lo, dst)
	}
	if isPrime128(v) {
		return false
	}
	d := findFactor128(v)
	if d.isZero() {
		return false
	}
	q := exactDiv128(v, d)
	return factorU128Into(d, dst) && factorU128Into(q, dst)
}

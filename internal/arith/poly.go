package arith

import (
	"math/big"
	"math/rand"
	"sort"
)

// Polynomials over F_p are represented as []uint64 with coeff of x^i at [i].
// Canonical polynomials have no trailing zero coefficients; the zero polynomial is nil/empty.

// PolyTrim strips trailing zero coefficients.
func PolyTrim(a []uint64) []uint64 {
	n := len(a)
	for n > 0 && a[n-1] == 0 {
		n--
	}
	return a[:n]
}

// PolyReduce reduces a polynomial's coefficients mod p and trims trailing zeros.
func PolyReduce(a []uint64, p uint64) []uint64 {
	out := make([]uint64, len(a))
	for i, c := range a {
		out[i] = c % p
	}
	return PolyTrim(out)
}

// PolyEqual reports whether two trimmed polynomials are equal.
func PolyEqual(a, b []uint64) bool {
	a, b = PolyTrim(a), PolyTrim(b)
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// PolyMonic returns the monic associate of a mod prime p.
func PolyMonic(a []uint64, p uint64) []uint64 {
	a = PolyTrim(a)
	if len(a) == 0 {
		return nil
	}
	inv := InvMod(a[len(a)-1], p)
	out := make([]uint64, len(a))
	for i, c := range a {
		out[i] = MulMod64(c, inv, p)
	}
	return out
}

// PolyAddModP returns (a + b) mod p.
func PolyAddModP(a, b []uint64, p uint64) []uint64 {
	n := len(a)
	if len(b) > n {
		n = len(b)
	}
	out := make([]uint64, n)
	for i := 0; i < n; i++ {
		var av, bv uint64
		if i < len(a) {
			av = a[i] % p
		}
		if i < len(b) {
			bv = b[i] % p
		}
		s := av + bv
		if s >= p || s < av {
			s -= p
		}
		out[i] = s
	}
	return PolyTrim(out)
}

// PolySubModP returns (a - b) mod p.
func PolySubModP(a, b []uint64, p uint64) []uint64 {
	n := len(a)
	if len(b) > n {
		n = len(b)
	}
	out := make([]uint64, n)
	for i := 0; i < n; i++ {
		var av, bv uint64
		if i < len(a) {
			av = a[i] % p
		}
		if i < len(b) {
			bv = b[i] % p
		}
		if av >= bv {
			out[i] = av - bv
		} else {
			out[i] = p - (bv - av)
		}
	}
	return PolyTrim(out)
}

// PolyMulModP returns (a * b) mod p.
func PolyMulModP(a, b []uint64, p uint64) []uint64 {
	a, b = PolyTrim(a), PolyTrim(b)
	if len(a) == 0 || len(b) == 0 {
		return nil
	}
	out := make([]uint64, len(a)+len(b)-1)
	if p <= 1<<32 {
		for i, ai := range a {
			if ai == 0 {
				continue
			}
			for j, bj := range b {
				out[i+j] = (out[i+j] + ai*bj) % p
			}
		}
		return PolyTrim(out)
	}
	for i, ai := range a {
		if ai == 0 {
			continue
		}
		for j, bj := range b {
			prod := MulMod64(ai, bj, p)
			s := out[i+j] + prod
			if s >= p || s < out[i+j] {
				s -= p
			}
			out[i+j] = s
		}
	}
	return PolyTrim(out)
}

// PolyDivRemModP divides a by b over F_p, returning quotient q and remainder r.
func PolyDivRemModP(a, b []uint64, p uint64) (q, r []uint64) {
	b = PolyTrim(b)
	if len(b) == 0 {
		panic("arith: division by zero polynomial")
	}
	r = append([]uint64(nil), PolyTrim(a)...)
	for i := range r {
		r[i] %= p
	}
	r = PolyTrim(r)
	if len(r) < len(b) {
		return nil, r
	}
	degB := len(b) - 1
	invLead := InvMod(b[degB], p)
	q = make([]uint64, len(r)-degB)
	for len(r) >= len(b) {
		degR := len(r) - 1
		shift := degR - degB
		coeff := MulMod64(r[degR], invLead, p)
		q[shift] = coeff
		for i := 0; i <= degB; i++ {
			sub := MulMod64(coeff, b[i], p)
			if r[shift+i] >= sub {
				r[shift+i] -= sub
			} else {
				r[shift+i] = p - (sub - r[shift+i])
			}
		}
		r = PolyTrim(r)
	}
	return PolyTrim(q), r
}

// PolyModModP returns a mod h over F_p.
func PolyModModP(a, h []uint64, p uint64) []uint64 {
	_, r := PolyDivRemModP(a, h, p)
	return r
}

// PolyGcdModP returns the monic greatest common divisor of a and b over F_p.
func PolyGcdModP(a, b []uint64, p uint64) []uint64 {
	a, b = PolyReduce(a, p), PolyReduce(b, p)
	for len(b) > 0 {
		_, r := PolyDivRemModP(a, b, p)
		a, b = b, r
	}
	return PolyMonic(a, p)
}

// PolyInvModP computes the inverse of a mod (h, p). Returns (nil, false) if not invertible.
func PolyInvModP(a, h []uint64, p uint64) ([]uint64, bool) {
	r0, r1 := PolyReduce(h, p), PolyReduce(a, p)
	var t0 []uint64
	t1 := []uint64{1}
	for len(r1) > 0 {
		q, r2 := PolyDivRemModP(r0, r1, p)
		qt := PolyMulModP(q, t1, p)
		t2 := PolySubModP(t0, qt, p)
		r0, r1 = r1, r2
		t0, t1 = t1, t2
	}
	if len(r0) != 1 {
		return nil, false
	}
	inv := InvMod(r0[0], p)
	for i := range t0 {
		t0[i] = MulMod64(t0[i], inv, p)
	}
	return PolyTrim(t0), true
}

// PolyDerivModP returns the formal derivative of a mod p.
func PolyDerivModP(a []uint64, p uint64) []uint64 {
	if len(a) <= 1 {
		return nil
	}
	out := make([]uint64, len(a)-1)
	for i := 1; i < len(a); i++ {
		out[i-1] = MulMod64(a[i], uint64(i)%p, p)
	}
	return PolyTrim(out)
}

// PolyEvalModP evaluates a(x) mod p using Horner's method.
func PolyEvalModP(a []uint64, x, p uint64) uint64 {
	var res uint64
	for i := len(a) - 1; i >= 0; i-- {
		prod := MulMod64(res, x, p)
		av := a[i] % p
		res = prod + av
		if res >= p || res < prod {
			res -= p
		}
	}
	return res
}

// PolyPowMod computes a^e mod (h, p).
func PolyPowMod(a []uint64, e uint64, h []uint64, p uint64) []uint64 {
	res := []uint64{1}
	if len(PolyTrim(h)) <= 1 {
		return nil
	}
	base := PolyModModP(a, h, p)
	for e > 0 {
		if e&1 != 0 {
			res = PolyModModP(PolyMulModP(res, base, p), h, p)
		}
		e >>= 1
		if e > 0 {
			base = PolyModModP(PolyMulModP(base, base, p), h, p)
		}
	}
	return res
}

// PolyPowBigMod computes a^e mod (h, p) for a big.Int exponent e >= 0.
func PolyPowBigMod(a []uint64, e *big.Int, h []uint64, p uint64) []uint64 {
	res := []uint64{1}
	if len(PolyTrim(h)) <= 1 {
		return nil
	}
	base := PolyModModP(a, h, p)
	bits := e.BitLen()
	for i := 0; i < bits; i++ {
		if e.Bit(i) != 0 {
			res = PolyModModP(PolyMulModP(res, base, p), h, p)
		}
		if i+1 < bits {
			base = PolyModModP(PolyMulModP(base, base, p), h, p)
		}
	}
	return res
}

// PolyPowXMod computes x^e mod (h, p).
func PolyPowXMod(e uint64, h []uint64, p uint64) []uint64 {
	return PolyPowMod([]uint64{0, 1}, e, h, p)
}

// PolyIsIrreducibleModP runs Rabin's irreducibility test on f over F_p.
// p must be prime.
func PolyIsIrreducibleModP(f []uint64, p uint64) bool {
	f = PolyReduce(f, p)
	d := len(f) - 1
	if d <= 0 {
		return false
	}
	if d == 1 {
		return true
	}
	f = PolyMonic(f, p)

	// Distinct prime factors of d.
	var primeDivs []int
	rem := d
	for q := 2; q*q <= rem; q++ {
		if rem%q == 0 {
			primeDivs = append(primeDivs, q)
			for rem%q == 0 {
				rem /= q
			}
		}
	}
	if rem > 1 {
		primeDivs = append(primeDivs, rem)
	}
	checkAt := make(map[int]bool, len(primeDivs))
	for _, q := range primeDivs {
		checkAt[d/q] = true
	}

	xPoly := []uint64{0, 1}
	g := []uint64{0, 1} // x^{p^k} mod f
	for k := 1; k <= d; k++ {
		g = PolyPowMod(g, p, f, p)
		if checkAt[k] {
			diff := PolySubModP(g, xPoly, p)
			if len(PolyGcdModP(f, diff, p)) > 1 {
				return false
			}
		}
	}
	return PolyEqual(g, xPoly)
}

// PolyRootsModP returns all distinct roots of f in F_p, sorted ascending.
// p must be prime and f must not be the zero polynomial mod p.
func PolyRootsModP(f []uint64, p uint64) []uint64 {
	f = PolyReduce(f, p)
	if len(f) <= 1 {
		return nil
	}
	if len(f) == 2 {
		// f = c1*x + c0
		c0 := (p - f[0]%p) % p
		return []uint64{MulMod64(c0, InvMod(f[1], p), p)}
	}

	var roots []uint64
	if f[0] == 0 {
		roots = append(roots, 0)
		i := 0
		for i < len(f) && f[i] == 0 {
			i++
		}
		f = f[i:]
	}
	if len(f) <= 1 {
		return roots
	}
	if p == 2 {
		if PolyEvalModP(f, 1, 2) == 0 {
			roots = append(roots, 1)
		}
		return roots
	}

	xpMinusX := PolySubModP(PolyPowXMod(p, f, p), []uint64{0, 1}, p)
	g := PolyGcdModP(f, xpMinusX, p)
	if len(g) <= 1 {
		return roots
	}

	rng := rand.New(rand.NewSource(1))
	var split func(h []uint64)
	split = func(h []uint64) {
		deg := len(h) - 1
		if deg <= 0 {
			return
		}
		if deg == 1 {
			r := (p - h[0]%p) % p
			r = MulMod64(r, InvMod(h[1], p), p)
			roots = append(roots, r)
			return
		}
		for {
			t := uint64(rng.Int63n(int64(p & ((1 << 62) - 1))))
			if p > (1 << 62) {
				t = rng.Uint64() % p
			}
			// (x + t)^((p-1)/2) - 1 mod h
			xt := []uint64{t, 1}
			pow := PolyPowMod(xt, (p-1)>>1, h, p)
			s := PolyGcdModP(h, PolySubModP(pow, []uint64{1}, p), p)
			if len(s) > 1 && len(s) < len(h) {
				q, _ := PolyDivRemModP(h, s, p)
				split(s)
				split(q)
				return
			}
		}
	}
	split(g)

	sort.Slice(roots, func(i, j int) bool { return roots[i] < roots[j] })
	return roots
}

// PolyFactorModP factors a squarefree polynomial f over F_p into monic
// irreducible factors using distinct-degree + Cantor-Zassenhaus equal-degree
// factorization. p must be prime and gcd(f, f') = 1 mod p.
func PolyFactorModP(f []uint64, p uint64) [][]uint64 {
	f = PolyMonic(PolyReduce(f, p), p)
	if len(f) <= 1 {
		return nil
	}
	if len(f) == 2 {
		return [][]uint64{f}
	}

	var factors [][]uint64
	rng := rand.New(rand.NewSource(1))

	// Distinct-degree factorization.
	h := append([]uint64(nil), f...)
	xPow := []uint64{0, 1} // x^{p^d} mod h
	xPoly := []uint64{0, 1}
	for d := 1; len(h)-1 >= 2*d; d++ {
		xPow = PolyPowMod(xPow, p, h, p)
		g := PolyGcdModP(h, PolySubModP(xPow, xPoly, p), p)
		if len(g) > 1 {
			factors = append(factors, edf(g, d, p, rng)...)
			h, _ = PolyDivRemModP(h, g, p)
			h = PolyMonic(h, p)
			xPow = PolyModModP(xPow, h, p)
		}
	}
	if len(h) > 1 {
		factors = append(factors, h)
	}
	return factors
}

// edf performs equal-degree factorization of g into factors of degree d over F_p.
func edf(g []uint64, d int, p uint64, rng *rand.Rand) [][]uint64 {
	deg := len(g) - 1
	if deg == d {
		return [][]uint64{PolyMonic(g, p)}
	}
	if d == 1 {
		rs := PolyRootsModP(g, p)
		out := make([][]uint64, len(rs))
		for i, r := range rs {
			out[i] = []uint64{(p - r) % p, 1}
		}
		return out
	}

	var res [][]uint64
	var split func(h []uint64)
	split = func(h []uint64) {
		if len(h)-1 == d {
			res = append(res, PolyMonic(h, p))
			return
		}
		for {
			// Random non-constant polynomial a of degree < len(h)-1.
			a := make([]uint64, len(h)-1)
			for i := range a {
				a[i] = rng.Uint64() % p
			}
			a = PolyTrim(a)
			if len(a) <= 1 {
				continue
			}
			if g1 := PolyGcdModP(a, h, p); len(g1) > 1 && len(g1) < len(h) {
				q, _ := PolyDivRemModP(h, g1, p)
				split(g1)
				split(q)
				return
			}
			var b []uint64
			if p == 2 {
				// Trace map: a + a^2 + a^4 + ... + a^{2^{d-1}} mod h
				b = append([]uint64(nil), a...)
				cur := append([]uint64(nil), a...)
				for i := 1; i < d; i++ {
					cur = PolyModModP(PolyMulModP(cur, cur, p), h, p)
					b = PolyAddModP(b, cur, p)
				}
			} else {
				// Compute (p^d - 1) / 2.
				exp := new(big.Int).Exp(new(big.Int).SetUint64(p), big.NewInt(int64(d)), nil)
				exp.Sub(exp, big.NewInt(1))
				exp.Rsh(exp, 1)
				b = PolySubModP(PolyPowBigMod(a, exp, h, p), []uint64{1}, p)
			}
			s := PolyGcdModP(h, b, p)
			if len(s) > 1 && len(s) < len(h) {
				q, _ := PolyDivRemModP(h, s, p)
				split(s)
				split(q)
				return
			}
		}
	}
	split(g)
	return res
}

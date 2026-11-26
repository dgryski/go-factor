// Package gnfs implements the General Number Field Sieve (GNFS) integer factorization algorithm.
package gnfs

import (
	"math"
	"math/big"
	"runtime"
	"sort"
	"sync"
	"sync/atomic"

	"github.com/dgryski/go-factor/internal/arith"
)

// Poly represents an algebraic/rational polynomial pair (f, g) with common
// root M mod N: f(M) ≡ 0 (mod N) and g(M) = A1*M - A0 ≡ 0 (mod N).
type Poly struct {
	Degree int
	Coeffs []*big.Int // c_0 .. c_d
	A1, A0 *big.Int   // g(x) = A1*x - A0, gcd(A1, A0) = 1
	M      *big.Int   // A0 * A1^-1 mod N
	Skew   float64
	Alpha  float64
}

// EvalF evaluates the homogeneous polynomial F(a, b) = b^d * f(a/b) = sum_{i=0}^d c_i * a^i * b^{d-i}.
func (p *Poly) EvalF(a, b int64) *big.Int {
	d := p.Degree
	bigA := big.NewInt(a)
	bigB := big.NewInt(b)

	aPows := make([]*big.Int, d+1)
	bPows := make([]*big.Int, d+1)
	aPows[0] = big.NewInt(1)
	bPows[0] = big.NewInt(1)
	for i := 1; i <= d; i++ {
		aPows[i] = new(big.Int).Mul(aPows[i-1], bigA)
		bPows[i] = new(big.Int).Mul(bPows[i-1], bigB)
	}

	res := new(big.Int)
	term := new(big.Int)
	for i := 0; i <= d; i++ {
		term.Mul(p.Coeffs[i], aPows[i])
		term.Mul(term, bPows[d-i])
		res.Add(res, term)
	}
	return res
}

// EvalG evaluates G(a, b) = A1*a - A0*b.
func (p *Poly) EvalG(a, b int64) *big.Int {
	t1 := new(big.Int).Mul(p.A1, big.NewInt(a))
	t2 := new(big.Int).Mul(p.A0, big.NewInt(b))
	return t1.Sub(t1, t2)
}

// Monic returns the monic polynomial f_hat(y) = c_d^(d-1) * f(y / c_d)
// and its derivative f_hat'(y).
func (p *Poly) Monic() (fHat, fHatDeriv []*big.Int) {
	d := p.Degree
	cd := p.Coeffs[d]
	fHat = make([]*big.Int, d+1)
	cdPow := big.NewInt(1)
	for i := d - 1; i >= 0; i-- {
		fHat[i] = new(big.Int).Mul(p.Coeffs[i], cdPow)
		cdPow = new(big.Int).Mul(cdPow, cd)
	}
	fHat[d] = big.NewInt(1)

	fHatDeriv = make([]*big.Int, d)
	for i := 1; i <= d; i++ {
		fHatDeriv[i-1] = new(big.Int).Mul(fHat[i], big.NewInt(int64(i)))
	}
	return fHat, fHatDeriv
}

// CoeffsModP returns f(x) mod p as []uint64.
func (p *Poly) CoeffsModP(mod uint64) []uint64 {
	out := make([]uint64, p.Degree+1)
	var bigP, r big.Int
	bigP.SetUint64(mod)
	for i, c := range p.Coeffs {
		r.Mod(c, &bigP)
		out[i] = r.Uint64()
	}
	return arith.PolyTrim(out)
}

// BaseMPoly constructs a degree-d polynomial with leading coefficient cd (or cd+1
// after balanced carry) using base-m expansion m = floor((n / cd)^(1/d)).
// If a non-trivial factor of n is discovered directly (e.g., via content or
// reducible f), it is returned as factor.
func BaseMPoly(n *big.Int, d int, cd uint64) (poly *Poly, factor *big.Int) {
	if cd == 0 {
		cd = 1
	}
	cdBig := new(big.Int).SetUint64(cd)
	q := new(big.Int).Quo(n, cdBig)
	if q.Sign() <= 0 {
		return nil, nil
	}
	m := arith.NthRoot(q, d)
	if m.Cmp(big.NewInt(2)) < 0 {
		return nil, nil
	}

	coeffs := make([]*big.Int, d+1)
	rem := new(big.Int).Set(n)
	for i := 0; i < d; i++ {
		quot, r := new(big.Int).QuoRem(rem, m, new(big.Int))
		coeffs[i] = r
		rem = quot
	}
	coeffs[d] = rem

	// Balance digits into (-m/2, m/2] by carrying.
	halfM := new(big.Int).Rsh(m, 1)
	one := big.NewInt(1)
	for i := 0; i < d; i++ {
		if coeffs[i].Cmp(halfM) > 0 {
			coeffs[i].Sub(coeffs[i], m)
			coeffs[i+1].Add(coeffs[i+1], one)
		}
	}

	// Remove content if > 1.
	content := new(big.Int).Abs(coeffs[0])
	for i := 1; i <= d; i++ {
		content.GCD(nil, nil, content, new(big.Int).Abs(coeffs[i]))
	}
	if content.Cmp(one) > 0 {
		var g big.Int
		g.GCD(nil, nil, content, n)
		if g.Cmp(one) > 0 && g.Cmp(n) < 0 {
			return nil, &g
		}
		for i := 0; i <= d; i++ {
			coeffs[i].Quo(coeffs[i], content)
		}
	}

	// If c_0 == 0, then x | f(x), so m | f(m) = n.
	if coeffs[0].Sign() == 0 {
		var g big.Int
		g.GCD(nil, nil, m, n)
		if g.Cmp(one) > 0 && g.Cmp(n) < 0 {
			return nil, &g
		}
		return nil, nil
	}

	p := &Poly{
		Degree: d,
		Coeffs: coeffs,
		A1:     big.NewInt(1),
		A0:     new(big.Int).Set(m),
		M:      new(big.Int).Mod(m, n),
	}

	// Check irreducibility; if reducible and f(m) = n, we may get a factor.
	if !IsIrreducible(p) {
		if fac := findFactorFromReducible(p, n); fac != nil {
			return nil, fac
		}
		return nil, nil
	}

	p.Skew = OptimizeSkewness(p)
	p.Alpha = ComputeAlpha(p, 200)
	return p, nil
}

// IsIrreducible reports whether p.Coeffs is irreducible over Z.
func IsIrreducible(p *Poly) bool {
	d := p.Degree
	if d <= 1 {
		return d == 1
	}
	if p.Coeffs[0].Sign() == 0 || p.Coeffs[d].Sign() == 0 {
		return false
	}

	// Fast check: if f mod prime is irreducible for any prime not dividing c_d, f is irreducible over Z.
	var bigP, r big.Int
	for _, prime := range arith.PrimesUpTo(250) {
		up := uint64(prime)
		bigP.SetUint64(up)
		if r.Mod(p.Coeffs[d], &bigP).Sign() == 0 {
			continue
		}
		fp := p.CoeffsModP(up)
		if len(fp)-1 == d && arith.PolyIsIrreducibleModP(fp, up) {
			return true
		}
	}

	// Fallback (especially for d = 4 when Gal(f) has no d-cycle): check rational roots
	// and, for d >= 4, check mod-p factorizations with Hensel-lifted recombination.
	return isIrreducibleFallback(p)
}

func isIrreducibleFallback(p *Poly) bool {
	fHat, _ := p.Monic()
	d := p.Degree
	// Find a prime p not dividing disc(fHat), i.e., gcd(fHat, fHat') == 1 mod p.
	var bigP, r big.Int
	var chosenP uint64
	var facs [][]uint64
	for prime := int64(3); prime < 10000; prime += 2 {
		if !arith.IsPrime64(uint64(prime)) {
			continue
		}
		up := uint64(prime)
		bigP.SetUint64(up)
		fp := make([]uint64, d+1)
		for i, c := range fHat {
			fp[i] = r.Mod(c, &bigP).Uint64()
		}
		dfp := arith.PolyDerivModP(fp, up)
		if len(arith.PolyGcdModP(fp, dfp, up)) != 1 {
			continue
		}
		facs = arith.PolyFactorModP(fp, up)
		if len(facs) == 1 {
			return true
		}
		chosenP = up
		break
	}
	if chosenP == 0 {
		return false
	}

	// Landau-Mignotte bound for factors of monic fHat in Z[y]:
	// ||h||_inf <= 2^d * ||fHat||_2.
	maxLog2 := 0.0
	for _, c := range fHat {
		if c.Sign() != 0 {
			l := arith.Log2Big(new(big.Int).Abs(c))
			if l > maxLog2 {
				maxLog2 = l
			}
		}
	}
	boundBits := int(math.Ceil(maxLog2)) + d + 8
	mod := big.NewInt(1)
	bigP.SetUint64(chosenP)
	for mod.BitLen() <= boundBits+1 {
		mod.Mul(mod, &bigP)
	}

	// Hensel lift each irreducible factor mod p up to mod, and test subsets of degree <= d/2.
	lifted := henselLiftFactors(fHat, facs, chosenP, mod)
	k := len(lifted)
	for mask := 1; mask < (1<<uint(k))-1; mask++ {
		deg := 0
		for i := 0; i < k; i++ {
			if mask&(1<<uint(i)) != 0 {
				deg += len(lifted[i]) - 1
			}
		}
		if deg == 0 || deg > d/2 {
			continue
		}
		cand := []*big.Int{big.NewInt(1)}
		for i := 0; i < k; i++ {
			if mask&(1<<uint(i)) != 0 {
				cand = polyMulBigMod(cand, lifted[i], mod)
			}
		}
		// Convert to symmetric residues in (-mod/2, mod/2].
		half := new(big.Int).Rsh(mod, 1)
		for i := range cand {
			if cand[i].Cmp(half) > 0 {
				cand[i].Sub(cand[i], mod)
			}
		}
		if polyDividesZ(cand, fHat) {
			return false
		}
	}
	return true
}

func henselLiftFactors(f []*big.Int, facs [][]uint64, p uint64, targetMod *big.Int) [][]*big.Int {
	k := len(facs)
	cur := make([][]*big.Int, k)
	for i, fac := range facs {
		cur[i] = make([]*big.Int, len(fac))
		for j, c := range fac {
			cur[i][j] = new(big.Int).SetUint64(c)
		}
	}
	if k == 1 {
		return cur
	}
	// Pairwise quadratic Hensel lifting via a binary split.
	var liftPair func(f []*big.Int, g, h []*big.Int) ([]*big.Int, []*big.Int)
	liftPair = func(f, g, h []*big.Int) ([]*big.Int, []*big.Int) {
		// Compute Bezout identity s*g + t*h = 1 mod p.
		gp := make([]uint64, len(g))
		hp := make([]uint64, len(h))
		bigP := new(big.Int).SetUint64(p)
		var rem big.Int
		for i, c := range g {
			gp[i] = rem.Mod(c, bigP).Uint64()
		}
		for i, c := range h {
			hp[i] = rem.Mod(c, bigP).Uint64()
		}
		s0, ok := arith.PolyInvModP(gp, hp, p)
		if !ok {
			return g, h
		}
		// t0 = (1 - s0*gp) / hp mod p
		oneMinusSG := arith.PolySubModP([]uint64{1}, arith.PolyMulModP(s0, gp, p), p)
		t0, _ := arith.PolyDivRemModP(oneMinusSG, hp, p)

		s := make([]*big.Int, len(s0))
		for i, c := range s0 {
			s[i] = new(big.Int).SetUint64(c)
		}
		t := make([]*big.Int, len(t0))
		for i, c := range t0 {
			t[i] = new(big.Int).SetUint64(c)
		}

		curMod := new(big.Int).SetUint64(p)
		for curMod.Cmp(targetMod) < 0 {
			nextMod := new(big.Int).Mul(curMod, curMod)
			if nextMod.Cmp(targetMod) > 0 {
				nextMod.Set(targetMod)
			}
			// Lift factors g, h from curMod to nextMod.
			e := polySubBigMod(f, polyMulBigMod(g, h, nextMod), nextMod)
			// q, r := (s*e) / h mod nextMod, where h is monic.
			se := polyMulBigMod(s, e, nextMod)
			q, r := polyDivRemMonicBigMod(se, h, nextMod)
			te := polyMulBigMod(t, e, nextMod)
			qh := polyMulBigMod(q, g, nextMod)
			gNew := polyAddBigMod(g, polyAddBigMod(te, qh, nextMod), nextMod)
			hNew := polyAddBigMod(h, r, nextMod)
			g, h = gNew, hNew

			// Lift Bezout cofactors s, t from curMod to nextMod.
			// d = s*g + t*h - 1 mod nextMod
			d := polySubBigMod(polyAddBigMod(polyMulBigMod(s, g, nextMod), polyMulBigMod(t, h, nextMod), nextMod), []*big.Int{big.NewInt(1)}, nextMod)
			sd := polyMulBigMod(s, d, nextMod)
			q2, r2 := polyDivRemMonicBigMod(sd, h, nextMod)
			td := polyMulBigMod(t, d, nextMod)
			qg := polyMulBigMod(q2, g, nextMod)
			s = polySubBigMod(s, r2, nextMod)
			t = polySubBigMod(t, polyAddBigMod(td, qg, nextMod), nextMod)

			curMod = nextMod
		}
		return g, h
	}

	var splitLift func(poly []*big.Int, parts [][]*big.Int) [][]*big.Int
	splitLift = func(poly []*big.Int, parts [][]*big.Int) [][]*big.Int {
		if len(parts) == 1 {
			res := make([]*big.Int, len(poly))
			for i, c := range poly {
				res[i] = new(big.Int).Mod(c, targetMod)
			}
			return [][]*big.Int{res}
		}
		mid := len(parts) / 2
		bigP := new(big.Int).SetUint64(p)
		left := []*big.Int{big.NewInt(1)}
		for _, pt := range parts[:mid] {
			left = polyMulBigMod(left, pt, bigP)
		}
		right := []*big.Int{big.NewInt(1)}
		for _, pt := range parts[mid:] {
			right = polyMulBigMod(right, pt, bigP)
		}
		lLifted, rLifted := liftPair(poly, left, right)
		out := splitLift(lLifted, parts[:mid])
		out = append(out, splitLift(rLifted, parts[mid:])...)
		return out
	}

	return splitLift(f, cur)
}

func polyTrimBig(a []*big.Int) []*big.Int {
	n := len(a)
	for n > 0 && a[n-1].Sign() == 0 {
		n--
	}
	return a[:n]
}

func polyAddBigMod(a, b []*big.Int, mod *big.Int) []*big.Int {
	n := len(a)
	if len(b) > n {
		n = len(b)
	}
	out := make([]*big.Int, n)
	for i := 0; i < n; i++ {
		s := new(big.Int)
		if i < len(a) {
			s.Add(s, a[i])
		}
		if i < len(b) {
			s.Add(s, b[i])
		}
		s.Mod(s, mod)
		out[i] = s
	}
	return polyTrimBig(out)
}

func polySubBigMod(a, b []*big.Int, mod *big.Int) []*big.Int {
	n := len(a)
	if len(b) > n {
		n = len(b)
	}
	out := make([]*big.Int, n)
	for i := 0; i < n; i++ {
		s := new(big.Int)
		if i < len(a) {
			s.Add(s, a[i])
		}
		if i < len(b) {
			s.Sub(s, b[i])
		}
		s.Mod(s, mod)
		out[i] = s
	}
	return polyTrimBig(out)
}

func polyMulBigMod(a, b []*big.Int, mod *big.Int) []*big.Int {
	a, b = polyTrimBig(a), polyTrimBig(b)
	if len(a) == 0 || len(b) == 0 {
		return nil
	}
	out := make([]*big.Int, len(a)+len(b)-1)
	for i := range out {
		out[i] = new(big.Int)
	}
	var prod big.Int
	for i, ai := range a {
		if ai.Sign() == 0 {
			continue
		}
		for j, bj := range b {
			if bj.Sign() == 0 {
				continue
			}
			prod.Mul(ai, bj)
			out[i+j].Add(out[i+j], &prod)
		}
	}
	for i := range out {
		out[i].Mod(out[i], mod)
	}
	return polyTrimBig(out)
}

func polyDivRemMonicBigMod(a, b []*big.Int, mod *big.Int) (q, r []*big.Int) {
	b = polyTrimBig(b)
	r = make([]*big.Int, len(a))
	for i, c := range a {
		r[i] = new(big.Int).Mod(c, mod)
	}
	r = polyTrimBig(r)
	if len(r) < len(b) {
		return nil, r
	}
	degB := len(b) - 1
	q = make([]*big.Int, len(r)-degB)
	for i := range q {
		q[i] = new(big.Int)
	}
	var prod big.Int
	for len(r) >= len(b) {
		degR := len(r) - 1
		shift := degR - degB
		coeff := new(big.Int).Set(r[degR])
		q[shift] = coeff
		for i := 0; i <= degB; i++ {
			prod.Mul(coeff, b[i])
			r[shift+i].Sub(r[shift+i], &prod)
			r[shift+i].Mod(r[shift+i], mod)
		}
		r = polyTrimBig(r)
	}
	return polyTrimBig(q), r
}

func polyDividesZ(h, f []*big.Int) bool {
	h = polyTrimBig(h)
	f = polyTrimBig(f)
	if len(h) == 0 || len(f) < len(h) {
		return false
	}
	r := make([]*big.Int, len(f))
	for i, c := range f {
		r[i] = new(big.Int).Set(c)
	}
	degH := len(h) - 1
	leadH := h[degH]
	var quot, rem, prod big.Int
	for len(r) >= len(h) {
		degR := len(r) - 1
		shift := degR - degH
		quot.QuoRem(r[degR], leadH, &rem)
		if rem.Sign() != 0 {
			return false
		}
		for i := 0; i <= degH; i++ {
			prod.Mul(&quot, h[i])
			r[shift+i].Sub(r[shift+i], &prod)
		}
		r = polyTrimBig(r)
	}
	return len(r) == 0
}

// findFactorFromReducible tries to extract a factor of n from a reducible
// monic-or-general polynomial f where f(m) ≡ 0 (mod n).
func findFactorFromReducible(p *Poly, n *big.Int) *big.Int {
	fHat, _ := p.Monic()
	d := p.Degree
	mHat := new(big.Int).Mul(p.M, p.Coeffs[d])
	mHat.Mod(mHat, n)

	// Try to factor mod a small squarefree prime and evaluate lifted factors at mHat.
	var bigP, r big.Int
	for prime := int64(3); prime < 2000; prime += 2 {
		if !arith.IsPrime64(uint64(prime)) {
			continue
		}
		up := uint64(prime)
		bigP.SetUint64(up)
		fp := make([]uint64, d+1)
		for i, c := range fHat {
			fp[i] = r.Mod(c, &bigP).Uint64()
		}
		dfp := arith.PolyDerivModP(fp, up)
		if len(arith.PolyGcdModP(fp, dfp, up)) != 1 {
			continue
		}
		facs := arith.PolyFactorModP(fp, up)
		if len(facs) <= 1 {
			return nil
		}
		maxLog2 := 0.0
		for _, c := range fHat {
			if c.Sign() != 0 {
				l := arith.Log2Big(new(big.Int).Abs(c))
				if l > maxLog2 {
					maxLog2 = l
				}
			}
		}
		boundBits := int(math.Ceil(maxLog2)) + d + 8
		mod := big.NewInt(1)
		for mod.BitLen() <= boundBits+1 {
			mod.Mul(mod, &bigP)
		}
		lifted := henselLiftFactors(fHat, facs, up, mod)
		k := len(lifted)
		for mask := 1; mask < (1<<uint(k))-1; mask++ {
			deg := 0
			for i := 0; i < k; i++ {
				if mask&(1<<uint(i)) != 0 {
					deg += len(lifted[i]) - 1
				}
			}
			if deg == 0 || deg > d/2 {
				continue
			}
			cand := []*big.Int{big.NewInt(1)}
			for i := 0; i < k; i++ {
				if mask&(1<<uint(i)) != 0 {
					cand = polyMulBigMod(cand, lifted[i], mod)
				}
			}
			half := new(big.Int).Rsh(mod, 1)
			for i := range cand {
				if cand[i].Cmp(half) > 0 {
					cand[i].Sub(cand[i], mod)
				}
			}
			if polyDividesZ(cand, fHat) {
				// Evaluate cand(mHat) mod n.
				val := big.NewInt(0)
				for i := len(cand) - 1; i >= 0; i-- {
					val.Mul(val, mHat)
					val.Add(val, cand[i])
					val.Mod(val, n)
				}
				var g big.Int
				g.GCD(nil, nil, val, n)
				if g.Cmp(big.NewInt(1)) > 0 && g.Cmp(n) < 0 {
					return &g
				}
			}
		}
		return nil
	}
	return nil
}

// OptimizeSkewness finds the skewness s >= 1 minimizing log(RMS(F)) + log(RMS(G)).
func OptimizeSkewness(p *Poly) float64 {
	return optimizeSkewnessBounded(p, 0)
}

func optimizeSkewnessBounded(p *Poly, maxSkew float64) float64 {
	d := p.Degree
	maxBits := 0
	for _, c := range p.Coeffs {
		if b := c.BitLen(); b > maxBits {
			maxBits = b
		}
	}
	shiftF := maxBits - 500
	if shiftF < 0 {
		shiftF = 0
	}
	cf := make([]float64, d+1)
	for i, c := range p.Coeffs {
		cf[i] = scaledFloat(c, shiftF)
	}

	maxGBits := max(p.A1.BitLen(), p.A0.BitLen())
	shiftG := maxGBits - 500
	if shiftG < 0 {
		shiftG = 0
	}
	a1f := scaledFloat(p.A1, shiftG)
	a0f := scaledFloat(p.A0, shiftG)

	s0 := 1.0
	if cf[d] != 0 && cf[0] != 0 {
		s0 = math.Pow(math.Abs(cf[0]/cf[d]), 1.0/float64(d))
	}
	if s0 < 1 {
		s0 = 1
	}

	evalLogNorm := func(logS float64) float64 {
		s := math.Exp(logS)
		var sumF float64
		for i := 0; i <= d; i++ {
			for j := 0; j <= d; j++ {
				if (i+j)&1 != 0 {
					continue
				}
				sumF += cf[i] * cf[j] * math.Pow(s, float64(i+j-d)) / float64((i+j+1)*(2*d-i-j+1))
			}
		}
		if sumF <= 0 {
			sumF = 1e-300
		}
		sumG := (a1f*a1f*s + a0f*a0f/s) / 3.0
		if sumG <= 0 {
			sumG = 1e-300
		}
		return 0.5*math.Log(sumF) + 0.5*math.Log(sumG)
	}

	center := math.Log(s0)
	lo := math.Max(0, center-8)
	hi := center + 12
	if maxSkew >= 1 {
		logMax := math.Log(maxSkew)
		if hi > logMax {
			hi = logMax
		}
		if lo > hi {
			lo = math.Max(0, hi-8)
		}
	}
	if hi <= lo {
		return math.Exp(lo)
	}

	const gridN = 32
	bestK := 0
	bestV := math.Inf(1)
	for k := 0; k <= gridN; k++ {
		ls := lo + (hi-lo)*float64(k)/gridN
		if v := evalLogNorm(ls); v < bestV {
			bestV = v
			bestK = k
		}
	}
	step := (hi - lo) / gridN
	lo = math.Max(lo, lo+float64(bestK-1)*step)
	hi = math.Min(hi, lo+2*step)

	phi := (math.Sqrt(5) - 1) / 2
	x1 := hi - phi*(hi-lo)
	x2 := lo + phi*(hi-lo)
	f1, f2 := evalLogNorm(x1), evalLogNorm(x2)
	for iter := 0; iter < 36; iter++ {
		if f1 < f2 {
			hi = x2
			x2 = x1
			f2 = f1
			x1 = hi - phi*(hi-lo)
			f1 = evalLogNorm(x1)
		} else {
			lo = x1
			x1 = x2
			f1 = f2
			x2 = lo + phi*(hi-lo)
			f2 = evalLogNorm(x2)
		}
	}
	return math.Exp(0.5 * (lo + hi))
}

func scaledFloat(x *big.Int, shift int) float64 {
	if x.Sign() == 0 {
		return 0
	}
	if shift == 0 {
		f, _ := new(big.Float).SetInt(x).Float64()
		return f
	}
	var y big.Int
	if x.Sign() < 0 {
		y.Neg(x)
		y.Rsh(&y, uint(shift))
		f, _ := new(big.Float).SetInt(&y).Float64()
		return -f
	}
	y.Rsh(x, uint(shift))
	f, _ := new(big.Float).SetInt(&y).Float64()
	return f
}

// ComputeAlpha computes Murphy's alpha(f) over primes <= bound.
func ComputeAlpha(p *Poly, bound int64) float64 {
	primes := arith.PrimesUpTo(bound)
	d := p.Degree
	alpha := 0.0
	for _, prime := range primes {
		up := uint64(prime)
		fp := make([]uint64, d+1)
		var bigP, r big.Int
		bigP.SetUint64(up)
		for i, c := range p.Coeffs {
			fp[i] = r.Mod(c, &bigP).Uint64()
		}
		ev := expectedValP(p.Coeffs, fp, up, d)
		alpha += (1.0/float64(prime-1) - ev) * math.Log(float64(prime))
	}
	return alpha
}

func expectedValP(coeffs []*big.Int, fp []uint64, p uint64, d int) float64 {
	fpTrim := arith.PolyTrim(fp)
	roots := arith.PolyRootsModP(fpTrim, p)
	dfp := arith.PolyDerivModP(fpTrim, p)

	allSimple := len(fpTrim)-1 == d || (fp[d] == 0 && d >= 1 && fp[d-1] != 0)
	if allSimple {
		for _, r := range roots {
			if arith.PolyEvalModP(dfp, r, p) == 0 {
				allSimple = false
				break
			}
		}
	}
	if allSimple {
		np := len(roots)
		if fp[d] == 0 {
			np++
		}
		return float64(np) * float64(p) / float64(p*p-1)
	}

	var total float64
	for _, r := range roots {
		total += liftExpectedVal(coeffs, r, 1, p, p)
	}
	if fp[d] == 0 {
		revCoeffs := make([]*big.Int, d+1)
		for i := 0; i <= d; i++ {
			revCoeffs[i] = coeffs[d-i]
		}
		total += liftExpectedVal(revCoeffs, 0, 1, p, p)
	}
	return total / float64(p+1)
}

func liftExpectedVal(coeffs []*big.Int, r uint64, level int, pk, p uint64) float64 {
	weight := math.Pow(float64(p), float64(-(level - 1)))
	dfModP := evalDerivModP(coeffs, r%p, p)
	if dfModP != 0 {
		return weight * float64(p) / float64(p-1)
	}
	if level >= 12 || pk > (^uint64(0))/p {
		return weight * float64(p) / float64(p-1)
	}
	nextPK := pk * p
	sum := weight
	for t := uint64(0); t < p; t++ {
		cand := r + t*pk
		if evalPolyBigMod(coeffs, cand, nextPK) == 0 {
			sum += liftExpectedVal(coeffs, cand, level+1, nextPK, p)
		}
	}
	return sum
}

func evalDerivModP(coeffs []*big.Int, x, p uint64) uint64 {
	var res uint64
	var bigP, r big.Int
	bigP.SetUint64(p)
	for i := len(coeffs) - 1; i >= 1; i-- {
		ci := r.Mod(coeffs[i], &bigP).Uint64()
		term := arith.MulMod64(ci, uint64(i)%p, p)
		res = (arith.MulMod64(res, x, p) + term) % p
	}
	return res
}

func evalPolyBigMod(coeffs []*big.Int, x, mod uint64) uint64 {
	var res uint64
	var bigM, r big.Int
	bigM.SetUint64(mod)
	for i := len(coeffs) - 1; i >= 0; i-- {
		ci := r.Mod(coeffs[i], &bigM).Uint64()
		prod := arith.MulMod64(res, x, mod)
		res = prod + ci
		if res >= mod || res < prod {
			res -= mod
		}
	}
	return res
}

// Cap skewness near the special-q lattice scale so |b| stays at least 1.
func maxSkewForParams(pr *Params) float64 {
	if !pr.UseLattice {
		return 1000
	}
	qMid := math.Sqrt(float64(pr.QMin) * float64(pr.QMax))
	if qMid < 10000 {
		qMid = 10000
	}
	return qMid * 1.5
}

// dickmanYieldScore returns u_F*ln(u_F) + u_G*ln(u_G) at typical sieve scale.
func dickmanYieldScore(p *Poly, pr *Params) float64 {
	logRMSF, logRMSG := logNormParts(p)
	var logA, logQ float64
	if pr.UseLattice {
		qTyp := math.Sqrt(float64(pr.QMin) * float64(pr.QMax))
		if qTyp < 2 {
			qTyp = 2
		}
		logQ = math.Log(qTyp)
		aTyp := math.Sqrt(float64(pr.IHalf) * float64(max(1, pr.JMax/2)) * qTyp)
		logA = math.Log(aTyp)
	} else {
		aTyp := math.Sqrt(float64(pr.LineA) * float64(max(1, pr.LineB/2)))
		logA = math.Log(aTyp)
	}

	lnBa := math.Log(float64(max(int64(100), pr.Ba)))
	lnBr := math.Log(float64(max(int64(100), pr.Br)))
	uF := math.Max(1.0, (logRMSF+float64(p.Degree)*logA-logQ+p.Alpha)/lnBa)
	uG := math.Max(1.0, (logRMSG+logA)/lnBr)
	return uF*math.Log(uF) + uG*math.Log(uG)
}

func logNormParts(p *Poly) (logRMSF, logRMSG float64) {
	d := p.Degree
	s := p.Skew
	if s <= 0 {
		s = 1
	}
	maxBits := 0
	for _, c := range p.Coeffs {
		if b := c.BitLen(); b > maxBits {
			maxBits = b
		}
	}
	shiftF := maxBits - 500
	if shiftF < 0 {
		shiftF = 0
	}
	var sumF float64
	for i := 0; i <= d; i++ {
		ci := scaledFloat(p.Coeffs[i], shiftF)
		for j := 0; j <= d; j++ {
			if (i+j)&1 != 0 {
				continue
			}
			cj := scaledFloat(p.Coeffs[j], shiftF)
			sumF += ci * cj * math.Pow(s, float64(i+j-d)) / float64((i+j+1)*(2*d-i-j+1))
		}
	}
	if sumF <= 0 {
		sumF = 1e-300
	}
	logRMSF = 0.5*math.Log(sumF) + float64(shiftF)*math.Ln2

	maxGBits := max(p.A1.BitLen(), p.A0.BitLen())
	shiftG := maxGBits - 500
	if shiftG < 0 {
		shiftG = 0
	}
	a1f := scaledFloat(p.A1, shiftG)
	a0f := scaledFloat(p.A0, shiftG)
	sumG := (a1f*a1f*s + a0f*a0f/s) / 3.0
	if sumG <= 0 {
		sumG = 1e-300
	}
	logRMSG = 0.5*math.Log(sumG) + float64(shiftG)*math.Ln2
	return logRMSF, logRMSG
}

func logNormScore(p *Poly) float64 {
	lf, lg := logNormParts(p)
	return lf + lg
}

// SelectPoly searches skewed base-m and Kleinjung polynomial pairs, size-optimizes
// them by translation and rotation, and runs a root sieve over (u, v) rotations.
func SelectPoly(n *big.Int, degree int) (*Poly, *big.Int) {
	pr := defaultParams(n)
	if degree <= 0 {
		degree = pr.Degree
	} else {
		cp := *pr
		cp.Degree = degree
		pr = &cp
	}
	bits := n.BitLen()
	maxSkew := maxSkewForParams(pr)

	if bits < 130 {
		return selectPolySmall(n, degree, pr, maxSkew)
	}

	var rawCands []*Poly
	for _, cd := range skewedCDPool(bits, degree) {
		p, fac := fixedCDBaseMPoly(n, degree, cd)
		if fac != nil {
			return nil, fac
		}
		if p != nil {
			rawCands = append(rawCands, p)
		}
	}

	if bits >= 170 && degree >= 3 {
		kjCands, fac := kleinjungCandidates(n, degree)
		if fac != nil {
			return nil, fac
		}
		rawCands = append(rawCands, kjCands...)
	}

	if len(rawCands) == 0 {
		return BaseMPoly(n, degree, 1)
	}

	type sizedCand struct {
		poly      *Poly
		sizeScore float64
	}
	sized := make([]sizedCand, 0, len(rawCands))
	for _, p := range rawCands {
		opt, fac := sizeOptimizePoly(p, n, maxSkew)
		if fac != nil {
			return nil, fac
		}
		if opt == nil {
			continue
		}
		opt.Alpha = projectiveAlphaBonus(opt.Coeffs[degree])
		sized = append(sized, sizedCand{poly: opt, sizeScore: dickmanYieldScore(opt, pr)})
	}
	if len(sized) == 0 {
		return BaseMPoly(n, degree, 1)
	}

	sort.Slice(sized, func(i, j int) bool { return sized[i].sizeScore < sized[j].sizeScore })
	topN := 24
	if len(sized) < topN {
		topN = len(sized)
	}
	sized = sized[:topN]

	nw := Workers
	if nw <= 0 {
		nw = runtime.GOMAXPROCS(0)
	}
	if nw > topN {
		nw = topN
	}
	if nw < 1 {
		nw = 1
	}

	results := make([]*Poly, topN)
	factors := make([]*big.Int, topN)
	var nextIdx atomic.Int64
	var wg sync.WaitGroup
	pks := buildRootSievePrimes(97)

	for w := 0; w < nw; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for {
				idx := int(nextIdx.Add(1) - 1)
				if idx >= topN {
					return
				}
				best, fac := rootSieveOptimize(sized[idx].poly, n, pr, maxSkew, pks)
				if fac != nil {
					factors[idx] = fac
					return
				}
				results[idx] = best
			}
		}()
	}
	wg.Wait()

	for _, fac := range factors {
		if fac != nil {
			return nil, fac
		}
	}

	var best *Poly
	bestScore := math.Inf(1)
	for _, cand := range results {
		if cand == nil {
			continue
		}
		sc := dickmanYieldScore(cand, pr)
		if sc < bestScore {
			bestScore = sc
			best = cand
		}
	}
	if best == nil {
		return BaseMPoly(n, degree, 1)
	}
	return best, nil
}

func selectPolySmall(n *big.Int, degree int, pr *Params, maxSkew float64) (*Poly, *big.Int) {
	var best *Poly
	bestScore := math.Inf(1)
	for _, cd := range []uint64{1, 2, 6, 12, 30, 60} {
		base, fac := BaseMPoly(n, degree, cd)
		if fac != nil {
			return nil, fac
		}
		if base == nil {
			continue
		}
		base.Skew = optimizeSkewnessBounded(base, maxSkew)
		sc := dickmanYieldScore(base, pr)
		if sc < bestScore {
			bestScore = sc
			best = base
		}
	}
	if best == nil {
		return BaseMPoly(n, degree, 1)
	}
	return best, nil
}

func chooseDegree(n *big.Int) int {
	return defaultParams(n).Degree
}

func skewedCDPool(bits, d int) []uint64 {
	maxLog2 := bits / (d + 3)
	if maxLog2 > 40 {
		maxLog2 = 40
	}
	if maxLog2 < 10 {
		maxLog2 = 10
	}
	limit := uint64(1) << uint(maxLog2)
	cands := []uint64{1}
	for _, p := range []uint64{2, 3, 5, 7, 11, 13} {
		n := len(cands)
		for i := 0; i < n; i++ {
			v := cands[i]
			for v <= limit/p {
				v *= p
				cands = append(cands, v)
			}
		}
	}
	sort.Slice(cands, func(i, j int) bool { return cands[i] < cands[j] })

	const targetCount = 64
	if len(cands) <= targetCount {
		return cands
	}
	out := make([]uint64, 0, targetCount)
	out = append(out, 1, 6, 12, 30, 60, 210, 420, 840, 2520)
	seen := make(map[uint64]bool, targetCount)
	for _, v := range out {
		seen[v] = true
	}
	maxL := math.Log(float64(cands[len(cands)-1]))
	minL := math.Log(60.0)
	steps := targetCount - len(out)
	for s := 0; s < steps; s++ {
		want := math.Exp(minL + (maxL-minL)*float64(s)/float64(max(1, steps-1)))
		idx := sort.Search(len(cands), func(i int) bool { return float64(cands[i]) >= want })
		for dIdx := 0; dIdx < len(cands); dIdx++ {
			if i := idx + dIdx; i < len(cands) && !seen[cands[i]] {
				seen[cands[i]] = true
				out = append(out, cands[i])
				break
			}
			if i := idx - 1 - dIdx; i >= 0 && !seen[cands[i]] {
				seen[cands[i]] = true
				out = append(out, cands[i])
				break
			}
		}
	}
	return out
}

// fixedCDBaseMPoly keeps c_d fixed and chooses nearest m = round((n/c_d)^(1/d))
// so c_d stays smooth and |c_{d-1}| <= (d/2)*c_d.
func fixedCDBaseMPoly(n *big.Int, d int, cd uint64) (*Poly, *big.Int) {
	if cd == 0 {
		cd = 1
	}
	cdBig := new(big.Int).SetUint64(cd)
	q := new(big.Int).Quo(n, cdBig)
	if q.Sign() <= 0 {
		return nil, nil
	}
	mFloor := arith.NthRoot(q, d)
	mCeil := new(big.Int).Add(mFloor, bigOne)
	dBig := big.NewInt(int64(d))

	rFloor := new(big.Int).Sub(n, new(big.Int).Mul(cdBig, new(big.Int).Exp(mFloor, dBig, nil)))
	rCeil := new(big.Int).Sub(n, new(big.Int).Mul(cdBig, new(big.Int).Exp(mCeil, dBig, nil)))
	m := mFloor
	rem := rFloor
	if new(big.Int).Abs(rCeil).Cmp(new(big.Int).Abs(rFloor)) < 0 {
		m = mCeil
		rem = rCeil
	}
	if m.Cmp(big.NewInt(4)) < 0 {
		return nil, nil
	}

	coeffs := make([]*big.Int, d+1)
	coeffs[d] = cdBig
	halfM := new(big.Int).Rsh(m, 1)
	negHalfM := new(big.Int).Neg(halfM)
	cur := new(big.Int).Set(rem)
	for i := 0; i < d-1; i++ {
		quot, r := new(big.Int).QuoRem(cur, m, new(big.Int))
		if r.Cmp(halfM) > 0 {
			r.Sub(r, m)
			quot.Add(quot, bigOne)
		} else if r.Cmp(negHalfM) < 0 {
			r.Add(r, m)
			quot.Sub(quot, bigOne)
		}
		coeffs[i] = r
		cur = quot
	}
	coeffs[d-1] = cur

	p := &Poly{
		Degree: d,
		Coeffs: coeffs,
		A1:     big.NewInt(1),
		A0:     new(big.Int).Set(m),
		M:      new(big.Int).Mod(m, n),
	}
	return normalizePoly(p, n)
}

// kleinjungCandidates builds non-monic g(x) = a_1*x - a_0 pairs with c_d*a_0^d ≡ N (mod a_1^2),
// following Kleinjung, Math. Comp. 75 (2006), Section 2.1.
func kleinjungCandidates(n *big.Int, d int) ([]*Poly, *big.Int) {
	bits := n.BitLen()
	maxA1Bits := (bits/d - 16) / 2
	if maxA1Bits < 12 {
		return nil, nil
	}
	if maxA1Bits > 29 {
		maxA1Bits = 29
	}

	var splitCandPrimes []uint64
	var bigP, rTmp big.Int
	for _, sp := range arith.PrimesUpTo(600) {
		if sp < 11 || int64(d)%sp == 0 {
			continue
		}
		up := uint64(sp)
		if rTmp.Mod(n, bigP.SetUint64(up)).Sign() == 0 {
			return nil, new(big.Int).SetUint64(up)
		}
		splitCandPrimes = append(splitCandPrimes, up)
	}

	cds := []uint64{
		60, 120, 210, 420, 840, 1260, 2520, 5040, 7560, 10080,
		15120, 25200, 27720, 50400, 60060, 110880, 180180, 360360, 720720,
	}

	type primeRoots struct {
		p  uint64
		p2 uint64
		rs []uint64
	}

	var out []*Poly
	dBig := big.NewInt(int64(d))
	for _, cd := range cds {
		var prs []primeRoots
		for _, p := range splitCandPrimes {
			if cd%p == 0 {
				continue
			}
			rs := rootsNthModP2(cd, n, d, p)
			if len(rs) == d {
				prs = append(prs, primeRoots{p: p, p2: p * p, rs: rs})
				if len(prs) >= 8 {
					break
				}
			}
		}
		numFactors := 3
		if (d == 4 && maxA1Bits >= 21 || d == 3 && maxA1Bits >= 20) && len(prs) >= 4 {
			numFactors = 4
		}
		if len(prs) < numFactors {
			continue
		}

		for offset := 0; offset+numFactors <= len(prs) && offset < 2; offset++ {
			tuple := prs[offset : offset+numFactors]
			a1 := uint64(1)
			tooBig := false
			for _, pr := range tuple {
				if a1 > (uint64(1)<<uint(maxA1Bits))/pr.p {
					tooBig = true
					break
				}
				a1 *= pr.p
			}
			if tooBig || a1 <= 1 {
				continue
			}
			a1Sq := a1 * a1

			rems := []uint64{0}
			mod := uint64(1)
			for _, pr := range tuple {
				invM, ok := arith.InvMod64(mod%pr.p2, pr.p2)
				if !ok {
					rems = nil
					break
				}
				next := make([]uint64, 0, len(rems)*len(pr.rs))
				for _, r0 := range rems {
					r0m := r0 % pr.p2
					for _, r1 := range pr.rs {
						diff := (r1 + pr.p2 - r0m) % pr.p2
						k := arith.MulMod64(diff, invM, pr.p2)
						next = append(next, r0+k*mod)
					}
				}
				rems = next
				mod *= pr.p2
			}
			if len(rems) == 0 {
				continue
			}

			cdBig := new(big.Int).SetUint64(cd)
			targetA0 := arith.NthRoot(new(big.Int).Quo(n, cdBig), d)
			if targetA0.BitLen() <= maxA1Bits*2+4 {
				continue
			}
			a1SqBig := new(big.Int).SetUint64(a1Sq)
			tMod := new(big.Int).Mod(targetA0, a1SqBig).Uint64()

			type rootCand struct {
				a0     *big.Int
				M      *big.Int
				td1    *big.Int
				remAbs *big.Int
			}
			var bestRoots []rootCand
			var a0, a0d, a0d1, num, M, remCheck, td1, rem3, remAbs big.Int

			for _, r := range rems {
				var delta int64
				if r >= tMod {
					dPos := r - tMod
					if dPos > a1Sq/2 {
						delta = -int64(a1Sq - dPos)
					} else {
						delta = int64(dPos)
					}
				} else {
					dNeg := tMod - r
					if dNeg > a1Sq/2 {
						delta = int64(a1Sq - dNeg)
					} else {
						delta = -int64(dNeg)
					}
				}
				a0.Add(targetA0, big.NewInt(delta))
				if a0.Sign() <= 0 {
					continue
				}
				a0ModA1 := rTmp.Mod(&a0, bigP.SetUint64(a1)).Uint64()
				if arith.Gcd64(a0ModA1, a1) != 1 {
					continue
				}

				a0d1.Exp(&a0, big.NewInt(int64(d-1)), nil)
				a0d.Mul(&a0d1, &a0)
				num.Mul(cdBig, &a0d)
				num.Sub(n, &num)
				M.QuoRem(&num, a1SqBig, &remCheck)
				if remCheck.Sign() != 0 {
					continue
				}

				roundQuoRem(&td1, &rem3, &M, &a0d1)
				remAbs.Abs(&rem3)

				if len(bestRoots) < 3 {
					bestRoots = append(bestRoots, rootCand{
						a0:     new(big.Int).Set(&a0),
						M:      new(big.Int).Set(&M),
						td1:    new(big.Int).Set(&td1),
						remAbs: new(big.Int).Set(&remAbs),
					})
					sort.Slice(bestRoots, func(i, j int) bool { return bestRoots[i].remAbs.Cmp(bestRoots[j].remAbs) < 0 })
				} else if remAbs.Cmp(bestRoots[2].remAbs) < 0 {
					bestRoots[2] = rootCand{
						a0:     new(big.Int).Set(&a0),
						M:      new(big.Int).Set(&M),
						td1:    new(big.Int).Set(&td1),
						remAbs: new(big.Int).Set(&remAbs),
					}
					sort.Slice(bestRoots, func(i, j int) bool { return bestRoots[i].remAbs.Cmp(bestRoots[j].remAbs) < 0 })
				}
			}

			a1Big := new(big.Int).SetUint64(a1)
			invA1ModN := new(big.Int).ModInverse(a1Big, n)
			if invA1ModN == nil {
				var g big.Int
				g.GCD(nil, nil, a1Big, n)
				if g.Cmp(bigOne) > 0 && g.Cmp(n) < 0 {
					return nil, &g
				}
				continue
			}

			for _, rc := range bestRoots {
				p := expandKleinjungRoot(n, d, cdBig, a1, a1Big, rc.a0, rc.M, rc.td1, invA1ModN, dBig)
				if p == nil {
					continue
				}
				normed, fac := normalizePoly(p, n)
				if fac != nil {
					return nil, fac
				}
				if normed != nil {
					out = append(out, normed)
				}
			}
		}
	}
	return out, nil
}

func rootsNthModP2(cd uint64, n *big.Int, d int, p uint64) []uint64 {
	var bigP, rem big.Int
	nModP := rem.Mod(n, bigP.SetUint64(p)).Uint64()
	if nModP == 0 {
		return nil
	}
	invCD := arith.InvMod(cd%p, p)
	rhs := arith.MulMod64(nModP, invCD, p)
	fp := make([]uint64, d+1)
	fp[0] = (p - rhs) % p
	fp[d] = 1
	roots := arith.PolyRootsModP(fp, p)
	if len(roots) != d {
		return nil
	}

	p2 := p * p
	nModP2 := rem.Mod(n, bigP.SetUint64(p2)).Uint64()
	cdP2 := cd % p2
	cdP := cd % p
	out := make([]uint64, 0, d)
	for _, x := range roots {
		if x == 0 {
			return nil
		}
		xd := arith.PowMod64(x, uint64(d), p2)
		val := (arith.MulMod64(cdP2, xd, p2) + p2 - nModP2) % p2
		deriv := arith.MulMod64(uint64(d)%p*cdP%p, arith.PowMod64(x, uint64(d-1), p), p)
		if deriv == 0 {
			return nil
		}
		k := arith.MulMod64((p-(val/p))%p, arith.InvMod(deriv, p), p)
		out = append(out, x+k*p)
	}
	return out
}

func roundQuoRem(q, r, a, b *big.Int) {
	q.QuoRem(a, b, r)
	var twoR big.Int
	twoR.Lsh(r, 1)
	if twoR.Sign() >= 0 {
		if twoR.Cmp(b) > 0 {
			q.Add(q, bigOne)
			r.Sub(r, b)
		}
	} else {
		var negB big.Int
		negB.Neg(b)
		if twoR.Cmp(&negB) < 0 {
			q.Sub(q, bigOne)
			r.Add(r, b)
		}
	}
}

func modBigUint64(x *big.Int, m uint64) uint64 {
	var bigM, r big.Int
	return r.Mod(x, bigM.SetUint64(m)).Uint64()
}

func expandKleinjungRoot(n *big.Int, d int, cdBig *big.Int, a1 uint64, a1Big, a0, M, td1, invA1ModN, _ *big.Int) *Poly {
	a0Pows := make([]*big.Int, d)
	a0Pows[0] = big.NewInt(1)
	for i := 1; i < d; i++ {
		a0Pows[i] = new(big.Int).Mul(a0Pows[i-1], a0)
	}

	coeffs := make([]*big.Int, d+1)
	coeffs[d] = new(big.Int).Set(cdBig)
	coeffs[d-1] = new(big.Int).Mul(a1Big, td1)

	Q := new(big.Int).Sub(M, new(big.Int).Mul(td1, a0Pows[d-1]))
	a0Mod := modBigUint64(a0, a1)
	invA0, ok := arith.InvMod64(a0Mod, a1)
	if !ok {
		return nil
	}

	var approx, remTmp, checkRem big.Int
	for i := d - 2; i >= 1; i-- {
		roundQuoRem(&approx, &remTmp, Q, a0Pows[i])
		qMod := modBigUint64(Q, a1)
		ciMod := arith.MulMod64(qMod, arith.PowMod64(invA0, uint64(i), a1), a1)
		approxMod := modBigUint64(&approx, a1)
		var shift int64
		if ciMod >= approxMod {
			dp := ciMod - approxMod
			if dp > a1/2 {
				shift = -int64(a1 - dp)
			} else {
				shift = int64(dp)
			}
		} else {
			dn := approxMod - ciMod
			if dn > a1/2 {
				shift = int64(a1 - dn)
			} else {
				shift = -int64(dn)
			}
		}
		ci := new(big.Int).Add(&approx, big.NewInt(shift))
		coeffs[i] = ci
		Q.Sub(Q, new(big.Int).Mul(ci, a0Pows[i]))
		Q.QuoRem(Q, a1Big, &checkRem)
		if checkRem.Sign() != 0 {
			return nil
		}
	}
	coeffs[0] = Q

	mMod := new(big.Int).Mul(a0, invA1ModN)
	mMod.Mod(mMod, n)
	return &Poly{
		Degree: d,
		Coeffs: coeffs,
		A1:     new(big.Int).Set(a1Big),
		A0:     new(big.Int).Set(a0),
		M:      mMod,
	}
}

func clonePoly(p *Poly) *Poly {
	coeffs := make([]*big.Int, len(p.Coeffs))
	for i, c := range p.Coeffs {
		coeffs[i] = new(big.Int).Set(c)
	}
	return &Poly{
		Degree: p.Degree,
		Coeffs: coeffs,
		A1:     new(big.Int).Set(p.A1),
		A0:     new(big.Int).Set(p.A0),
		M:      new(big.Int).Set(p.M),
		Skew:   p.Skew,
		Alpha:  p.Alpha,
	}
}

// translatePoly shifts x -> x + k, setting f_new(x) = f(x + k) and A0_new = A0 - k*A1.
func translatePoly(p *Poly, k *big.Int, n *big.Int) *Poly {
	if k.Sign() == 0 {
		return clonePoly(p)
	}
	out := clonePoly(p)
	d := out.Degree
	var prod big.Int
	for i := 0; i < d; i++ {
		for j := d - 1; j >= i; j-- {
			prod.Mul(k, out.Coeffs[j+1])
			out.Coeffs[j].Add(out.Coeffs[j], &prod)
		}
	}
	out.A0.Sub(out.A0, new(big.Int).Mul(k, out.A1))
	out.M.Sub(out.M, k)
	out.M.Mod(out.M, n)
	return out
}

func rotatePolyBig(base *Poly, u, v *big.Int) *Poly {
	out := clonePoly(base)
	if u.Sign() == 0 && v.Sign() == 0 {
		return out
	}
	t2 := new(big.Int).Mul(u, out.A1)
	t1 := new(big.Int).Sub(new(big.Int).Mul(v, out.A1), new(big.Int).Mul(u, out.A0))
	t0 := new(big.Int).Neg(new(big.Int).Mul(v, out.A0))
	out.Coeffs[2].Add(out.Coeffs[2], t2)
	out.Coeffs[1].Add(out.Coeffs[1], t1)
	out.Coeffs[0].Add(out.Coeffs[0], t0)
	return out
}

// sizeOptimizePoly minimizes the skewed norm via translation and 2x2 continuous rotation.
func sizeOptimizePoly(base *Poly, n *big.Int, maxSkew float64) (*Poly, *big.Int) {
	cur := clonePoly(base)
	cur.Skew = optimizeSkewnessBounded(cur, maxSkew)
	curScore := logNormScore(cur)
	d := cur.Degree

	for iter := 0; iter < 6; iter++ {
		improved := false

		denom := new(big.Int).Mul(big.NewInt(int64(d)), cur.Coeffs[d])
		if denom.Sign() != 0 {
			var kt, rem big.Int
			negC := new(big.Int).Neg(cur.Coeffs[d-1])
			if denom.Sign() < 0 {
				denom.Neg(denom)
				negC.Neg(negC)
			}
			roundQuoRem(&kt, &rem, negC, denom)
			if kt.Sign() != 0 {
				cand := translatePoly(cur, &kt, n)
				if cand.A0.Sign() > 0 {
					cand.Skew = optimizeSkewnessBounded(cand, maxSkew)
					if sc := logNormScore(cand); sc < curScore-1e-6 {
						cur, curScore, improved = cand, sc, true
					}
				}
			}
		}

		if d >= 3 {
			uBig, vBig, ok := solveBestUVBig(cur)
			if ok && (uBig.Sign() != 0 || vBig.Sign() != 0) {
				cand := rotatePolyBig(cur, uBig, vBig)
				cand.Skew = optimizeSkewnessBounded(cand, maxSkew)
				if sc := logNormScore(cand); sc < curScore-1e-6 {
					cur, curScore, improved = cand, sc, true
				}
			}
		}

		if !improved {
			break
		}
	}

	return normalizePoly(cur, n)
}

func solveBestUVBig(p *Poly) (uOut, vOut *big.Int, ok bool) {
	d := p.Degree
	s := p.Skew
	if s <= 0 {
		s = 1
	}
	maxBits := max(p.A0.BitLen(), p.A1.BitLen())
	for _, c := range p.Coeffs {
		if b := c.BitLen(); b > maxBits {
			maxBits = b
		}
	}
	shift := maxBits - 450
	if shift < 0 {
		shift = 0
	}
	cf := make([]float64, d+1)
	for i, c := range p.Coeffs {
		cf[i] = scaledFloat(c, shift)
	}
	a1f := scaledFloat(p.A1, shift)
	a0f := scaledFloat(p.A0, shift)

	du := make([]float64, d+1)
	dv := make([]float64, d+1)
	du[1], du[2] = -a0f, a1f
	dv[0], dv[1] = -a0f, a1f

	inner := func(x, y []float64) float64 {
		var tot float64
		for i := 0; i <= d; i++ {
			if x[i] == 0 {
				continue
			}
			for j := 0; j <= d; j++ {
				if (i+j)&1 != 0 || y[j] == 0 {
					continue
				}
				tot += x[i] * y[j] * math.Pow(s, float64(i+j-d)) / float64((i+j+1)*(2*d-i-j+1))
			}
		}
		return tot
	}

	a11 := inner(du, du)
	a12 := inner(du, dv)
	a22 := inner(dv, dv)
	det := a11*a22 - a12*a12
	if !(det > 0) || math.IsInf(det, 0) || math.IsNaN(det) {
		return nil, nil, false
	}
	b1 := -inner(cf, du)
	b2 := -inner(cf, dv)
	uf := (b1*a22 - b2*a12) / det
	vf := (a11*b2 - a12*b1) / det
	if math.IsNaN(uf) || math.IsInf(uf, 0) || math.IsNaN(vf) || math.IsInf(vf, 0) {
		return nil, nil, false
	}
	return floatToRoundedBig(uf), floatToRoundedBig(vf), true
}

func floatToRoundedBig(x float64) *big.Int {
	if math.Abs(x) < 0.5 {
		return new(big.Int)
	}
	if math.Abs(x) < 1e15 {
		return big.NewInt(int64(math.Round(x)))
	}
	bf := new(big.Float).SetPrec(128).SetFloat64(x)
	if x >= 0 {
		bf.Add(bf, big.NewFloat(0.5))
	} else {
		bf.Sub(bf, big.NewFloat(0.5))
	}
	out, _ := bf.Int(nil)
	return out
}

func projectiveAlphaBonus(cd *big.Int) float64 {
	var bonus float64
	var bigP, r big.Int
	for _, p := range []int64{2, 3, 5, 7, 11, 13, 17, 19, 23} {
		if r.Mod(cd, bigP.SetInt64(p)).Sign() == 0 {
			pf := float64(p)
			bonus -= pf / (pf*pf - 1) * math.Log(pf)
		}
	}
	return bonus
}

type rootSievePK struct {
	p      uint64
	pk     uint64
	weight float32
}

func buildRootSievePrimes(bound int64) []rootSievePK {
	var out []rootSievePK
	for _, p64 := range arith.PrimesUpTo(bound) {
		p := uint64(p64)
		pk := p
		level := 1
		for pk <= uint64(bound) {
			w := math.Log(float64(p)) / (float64(p+1) * math.Pow(float64(p), float64(level-1)))
			out = append(out, rootSievePK{p: p, pk: pk, weight: float32(w)})
			pk *= p
			level++
		}
	}
	return out
}

// rootSieveOptimize sieves (u, v) in [-U, U] x [-V, V] around base for small-prime root density,
// following Murphy (1999), Chapter 5.2.
func rootSieveOptimize(base *Poly, n *big.Int, pr *Params, maxSkew float64, pks []rootSievePK) (*Poly, *big.Int) {
	rmsF, _ := logNormParts(base)
	s := base.Skew
	if s < 1 {
		s = 1
	}
	logA0 := arith.Log2Big(new(big.Int).Abs(base.A0)) * math.Ln2
	maxU := int(math.Exp(rmsF + 0.4 + math.Log(s) - logA0))
	if maxU < 2 {
		maxU = 2
	}
	if maxU > 12 {
		maxU = 12
	}
	maxV := int(math.Exp(rmsF + 0.4 + 2*math.Log(s) - logA0))
	if maxV < 16 {
		maxV = 16
	}
	if maxV > 768 {
		maxV = 768
	}

	lenV := 2*maxV + 1
	gains := make([]float32, lenV)

	d := base.Degree
	numPK := len(pks)
	baseMod := make([][]uint64, numPK)
	a1Mod := make([]uint64, numPK)
	a0Mod := make([]uint64, numPK)
	var bigM, rTmp big.Int
	for idx, entry := range pks {
		bigM.SetUint64(entry.pk)
		cm := make([]uint64, d+1)
		for i, c := range base.Coeffs {
			cm[i] = rTmp.Mod(c, &bigM).Uint64()
		}
		baseMod[idx] = cm
		a1Mod[idx] = rTmp.Mod(base.A1, &bigM).Uint64()
		a0Mod[idx] = rTmp.Mod(base.A0, &bigM).Uint64()
	}

	type uvCand struct {
		u, v int64
		gain float32
	}
	const keepTop = 8
	top := make([]uvCand, 0, keepTop+1)
	top = append(top, uvCand{0, 0, 0})

	cu := make([]uint64, d+1)
	for u := int64(-maxU); u <= int64(maxU); u++ {
		clear(gains)
		for idx, entry := range pks {
			pk := entry.pk
			p := entry.p
			w := entry.weight
			copy(cu, baseMod[idx])
			um := modInt64(u, pk)
			a1m := a1Mod[idx]
			a0m := a0Mod[idx]
			cu[2] = (cu[2] + um*a1m) % pk
			sub1 := (um * a0m) % pk
			if cu[1] >= sub1 {
				cu[1] -= sub1
			} else {
				cu[1] = cu[1] + pk - sub1
			}

			pkInt := int(pk)
			vBias := maxV % pkInt
			for x := uint64(0); x < pk; x++ {
				gx := (a1m*x + pk - a0m) % pk
				var fx uint64
				for i := d; i >= 0; i-- {
					fx = (fx*x + cu[i]) % pk
				}
				if gx%p != 0 {
					invG, ok := arith.InvMod64(gx, pk)
					if !ok {
						continue
					}
					vMod := ((pk - fx) % pk) * invG % pk
					off := (int(vMod) + vBias) % pkInt
					for j := off; j < lenV; j += pkInt {
						gains[j] += w
					}
				} else {
					g := arith.Gcd64(gx, pk)
					if fx%g != 0 {
						continue
					}
					step := int(pk / g)
					if step == 1 {
						for j := 0; j < lenV; j++ {
							gains[j] += w
						}
					} else {
						stepU := uint64(step)
						invG, ok := arith.InvMod64((gx/g)%stepU, stepU)
						if !ok {
							continue
						}
						fxRed := (fx / g) % stepU
						vMod := ((stepU - fxRed) % stepU) * invG % stepU
						off := (int(vMod) + maxV%step) % step
						for j := off; j < lenV; j += step {
							gains[j] += w
						}
					}
				}
			}
		}

		for j, gVal := range gains {
			v := int64(j - maxV)
			if u == 0 && v == 0 {
				top[0].gain = gVal
				continue
			}
			if len(top) < keepTop || gVal > top[len(top)-1].gain {
				cand := uvCand{u: u, v: v, gain: gVal}
				if len(top) < keepTop {
					top = append(top, cand)
				} else {
					top[keepTop-1] = cand
				}
				for k := len(top) - 1; k > 1 && top[k].gain > top[k-1].gain; k-- {
					top[k], top[k-1] = top[k-1], top[k]
				}
			}
		}
	}

	var best *Poly
	bestScore := math.Inf(1)
	for _, c := range top {
		rot, fac := normalizePoly(rotatePoly(base, c.u, c.v), n)
		if fac != nil {
			return nil, fac
		}
		if rot == nil {
			continue
		}
		if !IsIrreducible(rot) {
			if fac := findFactorFromReducible(rot, n); fac != nil {
				return nil, fac
			}
			continue
		}
		rot.Skew = optimizeSkewnessBounded(rot, maxSkew)
		rot.Alpha = ComputeAlpha(rot, 200)
		sc := dickmanYieldScore(rot, pr)
		if sc < bestScore {
			bestScore = sc
			best = rot
		}
	}
	return best, nil
}

func normalizePoly(p *Poly, n *big.Int) (*Poly, *big.Int) {
	if p == nil || p.Coeffs[p.Degree].Sign() == 0 {
		return nil, nil
	}
	if p.Coeffs[p.Degree].Sign() < 0 {
		for i := range p.Coeffs {
			p.Coeffs[i].Neg(p.Coeffs[i])
		}
	}
	content := new(big.Int).Abs(p.Coeffs[0])
	for i := 1; i <= p.Degree; i++ {
		content.GCD(nil, nil, content, new(big.Int).Abs(p.Coeffs[i]))
	}
	if content.Sign() == 0 {
		return nil, nil
	}
	if content.Cmp(bigOne) > 0 {
		var g big.Int
		g.GCD(nil, nil, content, n)
		if g.Cmp(bigOne) > 0 && g.Cmp(n) < 0 {
			return nil, &g
		}
		for i := 0; i <= p.Degree; i++ {
			p.Coeffs[i].Quo(p.Coeffs[i], content)
		}
	}
	var gLin big.Int
	gLin.GCD(nil, nil, new(big.Int).Abs(p.A1), new(big.Int).Abs(p.A0))
	if gLin.Cmp(bigOne) > 0 {
		p.A1.Quo(p.A1, &gLin)
		p.A0.Quo(p.A0, &gLin)
	}
	if p.Coeffs[0].Sign() == 0 {
		var g big.Int
		g.GCD(nil, nil, p.A0, n)
		if g.Cmp(bigOne) > 0 && g.Cmp(n) < 0 {
			return nil, &g
		}
		return nil, nil
	}
	return p, nil
}

func rotatePoly(base *Poly, u, v int64) *Poly {
	return rotatePolyBig(base, big.NewInt(u), big.NewInt(v))
}

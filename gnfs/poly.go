// Package gnfs implements the General Number Field Sieve (GNFS) integer factorization algorithm.
package gnfs

import (
	"math"
	"math/big"
	"sort"

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

// OptimizeSkewness finds the skewness s >= 1 minimizing the combined skewed
// L2 norm log(RMS(F)) + log(RMS(G)) over the unit rectangle.
func OptimizeSkewness(p *Poly) float64 {
	d := p.Degree
	// Represent coefficients as float64 scaled by 2^(-shift) so nothing overflows.
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

	// Starting estimate |c_0 / c_d|^(1/d).
	s0 := 1.0
	if cf[d] != 0 && cf[0] != 0 {
		s0 = math.Pow(math.Abs(cf[0]/cf[d]), 1.0/float64(d))
	}
	if s0 < 1 {
		s0 = 1
	}

	evalLogNorm := func(logS float64) float64 {
		s := math.Exp(logS)
		// Compute I_F(s) / 4 = sum_{(i+j) even} c_i*c_j*s^(i+j-d) / ((i+j+1)*(2d-i-j+1))
		var sumF float64
		for i := 0; i <= d; i++ {
			for j := 0; j <= d; j++ {
				if (i+j)&1 != 0 {
					continue
				}
				term := cf[i] * cf[j] * math.Pow(s, float64(i+j-d)) / float64((i+j+1)*(2*d-i-j+1))
				sumF += term
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

	// Golden section search on [max(0, logS0 - 8), logS0 + 12].
	center := math.Log(s0)
	lo := math.Max(0, center-8)
	hi := center + 12
	phi := (math.Sqrt(5) - 1) / 2
	x1 := hi - phi*(hi-lo)
	x2 := lo + phi*(hi-lo)
	f1, f2 := evalLogNorm(x1), evalLogNorm(x2)
	for iter := 0; iter < 45; iter++ {
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

// ComputeAlpha computes Murphy's alpha(f) over primes <= bound:
//
//	alpha(f) = sum_{p <= bound} (1/(p-1) - E[v_p(F(a, b))]) * log(p)
//
// Negative alpha means F has more small prime factors than a random integer.
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

// expectedValP computes E[v_p(F(a, b))] over coprime (a, b) for prime p.
func expectedValP(coeffs []*big.Int, fp []uint64, p uint64, d int) float64 {
	// Fast check: if f mod p is squarefree and p does not divide c_d,
	// or if all roots (including projective) are simple:
	// Every simple affine root contributes p/(p^2-1); projective root (p | c_d, p !| c_{d-1}) contributes p/(p^2-1).
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
			np++ // simple projective root
		}
		return float64(np) * float64(p) / float64(p*p-1)
	}

	// General case (multiple roots or multiple projective root): lift mod p^k.
	var total float64
	// Affine classes: roots mod p.
	for _, r := range roots {
		total += liftExpectedVal(coeffs, r, 1, p, p)
	}
	// Projective class (1 : 0): evaluate reversed polynomial at y ≡ 0 (mod p).
	if fp[d] == 0 {
		revCoeffs := make([]*big.Int, d+1)
		for i := 0; i <= d; i++ {
			revCoeffs[i] = coeffs[d-i]
		}
		total += liftExpectedVal(revCoeffs, 0, 1, p, p)
	}
	return total / float64(p+1)
}

// liftExpectedVal computes sum_{k >= level} Pr(p^k | f(x) | x ≡ r mod p^level)
// weighted so that level=1 corresponds to conditioning on x ≡ r mod p.
func liftExpectedVal(coeffs []*big.Int, r uint64, level int, pk, p uint64) float64 {
	// Weight of matching r mod p^level within x ≡ r_1 mod p is p^(-(level-1)).
	weight := math.Pow(float64(p), float64(-(level - 1)))

	// Check Hensel condition: if p !| f'(r), every lift is unique at all higher levels.
	dfModP := evalDerivModP(coeffs, r%p, p)
	if dfModP != 0 {
		return weight * float64(p) / float64(p-1)
	}

	// Stop if p^(level+1) overflows uint64 or weight is negligible.
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

// SelectPoly chooses a degree and searches for a good polynomial pair (f, g)
// for n using chosen smooth leading coefficients c_d and rotations f + (ux+v)g.
// If a non-trivial factor of n is found during selection, it returns nil, factor.
func SelectPoly(n *big.Int, degree int) (*Poly, *big.Int) {
	if degree <= 0 {
		degree = chooseDegree(n)
	}

	// Candidate smooth leading coefficients.
	cds := smoothCDs(n.BitLen())

	var best *Poly
	bestScore := math.Inf(1)

	scorePoly := func(p *Poly) float64 {
		// Combined skewed log-norm of F and G plus alpha(f).
		return logNormScore(p) + p.Alpha
	}

	for _, cd := range cds {
		base, fac := BaseMPoly(n, degree, cd)
		if fac != nil {
			return nil, fac
		}
		if base == nil {
			continue
		}
		sc := scorePoly(base)
		if sc < bestScore {
			bestScore = sc
			best = base
		}

		// Try small rotations f + (u*x + v) * (x - m).
		if degree >= 3 && n.BitLen() >= 140 {
			rotated, fac := bestRotation(base, n)
			if fac != nil {
				return nil, fac
			}
			if rotated != nil {
				rsc := scorePoly(rotated)
				if rsc < bestScore {
					bestScore = rsc
					best = rotated
				}
			}
		}
	}
	if best == nil {
		// Fallback to cd = 1.
		return BaseMPoly(n, degree, 1)
	}
	return best, nil
}

// chooseDegree returns the degree from the parameter table, so the size
// thresholds live in one place.
func chooseDegree(n *big.Int) int {
	return defaultParams(n).Degree
}

func smoothCDs(bits int) []uint64 {
	if bits < 150 {
		return []uint64{1, 2, 6, 12, 30, 60}
	}
	// Smooth products of 2, 3, 5, 7, 11.
	cands := []uint64{1}
	primes := []uint64{2, 3, 5, 7, 11}
	limit := uint64(5000)
	if bits >= 280 {
		limit = 50000
	}
	for _, p := range primes {
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
	if len(cands) > 48 {
		// Keep a spread of sizes (first 24 and 24 evenly spaced).
		picked := make([]uint64, 0, 48)
		picked = append(picked, cands[:24]...)
		step := (len(cands) - 24) / 24
		for i := 24; i < len(cands) && len(picked) < 48; i += step {
			picked = append(picked, cands[i])
		}
		return picked
	}
	return cands
}

func logNormScore(p *Poly) float64 {
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
	logRMSF := 0.5*math.Log(sumF) + float64(shiftF)*math.Ln2

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
	logRMSG := 0.5*math.Log(sumG) + float64(shiftG)*math.Ln2
	return logRMSF + logRMSG
}

func bestRotation(base *Poly, n *big.Int) (*Poly, *big.Int) {
	m := base.A0 // A1 == 1 here
	var best *Poly
	bestScore := logNormScore(base) + base.Alpha

	type cand struct {
		u, v  int64
		alpha float64
	}
	var top []cand
	for u := int64(-2); u <= 2; u++ {
		for v := int64(-8); v <= 8; v++ {
			if u == 0 && v == 0 {
				continue
			}
			rot, fac := normalizePoly(rotatePoly(base, u, v), n)
			if fac != nil {
				return nil, fac
			}
			if rot == nil {
				continue
			}
			a := ComputeAlpha(rot, 80)
			top = append(top, cand{u, v, a})
		}
	}
	sort.Slice(top, func(i, j int) bool { return top[i].alpha < top[j].alpha })
	if len(top) > 6 {
		top = top[:6]
	}

	for _, c := range top {
		rot, fac := normalizePoly(rotatePoly(base, c.u, c.v), n)
		if fac != nil {
			return nil, fac
		}
		if rot == nil || rot.Coeffs[0].Sign() == 0 {
			if rot == nil {
				continue
			}
			var g big.Int
			g.GCD(nil, nil, m, n)
			if g.Cmp(big.NewInt(1)) > 0 && g.Cmp(n) < 0 {
				return nil, &g
			}
			continue
		}
		if !IsIrreducible(rot) {
			if fac := findFactorFromReducible(rot, n); fac != nil {
				return nil, fac
			}
			continue
		}
		rot.Skew = OptimizeSkewness(rot)
		rot.Alpha = ComputeAlpha(rot, 200)
		sc := logNormScore(rot) + rot.Alpha
		if sc < bestScore {
			bestScore = sc
			best = rot
		}
	}
	return best, nil
}

func normalizePoly(p *Poly, n *big.Int) (*Poly, *big.Int) {
	one := big.NewInt(1)
	content := new(big.Int).Abs(p.Coeffs[0])
	for i := 1; i <= p.Degree; i++ {
		content.GCD(nil, nil, content, new(big.Int).Abs(p.Coeffs[i]))
	}
	if content.Sign() == 0 {
		return nil, nil
	}
	if content.Cmp(one) > 0 {
		var g big.Int
		g.GCD(nil, nil, content, n)
		if g.Cmp(one) > 0 && g.Cmp(n) < 0 {
			return nil, &g
		}
		for i := 0; i <= p.Degree; i++ {
			p.Coeffs[i].Quo(p.Coeffs[i], content)
		}
	}
	return p, nil
}

// rotatePoly returns base + (u*x + v)*(A1*x - A0).
func rotatePoly(base *Poly, u, v int64) *Poly {
	d := base.Degree
	coeffs := make([]*big.Int, d+1)
	for i, c := range base.Coeffs {
		coeffs[i] = new(big.Int).Set(c)
	}
	// (u*x + v)*(A1*x - A0) = u*A1*x^2 + (v*A1 - u*A0)*x - v*A0
	bu := big.NewInt(u)
	bv := big.NewInt(v)
	t2 := new(big.Int).Mul(bu, base.A1)
	t1 := new(big.Int).Sub(new(big.Int).Mul(bv, base.A1), new(big.Int).Mul(bu, base.A0))
	t0 := new(big.Int).Neg(new(big.Int).Mul(bv, base.A0))
	coeffs[2].Add(coeffs[2], t2)
	coeffs[1].Add(coeffs[1], t1)
	coeffs[0].Add(coeffs[0], t0)
	return &Poly{
		Degree: d,
		Coeffs: coeffs,
		A1:     new(big.Int).Set(base.A1),
		A0:     new(big.Int).Set(base.A0),
		M:      new(big.Int).Set(base.M),
	}
}

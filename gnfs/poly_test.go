package gnfs

import (
	"math"
	"math/big"
	"math/rand"
	"testing"

	"github.com/dgryski/go-factor/internal/arith"
	"github.com/dgryski/go-factor/internal/testutil"
)

func evalMonicBig(f []*big.Int, x *big.Int) *big.Int {
	res := big.NewInt(0)
	for i := len(f) - 1; i >= 0; i-- {
		res.Mul(res, x)
		res.Add(res, f[i])
	}
	return res
}

func normInAlphaHat(fHat []*big.Int, a, b, cd int64) *big.Int {
	// Compute resultant / characteristic norm of (cd*a - b*alpha_hat) mod fHat.
	// Using Res_y(fHat(y), cd*a - b*y): since fHat is monic of degree d,
	// N(cd*a - b*alpha_hat) = b^d * fHat(cd*a / b) = sum_i cHat_i * (cd*a)^i * b^{d-i}.
	d := len(fHat) - 1
	cda := big.NewInt(cd * a)
	bigB := big.NewInt(b)
	res := big.NewInt(0)
	aPow := big.NewInt(1)
	bPows := make([]*big.Int, d+1)
	bPows[0] = big.NewInt(1)
	for i := 1; i <= d; i++ {
		bPows[i] = new(big.Int).Mul(bPows[i-1], bigB)
	}
	for i := 0; i <= d; i++ {
		term := new(big.Int).Mul(fHat[i], aPow)
		term.Mul(term, bPows[d-i])
		res.Add(res, term)
		aPow.Mul(aPow, cda)
	}
	return res
}

func TestPolyIdentities(t *testing.T) {
	rnd := rand.New(rand.NewSource(1))
	for _, bits := range []int{50, 80, 120} {
		n := new(big.Int).Mul(testutil.RandomPrime(rnd, bits), testutil.RandomPrime(rnd, bits))
		for _, deg := range []int{3, 4, 5} {
			for _, cd := range []uint64{1, 6, 30} {
				p, fac := BaseMPoly(n, deg, cd)
				if fac != nil {
					continue
				}
				if p == nil {
					t.Fatalf("BaseMPoly(%v, %d, %d) returned nil", n, deg, cd)
				}
				// Check f(m) == 0 (mod N)
				fm := evalMonicBig(p.Coeffs, p.M)
				if new(big.Int).Mod(fm, n).Sign() != 0 {
					t.Fatalf("f(m) != 0 mod N")
				}
				// Check f_hat(m_hat) == 0 (mod N)
				fHat, _ := p.Monic()
				cdBig := p.Coeffs[p.Degree]
				mHat := new(big.Int).Mul(cdBig, p.M)
				mHat.Mod(mHat, n)
				if new(big.Int).Mod(evalMonicBig(fHat, mHat), n).Sign() != 0 {
					t.Fatalf("f_hat(m_hat) != 0 mod N")
				}
				// Check c_d^(d-1) * F(a, b) == N(c_d * a - b * alpha_hat)
				cdInt := cdBig.Int64()
				cdPow := new(big.Int).Exp(cdBig, big.NewInt(int64(deg-1)), nil)
				for k := 0; k < 10; k++ {
					a := int64(rnd.Intn(201) - 100)
					b := int64(1 + rnd.Intn(100))
					lhs := new(big.Int).Mul(cdPow, p.EvalF(a, b))
					rhs := normInAlphaHat(fHat, a, b, cdInt)
					if lhs.Cmp(rhs) != 0 {
						t.Fatalf("norm mismatch for (a,b)=(%d,%d): lhs=%v rhs=%v", a, b, lhs, rhs)
					}
				}
			}
		}
	}
}

func TestIsIrreducibleSwinnertonDyer(t *testing.T) {
	// x^4 + 1 is irreducible over Z, but reducible mod EVERY prime!
	// x^4 - 10x^2 + 1 is Swinnerton-Dyer S_3, also reducible mod every prime!
	for _, coeffs := range [][]*big.Int{
		{big.NewInt(1), big.NewInt(0), big.NewInt(0), big.NewInt(0), big.NewInt(1)},
		{big.NewInt(1), big.NewInt(0), big.NewInt(-10), big.NewInt(0), big.NewInt(1)},
	} {
		p := &Poly{Degree: 4, Coeffs: coeffs, A1: big.NewInt(1), A0: big.NewInt(1)}
		if !IsIrreducible(p) {
			t.Fatalf("expected %v to be irreducible", coeffs)
		}
	}
	// Reducible degree 4: (x^2 + 1)(x^2 + 2) = x^4 + 3x^2 + 2
	red := &Poly{
		Degree: 4,
		Coeffs: []*big.Int{big.NewInt(2), big.NewInt(0), big.NewInt(3), big.NewInt(0), big.NewInt(1)},
		A1:     big.NewInt(1),
		A0:     big.NewInt(1),
	}
	if IsIrreducible(red) {
		t.Fatal("(x^2+1)(x^2+2) reported irreducible")
	}
}

func TestAlphaBruteForce(t *testing.T) {
	// Test expectedValP and ComputeAlpha against brute force over coprime (a, b).
	// Include a polynomial with a repeated root mod 2 and mod 3, and projective root at 2.
	// f(x) = 2*x^3 + 3*x^2 - 6*x + 9 (c_3=2, mod 2 has projective root, mod 3 is 2x^3 so triple root at 0).
	p := &Poly{
		Degree: 3,
		Coeffs: []*big.Int{big.NewInt(9), big.NewInt(-6), big.NewInt(3), big.NewInt(2)},
		A1:     big.NewInt(1),
		A0:     big.NewInt(5),
	}

	for _, prime := range []int64{2, 3, 5, 7, 11} {
		up := uint64(prime)
		fp := p.CoeffsModP(up)
		for len(fp) <= p.Degree {
			fp = append(fp, 0)
		}
		got := expectedValP(p.Coeffs, fp, up, p.Degree)

		// Empirical average over random coprime (a, b) with large modulus p^k.
		rng := rand.New(rand.NewSource(prime))
		const trials = 40000
		sumVp := 0
		for i := 0; i < trials; i++ {
			var a, b int64
			for {
				a = int64(rng.Int63n(1000000)) - 500000
				b = int64(rng.Int63n(500000)) + 1
				ua := uint64(a)
				if a < 0 {
					ua = uint64(-a)
				}
				if arith.Gcd64(ua, uint64(b)) == 1 {
					break
				}
			}
			val := p.EvalF(a, b)
			val.Abs(val)
			bigP := big.NewInt(prime)
			var rem big.Int
			for val.Sign() > 0 {
				val.QuoRem(val, bigP, &rem)
				if rem.Sign() != 0 {
					break
				}
				sumVp++
			}
		}
		emp := float64(sumVp) / trials
		if math.Abs(got-emp) > 0.03 {
			t.Errorf("p=%d: expectedValP=%.4f, empirical=%.4f", prime, got, emp)
		}
	}

	// Polynomial with no roots mod any p <= 13: x^2 + 1? Let's test x^3 + 7 (mod 2: x^3+1 has root 1;
	// check a polynomial with no roots mod 2, 3, 5: e.g. 2x^2 + 2x + 1? no, c_d=2 has projective root!).
	// Check alpha sign on x^4 + x^3 + x^2 + x + 1 (no roots mod 2, 3; only p ≡ 1 mod 5!).
	cyclo := &Poly{
		Degree: 4,
		Coeffs: []*big.Int{big.NewInt(1), big.NewInt(1), big.NewInt(1), big.NewInt(1), big.NewInt(1)},
		A1:     big.NewInt(1),
		A0:     big.NewInt(1),
	}
	if a := ComputeAlpha(cyclo, 10); a <= 0 {
		t.Errorf("x^4+x^3+x^2+x+1 has no roots mod <= 10, expected alpha > 0, got %v", a)
	}
}

func TestSelectPolyKleinjungAndRootSieve(t *testing.T) {
	rnd := rand.New(rand.NewSource(19))
	for _, halfBits := range []int{95, 132} {
		n := new(big.Int).Mul(testutil.RandomPrime(rnd, halfBits), testutil.RandomPrime(rnd, halfBits))
		p, fac := SelectPoly(n, 4)
		if fac != nil {
			continue
		}
		if p == nil {
			t.Fatalf("SelectPoly(%d bits) returned nil", n.BitLen())
		}
		if !IsIrreducible(p) {
			t.Fatalf("SelectPoly(%d bits) returned reducible f", n.BitLen())
		}
		var g big.Int
		g.GCD(nil, nil, new(big.Int).Abs(p.A1), new(big.Int).Abs(p.A0))
		if g.Cmp(bigOne) != 0 {
			t.Fatalf("gcd(A1, A0) = %v, want 1", &g)
		}
		linMod := new(big.Int).Sub(new(big.Int).Mul(p.A1, p.M), p.A0)
		if linMod.Mod(linMod, n).Sign() != 0 {
			t.Fatalf("A1*M != A0 mod N")
		}
		if fm := evalMonicBig(p.Coeffs, p.M); new(big.Int).Mod(fm, n).Sign() != 0 {
			t.Fatalf("f(M) != 0 mod N")
		}

		base, bfac := BaseMPoly(n, 4, 1)
		if bfac == nil && base != nil {
			pr := defaultParams(n)
			base.Skew = optimizeSkewnessBounded(base, maxSkewForParams(pr))
			gotScore := murphyYieldScore(p, pr)
			baseScore := murphyYieldScore(base, pr)
			if gotScore >= baseScore {
				t.Errorf("%d bits got Murphy score %.2f >= monic base-m %.2f", n.BitLen(), gotScore, baseScore)
			}
		}
	}
}

func TestDickmanRhoAndKleinjungCollision(t *testing.T) {
	ref := [][2]float64{
		{1.0, 1.0},
		{2.0, 0.3068528194},
		{3.0, 0.0486083883},
		{4.0, 0.0049109256},
		{5.0, 0.0003547247},
		{6.0, 1.9649696e-5},
		{7.0, 8.7456699e-7},
		{8.0, 3.2320693e-8},
	}
	for _, tc := range ref {
		got := dickmanRho(tc[0])
		relErr := math.Abs(got-tc[1]) / tc[1]
		if relErr > 2e-5 {
			t.Fatalf("dickmanRho(%.0f) = %.9g, want %.9g (relErr %.2g)", tc[0], got, tc[1], relErr)
		}
	}

	// Verify fractional collision finds the exact global minimum |M - t_{d-1}*a_0^{d-1}|.
	rnd := rand.New(rand.NewSource(53))
	n := new(big.Int).Mul(testutil.RandomPrime(rnd, 166), testutil.RandomPrime(rnd, 166))
	const d = 4
	const cd = 420
	var prs []kjPrimeRoots
	for _, p64 := range arith.PrimesUpTo(400) {
		p := uint64(p64)
		if p < 5 || cd%p == 0 || d%p64 == 0 {
			continue
		}
		if rs := rootsNthModP2(cd, n, d, p); len(rs) == d {
			prs = append(prs, kjPrimeRoots{p: p, p2: p * p, rs: rs})
			if len(prs) == 4 {
				break
			}
		}
	}
	if len(prs) < 4 {
		t.Skip("fewer than 4 splitting primes")
	}

	cdBig := big.NewInt(cd)
	m0 := arith.NthRoot(new(big.Int).Quo(n, cdBig), d)
	a1Check := uint64(1)
	for _, pr := range prs {
		a1Check *= pr.p
	}
	bound := new(big.Int).Mul(new(big.Int).SetUint64(a1Check*a1Check), big.NewInt(cd<<8))
	if bound.Cmp(m0) > 0 {
		t.Skip("a_1 too large for linearized collision bound")
	}
	dBig := big.NewInt(d)
	gotPolys, fac := searchKleinjungTuple(n, d, cd, cdBig, m0, dBig, prs)
	if fac != nil {
		t.Skip("lucky factor")
	}
	if len(gotPolys) == 0 {
		t.Fatal("searchKleinjungTuple returned 0 polynomials")
	}

	// Check all 4^4 = 256 CRT roots directly against gotPolys[0].
	_, a1Sq, allRems := buildHalfCRT(prs)
	a1 := uint64(1)
	for _, pr := range prs {
		a1 *= pr.p
	}
	a1SqBig := new(big.Int).SetUint64(a1Sq)
	m0Mod := new(big.Int).Mod(m0, a1SqBig)
	var minRem *big.Int
	for _, r := range allRems {
		delta := new(big.Int).Sub(new(big.Int).SetUint64(r), m0Mod)
		if delta.Cmp(new(big.Int).Rsh(a1SqBig, 1)) > 0 {
			delta.Sub(delta, a1SqBig)
		} else if delta.Cmp(new(big.Int).Neg(new(big.Int).Rsh(a1SqBig, 1))) <= 0 {
			delta.Add(delta, a1SqBig)
		}
		a0 := new(big.Int).Add(m0, delta)
		if arith.Gcd64(new(big.Int).Mod(a0, new(big.Int).SetUint64(a1)).Uint64(), a1) != 1 {
			continue
		}
		a0d1 := new(big.Int).Exp(a0, big.NewInt(d-1), nil)
		a0d := new(big.Int).Mul(a0d1, a0)
		M := new(big.Int).Quo(new(big.Int).Sub(n, new(big.Int).Mul(cdBig, a0d)), a1SqBig)
		var td1, rem3 big.Int
		roundQuoRem(&td1, &rem3, M, a0d1)
		rem3.Abs(&rem3)
		if minRem == nil || rem3.Cmp(minRem) < 0 {
			minRem = new(big.Int).Set(&rem3)
		}
	}

	// Recompute |M - td1*a0^(d-1)| for gotPolys[0]'s (A1, A0).
	p0 := gotPolys[0]
	a0d1 := new(big.Int).Exp(p0.A0, big.NewInt(d-1), nil)
	a0d := new(big.Int).Mul(a0d1, p0.A0)
	M0 := new(big.Int).Quo(new(big.Int).Sub(n, new(big.Int).Mul(p0.Coeffs[d], a0d)), a1SqBig)
	var td1, gotRem big.Int
	roundQuoRem(&td1, &gotRem, M0, a0d1)
	gotRem.Abs(&gotRem)
	if gotRem.Cmp(minRem) != 0 {
		t.Fatalf("collision got |rem3|=%v, brute-force min=%v", &gotRem, minRem)
	}
}

func TestSemiSmoothCurveAndLatticeProbe(t *testing.T) {
	curve := buildSemiSmoothCurve(250000, 1<<21, 42, 2)
	if got := curve.eval(0.8); got != 1.0 {
		t.Fatalf("curve.eval(0.8) = %v, want 1.0", got)
	}
	prev := curve.eval(1.0)
	for step := 1; step <= 40; step++ {
		u := 1.0 + 0.25*float64(step)
		got := curve.eval(u)
		if got > prev+1e-12 {
			t.Fatalf("curve not monotone at u=%.2f: %g > %g", u, got, prev)
		}
		rho := dickmanRho(u)
		if got+1e-12 < rho {
			t.Fatalf("curve(%.2f) = %g < dickmanRho = %g", u, got, rho)
		}
		if u >= 2.0 && u <= 6.0 && got <= rho*1.2 {
			t.Fatalf("curve(%.2f) = %g lacks 1LP/2LP contribution above rho = %g", u, got, rho)
		}
		prev = got
	}

	rnd := rand.New(rand.NewSource(71))
	n := new(big.Int).Mul(testutil.RandomPrime(rnd, 115), testutil.RandomPrime(rnd, 115))
	pr := defaultParams(n)
	poly, fac := BaseMPoly(n, pr.Degree, 60)
	if fac != nil {
		t.Skip("lucky factor")
	}
	poly.Skew = optimizeSkewnessBounded(poly, maxSkewForParams(pr))
	y1 := probeLatticeYield(poly, pr)
	y2 := probeLatticeYield(poly, pr)
	if y1 <= 0 {
		t.Fatalf("probeLatticeYield = %v, want > 0", y1)
	}
	if y1 != y2 {
		t.Fatalf("probeLatticeYield non-deterministic: %v vs %v", y1, y2)
	}
}

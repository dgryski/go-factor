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
			gotScore := dickmanYieldScore(p, pr)
			baseScore := dickmanYieldScore(base, pr)
			if gotScore >= baseScore {
				t.Errorf("%d bits got Dickman score %.2f >= monic base-m %.2f", n.BitLen(), gotScore, baseScore)
			}
		}
	}
}

package arith

import (
	"fmt"
	"math"
	"math/big"
	"math/rand"
	"testing"
)

func TestIsPrime64(t *testing.T) {
	limit := int64(200000)
	isComposite := make([]bool, limit)
	for _, p := range PrimesUpTo(limit - 1) {
		for m := p * p; m < limit; m += p {
			isComposite[m] = true
		}
	}
	for n := int64(0); n < limit; n++ {
		if want := n >= 2 && !isComposite[n]; IsPrime64(uint64(n)) != want {
			t.Fatalf("IsPrime64(%d) = %v", n, !want)
		}
	}
	for _, n := range []uint64{3215031751, 2152302898747, 3474749660383, 341550071728321, 3825123056546413051} {
		if IsPrime64(n) {
			t.Errorf("IsPrime64(%d) = true", n)
		}
	}
	if !IsPrime64(18446744073709551557) {
		t.Error("IsPrime64(2^64−59) = false")
	}
}

func TestSplitCofactor(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	const pMax, lpBound = 1 << 20, 64 << 20
	randPrime := func(lo, hi uint64) uint64 {
		for {
			if p := lo + uint64(rng.Int63n(int64(hi-lo))); IsPrime64(p) {
				return p
			}
		}
	}
	for i := 0; i < 200; i++ {
		a, b := randPrime(pMax, lpBound), randPrime(pMax, lpBound)
		lp1, lp2, ok := SplitCofactor(a*b, pMax, lpBound)
		if !ok || lp1*lp2 != a*b || lp1 > lp2 {
			t.Fatalf("SplitCofactor(%d·%d) = %d, %d, %v", a, b, lp1, lp2, ok)
		}
	}
	l := randPrime(lpBound-1<<16, lpBound)
	if lp1, lp2, ok := SplitCofactor(l*l, pMax, lpBound); !ok || lp1 != l || lp2 != l {
		t.Errorf("SplitCofactor(%d²) = %d, %d, %v", l, lp1, lp2, ok)
	}
	if _, _, ok := SplitCofactor(randPrime(pMax, 2*pMax)*randPrime(lpBound, 2*lpBound), pMax, lpBound); ok {
		t.Error("accepted a factor above lpBound")
	}
	if _, _, ok := SplitCofactor(lpBound*lpBound+1, pMax, lpBound); ok {
		t.Error("accepted a cofactor above lpBound²")
	}
	if _, _, ok := SplitCofactor(randPrime(pMax, 2*pMax)*randPrime(pMax, 2*pMax)*randPrime(pMax, 2*pMax), pMax, lpBound); ok {
		t.Error("accepted three primes")
	}
	if _, _, ok := SplitCofactor(randPrime(pMax*pMax, 2*pMax*pMax), pMax, lpBound); ok {
		t.Error("accepted a prime cofactor")
	}
}

func TestRho64(t *testing.T) {
	for _, n := range []uint64{0, 1} {
		if d := Rho64(n); d != 0 {
			t.Errorf("Rho64(%d) = %d", n, d)
		}
	}
	if d := Rho64(18); d != 2 {
		t.Errorf("Rho64(18) = %d", d)
	}
	if d := Rho64(1000000007); d != 0 {
		t.Errorf("Rho64(1000000007) = %d", d)
	}
	rng := rand.New(rand.NewSource(2))
	for i := 0; i < 200; i++ {
		var p, q uint64
		for p = uint64(rng.Int63n(1e8)); !IsPrime64(p) || p < 3; p++ {
		}
		for q = uint64(rng.Int63n(1e8)); !IsPrime64(q) || q < 3; q++ {
		}
		n := p * q
		if d := Rho64(n); d == 0 || d == 1 || d == n || n%d != 0 {
			t.Fatalf("Rho64(%d·%d) = %d", p, q, d)
		}
	}
}

func TestSqrtModPrime(t *testing.T) {
	for _, p := range PrimesUpTo(1000) {
		for n := int64(0); n < p; n++ {
			if n != 0 && p > 2 && PowMod64(uint64(n), uint64(p-1)/2, uint64(p)) != 1 {
				continue
			}
			if r := int64(SqrtModPrime(uint32(n), uint32(p))); r*r%p != n {
				t.Fatalf("SqrtModPrime(%d, %d) = %d", n, p, r)
			}
		}
	}
}

func TestLog2Big(t *testing.T) {
	one := big.NewInt(1)
	rng := rand.New(rand.NewSource(3))
	lengths := []int{52, 53, 54, 63, 64, 65, 117, 128, 266}
	for n := 1; n <= 200; n++ {
		lengths = append(lengths, n)
	}
	for _, n := range lengths {
		for i := 0; i < 10; i++ {
			v := new(big.Int).Rand(rng, new(big.Int).Lsh(one, uint(n-1)))
			v.SetBit(v, n-1, 1)
			if i == 0 {
				v.Lsh(one, uint(n-1))
			}
			want := float64(n - 1)
			if n > 1 {
				shift := n - 53
				if shift < 0 {
					shift = 0
				}
				top := new(big.Int).Rsh(v, uint(shift))
				want = math.Log2(float64(top.Uint64())) + float64(shift)
			}
			if got := Log2Big(v); got != want {
				t.Fatalf("Log2Big(%v) = %v, want %v (%d bits)", v, got, want, n)
			}
		}
	}
}

func TestPolyRootsModP(t *testing.T) {
	rng := rand.New(rand.NewSource(42))
	// Brute force check for all primes <= 500 and several primes up to 5000.
	primes := PrimesUpTo(5000)
	for idx, p64 := range primes {
		p := uint64(p64)
		if p > 500 && idx%20 != 0 {
			continue
		}
		for trial := 0; trial < 3; trial++ {
			deg := 3 + trial%3 // 3, 4, 5
			f := make([]uint64, deg+1)
			for i := range f {
				f[i] = uint64(rng.Int63n(p64))
			}
			if f[deg] == 0 {
				f[deg] = 1
			}
			got := PolyRootsModP(f, p)
			var want []uint64
			for x := uint64(0); x < p; x++ {
				if PolyEvalModP(f, x, p) == 0 {
					want = append(want, x)
				}
			}
			if len(got) != len(want) {
				t.Fatalf("p=%d f=%v: got %v, want %v", p, f, got, want)
			}
			for i := range got {
				if got[i] != want[i] {
					t.Fatalf("p=%d f=%v: got %v, want %v", p, f, got, want)
				}
			}
		}
	}

	// 64-bit prime test.
	bigP := uint64(18446744073709551557) // 2^64 - 59
	// Build (x - r1)(x - r2)(x^2 + 1) mod bigP (note bigP ≡ 1 mod 4 so x^2+1 splits too, let's check!)
	r1 := uint64(12345678901234567)
	r2 := uint64(98765432109876543)
	f := PolyMulModP([]uint64{bigP - r1, 1}, []uint64{bigP - r2, 1}, bigP)
	f = PolyMulModP(f, []uint64{3, 1, 1}, bigP) // x^2 + x + 3
	roots := PolyRootsModP(f, bigP)
	for _, r := range []uint64{r1, r2} {
		found := false
		for _, g := range roots {
			if g == r {
				found = true
			}
			if PolyEvalModP(f, g, bigP) != 0 {
				t.Fatalf("root %d does not evaluate to 0 mod %d", g, bigP)
			}
		}
		if !found {
			t.Fatalf("missing root %d in %v", r, roots)
		}
	}
}

func TestPolyIsIrreducibleAndFactorModP(t *testing.T) {
	rng := rand.New(rand.NewSource(7))
	for _, p := range []uint64{2, 3, 5, 17, 1000003, 18446744073709551557} {
		// Check x^2 + 1 mod 3, mod 2, etc.
		if p == 3 && !PolyIsIrreducibleModP([]uint64{1, 0, 1}, p) {
			t.Fatal("x^2+1 should be irreducible mod 3")
		}
		if p == 5 && PolyIsIrreducibleModP([]uint64{1, 0, 1}, p) {
			t.Fatal("x^2+1 should be reducible mod 5")
		}
		// Factor random squarefree degree 4 and 5 polynomials and verify product.
		for _, deg := range []int{3, 4, 5} {
			for attempt := 0; attempt < 5; attempt++ {
				f := make([]uint64, deg+1)
				for i := 0; i < deg; i++ {
					f[i] = rng.Uint64() % p
				}
				f[deg] = 1
				df := PolyDerivModP(f, p)
				if len(PolyGcdModP(f, df, p)) > 1 {
					continue // not squarefree
				}
				irred := PolyIsIrreducibleModP(f, p)
				facs := PolyFactorModP(f, p)
				if irred != (len(facs) == 1) {
					t.Fatalf("p=%d f=%v: irred=%v but got %d factors: %v", p, f, irred, len(facs), facs)
				}
				prod := []uint64{1}
				for _, h := range facs {
					if !PolyIsIrreducibleModP(h, p) {
						t.Fatalf("factor %v mod %d is not irreducible", h, p)
					}
					prod = PolyMulModP(prod, h, p)
				}
				if !PolyEqual(prod, f) {
					t.Fatalf("product %v != %v mod %d", prod, f, p)
				}
			}
		}
	}
}

func TestNthRoot(t *testing.T) {
	rng := rand.New(rand.NewSource(9))
	for d := 2; d <= 6; d++ {
		for i := 0; i < 50; i++ {
			bits := 20 + i*5
			n := new(big.Int).Rand(rng, new(big.Int).Lsh(big.NewInt(1), uint(bits)))
			r := NthRoot(n, d)
			rPow := new(big.Int).Exp(r, big.NewInt(int64(d)), nil)
			r1Pow := new(big.Int).Exp(new(big.Int).Add(r, big.NewInt(1)), big.NewInt(int64(d)), nil)
			if rPow.Cmp(n) > 0 || r1Pow.Cmp(n) <= 0 {
				t.Fatalf("NthRoot(%v, %d) = %v wrong", n, d, r)
			}
		}
	}
}

func BenchmarkBigMul(b *testing.B) {
	rng := rand.New(rand.NewSource(1))
	one := big.NewInt(1)
	for _, mbits := range []int{1, 2, 5, 10} {
		nbits := mbits * 1_000_000
		lim := new(big.Int).Lsh(one, uint(nbits))
		x := new(big.Int).Rand(rng, lim)
		y := new(big.Int).Rand(rng, lim)
		var z big.Int
		b.Run(fmt.Sprintf("Mbit=%d", mbits), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				z.Mul(x, y)
			}
		})
	}
}

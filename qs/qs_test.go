package qs

import (
	"fmt"
	"math"
	"math/big"
	"math/rand"
	"testing"

	"github.com/dgryski/go-factor/internal/testutil"
)

func TestQuadraticSieve(t *testing.T) {
	n := big.NewInt(1000036000099) // 1000003 * 1000033
	orig := new(big.Int).Set(n)
	p, c := QuadraticSieve(n)
	testutil.VerifyFactoring(t, orig, p, c)
}

func TestQuadraticSieveEdgeCases(t *testing.T) {
	for _, tt := range []struct {
		n      int64
		primes int
	}{
		{4, 2},                   // smallest perfect square
		{10201, 2},               // 101^2: would give Y = 0
		{1 << 40, 40},            // even, repeated squares
		{15, 2},                  // factors found by trial division
		{3 * 5 * 7 * 1000003, 4}, // small factors plus a large one
		{1000003, 1},             // prime
	} {
		n := big.NewInt(tt.n)
		p, c := QuadraticSieve(n)
		testutil.VerifyFactoring(t, big.NewInt(tt.n), p, c)
		if len(p) != tt.primes || len(c) != 0 {
			t.Errorf("QuadraticSieve(%d) = %v, %v; want %d primes", tt.n, p, c, tt.primes)
		}
	}
}

func testSemiprimes(t *testing.T, minBits, maxBits, step, count int) {
	rnd := rand.New(rand.NewSource(1))
	for bits := minBits; bits <= maxBits; bits += step {
		for i := 0; i < count; i++ {
			a, b := testutil.RandomPrime(rnd, bits), testutil.RandomPrime(rnd, bits)
			n := new(big.Int).Mul(a, b)
			p, c := QuadraticSieve(new(big.Int).Set(n))
			testutil.VerifyFactoring(t, n, p, c)
			if len(c) != 0 {
				t.Errorf("QuadraticSieve(%v = %v * %v) left composites %v", n, a, b, c)
			}
		}
	}
}

func TestQuadraticSieveRandomSemiprimes(t *testing.T) {
	testSemiprimes(t, 10, 38, 4, 20)
}

func TestQuadraticSieveSIQS(t *testing.T) {
	if testing.Short() {
		t.Skip("slow")
	}
	testSemiprimes(t, siqsMinBits/2, 75, 5, 3)
}

func TestLog2Big(t *testing.T) {
	rng := rand.New(rand.NewSource(3))
	lengths := []int{52, 53, 54, 63, 64, 65, 117, 128, 266}
	for n := 1; n <= 512; n++ {
		lengths = append(lengths, n)
	}
	for _, n := range lengths {
		for i := 0; i < 20; i++ {
			v := new(big.Int).Rand(rng, new(big.Int).Lsh(one, uint(n-1)))
			v.SetBit(v, n-1, 1)
			if i == 0 {
				v.Lsh(one, uint(n-1)) // a power of two
			}
			// Exact for the top 53 bits, then the logarithm of the rest.
			want := float64(n - 1)
			if n > 1 {
				shift := n - 53
				if shift < 0 {
					shift = 0
				}
				top := new(big.Int).Rsh(v, uint(shift))
				want = math.Log2(float64(top.Uint64())) + float64(shift)
			}
			if got := log2Big(v); got != want {
				t.Fatalf("log2Big(%v) = %v, want %v (%d bits)", v, got, want, n)
			}
		}
	}
}

func TestSqrtModPrime(t *testing.T) {
	for _, p := range primesUpTo(3000) {
		for n := int64(0); n < p; n++ {
			if n != 0 && p > 2 && powMod64(uint64(n), uint64(p-1)/2, uint64(p)) != 1 {
				continue
			}
			if r := int64(sqrtModPrime(uint32(n), uint32(p))); r*r%p != n {
				t.Fatalf("sqrtModPrime(%d, %d) = %d", n, p, r)
			}
		}
	}
}

// TestRelations checks that every relation gathered by both sieves is
// correct: X² ≡ Y (mod n), with Y rebuilt from its factors. The factors of
// these semiprimes are all far above the factor base.
func TestRelations(t *testing.T) {
	rnd := rand.New(rand.NewSource(2))
	for _, bits := range []int{30, 45, 60, 80} {
		n := new(big.Int).Mul(testutil.RandomPrime(rnd, bits), testutil.RandomPrime(rnd, bits))
		b, siqs, sieveRange := params(n)
		kn := new(big.Int).Mul(n, big.NewInt(knuthSchroeppel(n)))
		factorBase := buildFactorBase(kn, primesUpTo(b))
		sv := newSieve(n, kn, factorBase, siqs)
		if siqs {
			sv.siqs(sieveRange)
		} else {
			sv.singlePoly(sieveRange)
		}
		if !sv.done() {
			t.Errorf("%v: found %d relations and %d partials, want %d", n, len(sv.relations), len(sv.partials), sv.want)
		}
		primes := sv.finish(factorBase)
		for _, rel := range sv.relations {
			y := big.NewInt(1)
			for _, c := range rel.Factors {
				if c == 0 {
					y.Neg(y)
				} else {
					y.Mul(y, big.NewInt(primes[c-1]))
				}
			}
			lhs := new(big.Int).Mul(rel.X, rel.X)
			lhs.Sub(lhs, y)
			if lhs.Mod(lhs, n).Sign() != 0 {
				t.Fatalf("%v: X=%v does not square to its Y=%v", n, rel.X, y)
			}
		}
	}
}

func BenchmarkQuadraticSieve(b *testing.B) {
	for _, bits := range []int{30, 40, 50, 60, 70, 80} {
		rnd := rand.New(rand.NewSource(1))
		n := new(big.Int).Mul(testutil.RandomPrime(rnd, bits), testutil.RandomPrime(rnd, bits))
		b.Run(fmt.Sprintf("bits=%d", 2*bits), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if p, _ := split(n); p == nil {
					b.Fatalf("split(%v) failed", n)
				}
			}
		})
	}
}

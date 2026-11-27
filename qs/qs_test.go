package qs

import (
	"fmt"
	"math/big"
	"math/rand"
	"testing"

	"github.com/dgryski/go-factor/internal/arith"
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

// TestRelations checks that every relation gathered by both sieves is
// correct: X² ≡ Y (mod n), with Y rebuilt from its factors. The factors of
// these semiprimes are all far above the factor base.
func TestRelations(t *testing.T) {
	rnd := rand.New(rand.NewSource(2))
	for _, bits := range []int{30, 45, 60, 80} {
		n := new(big.Int).Mul(testutil.RandomPrime(rnd, bits), testutil.RandomPrime(rnd, bits))
		b, siqs, sieveRange := params(n)
		kn := new(big.Int).Mul(n, big.NewInt(knuthSchroeppel(n)))
		factorBase := buildFactorBase(kn, arith.PrimesUpTo(b))
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

// TestWorkers checks that SIQS collects the same relations and partials, in
// the same order, with one worker and with several, and with double large
// primes on and off.
func TestWorkers(t *testing.T) {
	defer func(w, d int) { Workers, dlpMinBits = w, d }(Workers, dlpMinBits)
	rnd := rand.New(rand.NewSource(5))
	for _, c := range []struct{ bits, dlpMinBits int }{{60, 1 << 30}, {80, 1 << 30}, {80, 0}} {
		dlpMinBits = c.dlpMinBits
		n := new(big.Int).Mul(testutil.RandomPrime(rnd, c.bits), testutil.RandomPrime(rnd, c.bits))
		b, _, sieveRange := params(n)
		kn := new(big.Int).Mul(n, big.NewInt(knuthSchroeppel(n)))
		factorBase := buildFactorBase(kn, arith.PrimesUpTo(b))
		var runs [2]*sieve
		for i, w := range []int{1, 4} {
			Workers = w
			runs[i] = newSieve(n, kn, factorBase, true)
			if f := runs[i].siqs(sieveRange); f != nil {
				t.Fatalf("%v: found factor %v while sieving", n, f)
			}
		}
		one, four := runs[0], runs[1]
		if four.workers() != 4 {
			t.Fatalf("%d bits: %d workers, want 4", n.BitLen(), four.workers())
		}
		if one.polys != four.polys || len(one.relations) != len(four.relations) || len(one.partials) != len(four.partials) {
			t.Fatalf("%d bits: %d polynomials, %d relations, %d partials with one worker; %d, %d, %d with four",
				n.BitLen(), one.polys, len(one.relations), len(one.partials), four.polys, len(four.relations), len(four.partials))
		}
		same := func(a, b relation) bool {
			return a.X.Cmp(b.X) == 0 && fmt.Sprint(a.Factors) == fmt.Sprint(b.Factors)
		}
		for i := range one.relations {
			if !same(one.relations[i], four.relations[i]) {
				t.Fatalf("%d bits: relation %d differs", n.BitLen(), i)
			}
		}
		for i := range one.partials {
			p, q := one.partials[i], four.partials[i]
			if !same(p.relation, q.relation) || p.lp1 != q.lp1 || p.lp2 != q.lp2 {
				t.Fatalf("%d bits: partial %d differs", n.BitLen(), i)
			}
		}
	}
}

func TestAppendGe(t *testing.T) {
	rng := rand.New(rand.NewSource(17))
	buf := make([]byte, 261)
	for i := range buf {
		buf[i] = byte(rng.Intn(256))
	}
	buf[0], buf[7], buf[8], buf[260] = 0, 128, 127, 255
	var got []int
	for th := -1; th <= 257; th += 3 {
		got = appendGe(got[:0], buf, th)
		var want []int
		for i, b := range buf {
			if int(b) >= th {
				want = append(want, i)
			}
		}
		if fmt.Sprint(got) != fmt.Sprint(want) {
			t.Fatalf("appendGe(th=%d) = %v, want %v", th, got, want)
		}
	}
}

func TestQ2Relations(t *testing.T) {
	rnd := rand.New(rand.NewSource(23))
	foundQ2 := false
	for trial := 0; trial < 20; trial++ {
		n := new(big.Int).Mul(testutil.RandomPrime(rnd, 55), testutil.RandomPrime(rnd, 55))
		b, siqs, sieveRange := params(n)
		kn := new(big.Int).Mul(n, big.NewInt(knuthSchroeppel(n)))
		factorBase := buildFactorBase(kn, arith.PrimesUpTo(b))
		sv := newSieve(n, kn, factorBase, siqs)
		if !sv.q2 {
			continue
		}
		foundQ2 = true
		sv.siqs(sieveRange)
		if !sv.done() {
			t.Fatalf("%v: Q2 sieve did not reach target", n)
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
			if lhs.Sub(lhs, y).Mod(lhs, n).Sign() != 0 {
				t.Fatalf("Q2 relation X=%v does not square to Y=%v mod N", rel.X, y)
			}
		}
		break
	}
	if !foundQ2 {
		t.Fatal("no kN ≡ 1 (mod 8) input encountered")
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

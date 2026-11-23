package factor

import (
	"math/big"
	"math/rand"
	"testing"
)

func TestQuadraticSieve(t *testing.T) {
	n := big.NewInt(1000036000099) // 1000003 * 1000033
	orig := (&big.Int{}).Set(n)
	p, c := QuadraticSieve(n)
	t.Log(p, c)
	verifyFactoring(t, orig, p, c)
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
		verifyFactoring(t, big.NewInt(tt.n), p, c)
		if len(p) != tt.primes || len(c) != 0 {
			t.Errorf("QuadraticSieve(%d) = %v, %v; want %d primes", tt.n, p, c, tt.primes)
		}
	}
}

func randomPrime(rnd *rand.Rand, bits int) *big.Int {
	min := new(big.Int).Lsh(one, uint(bits-1))
	for {
		p := new(big.Int).Rand(rnd, min)
		p.Add(p, min)
		if p.ProbablyPrime(20) {
			return p
		}
	}
}

func TestQuadraticSieveRandomSemiprimes(t *testing.T) {
	rnd := rand.New(rand.NewSource(1))
	for bits := 10; bits <= 26; bits += 4 {
		for i := 0; i < 20; i++ {
			a, b := randomPrime(rnd, bits), randomPrime(rnd, bits)
			n := new(big.Int).Mul(a, b)
			p, c := QuadraticSieve(new(big.Int).Set(n))
			verifyFactoring(t, n, p, c)
			if len(c) != 0 {
				t.Errorf("QuadraticSieve(%v = %v * %v) left composites %v", n, a, b, c)
			}
		}
	}
}

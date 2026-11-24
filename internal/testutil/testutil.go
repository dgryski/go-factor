// Package testutil has helpers shared by the factoring tests.
package testutil

import (
	"math/big"
	"math/rand"
	"testing"
	"time"
)

// VerifyFactoring checks that p holds at least one prime and only primes, c
// only composites, and that together they multiply back to n.
func VerifyFactoring(t *testing.T, n *big.Int, p []*big.Int, c []*big.Int) {
	t.Helper()

	if len(p) == 0 {
		t.Errorf("no prime factors found for %v", n)
	}

	for _, pp := range p {
		if !pp.ProbablyPrime(20) {
			t.Errorf("pp=%v is not prime", pp)
		}
	}

	for _, cc := range c {
		if cc.ProbablyPrime(20) {
			t.Errorf("cc=%v is probably prime", cc)
		}
	}

	nn := big.NewInt(1)

	for _, pp := range p {
		nn = nn.Mul(nn, pp)
	}

	for _, cc := range c {
		nn = nn.Mul(nn, cc)
	}

	if n.Cmp(nn) != 0 {
		t.Errorf("Product(p, c)=%v, want %v", nn, n)
	}
}

// VerifyProduct checks that the factors multiply back to n.
func VerifyProduct(t *testing.T, n *big.Int, p []*big.Int, c []*big.Int) {
	t.Helper()

	nn := big.NewInt(1)
	for _, pp := range p {
		if !pp.ProbablyPrime(20) {
			t.Errorf("pp=%v is not prime", pp)
		}
		nn.Mul(nn, pp)
	}
	for _, cc := range c {
		nn.Mul(nn, cc)
	}
	if n.Cmp(nn) != 0 {
		t.Errorf("Product(p, c)=%v, want %v", nn, n)
	}
}

// WithTimeout fails the test if f does not return within d.
func WithTimeout(t *testing.T, d time.Duration, name string, f func()) {
	t.Helper()
	done := make(chan struct{})
	go func() {
		defer close(done)
		f()
	}()
	select {
	case <-done:
	case <-time.After(d):
		t.Fatalf("%s did not return within %v", name, d)
	}
}

// RandomPrime returns a random prime with exactly bits bits.
func RandomPrime(rnd *rand.Rand, bits int) *big.Int {
	min := new(big.Int).Lsh(big.NewInt(1), uint(bits-1))
	for {
		p := new(big.Int).Rand(rnd, min)
		p.Add(p, min)
		if p.ProbablyPrime(20) {
			return p
		}
	}
}

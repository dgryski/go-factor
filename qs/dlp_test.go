package qs

import (
	"math/big"
	"math/rand"
	"testing"

	"github.com/dgryski/go-factor/internal/arith"
	"github.com/dgryski/go-factor/internal/testutil"
)

// TestDLPRelations checks every relation from a sieve with double large
// primes on, large prime columns included, and that the run finds a factor.
func TestDLPRelations(t *testing.T) {
	defer func(old int) { dlpMinBits = old }(dlpMinBits)
	dlpMinBits = 0
	rnd := rand.New(rand.NewSource(4))
	n := new(big.Int).Mul(testutil.RandomPrime(rnd, 80), testutil.RandomPrime(rnd, 80))
	b, _, sieveRange := params(n)
	kn := new(big.Int).Mul(n, big.NewInt(knuthSchroeppel(n)))
	factorBase := buildFactorBase(kn, arith.PrimesUpTo(b))
	s := newSieve(n, kn, factorBase, true)
	if s.dlpCap == 0 {
		t.Fatal("double large primes are off")
	}
	s.siqs(sieveRange)
	two := 0
	for _, p := range s.partials {
		if p.lp2 != 0 {
			two++
		}
	}
	if two == 0 {
		t.Fatal("no partials with two large primes")
	}
	primes := s.finish(factorBase)
	for _, rel := range s.relations {
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
			t.Fatalf("X=%v does not square to its Y=%v", rel.X, y)
		}
	}
	for _, dep := range findDependencies(s.relations, len(primes)+1) {
		if tryDependency(n, primes, s.relations, dep) != nil {
			t.Logf("%d partials with two large primes, %d large prime columns", two, len(primes)-len(factorBase))
			return
		}
	}
	t.Fatal("no factor found")
}

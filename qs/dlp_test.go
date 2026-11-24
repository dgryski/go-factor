package qs

import (
	"math/big"
	"math/rand"
	"testing"

	"github.com/dgryski/go-factor/internal/testutil"
)

func TestIsPrime64(t *testing.T) {
	limit := int64(200000)
	isComposite := make([]bool, limit)
	for _, p := range primesUpTo(limit - 1) {
		for m := p * p; m < limit; m += p {
			isComposite[m] = true
		}
	}
	for n := int64(0); n < limit; n++ {
		if want := n >= 2 && !isComposite[n]; isPrime64(uint64(n)) != want {
			t.Fatalf("isPrime64(%d) = %v", n, !want)
		}
	}
	// Strong pseudoprimes to several small bases.
	for _, n := range []uint64{3215031751, 2152302898747, 3474749660383, 341550071728321, 3825123056546413051} {
		if isPrime64(n) {
			t.Errorf("isPrime64(%d) = true", n)
		}
	}
	if !isPrime64(18446744073709551557) {
		t.Error("isPrime64(2^64−59) = false")
	}
}

func TestSplitCofactor(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	const pMax, lpBound = 1 << 20, 64 << 20
	randPrime := func(lo, hi uint64) uint64 {
		for {
			if p := lo + uint64(rng.Int63n(int64(hi-lo))); isPrime64(p) {
				return p
			}
		}
	}
	for i := 0; i < 200; i++ {
		a, b := randPrime(pMax, lpBound), randPrime(pMax, lpBound)
		lp1, lp2, ok := splitCofactor(a*b, pMax, lpBound)
		if !ok || lp1*lp2 != a*b || lp1 > lp2 {
			t.Fatalf("splitCofactor(%d·%d) = %d, %d, %v", a, b, lp1, lp2, ok)
		}
	}
	// The square of the bound, and the largest square below it.
	l := randPrime(lpBound-1<<16, lpBound)
	if lp1, lp2, ok := splitCofactor(l*l, pMax, lpBound); !ok || lp1 != l || lp2 != l {
		t.Errorf("splitCofactor(%d²) = %d, %d, %v", l, lp1, lp2, ok)
	}
	// One factor too large, a cofactor above lpBound², three primes, and a
	// prime cofactor.
	if _, _, ok := splitCofactor(randPrime(pMax, 2*pMax)*randPrime(lpBound, 2*lpBound), pMax, lpBound); ok {
		t.Error("accepted a factor above lpBound")
	}
	if _, _, ok := splitCofactor(lpBound*lpBound+1, pMax, lpBound); ok {
		t.Error("accepted a cofactor above lpBound²")
	}
	if _, _, ok := splitCofactor(randPrime(pMax, 2*pMax)*randPrime(pMax, 2*pMax)*randPrime(pMax, 2*pMax), pMax, lpBound); ok {
		t.Error("accepted three primes")
	}
	if _, _, ok := splitCofactor(randPrime(pMax*pMax, 2*pMax*pMax), pMax, lpBound); ok {
		t.Error("accepted a prime cofactor")
	}
}

func TestRho64(t *testing.T) {
	for _, n := range []uint64{0, 1} {
		if d := rho64(n); d != 0 {
			t.Errorf("rho64(%d) = %d", n, d)
		}
	}
	if d := rho64(18); d != 2 {
		t.Errorf("rho64(18) = %d", d)
	}
	// Small odd composites may defeat rho, but any answer must divide n.
	for _, n := range []uint64{9, 15, 21, 25, 49, 91} {
		if d := rho64(n); d != 0 && (d == 1 || d == n || n%d != 0) {
			t.Errorf("rho64(%d) = %d", n, d)
		}
	}
	// A prime gives up rather than spinning.
	if d := rho64(1000000007); d != 0 {
		t.Errorf("rho64(1000000007) = %d", d)
	}
	rng := rand.New(rand.NewSource(2))
	for i := 0; i < 200; i++ {
		var p, q uint64
		for p = uint64(rng.Int63n(1e8)); !isPrime64(p) || p < 3; p++ {
		}
		for q = uint64(rng.Int63n(1e8)); !isPrime64(q) || q < 3; q++ {
		}
		n := p * q
		if d := rho64(n); d == 0 || d == 1 || d == n || n%d != 0 {
			t.Fatalf("rho64(%d·%d) = %d", p, q, d)
		}
	}
}

// TestDLPRelations checks every relation from a sieve with double large
// primes on, large prime columns included, and that the run finds a factor.
func TestDLPRelations(t *testing.T) {
	defer func(old int) { dlpMinBits = old }(dlpMinBits)
	dlpMinBits = 0
	rnd := rand.New(rand.NewSource(4))
	n := new(big.Int).Mul(testutil.RandomPrime(rnd, 80), testutil.RandomPrime(rnd, 80))
	b, _, sieveRange := params(n)
	kn := new(big.Int).Mul(n, big.NewInt(knuthSchroeppel(n)))
	factorBase := buildFactorBase(kn, primesUpTo(b))
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

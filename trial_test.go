package factor

import (
	"math/big"
	"testing"
)

func TestTrial(t *testing.T) {
	n := big.NewInt(275)
	orig := (&big.Int{}).Set(n)
	p, c := Trial(n)
	t.Log(p, c)
	verifyFactoring(t, orig, p, c)
}

func verifyFactoring(t *testing.T, n *big.Int, p []*big.Int, c []*big.Int) {
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

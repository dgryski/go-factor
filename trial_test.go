package factor

import (
	"math/big"
	"testing"
	"time"
)

func TestTrial(t *testing.T) {
	n := big.NewInt(275)
	orig := (&big.Int{}).Set(n)
	p, c := Trial(n)
	t.Log(p, c)
	verifyFactoring(t, orig, p, c)
}

func TestInputUnchanged(t *testing.T) {
	for _, tt := range []struct {
		name string
		n    string
		f    func(*big.Int) ([]*big.Int, []*big.Int)
	}{
		{"Trial", "275", Trial},
		{"PMinus1", "3928471", PMinus1},
		{"Brent", "3928471", func(n *big.Int) ([]*big.Int, []*big.Int) { return Brent(n, 19, 1) }},
		{"Squfof", "29384129384129384712938471", func(n *big.Int) ([]*big.Int, []*big.Int) { return Squfof(n, 5) }},
	} {
		n, _ := new(big.Int).SetString(tt.n, 10)
		tt.f(n)
		if n.String() != tt.n {
			t.Errorf("%s changed its input %v to %v", tt.name, tt.n, n)
		}
	}
}

func TestTrialNonPositive(t *testing.T) {
	for _, v := range []int64{0, -1, -12} {
		n := big.NewInt(v)
		var p, c []*big.Int
		withTimeout(t, 5*time.Second, "Trial("+n.String()+")", func() {
			p, c = Trial(n)
		})
		verifyProduct(t, big.NewInt(v), p, c)
	}
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

// verifyProduct checks that the factors multiply back to n.
func verifyProduct(t *testing.T, n *big.Int, p []*big.Int, c []*big.Int) {
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

// withTimeout fails the test if f does not return within d.
func withTimeout(t *testing.T, d time.Duration, name string, f func()) {
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

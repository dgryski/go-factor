package factor

import (
	"math/big"
	"testing"

	"github.com/dgryski/go-factor/internal/testutil"
)

func nextPrime(s string) *big.Int {
	n, _ := new(big.Int).SetString(s, 10)
	for !n.ProbablyPrime(20) {
		n.Add(n, one)
	}
	return n
}

func product(ps ...*big.Int) *big.Int {
	n := big.NewInt(1)
	for _, p := range ps {
		n.Mul(n, p)
	}
	return n
}

func TestECM(t *testing.T) {
	tests := []struct {
		name    string
		factors []*big.Int
	}{
		{"13+27 digits", []*big.Int{nextPrime("1000000000039"), nextPrime("314159265358979323846264338")}},
		{"three 12 digit", []*big.Int{nextPrime("100000000000"), nextPrime("271828182845"), nextPrime("577215664901")}},
		{"repeated", []*big.Int{nextPrime("1000003"), nextPrime("1000003"), nextPrime("99999999977")}},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			n := product(tt.factors...)
			orig := new(big.Int).Set(n)
			p, c := ECM(n, 2000, 50)
			testutil.VerifyFactoring(t, orig, p, c)
			if len(c) != 0 {
				t.Errorf("composites left: %v", c)
			}
			if len(p) != len(tt.factors) {
				t.Errorf("got %d primes %v, want %d", len(p), p, len(tt.factors))
			}
		})
	}
}

func TestECMPrime(t *testing.T) {
	n := nextPrime("1000000000039")
	p, c := ECM(n, 2000, 10)
	if len(p) != 1 || p[0].Cmp(n) != 0 || len(c) != 0 {
		t.Errorf("ECM(prime) = %v, %v", p, c)
	}
	p, c = ECM(big.NewInt(1), 2000, 10)
	if len(p) != 0 || len(c) != 0 {
		t.Errorf("ECM(1) = %v, %v", p, c)
	}
}

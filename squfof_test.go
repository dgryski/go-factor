package factor

import (
	"math/big"
	"testing"

	"github.com/dgryski/go-factor/internal/testutil"
)

func TestSqufof(t *testing.T) {
	n, _ := big.NewInt(0).SetString("29384129384129384712938471", 10)
	orig := (&big.Int{}).Set(n)
	p, c := Squfof(n, 5)
	t.Log(p, c)
	testutil.VerifyFactoring(t, orig, p, c)
}

func TestSqufofSquare(t *testing.T) {
	p := big.NewInt(1000003)
	for _, tt := range []struct {
		n *big.Int
		k int64
	}{
		{new(big.Int).Mul(p, p), 1},
		{new(big.Int).Mul(p, p), 9},
		{new(big.Int).Mul(big.NewInt(3), new(big.Int).Mul(p, p)), 3},
	} {
		orig := new(big.Int).Set(tt.n)
		pr, co := Squfof(tt.n, tt.k)
		testutil.VerifyFactoring(t, orig, pr, co)
	}
}

func TestSqufofMultiplier(t *testing.T) {
	// These fail with k=3 if the queue ignores the multiplier.
	for _, s := range []string{"3017139931213241", "3027440398119283", "3034316960173133"} {
		n, _ := big.NewInt(0).SetString(s, 10)
		orig := new(big.Int).Set(n)
		p, c := Squfof(n, 3)
		testutil.VerifyFactoring(t, orig, p, c)
	}
}

func TestIsSquare(t *testing.T) {
	r, _ := big.NewInt(0).SetString("123456789012345678901", 10)
	sq := new(big.Int).Mul(r, r)
	if got, ok := isSquare(sq); !ok || got.Cmp(r) != 0 {
		t.Errorf("isSquare(%v) = %v, %v", sq, got, ok)
	}
	sq.Add(sq, one)
	if _, ok := isSquare(sq); ok {
		t.Errorf("isSquare(%v) = true", sq)
	}
	if _, ok := isSquare(big.NewInt(0)); !ok {
		t.Errorf("isSquare(0) = false")
	}
}

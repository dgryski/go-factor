package factor

import (
	"math/big"
	"testing"
)

func TestSqufof(t *testing.T) {
	n, _ := big.NewInt(0).SetString("29384129384129384712938471", 10)
	orig := (&big.Int{}).Set(n)
	p, c := Squfof(n, 5)
	t.Log(p, c)
	verifyFactoring(t, orig, p, c)
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
		verifyFactoring(t, orig, pr, co)
	}
}

package factor

import (
	"math/big"
	"testing"
	"time"
)

func TestBrent(t *testing.T) {
	n := big.NewInt(3928471)
	orig := (&big.Int{}).Set(n)
	p, c := Brent(n, 19, 1)
	t.Log(p, c)
	verifyFactoring(t, orig, p, c)
}

func TestBrentCycle(t *testing.T) {
	// x^2+1 starting at 2 cycles mod n before a factor is found.
	for _, v := range []int64{15, 21, 35, 55, 95, 133, 589, 667} {
		n := big.NewInt(v)
		orig := new(big.Int).Set(n)
		var p, c []*big.Int
		withTimeout(t, 5*time.Second, "Brent("+orig.String()+")", func() {
			p, c = Brent(n, 2, 1)
		})
		verifyProduct(t, orig, p, c)
	}
}

func TestBrentPartialBatch(t *testing.T) {
	// These are found in a batch of fewer than 16 terms.
	for _, v := range []int64{21, 35, 55, 91, 323} {
		n := big.NewInt(v)
		p, c := Brent(n, 2, 1)
		verifyFactoring(t, big.NewInt(v), p, c)
	}
}

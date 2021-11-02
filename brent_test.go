package factor

import (
	"math/big"
	"testing"
)

func TestBrent(t *testing.T) {
	n := big.NewInt(3928471)
	orig := (&big.Int{}).Set(n)
	p, c := Brent(n, 19, 1)
	t.Log(p, c)
	verifyFactoring(t, orig, p, c)
}

package factor

import (
	"math/big"
	"testing"
)

func TestBrent(t *testing.T) {
	p, c := Brent(big.NewInt(3928471), 19, 1)
	t.Log(p, c)
}

package factor

import (
	"math/big"
	"testing"
)

func TestSqufof(t *testing.T) {
	n, _ := big.NewInt(0).SetString("29384129384129384712938471", 10)
	p, c := Squfof(n, 5)
	t.Log(p, c)
}

package factor

import (
	"math/big"
	"testing"
)

func TestSqufof(t *testing.T) {
	p, c := Squfof(big.NewInt(1467937633499), 3)
	t.Log(p, c)
}

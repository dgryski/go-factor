package factor

import (
	"math/big"
	"testing"
)

func TestTrial(t *testing.T) {
	p, c := Trial(big.NewInt(725))
	t.Log(p, c)
}

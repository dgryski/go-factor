package factor

import (
	"math/big"
	"testing"
)

func TestPMinus1(t *testing.T) {
	// 2500069784293916113216294633942411
	// 923874293742934712938741239847123984712394871
	p, c := PMinus1(big.NewInt(3928471))
	t.Log(p, c)
}

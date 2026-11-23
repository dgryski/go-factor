package factor

import (
	"math/big"
	"testing"
)

func TestPMinus1(t *testing.T) {
	// 2500069784293916113216294633942411
	// 923874293742934712938741239847123984712394871
	n := big.NewInt(3928471)
	orig := (&big.Int{}).Set(n)
	p, c := PMinus1(n)
	t.Log(p, c)
	verifyFactoring(t, orig, p, c)
}

func TestPMinus1SmallFactor(t *testing.T) {
	for _, v := range []int64{7 * 1000003, 11 * 1000003, 2 * 3928471, 4 * 3928471, 2, 4} {
		n := big.NewInt(v)
		orig := new(big.Int).Set(n)
		p, c := PMinus1(n)
		verifyFactoring(t, orig, p, c)
		if len(c) != 0 {
			t.Errorf("PMinus1(%v) left composites %v", orig, c)
		}
	}
}

func TestPMinus1LastBatch(t *testing.T) {
	// p-1 has the largest small prime as a factor.
	last := smallPrimes[len(smallPrimes)-1]
	p := new(big.Int)
	for k := int64(1); ; k++ {
		p.SetInt64(2*k*last + 1)
		if p.ProbablyPrime(20) {
			break
		}
	}
	n := new(big.Int).Mul(p, big.NewInt(1000003))
	orig := new(big.Int).Set(n)
	pr, c := PMinus1(n)
	verifyFactoring(t, orig, pr, c)
}

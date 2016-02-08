package factor

import (
	"math/big"
)

func PMinus1(n *big.Int) (primes, composites []*big.Int) {

	if n.ProbablyPrime(10) {
		return []*big.Int{n}, nil
	}

	a := big.NewInt(2 * 3 * 5 * 7 * 11)

	var one big.Int
	one.SetInt64(1)

	var bigp big.Int
	var m big.Int
	for i, p := range smallPrimes {

		bigp.SetInt64(p)

		m.SetInt64(1)
		for n.Cmp(&m) == 1 {
			m.Mul(&m, &bigp)
		}
		// went one too far
		m.Div(&m, &bigp)

		a.Exp(a, &m, n)

		if i > 0 && i%16 == 0 {
			var aMinus1 big.Int
			aMinus1.Sub(a, &one)
			if newN, newG, ok := checkGCD(n, &aMinus1); ok {
				pr, co := PMinus1(newN)
				primes, composites = append(primes, pr...), append(composites, co...)
				pr, co = PMinus1(newG)
				primes, composites = append(primes, pr...), append(composites, co...)
				return primes, composites
			}
		}
	}

	return nil, []*big.Int{n}
}

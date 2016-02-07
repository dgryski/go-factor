package factor

import (
	"math/big"
)

func Trial(n *big.Int) (primes []*big.Int, composites []*big.Int) {

	var z big.Int
	var m big.Int

	one := big.NewInt(1)

	for _, p := range smallPrimes {
		for {
			bigp := big.NewInt(p)
			z.DivMod(n, bigp, &m)
			if m.Sign() != 0 {
				break
			}
			primes = append(primes, bigp)
			n.Set(&z)
		}

		if n.Cmp(one) == 0 {
			break
		}
	}

	if n.Cmp(one) > 0 {
		if n.ProbablyPrime(10) {
			primes = append(primes, n)
		} else {
			composites = append(composites, n)
		}
	}

	return primes, composites
}

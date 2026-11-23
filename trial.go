package factor

import (
	"math/big"
)

// Trial divides n by the small primes. Values below 1 are returned unchanged.
func Trial(n *big.Int) (primes []*big.Int, composites []*big.Int) {
	if n.Sign() <= 0 {
		return nil, []*big.Int{n}
	}
	n = new(big.Int).Set(n)

	var z big.Int
	var m big.Int

	bigp := big.NewInt(0)

	for _, p := range smallPrimes {
		for {
			// reuse the existing bigint
			bigp.SetInt64(p)
			z.DivMod(n, bigp, &m)
			if m.Sign() != 0 {
				break
			}
			n.Set(&z)
			primes = append(primes, bigp)
			// create a new bigint, since we just stored the one we were using
			bigp = big.NewInt(0)
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

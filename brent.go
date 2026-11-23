package factor

import (
	"math/big"
)

func checkGCD(n, g *big.Int) (newN, newG *big.Int, ok bool) {

	var z big.Int

	g.Abs(g)

	z.GCD(nil, nil, n, g)

	if z.Cmp(n) == 0 {
		return nil, nil, false
	}

	if z.Cmp(one) == 0 {
		return nil, nil, false
	}

	n.Div(n, &z)

	return n, &z, true
}

func Brent(n *big.Int, start, c int64) (primes, composites []*big.Int) {

	if n.ProbablyPrime(10) {
		return []*big.Int{n}, nil
	}

	bigc := big.NewInt(c)
	f := func(x *big.Int) {
		x.Mul(x, x)
		x.Add(x, bigc)
		x.Mod(x, n)
	}

	x1 := big.NewInt(start)
	x1.Mod(x1, n)
	x2 := new(big.Int).Set(x1)
	f(x2)

	// try to get off the tail
	for j := 0; j < 1000; j++ {
		f(x2)
	}

	split := func(newN, newG *big.Int) (primes, composites []*big.Int) {
		pr, co := Brent(newN, start, c)
		primes, composites = append(primes, pr...), append(composites, co...)
		pr, co = Brent(newG, start, c)
		primes, composites = append(primes, pr...), append(composites, co...)
		return primes, composites
	}

	limit := 1
	product := big.NewInt(1)
	var tmp big.Int

	// A batch never spans two rounds, since x1 changes between them.
	terms, batch := 0, 0
loop:
	for terms < (1 << 16) {
		for j := 0; j < limit; j++ {
			f(x2)

			if x1.Cmp(x2) == 0 {
				// The sequence cycled mod n without finding a factor.
				break loop
			}

			tmp.Sub(x1, x2)
			product.Mul(product, &tmp)
			product.Mod(product, n)

			terms++
			batch++
			if batch == 16 || j == limit-1 {
				if newN, newG, ok := checkGCD(n, product); ok {
					return split(newN, newG)
				}
				product.SetInt64(1)
				batch = 0
			}
		}

		x1.Set(x2)
		limit *= 2
		for j := 0; j < limit; j++ {
			f(x2)
		}
	}

	if batch > 0 {
		if newN, newG, ok := checkGCD(n, product); ok {
			return split(newN, newG)
		}
	}

	return nil, []*big.Int{n}
}

package factor

import (
	"context"
	"math/big"
)

func checkGCD(n, g *big.Int) (newN, newG *big.Int, ok bool) {
	z := new(big.Int).GCD(nil, nil, n, g)

	if z.Cmp(n) == 0 || z.Cmp(one) == 0 {
		return nil, nil, false
	}

	return new(big.Int).Div(n, z), z, true
}

// Brent factors n with Brent's variant of Pollard's rho, iterating
// x^2 + c from start for up to 65536 terms.
func Brent(n *big.Int, start, c int64) (primes, composites []*big.Int) {
	return BrentContext(context.Background(), n, start, c)
}

// BrentContext is Brent that stops early when ctx is done, after a last gcd,
// and returns what it has found so far.
func BrentContext(ctx context.Context, n *big.Int, start, c int64) (primes, composites []*big.Int) {

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
		pr, co := BrentContext(ctx, newN, start, c)
		primes, composites = append(primes, pr...), append(composites, co...)
		pr, co = BrentContext(ctx, newG, start, c)
		primes, composites = append(primes, pr...), append(composites, co...)
		return primes, composites
	}

	limit := 1
	product := big.NewInt(1)
	var tmp big.Int

	// x2 at the start of the batch, for backtracking
	ys := new(big.Int)

	// batchGCD returns gcd(product, n), or nil if it is 1. If it is n the
	// batch is redone one term at a time to try to separate the factors.
	batchGCD := func(count int) *big.Int {
		g := new(big.Int).GCD(nil, nil, product, n)
		if g.Cmp(n) == 0 {
			y := new(big.Int).Set(ys)
			for k := 0; k < count; k++ {
				f(y)
				tmp.Sub(x1, y)
				tmp.Mod(&tmp, n)
				if g.GCD(nil, nil, &tmp, n); g.Cmp(one) != 0 {
					break
				}
			}
		}
		if g.Cmp(one) == 0 {
			return nil
		}
		return g
	}

	// found splits n by g, or gives up if g is n.
	found := func(g *big.Int) (primes, composites []*big.Int) {
		if g.Cmp(n) == 0 {
			return nil, []*big.Int{n}
		}
		return split(new(big.Int).Div(n, g), g)
	}

	// A batch never spans two rounds, since x1 changes between them.
	terms, batch := 0, 0
loop:
	for terms < (1 << 16) {
		for j := 0; j < limit; j++ {
			if batch == 0 {
				ys.Set(x2)
			}
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
				if g := batchGCD(batch); g != nil {
					return found(g)
				}
				if ctx.Err() != nil {
					return nil, []*big.Int{n}
				}
				product.SetInt64(1)
				batch = 0
			}
		}

		x1.Set(x2)
		limit *= 2
		for j := 0; j < limit; j++ {
			if batch == 0 {
				ys.Set(x2)
			}
			f(x2)
		}
	}

	if batch > 0 {
		if g := batchGCD(batch); g != nil {
			return found(g)
		}
	}

	return nil, []*big.Int{n}
}

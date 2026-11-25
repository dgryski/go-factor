package factor

import (
	"context"
	"math/big"
	"time"
)

// PMinus1 factors n with Pollard's p-1 method, using the primes below 50000
// each raised to the largest power below n.
func PMinus1(n *big.Int) (primes, composites []*big.Int) {
	return pminus1(context.Background(), n, smallPrimes, nil)
}

// PMinus1Context factors n with Pollard's p-1 method using the prime powers
// up to b1, and powers of 2 up to n. It stops early when ctx is done, after
// a last gcd, and returns what it has found so far.
func PMinus1Context(ctx context.Context, n *big.Int, b1 int) (primes, composites []*big.Int) {
	ps := sievePrimes(b1)
	primes64 := make([]int64, len(ps))
	for i, p := range ps {
		primes64[i] = int64(p)
	}
	return pminus1(ctx, n, primes64, big.NewInt(int64(b1)))
}

// pminus1 raises 2 to each prime in ps, to the largest power at most bound,
// or below n for 2 or when bound is nil.
func pminus1(ctx context.Context, n *big.Int, ps []int64, bound *big.Int) (primes, composites []*big.Int) {
	if n.Cmp(one) <= 0 {
		return nil, nil
	}

	if n.ProbablyPrime(10) {
		return []*big.Int{n}, nil
	}

	if n.Bit(0) == 0 {
		var half big.Int
		half.Rsh(n, 1)
		pr, co := pminus1(ctx, &half, ps, bound)
		return append([]*big.Int{big.NewInt(2)}, pr...), co
	}

	t0 := time.Now()
	if Logf != nil {
		b1 := "n"
		if bound != nil {
			b1 = bound.String()
		}
		logf("p-1: stage 1 with B1=%s (%d primes) on %d bits", b1, len(ps), n.BitLen())
	}

	// The base must be coprime to n or its common primes are never found.
	a := big.NewInt(2)

	// a at the start of the current batch, for backtracking
	saved := new(big.Int).Set(a)
	start := 0

	var g, aMinus1 big.Int
	gcd := func() {
		aMinus1.Sub(a, one)
		g.GCD(nil, nil, &aMinus1, n)
	}

	// within reports whether the power m of p is in the exponent: below n,
	// or at most bound.
	within := func(m *big.Int, p int64) bool {
		if bound == nil || p == 2 {
			return m.Cmp(n) < 0
		}
		return m.Cmp(bound) <= 0
	}

	// The batch's prime powers go into one exponent. Exp only uses Montgomery
	// multiplication for exponents of more than one word, and division for
	// each step is about twice as slow.
	var bigp, m, next big.Int
	exp := big.NewInt(1)
	for i, p := range ps {

		bigp.SetInt64(p)
		m.Set(&bigp)
		for within(next.Mul(&m, &bigp), p) {
			m.Set(&next)
		}
		exp.Mul(exp, &m)

		if (i > 0 && i%16 == 0) || i == len(ps)-1 {
			a.Exp(a, exp, n)
			exp.SetInt64(1)
			gcd()
			if g.Cmp(n) == 0 {
				// Every prime of n was found in this batch. Redo it one
				// prime power at a time to separate them.
				a.Set(saved)
				g.SetInt64(1)
			replay:
				for _, q := range ps[start : i+1] {
					bigp.SetInt64(q)
					for m.Set(&bigp); m.Cmp(&bigp) == 0 || within(&m, q); m.Mul(&m, &bigp) {
						a.Exp(a, &bigp, n)
						if gcd(); g.Cmp(one) != 0 {
							break replay
						}
					}
				}
				if g.Cmp(n) == 0 {
					logf("p-1: found only n after %d of %d primes in %v", i+1, len(ps), since(t0))
					return nil, []*big.Int{n}
				}
			}
			if newN, newG, ok := checkGCD(n, &aMinus1); ok {
				logf("p-1: found %v after %d of %d primes in %v", newG, i+1, len(ps), since(t0))
				pr, co := pminus1(ctx, newN, ps, bound)
				primes, composites = append(primes, pr...), append(composites, co...)
				pr, co = pminus1(ctx, newG, ps, bound)
				primes, composites = append(primes, pr...), append(composites, co...)
				return primes, composites
			}
			if ctx.Err() != nil {
				logf("p-1: stopped after %d of %d primes in %v", i+1, len(ps), since(t0))
				return nil, []*big.Int{n}
			}
			saved.Set(a)
			start = i + 1
		}
	}

	logf("p-1: found nothing in %v", since(t0))
	return nil, []*big.Int{n}
}

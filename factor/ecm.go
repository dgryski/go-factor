package main

import (
	"log"
	"math/big"
	"time"

	"github.com/dgryski/go-factor"
)

// ecmLevels gives B1, B2 and the expected number of curves to find a factor
// of the given size, from Table 1 of the GMP-ECM README. The 15 digit row
// uses GMP-ECM's default B2 for B1=2000.
var ecmLevels = []struct{ digits, b1, b2, curves int }{
	{15, 2000, 147396, 25},
	{20, 11000, 1900000, 74},
	{25, 50000, 13000000, 221},
	{30, 250000, 130000000, 453},
}

// With no quadratic sieve to follow, ECM runs every level up to here.
const ecmNoQSDigits = 25

// ecmPretest runs ECM on n with rising bounds until b's deadline, since more
// curves after that are likely to cost more than they save.
func ecmPretest(n *big.Int, b *budget) (primes, composites []*big.Int) {
	todo := []*big.Int{n}
	_, _, hasSieve := b.sieveEstimate(n)
	for _, lv := range ecmLevels {
		if len(todo) == 0 || !hasSieve && lv.digits > ecmNoQSDigits {
			break
		}
		var next []*big.Int
		for _, c := range todo {
			ctx, cancel, sieveName, ok := b.context(c, shareECM)
			if !ok {
				next = append(next, c)
				continue
			}
			log.Printf("starting ecm(b1=%v, b2=%v, curves=%v) on %d digits", lv.b1, lv.b2, lv.curves, len(c.String()))
			p, co := factor.ECMContext(ctx, c, lv.b1, lv.b2, lv.curves)
			if ctx.Err() != nil {
				log.Printf("ecm stopped %v after brent started, leaving the rest to %s", time.Since(b.start).Round(time.Millisecond), sieveName)
			}
			cancel()
			primes = append(primes, p...)
			next = append(next, co...)
		}
		todo = next
	}
	return primes, todo
}

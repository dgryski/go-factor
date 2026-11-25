package main

import (
	"context"
	"log"
	"math"
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

const (
	// ECM gets this share of the quadratic sieve's estimated time.
	ecmQSFraction = 1.0 / 3

	// With no quadratic sieve to follow, ECM runs every level up to here.
	ecmNoQSDigits = 25
)

// qsSingle is the quadratic sieve's time with one worker for balanced
// semiprimes, measured on an Apple M-series machine.
var qsSingle = []struct {
	digits  int
	seconds float64
}{
	{30, 0.004}, {39, 0.019}, {45, 0.081}, {50, 0.24}, {55, 0.60},
	{60, 1.54}, {64, 4.9}, {67, 9.5}, {71, 20.5}, {74, 52}, {79, 159},
}

// qsTime estimates the quadratic sieve's time for n with the given number
// of workers. It interpolates qsSingle in log space, extending the end
// slopes, and assumes each extra worker adds 0.55 of one worker's speed.
func qsTime(n *big.Int, workers int) time.Duration {
	d := float64(len(n.String()))
	i := 1
	for i < len(qsSingle)-1 && d > float64(qsSingle[i].digits) {
		i++
	}
	lo, hi := qsSingle[i-1], qsSingle[i]
	slope := (math.Log(hi.seconds) - math.Log(lo.seconds)) / float64(hi.digits-lo.digits)
	seconds := math.Exp(math.Log(lo.seconds) + slope*(d-float64(lo.digits)))
	speedup := 1 + 0.55*float64(workers-1)
	return time.Duration(seconds / speedup * float64(time.Second))
}

// ecmPretest runs ECM on n with rising bounds. When the quadratic sieve will
// run afterwards it stops once ECM has used ecmQSFraction of the sieve's
// estimated time for what is left, since more curves are then likely to cost
// more than they save.
func ecmPretest(n *big.Int, qsRuns bool, workers int) (primes, composites []*big.Int) {
	start := time.Now()
	todo := []*big.Int{n}
	for _, lv := range ecmLevels {
		if len(todo) == 0 || !qsRuns && lv.digits > ecmNoQSDigits {
			break
		}
		var next []*big.Int
		for _, c := range todo {
			ctx, cancel := context.Background(), context.CancelFunc(func() {})
			if qsRuns {
				deadline := start.Add(time.Duration(ecmQSFraction * float64(qsTime(c, workers))))
				if !time.Now().Before(deadline) {
					next = append(next, c)
					continue
				}
				ctx, cancel = context.WithDeadline(ctx, deadline)
			}
			log.Printf("starting ecm(b1=%v, b2=%v, curves=%v) on %d digits", lv.b1, lv.b2, lv.curves, len(c.String()))
			p, co := factor.ECMContext(ctx, c, lv.b1, lv.b2, lv.curves)
			if ctx.Err() != nil {
				log.Printf("ecm stopped after %v, leaving the rest to qs", time.Since(start).Round(time.Millisecond))
			}
			cancel()
			primes = append(primes, p...)
			next = append(next, co...)
		}
		todo = next
	}
	return primes, todo
}

package main

import (
	"context"
	"math"
	"math/big"
	"time"
)

// The methods before the quadratic sieve share pretestQSFraction of its
// estimated time, counted from the start of brent. Each may run until the
// cumulative share below is used, so time one method leaves goes to the next.
const (
	pretestQSFraction = 1.0 / 3

	shareBrent  = 0.10
	sharePMinus = 0.20
	sharePPlus  = 0.35
	shareECM    = 1.0
)

// budget hands out deadlines to the methods before the quadratic sieve.
type budget struct {
	start   time.Time
	workers int                 // quadratic sieve workers
	qsRuns  func(*big.Int) bool // whether the sieve will run on a composite
}

// context returns a context for working on c until share of its budget is
// used. ok is false if that time has passed. If the quadratic sieve will not
// run on c there is no deadline.
func (b *budget) context(c *big.Int, share float64) (ctx context.Context, cancel context.CancelFunc, ok bool) {
	if !b.qsRuns(c) {
		ctx, cancel = context.WithCancel(context.Background())
		return ctx, cancel, true
	}
	d := time.Duration(share * pretestQSFraction * float64(qsTime(c, b.workers)))
	deadline := b.start.Add(d)
	if !time.Now().Before(deadline) {
		return nil, nil, false
	}
	ctx, cancel = context.WithDeadline(context.Background(), deadline)
	return ctx, cancel, true
}

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

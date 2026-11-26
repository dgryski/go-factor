package main

import (
	"context"
	"math"
	"math/big"
	"time"
)

// Methods before the sieve share pretestFraction of the estimated sieve time.
// Unused time from earlier methods carries forward to later ones.
const (
	pretestFraction = 1.0 / 3

	shareBrent  = 0.10
	sharePMinus = 0.20
	sharePPlus  = 0.35
	shareECM    = 1.0
)

type timePoint struct {
	digits  int
	seconds float64
}

// budget hands out deadlines to the methods that run before qs or gnfs.
type budget struct {
	start    time.Time
	workers  int
	qsRuns   func(*big.Int) bool
	gnfsRuns func(*big.Int) bool
}

func (b *budget) sieveEstimate(c *big.Int) (d time.Duration, sieveName string, ok bool) {
	if b.qsRuns != nil && b.qsRuns(c) {
		return qsTime(c, b.workers), "qs", true
	}
	if b.gnfsRuns != nil && b.gnfsRuns(c) {
		return gnfsTime(c, b.workers), "gnfs", true
	}
	return 0, "", false
}

// context returns a deadline context for share of the pretest budget on c.
// If neither qs nor gnfs will run on c, no deadline is set.
func (b *budget) context(c *big.Int, share float64) (ctx context.Context, cancel context.CancelFunc, sieveName string, ok bool) {
	est, name, hasSieve := b.sieveEstimate(c)
	if !hasSieve {
		ctx, cancel = context.WithCancel(context.Background())
		return ctx, cancel, "", true
	}
	d := time.Duration(share * pretestFraction * float64(est))
	deadline := b.start.Add(d)
	if !time.Now().Before(deadline) {
		return nil, nil, name, false
	}
	ctx, cancel = context.WithDeadline(context.Background(), deadline)
	return ctx, cancel, name, true
}

// qsSingle holds 1-worker quadratic sieve times on an Apple M-series machine.
var qsSingle = []timePoint{
	{30, 0.004}, {39, 0.019}, {45, 0.081}, {50, 0.24}, {55, 0.60},
	{60, 1.54}, {64, 4.9}, {67, 9.5}, {71, 20.5}, {74, 52},
	{79, 159}, {90, 1140}, {95, 4750},
}

// gnfs15Cores holds 15-worker GNFS times on an Apple M5 Pro.
var gnfs15Cores = []timePoint{
	{35, 0.25}, {45, 0.60}, {59, 3.0}, {70, 17}, {75, 95}, {80, 274}, {89, 1311},
}

func interpLogTime(table []timePoint, digits int) float64 {
	d := float64(digits)
	i := 1
	for i < len(table)-1 && d > float64(table[i].digits) {
		i++
	}
	lo, hi := table[i-1], table[i]
	slope := (math.Log(hi.seconds) - math.Log(lo.seconds)) / float64(hi.digits-lo.digits)
	return math.Exp(math.Log(lo.seconds) + slope*(d-float64(lo.digits)))
}

func qsTime(n *big.Int, workers int) time.Duration {
	if workers < 1 {
		workers = 1
	}
	seconds := interpLogTime(qsSingle, len(n.String()))
	speedup := 1 + 0.55*float64(workers-1)
	return time.Duration(seconds / speedup * float64(time.Second))
}

func gnfsTime(n *big.Int, workers int) time.Duration {
	if workers < 1 {
		workers = 1
	}
	seconds15 := interpLogTime(gnfs15Cores, len(n.String()))
	scale := (1 + 0.65*14) / (1 + 0.65*float64(workers-1))
	return time.Duration(seconds15 * scale * float64(time.Second))
}

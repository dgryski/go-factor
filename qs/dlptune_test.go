package qs

import (
	"flag"
	"fmt"
	"math/big"
	"math/rand"
	"testing"
	"time"

	"github.com/dgryski/go-factor/internal/testutil"
)

var tune = flag.Bool("tune", false, "run the slow tuning sweeps")

// TestDLPTune times split with double large primes off and at several
// threshold allowances.
func TestDLPTune(t *testing.T) {
	if !*tune {
		t.Skip("slow; run with -tune")
	}
	for _, bits := range []int{110, 120} {
		rnd := rand.New(rand.NewSource(7))
		n := new(big.Int).Mul(testutil.RandomPrime(rnd, bits), testutil.RandomPrime(rnd, bits))
		for _, c := range []struct {
			minBits int
			frac    float64
		}{{9999, 0}, {0, 0.3}, {0, 0.35}, {0, 0.4}, {0, 0.45}, {0, 0.55}} {
			dlpMinBits, dlpSlackFrac = c.minBits, c.frac
			var last string
			Logf = func(f string, a ...interface{}) {
				if len(f) > 8 && f[:8] == "qs: %d/%" {
					last = fmt.Sprintf(f, a...)
				}
			}
			st := time.Now()
			p, _ := split(n)
			Logf = nil
			t.Logf("%d bits dlp=%v frac=%v: ok=%v %v | %s", n.BitLen(), c.minBits == 0, c.frac, p != nil, time.Since(st).Round(time.Millisecond), last)
		}
	}
}

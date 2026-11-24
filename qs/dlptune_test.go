package qs

import (
	"flag"
	"fmt"
	"math/big"
	"math/rand"
	"os"
	"strconv"
	"strings"
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

// TestSweep times split on QS_N for each QS_RANGE, QS_DLP slack fraction
// (0 means off, others start at 250 bits) and QS_WORKERS count.
func TestSweep(t *testing.T) {
	if !*tune || os.Getenv("QS_N") == "" {
		t.Skip("slow; run with -tune and QS_N, QS_RANGE, QS_DLP, QS_WORKERS")
	}
	n, _ := new(big.Int).SetString(os.Getenv("QS_N"), 10)
	list := func(k string) []string { return strings.Split(os.Getenv(k), ",") }
	for _, rs := range list("QS_RANGE") {
		largeSieveRange, _ = strconv.ParseInt(rs, 10, 64)
		for _, fs := range list("QS_DLP") {
			frac, _ := strconv.ParseFloat(fs, 64)
			if frac == 0 {
				dlpMinBits = 1 << 30
			} else {
				dlpMinBits, dlpSlackFrac = 250, frac
			}
			for _, ws := range list("QS_WORKERS") {
				Workers, _ = strconv.Atoi(ws)
				var prog, mat string
				var sieved time.Duration
				st := time.Now()
				Logf = func(f string, a ...interface{}) {
					switch {
					case strings.HasPrefix(f, "qs: %d/%"):
						prog = fmt.Sprintf(f, a...)
						sieved = time.Since(st)
					case strings.HasPrefix(f, "qs: matrix"):
						mat = fmt.Sprintf(f, a...)
					}
				}
				p, _ := split(n)
				Logf = nil
				t.Logf("M=%dK frac=%v workers=%d: ok=%v total %v sieve %v | %s | %s",
					largeSieveRange/1024, frac, Workers, p != nil,
					time.Since(st).Round(time.Millisecond), sieved.Round(time.Millisecond), prog, mat)
			}
		}
	}
}

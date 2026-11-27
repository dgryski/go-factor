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
	defer func(dm int, df float64) { dlpMinBits, dlpSlackFrac = dm, df }(dlpMinBits, dlpSlackFrac)
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

// TestSweep times split on QS_N (or seed-7 semiprimes of QS_BITS bits) across
// QS_B, QS_LPM, QS_SMALL, QS_RANGE, QS_DLP, and QS_WORKERS.
func TestSweep(t *testing.T) {
	if !*tune || (os.Getenv("QS_N") == "" && os.Getenv("QS_BITS") == "") {
		t.Skip("slow; run with -tune and QS_N or QS_BITS")
	}
	list := func(k, def string) []string {
		v := os.Getenv(k)
		if v == "" {
			v = def
		}
		return strings.Split(v, ",")
	}
	var inputs []*big.Int
	if os.Getenv("QS_N") != "" {
		for _, s := range list("QS_N", "") {
			n, _ := new(big.Int).SetString(s, 10)
			inputs = append(inputs, n)
		}
	} else {
		seeds, _ := strconv.Atoi(os.Getenv("QS_SEEDS"))
		if seeds <= 0 {
			seeds = 1
		}
		for _, bStr := range list("QS_BITS", "") {
			bits, _ := strconv.Atoi(bStr)
			for seed := 1; seed <= seeds; seed++ {
				src := int64(7)
				if seeds > 1 {
					src = int64(seed*1000 + bits)
				}
				rnd := rand.New(rand.NewSource(src))
				h1 := bits / 2
				h2 := bits - h1
				inputs = append(inputs, new(big.Int).Mul(testutil.RandomPrime(rnd, h1), testutil.RandomPrime(rnd, h2)))
			}
		}
	}
	defer func(b int64, l uint64, s uint32, r int64, dm int, df float64) {
		fbBoundOverride, largePrimeMultOverride, minSievePrimeOverride, largeSieveRange = b, l, s, r
		dlpMinBits, dlpSlackFrac = dm, df
	}(fbBoundOverride, largePrimeMultOverride, minSievePrimeOverride, largeSieveRange, dlpMinBits, dlpSlackFrac)
	defRange, defDM, defDF := largeSieveRange, dlpMinBits, dlpSlackFrac
	for _, n := range inputs {
		for _, bs := range list("QS_B", "0") {
			fbBoundOverride, _ = strconv.ParseInt(bs, 10, 64)
			for _, ls := range list("QS_LPM", "0") {
				largePrimeMultOverride, _ = strconv.ParseUint(ls, 10, 64)
				for _, ss := range list("QS_SMALL", "0") {
					v, _ := strconv.ParseUint(ss, 10, 32)
					minSievePrimeOverride = uint32(v)
					for _, rs := range list("QS_RANGE", "0") {
						largeSieveRange = defRange
						if rv, _ := strconv.ParseInt(rs, 10, 64); rv > 0 {
							largeSieveRange = rv
						}
						for _, fs := range list("QS_DLP", "default") {
							dlpMinBits, dlpSlackFrac = defDM, defDF
							if fs != "default" {
								frac, _ := strconv.ParseFloat(fs, 64)
								if frac == 0 {
									dlpMinBits = 1 << 30
								} else {
									dlpMinBits, dlpSlackFrac = 0, frac
								}
							}
							for _, ws := range list("QS_WORKERS", "0") {
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
								fmt.Printf("  [%db] B=%s lpm=%s small=%s M=%s frac=%s: ok=%v total %v sieve %v | %s | %s\n",
									n.BitLen(), bs, ls, ss, rs, fs, p != nil,
									time.Since(st).Round(time.Millisecond), sieved.Round(time.Millisecond), prog, mat)
							}
						}
					}
				}
			}
		}
	}
}

package gnfs

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

	"github.com/dgryski/go-factor/internal/arith"
	"github.com/dgryski/go-factor/internal/testutil"
)

var tune = flag.Bool("tune", false, "run the slow tuning sweeps")

// tuneCands are parameter sets to compare in TestSweep. "default" uses
// defaultParams for each number.
var tuneCands = map[string]Params{
	"line_1200": {Degree: 3, Br: 1200, Ba: 1800, LpBoundRat: 32000, LpBoundAlg: 48000,
		MfbRat: 30, MfbAlg: 32, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: false, LineA: 65536, LineB: 1000,
		SlackRat: 34, SlackAlg: 36, NumQCs: 64},
	"line_2000": {Degree: 3, Br: 2000, Ba: 2500, LpBoundRat: 48000, LpBoundAlg: 64000,
		MfbRat: 32, MfbAlg: 32, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: false, LineA: 65536, LineB: 1000,
		SlackRat: 34, SlackAlg: 36, NumQCs: 64},
	"lat_110": {Degree: 3, Br: 3000, Ba: 3000, LpBoundRat: 1 << 16, LpBoundAlg: 1 << 16,
		MfbRat: 32, MfbAlg: 32, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 3000, QMax: 24000, IHalf: 512, JMax: 256,
		SlackRat: 34, SlackAlg: 34, NumQCs: 64},
	// The 200-250 bit row before 2026-11-25.
	"old250": {Degree: 4, Br: 10000, Ba: 20000, LpBoundRat: 1 << 21, LpBoundAlg: 1 << 21,
		MfbRat: 42, MfbAlg: 42, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 20000, QMax: 250000, IHalf: 1024, JMax: 1024,
		SlackRat: 48, SlackAlg: 48, NumQCs: 64},
	"60k": {Degree: 4, Br: 60000, Ba: 60000, LpBoundRat: 1 << 19, LpBoundAlg: 1 << 19,
		MfbRat: 38, MfbAlg: 38, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 60000, QMax: 240000, IHalf: 1024, JMax: 512,
		SlackRat: 44, SlackAlg: 44, NumQCs: 64},
	"100k": {Degree: 4, Br: 100000, Ba: 100000, LpBoundRat: 1 << 20, LpBoundAlg: 1 << 20,
		MfbRat: 40, MfbAlg: 40, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 100000, QMax: 400000, IHalf: 1024, JMax: 1024,
		SlackRat: 46, SlackAlg: 46, NumQCs: 64},
	"160k_21": {Degree: 4, Br: 160000, Ba: 160000, LpBoundRat: 1 << 21, LpBoundAlg: 1 << 21,
		MfbRat: 42, MfbAlg: 42, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 160000, QMax: 800000, IHalf: 2048, JMax: 1024,
		SlackRat: 48, SlackAlg: 48, NumQCs: 64},
	"160k_21_1k": {Degree: 4, Br: 160000, Ba: 160000, LpBoundRat: 1 << 21, LpBoundAlg: 1 << 21,
		MfbRat: 42, MfbAlg: 42, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 160000, QMax: 800000, IHalf: 1024, JMax: 1024,
		SlackRat: 48, SlackAlg: 48, NumQCs: 64},
	"200k_21_1k": {Degree: 4, Br: 200000, Ba: 200000, LpBoundRat: 1 << 21, LpBoundAlg: 1 << 21,
		MfbRat: 42, MfbAlg: 42, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 200000, QMax: 1000000, IHalf: 1024, JMax: 1024,
		SlackRat: 48, SlackAlg: 48, NumQCs: 64},
	"200k_21": {Degree: 4, Br: 200000, Ba: 200000, LpBoundRat: 1 << 21, LpBoundAlg: 1 << 21,
		MfbRat: 42, MfbAlg: 42, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 200000, QMax: 1000000, IHalf: 2048, JMax: 1024,
		SlackRat: 48, SlackAlg: 48, NumQCs: 64},
	"200k_q80k": {Degree: 4, Br: 200000, Ba: 200000, LpBoundRat: 1 << 21, LpBoundAlg: 1 << 21,
		MfbRat: 42, MfbAlg: 42, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 80000, QMax: 800000, IHalf: 2048, JMax: 1024,
		SlackRat: 48, SlackAlg: 48, NumQCs: 64},
	"250k_q80k": {Degree: 4, Br: 250000, Ba: 250000, LpBoundRat: 1 << 21, LpBoundAlg: 1 << 21,
		MfbRat: 42, MfbAlg: 42, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 80000, QMax: 800000, IHalf: 2048, JMax: 1024,
		SlackRat: 43, SlackAlg: 43, NumQCs: 64},
	"300k_22_q120k": {Degree: 4, Br: 300000, Ba: 300000, LpBoundRat: 1 << 22, LpBoundAlg: 1 << 22,
		MfbRat: 44, MfbAlg: 44, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 120000, QMax: 1500000, IHalf: 2048, JMax: 1024,
		SlackRat: 50, SlackAlg: 50, NumQCs: 64},
	"400k_22_q150k": {Degree: 4, Br: 400000, Ba: 400000, LpBoundRat: 1 << 22, LpBoundAlg: 1 << 22,
		MfbRat: 44, MfbAlg: 44, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 150000, QMax: 2000000, IHalf: 2048, JMax: 2048,
		SlackRat: 50, SlackAlg: 50, NumQCs: 64},
	"250k_22_1k": {Degree: 4, Br: 250000, Ba: 250000, LpBoundRat: 1 << 22, LpBoundAlg: 1 << 22,
		MfbRat: 44, MfbAlg: 44, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 250000, QMax: 1500000, IHalf: 1024, JMax: 1024,
		SlackRat: 50, SlackAlg: 50, NumQCs: 64},
	"200k_22": {Degree: 4, Br: 200000, Ba: 250000, LpBoundRat: 1 << 22, LpBoundAlg: 1 << 22,
		MfbRat: 44, MfbAlg: 44, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 250000, QMax: 1200000, IHalf: 2048, JMax: 1024,
		SlackRat: 50, SlackAlg: 50, NumQCs: 64},
	"300k_22": {Degree: 4, Br: 300000, Ba: 300000, LpBoundRat: 1 << 22, LpBoundAlg: 1 << 22,
		MfbRat: 44, MfbAlg: 44, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 300000, QMax: 2000000, IHalf: 2048, JMax: 1024,
		SlackRat: 50, SlackAlg: 50, NumQCs: 64},
	"350k_22": {Degree: 4, Br: 350000, Ba: 400000, LpBoundRat: 1 << 22, LpBoundAlg: 1 << 22,
		MfbRat: 44, MfbAlg: 44, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 400000, QMax: 2500000, IHalf: 2048, JMax: 1024,
		SlackRat: 50, SlackAlg: 50, NumQCs: 64},
	"300k_23": {Degree: 4, Br: 300000, Ba: 350000, LpBoundRat: 1 << 23, LpBoundAlg: 1 << 23,
		MfbRat: 46, MfbAlg: 46, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 350000, QMax: 2000000, IHalf: 2048, JMax: 2048,
		SlackRat: 52, SlackAlg: 52, NumQCs: 64},
}

// TestSweep times splitWithParams on random semiprimes of each GNFS_BITS
// size, GNFS_SEEDS numbers per size, for each candidate in GNFS_CANDS. It
// logs the relations collected against π(LpBoundRat) + π(LpBoundAlg), which
// calibrates relTargetRatio.
func TestSweep(t *testing.T) {
	if !*tune || os.Getenv("GNFS_BITS") == "" {
		t.Skip("slow; run with -tune and GNFS_BITS, GNFS_CANDS (default " + candNames() + "), GNFS_SEEDS")
	}
	list := func(k, def string) []string {
		if v := os.Getenv(k); v != "" {
			return strings.Split(v, ",")
		}
		return strings.Split(def, ",")
	}
	seeds, _ := strconv.Atoi(os.Getenv("GNFS_SEEDS"))
	if seeds <= 0 {
		seeds = 1
	}
	for _, bs := range list("GNFS_BITS", "") {
		bits, err := strconv.Atoi(bs)
		if err != nil || bits < 8 {
			t.Fatalf("bad GNFS_BITS entry %q", bs)
		}
		for seed := 1; seed <= seeds; seed++ {
			rnd := rand.New(rand.NewSource(int64(seed*1000 + bits)))
			n := new(big.Int).Mul(testutil.RandomPrime(rnd, bits/2), testutil.RandomPrime(rnd, bits-bits/2))
			for _, name := range list("GNFS_CANDS", "default") {
				var p Params
				if name == "default" {
					p = *defaultParams(n)
				} else if c, ok := tuneCands[name]; ok {
					p = c
				} else {
					t.Fatalf("unknown candidate %q: want default or one of %s", name, candNames())
				}
				var fbLine, sieved, laLine, sqrtLine string
				lastPrint := time.Now()
				Logf = func(f string, a ...interface{}) {
					s := fmt.Sprintf(f, a...)
					switch {
					case strings.Contains(s, "cols="):
						fbLine = s
					case strings.Contains(s, "sieved"):
						sieved = s
					case strings.Contains(s, "dependencies in"):
						laLine = s
					case strings.Contains(s, "sqrt time"):
						sqrtLine = s
					case strings.Contains(s, "relations (excess"):
						if time.Since(lastPrint) >= 15*time.Second {
							lastPrint = time.Now()
							fmt.Printf("  [%d bits seed %d %s] %s\n", n.BitLen(), seed, name, s)
						}
					}
				}
				st := time.Now()
				f1, _ := splitWithParams(n, &p)
				Logf = nil
				ideals := len(arith.PrimesUpTo(int64(p.LpBoundRat))) + len(arith.PrimesUpTo(int64(p.LpBoundAlg)))
				t.Logf("%d bits (%d digits) seed %d %s: ok=%v %v, LP ideals %d | %s | %s | %s | %s",
					n.BitLen(), len(n.String()), seed, name, f1 != nil, time.Since(st).Round(time.Millisecond), ideals, fbLine, sieved, laLine, sqrtLine)
			}
		}
	}
}

func candNames() string {
	names := []string{"default"}
	for k := range tuneCands {
		names = append(names, k)
	}
	return strings.Join(names, ",")
}

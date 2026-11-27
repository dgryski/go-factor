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
	"m3_295_j1k": {Degree: 4, Br: 400000, Ba: 400000, LpBoundRat: 1 << 22, LpBoundAlg: 1 << 22,
		MfbRat: 44, MfbAlg: 44, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 120000, QMax: 2000000, IHalf: 2048, JMax: 1024,
		SlackRat: 45, SlackAlg: 45, NumQCs: 64, RelTargetRatio: 0.84},
	"m3_295_320k": {Degree: 4, Br: 320000, Ba: 350000, LpBoundRat: 3200000, LpBoundAlg: 3200000,
		MfbRat: 43, MfbAlg: 43, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 120000, QMax: 1600000, IHalf: 2048, JMax: 1024,
		SlackRat: 44, SlackAlg: 44, NumQCs: 64, RelTargetRatio: 0.84},
	"m3_315_4m_j2k": {Degree: 4, Br: 450000, Ba: 500000, LpBoundRat: 1 << 22, LpBoundAlg: 1 << 22,
		MfbRat: 44, MfbAlg: 44, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 180000, QMax: 2500000, IHalf: 2048, JMax: 2048,
		SlackRat: 45, SlackAlg: 45, NumQCs: 64, RelTargetRatio: 0.86},
	"m3_315_6m_j1k": {Degree: 4, Br: 480000, Ba: 520000, LpBoundRat: 6000000, LpBoundAlg: 6000000,
		MfbRat: 45, MfbAlg: 45, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 180000, QMax: 3000000, IHalf: 2048, JMax: 1024,
		SlackRat: 46, SlackAlg: 46, NumQCs: 64, RelTargetRatio: 0.85},
	"m3_315_6m_j2k": {Degree: 4, Br: 480000, Ba: 520000, LpBoundRat: 6000000, LpBoundAlg: 6000000,
		MfbRat: 45, MfbAlg: 45, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 200000, QMax: 3000000, IHalf: 2048, JMax: 2048,
		SlackRat: 46, SlackAlg: 46, NumQCs: 64, RelTargetRatio: 0.85},
	"m3_315_550k_6m": {Degree: 4, Br: 550000, Ba: 650000, LpBoundRat: 6000000, LpBoundAlg: 6000000,
		MfbRat: 45, MfbAlg: 45, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 180000, QMax: 3000000, IHalf: 2048, JMax: 2048,
		SlackRat: 46, SlackAlg: 46, NumQCs: 64, RelTargetRatio: 0.85},
	"m3_332_700k_23": {Degree: 4, Br: 650000, Ba: 750000, LpBoundRat: 1 << 23, LpBoundAlg: 1 << 23,
		MfbRat: 46, MfbAlg: 46, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 250000, QMax: 4500000, IHalf: 2048, JMax: 2048,
		SlackRat: 47, SlackAlg: 47, NumQCs: 64, RelTargetRatio: 0.78},
	"m3_332_750k_12m": {Degree: 4, Br: 700000, Ba: 800000, LpBoundRat: 12000000, LpBoundAlg: 12000000,
		MfbRat: 47, MfbAlg: 47, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 280000, QMax: 5000000, IHalf: 2048, JMax: 2048,
		SlackRat: 48, SlackAlg: 48, NumQCs: 64, RelTargetRatio: 0.78},
	"m3_332_700k_j1k": {Degree: 4, Br: 650000, Ba: 750000, LpBoundRat: 1 << 23, LpBoundAlg: 1 << 23,
		MfbRat: 46, MfbAlg: 46, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 250000, QMax: 4500000, IHalf: 2048, JMax: 1024,
		SlackRat: 47, SlackAlg: 47, NumQCs: 64, RelTargetRatio: 0.86},
	"m3_332_800k_11m": {Degree: 4, Br: 750000, Ba: 850000, LpBoundRat: 11000000, LpBoundAlg: 11000000,
		MfbRat: 47, MfbAlg: 47, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 350000, QMax: 5000000, IHalf: 2048, JMax: 2048,
		SlackRat: 48, SlackAlg: 48, NumQCs: 64, RelTargetRatio: 0.86},
	"m5_332_3lp": {Degree: 4, Br: 700000, Ba: 800000, LpBoundRat: 12000000, LpBoundAlg: 12000000,
		MfbRat: 62, MfbAlg: 64, MaxLPRat: 3, MaxLPAlg: 3,
		UseLattice: true, QMin: 280000, QMax: 5000000, IHalf: 4096, JMax: 2048,
		SlackRat: 50, SlackAlg: 53, NumQCs: 64, RelTargetRatio: 0.78},
	"m5_d4_2lp_4k": {Degree: 4, Br: 850000, Ba: 1000000, LpBoundRat: 1 << 24, LpBoundAlg: 1 << 24,
		MfbRat: 48, MfbAlg: 48, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 350000, QMax: 5000000, IHalf: 4096, JMax: 2048,
		SlackRat: 49, SlackAlg: 49, NumQCs: 64, RelTargetRatio: 0.80},
	"m5_d4_3lp_4k": {Degree: 4, Br: 850000, Ba: 1000000, LpBoundRat: 1 << 24, LpBoundAlg: 1 << 24,
		MfbRat: 64, MfbAlg: 66, MaxLPRat: 3, MaxLPAlg: 3,
		UseLattice: true, QMin: 350000, QMax: 5000000, IHalf: 4096, JMax: 2048,
		SlackRat: 52, SlackAlg: 55, NumQCs: 64, RelTargetRatio: 0.78},
	"m5_d4_3lp_hi": {Degree: 4, Br: 850000, Ba: 1000000, LpBoundRat: 1 << 24, LpBoundAlg: 1 << 24,
		MfbRat: 66, MfbAlg: 68, MaxLPRat: 3, MaxLPAlg: 3,
		UseLattice: true, QMin: 350000, QMax: 5000000, IHalf: 4096, JMax: 2048,
		SlackRat: 55, SlackAlg: 58, NumQCs: 64, RelTargetRatio: 0.78},
	"m5_d5_2lp_4k": {Degree: 5, Br: 850000, Ba: 1000000, LpBoundRat: 1 << 24, LpBoundAlg: 1 << 24,
		MfbRat: 48, MfbAlg: 48, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 350000, QMax: 5000000, IHalf: 4096, JMax: 2048,
		SlackRat: 49, SlackAlg: 49, NumQCs: 64, RelTargetRatio: 0.80},
	"m5_d5_3lp_4k": {Degree: 5, Br: 850000, Ba: 1000000, LpBoundRat: 1 << 24, LpBoundAlg: 1 << 24,
		MfbRat: 64, MfbAlg: 66, MaxLPRat: 3, MaxLPAlg: 3,
		UseLattice: true, QMin: 350000, QMax: 5000000, IHalf: 4096, JMax: 2048,
		SlackRat: 52, SlackAlg: 55, NumQCs: 64, RelTargetRatio: 0.78},
	"m5_d5_3lp_24m": {Degree: 5, Br: 1000000, Ba: 1200000, LpBoundRat: 24000000, LpBoundAlg: 24000000,
		MfbRat: 65, MfbAlg: 68, MaxLPRat: 3, MaxLPAlg: 3,
		UseLattice: true, QMin: 400000, QMax: 6000000, IHalf: 4096, JMax: 2048,
		SlackRat: 54, SlackAlg: 57, NumQCs: 64, RelTargetRatio: 0.78},
	"m5_105_d4_asym": {Degree: 4, Br: 800000, Ba: 1100000, LpBoundRat: 14000000, LpBoundAlg: 20000000,
		MfbRat: 48, MfbAlg: 66, MaxLPRat: 2, MaxLPAlg: 3,
		UseLattice: true, QMin: 350000, QMax: 5000000, IHalf: 4096, JMax: 2048,
		SlackRat: 49, SlackAlg: 54, NumQCs: 64, RelTargetRatio: 0.78},
	"m5_105_d5_asym": {Degree: 5, Br: 750000, Ba: 1100000, LpBoundRat: 14000000, LpBoundAlg: 20000000,
		MfbRat: 48, MfbAlg: 66, MaxLPRat: 2, MaxLPAlg: 3,
		UseLattice: true, QMin: 350000, QMax: 5000000, IHalf: 4096, JMax: 2048,
		SlackRat: 49, SlackAlg: 54, NumQCs: 64, RelTargetRatio: 0.78},
	"m5_110_d5_asym": {Degree: 5, Br: 900000, Ba: 1400000, LpBoundRat: 16000000, LpBoundAlg: 28000000,
		MfbRat: 48, MfbAlg: 68, MaxLPRat: 2, MaxLPAlg: 3,
		UseLattice: true, QMin: 400000, QMax: 6000000, IHalf: 4096, JMax: 2048,
		SlackRat: 49, SlackAlg: 56, NumQCs: 64, RelTargetRatio: 0.78},
	"m5_110_d5_asym3": {Degree: 5, Br: 1000000, Ba: 1500000, LpBoundRat: 20000000, LpBoundAlg: 32000000,
		MfbRat: 64, MfbAlg: 70, MaxLPRat: 3, MaxLPAlg: 3,
		UseLattice: true, QMin: 450000, QMax: 7000000, IHalf: 4096, JMax: 2048,
		SlackRat: 52, SlackAlg: 57, NumQCs: 64, RelTargetRatio: 0.78},
	"m5_110_d5_fast": {Degree: 5, Br: 1100000, Ba: 1600000, LpBoundRat: 16000000, LpBoundAlg: 28000000,
		MfbRat: 49, MfbAlg: 69, MaxLPRat: 2, MaxLPAlg: 3,
		UseLattice: true, QMin: 450000, QMax: 7000000, IHalf: 4096, JMax: 2048,
		SlackRat: 50, SlackAlg: 56, NumQCs: 64, RelTargetRatio: 0.78},
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
					}
					if strings.Contains(s, "relations (excess") {
						if time.Since(lastPrint) >= 15*time.Second {
							lastPrint = time.Now()
							fmt.Printf("  [%d bits seed %d %s] %s\n", n.BitLen(), seed, name, s)
						}
					} else {
						fmt.Printf("  [%d bits seed %d %s] %s\n", n.BitLen(), seed, name, s)
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

func TestYieldProbe(t *testing.T) {
	if !*tune || os.Getenv("GNFS_PROBE") == "" || os.Getenv("GNFS_BITS") == "" {
		t.Skip("run with -tune and GNFS_PROBE=1 GNFS_BITS=... GNFS_CANDS=...")
	}
	for _, bs := range strings.Split(os.Getenv("GNFS_BITS"), ",") {
		bits, err := strconv.Atoi(bs)
		if err != nil {
			t.Fatal(err)
		}
		rnd := rand.New(rand.NewSource(int64(1000 + bits)))
		n := new(big.Int).Mul(testutil.RandomPrime(rnd, bits/2), testutil.RandomPrime(rnd, bits-bits/2))
		for _, name := range strings.Split(os.Getenv("GNFS_CANDS"), ",") {
			var p Params
			if name == "default" {
				p = *defaultParams(n)
			} else {
				c, ok := tuneCands[name]
				if !ok {
					t.Fatalf("unknown cand %q", name)
				}
				p = c
			}
			t0 := time.Now()
			poly, fac := selectPolyWithParams(n, &p)
			if fac != nil || poly == nil {
				t.Fatalf("selectPoly failed")
			}
			polyDur := time.Since(t0).Round(time.Millisecond)
			qcMin := computeQCMin(&p)
			fb := BuildFactorBases(poly, p.Br, p.Ba, qcMin, -1)

			var sample []LatticeQ
			for _, off := range []uint64{10000, 400000, 1200000} {
				sub := p
				sub.QMin = p.QMin + off
				sub.QMax = sub.QMin + 4000
				sqs := buildSpecialQs(poly, fb, &sub)
				if len(sqs) > 30 {
					sqs = sqs[:30]
				}
				sample = append(sample, sqs...)
			}
			seen := make(map[[2]int64]struct{}, len(sample)*16)
			threeLP := 0
			ts := time.Now()
			sieveLatticeQs(poly, fb, &p, sample, func(r Relation) {
				key := [2]int64{r.A, r.B}
				if _, ok := seen[key]; !ok {
					seen[key] = struct{}{}
					if r.RatLPs[2] != 0 || len(r.AlgLPs) >= 3 {
						threeLP++
					}
				}
			}, func() bool { return false })
			sieveDur := time.Since(ts)
			relPerSec := float64(len(seen)) / sieveDur.Seconds()
			relPerQ := float64(len(seen)) / float64(len(sample))
			ideals := len(arith.PrimesUpTo(int64(p.LpBoundRat))) + len(arith.PrimesUpTo(int64(p.LpBoundAlg)))
			estMin := (0.80 * float64(ideals) / max(1.0, relPerSec)) / 60.0
			t.Logf("%d bits (%dd) %-15s deg=%d c_d=%v a_1=%v skew=%.0f alpha=%.2f poly=%v | %d qs -> %d rels (%.1f/q, %d 3LP+) in %v = %.0f rel/s | LP=%d est=%.1fm",
				n.BitLen(), len(n.String()), name, poly.Degree, poly.Coeffs[poly.Degree], poly.A1, poly.Skew, poly.Alpha, polyDur,
				len(sample), len(seen), relPerQ, threeLP, sieveDur.Round(time.Millisecond), relPerSec, ideals, estMin)
		}
	}
}

package main

import (
	"context"
	"flag"
	"fmt"
	"log"
	"math/big"
	"runtime"
	"strings"
	"time"

	"github.com/dgryski/go-factor"
	"github.com/dgryski/go-factor/gnfs"
	"github.com/dgryski/go-factor/qs"
)

// allMethods lists the factoring methods in the order they run.
var allMethods = []string{"trial", "brent", "p-1", "p+1", "squfof", "ecm", "qs", "gnfs"}

// Sieve bit limits when -methods is not explicitly set. QS hands off to GNFS
// above 342 bits (103 digits), near their measured 2LP crossover.
const (
	qsMaxBits   = 342
	gnfsMaxBits = 450
)

func main() {
	// 21251589199163473368424683523876245871651712859289

	verbose := flag.Bool("v", false, "log the progress of each method's stages")
	qsSave := flag.String("qs-save", "", "directory to save the quadratic sieve's relations in, as <n>.rel")
	gnfsSave := flag.String("gnfs-save", "", "directory to save GNFS relations in, as <n>.gnfs.rel")
	workers := flag.Int("workers", 0, "goroutines for ECM curves, quadratic sieving, and GNFS (0 means GOMAXPROCS)")
	methods := flag.String("methods", strings.Join(allMethods, ","), "comma-separated methods to run, always in the order "+strings.Join(allMethods, ","))
	flag.Parse()

	if *verbose {
		qs.Logf = log.Printf
		gnfs.Logf = log.Printf
		factor.Logf = log.Printf
	}
	qs.SaveDir = *qsSave
	gnfs.SaveDir = *gnfsSave
	qs.Workers = *workers
	gnfs.Workers = *workers
	factor.ECMWorkers = *workers

	run := make(map[string]bool)
	for _, m := range strings.Split(*methods, ",") {
		m = strings.TrimSpace(m)
		known := false
		for _, a := range allMethods {
			known = known || m == a
		}
		if !known {
			log.Fatalf("unknown method %q: want some of %s", m, strings.Join(allMethods, ","))
		}
		run[m] = true
	}
	methodsSet := false
	flag.Visit(func(f *flag.Flag) { methodsSet = methodsSet || f.Name == "methods" })

	var n big.Int

	log.Printf("factoring %+v\n", flag.Arg(0))

	if _, ok := n.SetString(flag.Arg(0), 10); !ok || n.Cmp(big.NewInt(1)) <= 0 {
		log.Fatalf("invalid number %q: want an integer greater than 1", flag.Arg(0))
	}

	if n.ProbablyPrime(20) {
		log.Printf("prime=%+v\n", &n)
		return
	}

	ctmp := []*big.Int{&n}

	if run["trial"] {
		log.Println("starting trial")
		var ptmp []*big.Int
		ptmp, ctmp = factor.Trial(&n)
		if len(ptmp) > 0 {
			log.Printf("trial=%+v\n", ptmp)
		}
	}

	qsRuns := func(n *big.Int) bool {
		return run["qs"] && (methodsSet || n.BitLen() <= qsMaxBits)
	}
	gnfsRuns := func(n *big.Int) bool {
		return run["gnfs"] && (methodsSet || n.BitLen() <= gnfsMaxBits)
	}

	qsWorkers := *workers
	if qsWorkers <= 0 {
		qsWorkers = runtime.GOMAXPROCS(0)
	}
	b := &budget{
		start:    time.Now(),
		workers:  qsWorkers,
		qsRuns:   qsRuns,
		gnfsRuns: gnfsRuns,
	}

	// each runs f on every composite left, until share of its budget is used.
	each := func(name string, share float64, f func(context.Context, *big.Int) ([]*big.Int, []*big.Int)) {
		log.Printf("starting %s", name)
		var comps []*big.Int
		for _, cc := range ctmp {
			ctx, cancel, sieveName, ok := b.context(cc, share)
			if !ok {
				log.Printf("skipping %s for %d digits: its share of the time before %s is used", name, len(cc.String()), sieveName)
				comps = append(comps, cc)
				continue
			}
			ptmp, c2 := f(ctx, cc)
			cancel()
			if len(ptmp) > 0 {
				log.Printf("%s=%+v\n", name, ptmp)
			}
			comps = append(comps, c2...)
		}
		ctmp = comps
	}

	if run["brent"] {
		for _, v := range []int64{1, -1, 3} {
			if len(ctmp) == 0 {
				break
			}
			each(fmt.Sprintf("brent(2,%v)", v), shareBrent, func(ctx context.Context, n *big.Int) ([]*big.Int, []*big.Int) {
				return factor.BrentContext(ctx, n, 2, v)
			})
		}
	}

	// GMP-ECM suggests p-1 with 10 times the B1 of the last ECM level.
	if run["p-1"] && len(ctmp) > 0 {
		each("p-1", sharePMinus, func(ctx context.Context, n *big.Int) ([]*big.Int, []*big.Int) {
			return factor.PMinus1Context(ctx, n, 1000000)
		})
	}

	// GMP-ECM suggests p+1 three times with 5 times the B1 of the last ECM
	// level.
	if run["p+1"] && len(ctmp) > 0 {
		each("p+1", sharePPlus, func(ctx context.Context, n *big.Int) ([]*big.Int, []*big.Int) {
			return factor.PPlus1Context(ctx, n, 250000, 3)
		})
	}

	if run["squfof"] && len(ctmp) > 0 {
		// Above 64 bits the quadratic sieve is faster, and from about 26
		// digits squfof runs out of iterations without finding anything.
		// Asking for squfof by name with -methods lifts the limit.
		const squfofMaxBits = 64
		multiplier := []int64{1, 3, 5, 7, 11, 3 * 5, 3 * 7, 3 * 11, 5 * 7, 5 * 11, 7 * 11, 3 * 5 * 7, 3 * 5 * 11, 3 * 7 * 11, 5 * 7 * 11, 3 * 5 * 7 * 11}
		log.Printf("starting squfof")
		for _, k := range multiplier {
			var comps []*big.Int
			for _, cc := range ctmp {
				if !methodsSet && cc.BitLen() > squfofMaxBits {
					if k == 1 {
						log.Printf("skipping squfof for %d-bit %v", cc.BitLen(), cc)
					}
					comps = append(comps, cc)
					continue
				}
				ptmp, c2 := factor.Squfof(cc, k)
				if len(ptmp) > 0 {
					log.Printf("squfof=%+v\n", ptmp)
				}
				comps = append(comps, c2...)
			}
			ctmp = comps
		}
	}

	if run["ecm"] && len(ctmp) > 0 {
		var comps []*big.Int
		for _, cc := range ctmp {
			ptmp, c2 := ecmPretest(cc, b)
			if len(ptmp) > 0 {
				log.Printf("ecm=%+v\n", ptmp)
			}
			comps = append(comps, c2...)
		}
		ctmp = comps
	}

	if run["qs"] && len(ctmp) > 0 {
		log.Println("starting qs")
		var comps []*big.Int
		for _, cc := range ctmp {
			if !qsRuns(cc) {
				log.Printf("skipping qs for %d-bit %v", cc.BitLen(), cc)
				comps = append(comps, cc)
				continue
			}
			ptmp, c2 := qs.QuadraticSieve(cc)
			if len(ptmp) > 0 {
				log.Printf("qs=%+v\n", ptmp)
			}
			comps = append(comps, c2...)
		}
		ctmp = comps
	}

	if run["gnfs"] && len(ctmp) > 0 {
		log.Println("starting gnfs")
		var comps []*big.Int
		for _, cc := range ctmp {
			if !gnfsRuns(cc) {
				log.Printf("skipping gnfs for %d-bit %v", cc.BitLen(), cc)
				comps = append(comps, cc)
				continue
			}
			ptmp, c2 := gnfs.Factor(cc)
			if len(ptmp) > 0 {
				log.Printf("gnfs=%+v\n", ptmp)
			}
			comps = append(comps, c2...)
		}
		ctmp = comps
	}

	if len(ctmp) > 0 {
		log.Printf("cfinal=%+v\n", ctmp)
	}

	runtime.GC()
}

package main

import (
	"flag"
	"log"
	"math/big"
	"runtime"
	"strings"

	"github.com/dgryski/go-factor"
	"github.com/dgryski/go-factor/qs"
)

// allMethods lists the factoring methods in the order they run.
var allMethods = []string{"trial", "brent", "p-1", "p+1", "squfof", "ecm", "qs"}

func main() {
	// 21251589199163473368424683523876245871651712859289

	verbose := flag.Bool("v", false, "log progress during the quadratic sieve")
	qsSave := flag.String("qs-save", "", "directory to save the quadratic sieve's relations in, as <n>.rel")
	workers := flag.Int("workers", 0, "goroutines sieving in the quadratic sieve (0 means GOMAXPROCS)")
	methods := flag.String("methods", strings.Join(allMethods, ","), "comma-separated methods to run, always in the order "+strings.Join(allMethods, ","))
	flag.Parse()

	if *verbose {
		qs.Logf = log.Printf
	}
	qs.SaveDir = *qsSave
	qs.Workers = *workers

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

	if run["brent"] && len(ctmp) > 0 {
		p := int64(2)
		for _, v := range []int64{1, -1, 3} {
			var comps []*big.Int
			for _, cc := range ctmp {
				log.Printf("starting brent(%v,%v)", p, v)
				ptmp, c2 := factor.Brent(cc, p, v)
				if len(ptmp) > 0 {
					log.Printf("brent(%v,%v)=%+v\n", p, v, ptmp)
				}
				comps = append(comps, c2...)
			}
			ctmp = comps
		}
	}

	if run["p-1"] && len(ctmp) > 0 {
		log.Println("starting p-1")
		var comps []*big.Int
		for _, cc := range ctmp {
			ptmp, c2 := factor.PMinus1(cc)
			if len(ptmp) > 0 {
				log.Printf("p-1=%+v\n", ptmp)
			}
			comps = append(comps, c2...)
		}
		ctmp = comps
	}

	if run["p+1"] && len(ctmp) > 0 {
		log.Println("starting p+1")
		var comps []*big.Int
		for _, cc := range ctmp {
			ptmp, c2 := factor.PPlus1(cc, 10000, 3)
			if len(ptmp) > 0 {
				log.Printf("p+1=%+v\n", ptmp)
			}
			comps = append(comps, c2...)
		}
		ctmp = comps
	}

	if run["squfof"] && len(ctmp) > 0 {
		multiplier := []int64{1, 3, 5, 7, 11, 3 * 5, 3 * 7, 3 * 11, 5 * 7, 5 * 11, 7 * 11, 3 * 5 * 7, 3 * 5 * 11, 3 * 7 * 11, 5 * 7 * 11, 3 * 5 * 7 * 11}
		log.Printf("starting squfof")
		for _, k := range multiplier {
			var comps []*big.Int
			for _, cc := range ctmp {
				ptmp, c2 := factor.Squfof(cc, k)
				if len(ptmp) > 0 {
					log.Printf("squfof=%+v\n", ptmp)
				}
				comps = append(comps, c2...)
			}
			ctmp = comps
		}
	}

	// B1 and curves for factors up to 15, 20 and 25 digits, from the GMP-ECM README.
	for _, s := range []struct{ b1, curves int }{{2000, 25}, {11000, 90}, {50000, 300}} {
		if !run["ecm"] || len(ctmp) == 0 {
			break
		}
		log.Printf("starting ecm(b1=%v, curves=%v)", s.b1, s.curves)
		var comps []*big.Int
		for _, cc := range ctmp {
			ptmp, c2 := factor.ECM(cc, s.b1, s.curves)
			if len(ptmp) > 0 {
				log.Printf("ecm=%+v\n", ptmp)
			}
			comps = append(comps, c2...)
		}
		ctmp = comps
	}

	if run["qs"] && len(ctmp) > 0 {
		// Above about 220 bits (66 digits) the quadratic sieve takes
		// minutes, and its dense matrix needs gigabytes. Asking for qs by
		// name with -methods lifts the limit.
		const qsMaxBits = 220
		log.Println("starting qs")
		var comps []*big.Int
		for _, cc := range ctmp {
			if !methodsSet && cc.BitLen() > qsMaxBits {
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

	if len(ctmp) > 0 {
		log.Printf("cfinal=%+v\n", ctmp)
	}

	runtime.GC()
}

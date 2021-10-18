package main

import (
	"flag"
	"log"
	"math/big"

	"github.com/dgryski/go-factor"
)

func main() {

	flag.Parse()

	var n big.Int

	log.Printf("factoring %+v\n", flag.Arg(0))

	n.SetString(flag.Arg(0), 10)

	log.Println("starting trial")
	ptmp, ctmp := factor.Trial(&n)

	if len(ptmp) > 0 {
		log.Printf("trial=%+v\n", ptmp)
	}

	if len(ctmp) > 0 {
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

	if len(ctmp) > 0 {
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

	if len(ctmp) > 0 {
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

	if len(ctmp) > 0 {
		log.Printf("cfinal=%+v\n", ctmp)
	}
}

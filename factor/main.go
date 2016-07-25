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

	var ptmp, ctmp []*big.Int

	log.Println("starting trial")
	ptmp, ctmp = factor.Trial(&n)

	if len(ptmp) > 0 {
		log.Printf("trial=%+v\n", ptmp)
	}

	var cfinal []*big.Int

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
		for _, cc := range ctmp {
			ptmp, c2 := factor.PMinus1(cc)
			if len(ptmp) > 0 {
				log.Printf("p-1=%+v\n", ptmp)
			}
			cfinal = append(cfinal, c2...)
		}
	}

	if len(cfinal) > 0 {
		log.Printf("cfinal=%+v\n", cfinal)
	}
}

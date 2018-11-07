// +build ignore

package main

import (
	"flag"
	"fmt"
	"math"
)

func main() {

	optmax := flag.Int("maxprime", 100000, "maximum prime to generate to")
	flag.Parse()

	maxprime := uint64(*optmax)
	sqrtmaxprime := uint64(math.Sqrt(float64(maxprime))) + 1

	bits := newbv(maxprime)

	for i := range bits {
		// even numbers are not prime
		bits[i] = 0xAAAAAAAAAAAAAAAA
	}

	bits.clear(0)
	bits.clear(1)
	bits.set(2)

	i := uint64(3)
	for i <= sqrtmaxprime {
		// clear all multiples of this prime */
		for j := 2 * i; j < maxprime; j += i {
			bits.clear(j)
		}

		i += 2
		for i < sqrtmaxprime && bits.get(i) == 0 {
			i += 2
		}
	}

	fmt.Println("package factor\nvar smallPrimes = []int64{")
	fmt.Println("\t2,")
	for i := uint64(3); i < maxprime; i += 2 {
		if bits.get(i) == 1 {
			fmt.Printf("\t%d,\n", i)
		}
	}
	fmt.Println("}")
}

type bitvector []uint64

func newbv(size uint64) bitvector {
	return make([]uint64, (size+63)/64)
}

// get bit 'bit' in the bitvector d
func (b bitvector) get(bit uint64) uint {
	shift := bit % 64
	bb := b[bit/64]
	bb &= (1 << shift)

	return uint(bb >> shift)
}

// set bit 'bit' in the bitvector d
func (b bitvector) set(bit uint64) {
	b[bit/64] |= (1 << (bit % 64))
}

// clear bit 'bit' in the bitvector d
func (b bitvector) clear(bit uint64) {
	b[bit/64] &^= (1 << (bit % 64))
}

func (b bitvector) size() uint64 {
	return uint64(len(b)) * 64
}

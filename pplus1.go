package factor

import (
	"math/big"
	"math/bits"
	"math/rand"
)

// PPlus1 factors n with Williams' p+1 method, trying up to tries starting values.
// See https://en.wikipedia.org/wiki/Williams%27s_p_%2B_1_algorithm
func PPlus1(n *big.Int, b1, tries int) (primes, composites []*big.Int) {
	e := &pplus1{
		rnd:    rand.New(rand.NewSource(1)),
		b1:     b1,
		primes: sievePrimes(pplus1B2Multiplier * b1),
		tries:  tries,
	}
	return e.factor(n)
}

const (
	pplus1B2Multiplier = 100
	pplus1Batch        = 16
)

type pplus1 struct {
	rnd    *rand.Rand
	b1     int
	primes []int
	tries  int // left to try, shared by recursive calls
}

func (e *pplus1) factor(n *big.Int) (primes, composites []*big.Int) {
	if n.Cmp(one) <= 0 {
		return nil, nil
	}
	if n.ProbablyPrime(10) {
		return []*big.Int{n}, nil
	}
	if n.Bit(0) == 0 {
		pr, co := e.factor(new(big.Int).Rsh(n, 1))
		return append([]*big.Int{big.NewInt(2)}, pr...), co
	}

	for e.tries > 0 {
		e.tries--
		g := e.try(n)
		if g == nil {
			continue
		}
		pr, co := e.factor(g)
		primes, composites = append(primes, pr...), append(composites, co...)
		pr, co = e.factor(new(big.Int).Div(n, g))
		primes, composites = append(primes, pr...), append(composites, co...)
		return primes, composites
	}

	return nil, []*big.Int{n}
}

// try runs both stages from a random start and returns a proper factor of n or nil.
func (e *pplus1) try(n *big.Int) *big.Int {
	var v, t, g big.Int
	v.Rand(e.rnd, n)

	// P^2 - 4 = 0 mod p makes the sequence degenerate mod p but exposes p.
	t.Mul(&v, &v)
	t.Sub(&t, big.NewInt(4))
	if g.GCD(nil, nil, t.Mod(&t, n), n); g.Cmp(one) != 0 {
		if g.Cmp(n) == 0 {
			return nil
		}
		return &g
	}

	// gcd(V - 2, n) for the current v.
	check := func(x *big.Int) *big.Int {
		t.Sub(x, two)
		t.Mod(&t, n)
		return g.GCD(nil, nil, &t, n)
	}

	// Stage 1 sets v = V_k(P) where k is the product of prime powers up to b1.
	saved := new(big.Int).Set(&v)
	start := 0
	for i, q := range e.primes {
		if q > e.b1 {
			break
		}
		k := uint64(q)
		for k*uint64(q) <= uint64(e.b1) {
			k *= uint64(q)
		}
		lucasV(&v, &v, k, n)

		last := i+1 == len(e.primes) || e.primes[i+1] > e.b1
		if (i+1)%pplus1Batch != 0 && !last {
			continue
		}
		if check(&v); g.Cmp(one) == 0 {
			saved.Set(&v)
			start = i + 1
			continue
		}
		if g.Cmp(n) == 0 {
			// Every prime of n was found in this batch. Redo it one prime at a time.
			v.Set(saved)
		replay:
			for _, q := range e.primes[start : i+1] {
				for k := uint64(q); k <= uint64(e.b1); k *= uint64(q) {
					lucasV(&v, &v, uint64(q), n)
					if check(&v); g.Cmp(one) != 0 {
						break replay
					}
				}
			}
			if g.Cmp(n) == 0 {
				return nil
			}
		}
		return new(big.Int).Set(&g)
	}

	// Stage 2 looks for one more prime q up to b2 using V_q(v) = V_kq(P).
	var w big.Int
	product := big.NewInt(1)
	batch := 0
	for i, q := range e.primes {
		if q <= e.b1 {
			continue
		}
		lucasV(&w, &v, uint64(q), n)
		t.Sub(&w, two)
		product.Mul(product, &t)
		product.Mod(product, n)
		batch++

		if batch < pplus1Batch && i+1 < len(e.primes) {
			continue
		}
		g.GCD(nil, nil, product, n)
		if g.Cmp(n) == 0 {
			for _, q := range e.primes[i+1-batch : i+1] {
				lucasV(&w, &v, uint64(q), n)
				if check(&w); g.Cmp(one) != 0 {
					break
				}
			}
		}
		if g.Cmp(n) == 0 {
			return nil
		}
		if g.Cmp(one) != 0 {
			return new(big.Int).Set(&g)
		}
		product.SetInt64(1)
		batch = 0
	}

	return nil
}

// lucasV sets r = V_k(p) mod n for the Lucas sequence with Q = 1, for k >= 1.
func lucasV(r, p *big.Int, k uint64, n *big.Int) {
	// x = V_i and y = V_i+1, starting from i = 1.
	x := new(big.Int).Set(p)
	y := new(big.Int).Mul(p, p)
	y.Sub(y, two)
	y.Mod(y, n)
	pp := new(big.Int).Set(p)

	var t big.Int
	for bit := 62 - bits.LeadingZeros64(k); bit >= 0; bit-- {
		t.Mul(x, y)
		t.Sub(&t, pp)
		t.Mod(&t, n)
		if k>>uint(bit)&1 == 1 {
			x.Set(&t)
			y.Mul(y, y)
			y.Sub(y, two)
			y.Mod(y, n)
		} else {
			y.Set(&t)
			x.Mul(x, x)
			x.Sub(x, two)
			x.Mod(x, n)
		}
	}
	r.Set(x)
}

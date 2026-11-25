package factor

import (
	"context"
	"math/big"
	"math/bits"
	"math/rand"
	"time"
)

// PPlus1 factors n with Williams' p+1 method, trying up to tries starting values.
// See https://en.wikipedia.org/wiki/Williams%27s_p_%2B_1_algorithm
func PPlus1(n *big.Int, b1, tries int) (primes, composites []*big.Int) {
	return PPlus1Context(context.Background(), n, b1, tries)
}

// PPlus1Context is PPlus1 that stops early when ctx is done, after a last
// gcd, and returns what it has found so far.
func PPlus1Context(ctx context.Context, n *big.Int, b1, tries int) (primes, composites []*big.Int) {
	e := &pplus1{
		stage2Plan: newStage2Plan(b1, pplus1B2Multiplier*b1),
		ctx:        ctx,
		rnd:        rand.New(rand.NewSource(1)),
		primes:     sievePrimes(b1),
		tries:      tries,
	}
	return e.factor(n)
}

const (
	pplus1B2Multiplier = 100

	// Stage 1 takes a gcd after this many primes, and stage 2 after this
	// many giant steps.
	pplus1Batch = 16
)

type pplus1 struct {
	*stage2Plan
	ctx    context.Context
	rnd    *rand.Rand
	primes []int // primes <= b1
	tries  int   // left to try, shared by recursive calls
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

	for e.tries > 0 && e.ctx.Err() == nil {
		e.tries--
		logf("p+1: try with B1=%d, B2=%d on %d bits, %d tries left", e.b1, e.b2, n.BitLen(), e.tries)
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
	t0 := time.Now()
	saved := new(big.Int).Set(&v)
	start := 0
	for i, q := range e.primes {
		k := uint64(q)
		for k*uint64(q) <= uint64(e.b1) {
			k *= uint64(q)
		}
		lucasV(&v, &v, k, n)

		if (i+1)%pplus1Batch != 0 && i+1 != len(e.primes) {
			continue
		}
		if check(&v); g.Cmp(one) == 0 {
			if e.ctx.Err() != nil {
				logf("p+1: stage 1 stopped after %d of %d primes in %v", i+1, len(e.primes), since(t0))
				return nil
			}
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
				logf("p+1: stage 1 found only n in %v", since(t0))
				return nil
			}
		}
		logf("p+1: stage 1 found %v in %v", &g, since(t0))
		return new(big.Int).Set(&g)
	}
	logf("p+1: stage 1 found nothing in %v", since(t0))

	return e.stage2(n, &v)
}

// stage2 looks for one more prime q in (b1, b2] with V_q(v) = 2 mod a prime
// of n. With v = a + 1/a,
//
//	V_c - V_m = (a^(c+m) - 1)(a^(c-m) - 1) / a^c,
//
// so the product of V_c - V_m over the baby steps m covers both c-m and c+m.
// Both sequences step with V_(j+d) = V_j V_d - V_(j-d).
func (e *pplus1) stage2(n, v *big.Int) *big.Int {
	start := time.Now()
	e.sieve()
	var t, q big.Int
	mulSub := func(z, x, y, w *big.Int) {
		// z = x*y - w mod n
		t.Mul(x, y)
		t.Sub(&t, w)
		q.QuoRem(&t, n, z)
		if z.Sign() < 0 {
			z.Add(z, n)
		}
	}

	// Baby steps: vm[m] = V_m(v) for odd m <= stage2D/2.
	vm := make([]big.Int, stage2D/2+1)
	var v2 big.Int
	mulSub(&v2, v, v, two)
	vm[1].Set(v)
	mulSub(&vm[3], v, &v2, v)
	for m := 5; m <= stage2D/2; m += 2 {
		mulSub(&vm[m], &vm[m-2], &v2, &vm[m-4])
	}

	// Giant steps: cur = V_(iv*D) and prev = V_((iv-1)*D), from V_ab = V_a(V_b).
	var vd, cur, prev, next big.Int
	lucasV(&vd, v, stage2D, n)
	iv := e.first()
	lucasV(&cur, &vd, uint64(iv), n)
	if iv == 1 {
		prev.Set(two)
	} else {
		lucasV(&prev, &vd, uint64(iv-1), n)
	}

	var g big.Int
	acc := big.NewInt(1)
	for steps := 1; ; steps++ {
		base := iv * stage2D
		for _, m := range e.babies {
			if !e.prime(base-m) && !e.prime(base+m) {
				continue
			}
			t.Sub(&cur, &vm[m])
			t.Mul(&t, acc)
			q.QuoRem(&t, n, acc)
			if acc.Sign() < 0 {
				acc.Add(acc, n)
			}
		}

		last := e.last(iv)
		if steps%pplus1Batch == 0 || last {
			g.GCD(nil, nil, acc, n)
			if g.Cmp(one) != 0 {
				if g.Cmp(n) == 0 {
					logf("p+1: stage 2 found only n in %v", since(start))
					return nil
				}
				logf("p+1: stage 2 found %v in %v", &g, since(start))
				return &g
			}
			if last {
				logf("p+1: stage 2 found nothing in %v", since(start))
				return nil
			}
			if e.ctx.Err() != nil {
				logf("p+1: stage 2 stopped at %d of %d in %v", base, e.b2, since(start))
				return nil
			}
		}

		mulSub(&next, &cur, &vd, &prev)
		prev.Set(&cur)
		cur.Set(&next)
		iv++
	}
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

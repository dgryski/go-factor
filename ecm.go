package factor

import (
	"math/big"
	"math/rand"
)

// ECM factors n with Lenstra's elliptic curve method, trying up to curves curves.
// See https://en.wikipedia.org/wiki/Lenstra_elliptic-curve_factorization
func ECM(n *big.Int, b1, curves int) (primes, composites []*big.Int) {
	if n.Cmp(one) <= 0 {
		return nil, nil
	}
	if n.ProbablyPrime(10) {
		return []*big.Int{n}, nil
	}

	e := &ecm{
		rnd:    rand.New(rand.NewSource(1)),
		b1:     b1,
		primes: sievePrimes(ecmB2Multiplier * b1),
		curves: curves,
	}
	return e.factor(n)
}

const (
	ecmB2Multiplier = 100

	// Stage 2 caches d*P for prime gaps d < 2*ecmCacheSize.
	ecmCacheSize = 128
)

type ecm struct {
	rnd    *rand.Rand
	b1     int
	primes []int
	curves int // left to try, shared by recursive calls
}

func (e *ecm) factor(n *big.Int) (primes, composites []*big.Int) {
	if n.Cmp(one) == 0 {
		return nil, nil
	}
	if n.ProbablyPrime(10) {
		return []*big.Int{n}, nil
	}

	for e.curves > 0 {
		e.curves--
		g := e.curve(n)
		if g == nil {
			continue
		}
		var q big.Int
		q.Div(n, g)
		pr, co := e.factor(g)
		primes, composites = append(primes, pr...), append(composites, co...)
		pr, co = e.factor(&q)
		primes, composites = append(primes, pr...), append(composites, co...)
		return primes, composites
	}

	return nil, []*big.Int{n}
}

// curve runs both stages on a random curve and returns a proper factor of n or nil.
func (e *ecm) curve(n *big.Int) *big.Int {
	c := &ecCurve{n: n}

	var p ecPoint
	p.x.Rand(e.rnd, n)
	p.y.Rand(e.rnd, n)
	c.a.Rand(e.rnd, n)

	// b = y^2 - x(x^2 + a)
	var b, t big.Int
	t.Mul(&p.x, &p.x)
	t.Add(&t, &c.a)
	t.Mul(&t, &p.x)
	b.Mul(&p.y, &p.y)
	b.Sub(&b, &t)
	b.Mod(&b, n)

	// disc = 4a^3 + 27b^2 is 0 for a singular curve but may share a factor with n.
	var disc big.Int
	disc.Mul(&c.a, &c.a)
	disc.Mul(&disc, &c.a)
	disc.Lsh(&disc, 2)
	t.Mul(&b, &b)
	t.Mul(&t, big.NewInt(27))
	disc.Add(&disc, &t)
	if g, ok := c.split(&disc); ok {
		return g
	} else if g != nil {
		return nil
	}

	// Stage 1 multiplies P by every prime power up to b1.
	for _, pr := range e.primes {
		if pr > e.b1 {
			break
		}
		k := uint64(pr)
		for k*uint64(pr) <= uint64(e.b1) {
			k *= uint64(pr)
		}
		if g := c.mul(&p, &p, k); g != nil {
			return c.proper(g)
		}
		if p.inf {
			// P has order dividing k mod every prime of n.
			return nil
		}
	}

	// Stage 2 steps Q = q*P over primes q above b1. When Q = O mod p the next
	// Q + d*P has no inverse mod p, which exposes p.
	start := 0
	for start < len(e.primes) && e.primes[start] <= e.b1 {
		start++
	}
	if start == len(e.primes) {
		return nil
	}

	var q ecPoint
	if g := c.mul(&q, &p, uint64(e.primes[start])); g != nil {
		return c.proper(g)
	}

	var cache [ecmCacheSize]*ecPoint
	var step ecPoint
	for i := start + 1; i < len(e.primes); i++ {
		if q.inf {
			return nil
		}
		d := e.primes[i] - e.primes[i-1]
		var dp *ecPoint
		if d/2 < ecmCacheSize {
			dp = cache[d/2]
			if dp == nil {
				dp = new(ecPoint)
				if g := c.mul(dp, &p, uint64(d)); g != nil {
					return c.proper(g)
				}
				cache[d/2] = dp
			}
		} else {
			if g := c.mul(&step, &p, uint64(d)); g != nil {
				return c.proper(g)
			}
			dp = &step
		}
		if g := c.add(&q, &q, dp); g != nil {
			return c.proper(g)
		}
	}

	return nil
}

type ecPoint struct {
	x, y big.Int
	inf  bool // the point at infinity
}

func (p *ecPoint) set(q *ecPoint) {
	p.x.Set(&q.x)
	p.y.Set(&q.y)
	p.inf = q.inf
}

type ecCurve struct {
	n *big.Int
	a big.Int

	l, t, u, inv big.Int
	r, s         ecPoint
}

// split returns (g, true) for a proper factor g = gcd(d, n). Otherwise g is
// nil if the gcd is 1 and n if it is n.
func (c *ecCurve) split(d *big.Int) (*big.Int, bool) {
	var dd, g big.Int
	dd.Mod(d, c.n)
	g.GCD(nil, nil, &dd, c.n)
	if g.Cmp(one) == 0 {
		return nil, false
	}
	if g.Cmp(c.n) == 0 {
		return c.n, false
	}
	return &g, true
}

// invert sets c.inv to 1/d mod n. If d is not invertible it returns gcd(d, n).
func (c *ecCurve) invert(d *big.Int) *big.Int {
	c.u.Mod(d, c.n)
	if c.inv.ModInverse(&c.u, c.n) != nil {
		return nil
	}
	g, _ := c.split(&c.u)
	return g
}

// add sets r = p + q. If a denominator is not invertible it returns its gcd
// with n, which is n itself if the curve failed mod every prime at once.
func (c *ecCurve) add(r, p, q *ecPoint) *big.Int {
	switch {
	case p.inf:
		r.set(q)
		return nil
	case q.inf:
		r.set(p)
		return nil
	}

	n := c.n
	if p.x.Cmp(&q.x) == 0 {
		c.t.Add(&p.y, &q.y)
		c.t.Mod(&c.t, n)
		if c.t.Sign() == 0 {
			// P + (-P) = O
			r.inf = true
			return nil
		}
		// l = (3x^2 + a) / 2y
		c.t.Lsh(&p.y, 1)
		if g := c.invert(&c.t); g != nil {
			return g
		}
		c.l.Mul(&p.x, &p.x)
		c.l.Mul(&c.l, big.NewInt(3))
		c.l.Add(&c.l, &c.a)
	} else {
		// l = (y2 - y1) / (x2 - x1)
		c.t.Sub(&q.x, &p.x)
		if g := c.invert(&c.t); g != nil {
			return g
		}
		c.l.Sub(&q.y, &p.y)
	}
	c.l.Mul(&c.l, &c.inv)
	c.l.Mod(&c.l, n)

	// x3 = l^2 - x1 - x2
	c.t.Mul(&c.l, &c.l)
	c.t.Sub(&c.t, &p.x)
	c.t.Sub(&c.t, &q.x)
	c.t.Mod(&c.t, n)

	// y3 = l(x1 - x3) - y1
	c.u.Sub(&p.x, &c.t)
	c.u.Mul(&c.u, &c.l)
	c.u.Sub(&c.u, &p.y)
	c.u.Mod(&c.u, n)

	r.x.Set(&c.t)
	r.y.Set(&c.u)
	r.inf = false
	return nil
}

// mul sets r = k*p. It returns the same values as add and leaves r unchanged
// when it fails.
func (c *ecCurve) mul(r, p *ecPoint, k uint64) *big.Int {
	acc, base := &c.r, &c.s
	acc.inf = true
	base.set(p)
	for k > 0 {
		if k&1 == 1 {
			if g := c.add(acc, acc, base); g != nil {
				return g
			}
		}
		k >>= 1
		if k > 0 {
			if g := c.add(base, base, base); g != nil {
				return g
			}
		}
	}
	r.set(acc)
	return nil
}

// proper returns g, or nil if g is n and the curve has no factor to give.
func (c *ecCurve) proper(g *big.Int) *big.Int {
	if g == c.n {
		return nil
	}
	return g
}

// sievePrimes returns the primes <= max.
func sievePrimes(max int) []int {
	if max < 2 {
		return nil
	}
	composite := make([]bool, max+1)
	primes := []int{2}
	for i := 3; i <= max; i += 2 {
		if composite[i] {
			continue
		}
		primes = append(primes, i)
		for j := i * i; j <= max; j += 2 * i {
			composite[j] = true
		}
	}
	return primes
}

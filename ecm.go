package factor

import (
	"context"
	"math/big"
	"math/rand"
)

// ECM factors n with Lenstra's elliptic curve method, trying up to curves
// curves with stage 2 bound 100*b1.
// See https://en.wikipedia.org/wiki/Lenstra_elliptic-curve_factorization
//
// It uses Montgomery curves with torsion Z/12 and x-only arithmetic, so stage
// 1 needs no inversions, and a baby-step giant-step stage 2 as in MIRACL's
// factor.c.
func ECM(n *big.Int, b1, curves int) (primes, composites []*big.Int) {
	return ECMContext(context.Background(), n, b1, ecmB2Multiplier*b1, curves)
}

// ECMContext is ECM with stage 2 bound b2. It stops early when ctx is done,
// abandoning the curve in progress, and returns what it has found so far.
func ECMContext(ctx context.Context, n *big.Int, b1, b2, curves int) (primes, composites []*big.Int) {
	if n.Cmp(one) <= 0 {
		return nil, nil
	}
	if n.ProbablyPrime(10) {
		return []*big.Int{n}, nil
	}

	e := &ecm{
		ctx:        ctx,
		rnd:        rand.New(rand.NewSource(1)),
		stage2Plan: newStage2Plan(b1, b2),
		primes:     sievePrimes(b1),
		curves:     curves,
	}
	return e.factor(n)
}

const (
	ecmB2Multiplier = 100

	// Stage 1 takes a gcd after this many primes.
	ecmStage1Batch = 50

	// Stage 2 takes a gcd after this many giant steps.
	ecmStage2Batch = 16
)

type ecm struct {
	*stage2Plan
	ctx    context.Context
	rnd    *rand.Rand
	primes []int // primes <= b1
	curves int   // left to try, shared by recursive calls
}

func (e *ecm) factor(n *big.Int) (primes, composites []*big.Int) {
	if n.Cmp(one) == 0 {
		return nil, nil
	}
	if n.ProbablyPrime(10) {
		return []*big.Int{n}, nil
	}

	for e.curves > 0 && e.ctx.Err() == nil {
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
	var p xzPoint
	if g := c.torsion12(&p, 2+e.rnd.Int63n(1<<32)); g != nil {
		return c.proper(g)
	}

	// Stage 1 multiplies P by every prime power up to b1. A prime of n shows
	// up in gcd(Z, n) once P = O mod that prime.
	saved := new(xzPoint)
	saved.set(&p)
	start := 0
	for i, q := range e.primes {
		k := uint64(q)
		for k*uint64(q) <= uint64(e.b1) {
			k *= uint64(q)
		}
		c.mul(&p, &p, k)

		if (i+1)%ecmStage1Batch != 0 && i+1 != len(e.primes) {
			continue
		}
		g, _ := c.split(&p.z)
		if g == nil {
			if e.ctx.Err() != nil {
				return nil
			}
			saved.set(&p)
			start = i + 1
			continue
		}
		if g != n {
			return g
		}
		// Every prime of n was found in this batch. Redo it one prime at a time.
		p.set(saved)
		for _, q := range e.primes[start : i+1] {
			for k := uint64(q); k <= uint64(e.b1); k *= uint64(q) {
				c.mul(&p, &p, uint64(q))
				if g, _ := c.split(&p.z); g != nil {
					return c.proper(g)
				}
			}
		}
		return nil
	}

	return e.stage2(c, &p)
}

// stage2 looks for one more prime q in (b1, b2] with q*P = O mod a prime of n.
// It writes q = c±m with c a multiple of stage2D and m a baby step. x(c*P) and
// x(m*P) agree mod p exactly when (c-m)*P or (c+m)*P is O mod p, so one
// product term covers both.
func (e *ecm) stage2(c *ecCurve, p *xzPoint) *big.Int {
	n := c.n
	e.sieve()

	// Baby steps: fu[m] = x(m*P) for m in e.babies.
	fu := make([]big.Int, stage2D/2+1)
	if g := c.normalise(&fu[1], p); g != nil {
		return c.proper(g)
	}
	var p2, prev, cur, next xzPoint
	c.double(&p2, p)
	prev.set(p)
	c.add(&cur, &p2, p, p) // 3P
	for m := 3; m <= stage2D/2; m += 2 {
		if gcdInt(m, stage2D) == 1 {
			if g := c.normalise(&fu[m], &cur); g != nil {
				return c.proper(g)
			}
		}
		// (m+2)P = mP + 2P, with difference (m-2)P.
		c.add(&next, &cur, &p2, &prev)
		prev.set(&cur)
		cur.set(&next)
	}

	// Giant steps: cur = iv*G and next = (iv+1)*G for G = stage2D*P, starting at
	// the multiple of stage2D nearest b1.
	var giant xzPoint
	c.mul(&giant, p, stage2D)
	iv := e.first()
	c.ladder(&cur, &next, &giant, uint64(iv))

	var x, t, g big.Int
	acc := big.NewInt(1)
	for steps := 1; ; steps++ {
		base := iv * stage2D
		if gg := c.normalise(&x, &cur); gg != nil {
			return c.proper(gg)
		}
		for _, m := range e.babies {
			if !e.prime(base-m) && !e.prime(base+m) {
				continue
			}
			c.sub(&t, &x, &fu[m])
			c.mulMod(acc, acc, &t)
		}

		last := e.last(iv)
		if steps%ecmStage2Batch == 0 || last {
			g.GCD(nil, nil, acc, n)
			if g.Cmp(one) != 0 {
				if g.Cmp(n) == 0 {
					return nil
				}
				return &g
			}
		}
		if last || steps%ecmStage2Batch == 0 && e.ctx.Err() != nil {
			return nil
		}

		// (iv+2)G = (iv+1)G + G, with difference iv*G.
		c.add(&prev, &next, &giant, &cur)
		cur.set(&next)
		next.set(&prev)
		iv++
	}
}

// xzPoint is a point (X:Z) on a Montgomery curve with y dropped.
// Z = 0 is the point at infinity.
type xzPoint struct {
	x, z big.Int
}

func (p *xzPoint) set(q *xzPoint) {
	p.x.Set(&q.x)
	p.z.Set(&q.z)
}

// ecCurve is the Montgomery curve By^2 = x^3 + Ax^2 + x mod n, stored as
// a24 = (A+2)/4.
type ecCurve struct {
	n   *big.Int
	a24 big.Int

	q, t           big.Int // mulMod scratch
	s1, d1, s2, d2 big.Int
	u, w           big.Int
	l1, l2         xzPoint // ladder scratch
}

// torsion12 sets p to the starting point of Montgomery's curve with torsion
// Z/12 for k >= 2. With (u, v) = k*(-2, 4) on v^2 = u^3 - 12u it uses
//
//	t^2 = (u^2-12)/4u, a = (t^2-1)/(t^2+3),
//	A = (-3a^4-6a^2+1)/4a^3, x = (3a^2+1)/4a.
//
// Compared with Suyama's curves the group order has more factors of 2 on
// average. If a denominator is not invertible it returns its gcd with n.
func (c *ecCurve) torsion12(p *xzPoint, k int64) *big.Int {
	n := c.n
	var w xyPoint
	if g := c.w12Mul(&w, k); g != nil {
		return g
	}

	var t, d, inv, a, num big.Int
	reduce := func(z *big.Int) *big.Int { return z.Mod(z, n) }

	// t^2 = (u^2-12)/4u
	t.Mul(&w.x, &w.x)
	reduce(t.Sub(&t, big.NewInt(12)))
	reduce(d.Lsh(&w.x, 2))
	if g := c.inverse(&inv, &d); g != nil {
		return g
	}
	c.mulMod(&t, &t, &inv)

	// a = (t^2-1)/(t^2+3)
	reduce(d.Add(&t, big.NewInt(3)))
	if g := c.inverse(&inv, &d); g != nil {
		return g
	}
	reduce(a.Sub(&t, one))
	c.mulMod(&a, &a, &inv)

	// x = (3a^2+1)/4a, z = 1
	reduce(d.Lsh(&a, 2))
	if g := c.inverse(&inv, &d); g != nil {
		return g
	}
	c.mulMod(&num, &a, &a)
	reduce(num.Add(num.Mul(&num, big.NewInt(3)), one))
	c.mulMod(&p.x, &num, &inv)
	p.z.SetInt64(1)

	// a24 = (A+2)/4 = -(a-1)^3 (3a+1) / 16a^3
	c.mulMod(&d, &a, &a)
	c.mulMod(&d, &d, &a)
	reduce(d.Lsh(&d, 4))
	if g := c.inverse(&inv, &d); g != nil {
		return g
	}
	c.sub(&t, &a, one)
	c.mulMod(&num, &t, &t)
	c.mulMod(&num, &num, &t)
	reduce(t.Add(t.Mul(&a, big.NewInt(3)), one))
	c.mulMod(&num, &num, &t)
	c.mulMod(&num, &num, &inv)
	c.sub(&c.a24, big.NewInt(0), &num)
	return nil
}

// xyPoint is an affine point on v^2 = u^3 - 12u, used to pick a curve.
type xyPoint struct {
	x, y big.Int
	inf  bool
}

// w12Mul sets r = k*(-2, 4) on v^2 = u^3 - 12u for k >= 1. If a
// denominator is not invertible it returns its gcd with n, and it returns n
// if r is the point at infinity.
func (c *ecCurve) w12Mul(r *xyPoint, k int64) *big.Int {
	var base xyPoint
	base.x.Sub(c.n, big.NewInt(2))
	base.y.SetInt64(4)
	r.x.Set(&base.x)
	r.y.Set(&base.y)
	top := 62
	for k>>uint(top)&1 == 0 {
		top--
	}
	for bit := top - 1; bit >= 0; bit-- {
		if g := c.w12Add(r, r, r); g != nil {
			return g
		}
		if k>>uint(bit)&1 == 1 {
			if g := c.w12Add(r, r, &base); g != nil {
				return g
			}
		}
	}
	if r.inf {
		return c.n
	}
	return nil
}

// w12Add sets r = p+q on v^2 = u^3 - 12u. r may alias p or q.
func (c *ecCurve) w12Add(r, p, q *xyPoint) *big.Int {
	switch {
	case p.inf:
		p, q = q, p
		fallthrough
	case q.inf:
		r.x.Set(&p.x)
		r.y.Set(&p.y)
		r.inf = p.inf
		return nil
	}
	var l, d, inv, x3 big.Int
	if p.x.Cmp(&q.x) == 0 {
		c.addMod(&d, &p.y, &q.y)
		if d.Sign() == 0 {
			r.inf = true
			return nil
		}
		// l = (3u^2 - 12) / 2v
		c.mulMod(&l, &p.x, &p.x)
		l.Mul(&l, big.NewInt(3))
		l.Sub(&l, big.NewInt(12))
		l.Mod(&l, c.n)
		c.addMod(&d, &p.y, &p.y)
	} else {
		// l = (v2 - v1) / (u2 - u1)
		c.sub(&l, &q.y, &p.y)
		c.sub(&d, &q.x, &p.x)
	}
	if g := c.inverse(&inv, &d); g != nil {
		return g
	}
	c.mulMod(&l, &l, &inv)

	// u3 = l^2 - u1 - u2, v3 = l(u1 - u3) - v1
	c.mulMod(&x3, &l, &l)
	c.sub(&x3, &x3, &p.x)
	c.sub(&x3, &x3, &q.x)
	c.sub(&d, &p.x, &x3)
	c.mulMod(&d, &d, &l)
	c.sub(&r.y, &d, &p.y)
	r.x.Set(&x3)
	r.inf = false
	return nil
}

// inverse sets z = 1/x mod n for x in [0, n). If x is not invertible it
// returns gcd(x, n).
func (c *ecCurve) inverse(z, x *big.Int) *big.Int {
	if z.ModInverse(x, c.n) != nil {
		return nil
	}
	if g, _ := c.split(x); g != nil {
		return g
	}
	return c.n // x = 0
}

// mulMod sets z = x*y mod n for x, y in [0, n).
func (c *ecCurve) mulMod(z, x, y *big.Int) {
	c.t.Mul(x, y)
	c.q.QuoRem(&c.t, c.n, z)
}

// addMod sets z = x+y mod n for x, y in [0, n).
func (c *ecCurve) addMod(z, x, y *big.Int) {
	z.Add(x, y)
	if z.Cmp(c.n) >= 0 {
		z.Sub(z, c.n)
	}
}

// sub sets z = x-y mod n for x, y in [0, n).
func (c *ecCurve) sub(z, x, y *big.Int) {
	z.Sub(x, y)
	if z.Sign() < 0 {
		z.Add(z, c.n)
	}
}

// double sets r = 2p. r may alias p.
func (c *ecCurve) double(r, p *xzPoint) {
	c.addMod(&c.s1, &p.x, &p.z)
	c.sub(&c.d1, &p.x, &p.z)
	c.mulMod(&c.s1, &c.s1, &c.s1) // (x+z)^2
	c.mulMod(&c.d1, &c.d1, &c.d1) // (x-z)^2
	c.sub(&c.u, &c.s1, &c.d1)     // 4xz
	c.mulMod(&r.x, &c.s1, &c.d1)
	c.mulMod(&c.w, &c.a24, &c.u)
	c.addMod(&c.w, &c.w, &c.d1)
	c.mulMod(&r.z, &c.w, &c.u)
}

// add sets r = p+q given diff = p-q. r may alias p or q but not diff.
func (c *ecCurve) add(r, p, q, diff *xzPoint) {
	c.addMod(&c.s1, &p.x, &p.z)
	c.sub(&c.d1, &p.x, &p.z)
	c.addMod(&c.s2, &q.x, &q.z)
	c.sub(&c.d2, &q.x, &q.z)
	c.mulMod(&c.u, &c.d1, &c.s2)
	c.mulMod(&c.w, &c.s1, &c.d2)
	c.addMod(&c.s1, &c.u, &c.w)
	c.sub(&c.d1, &c.u, &c.w)
	c.mulMod(&c.s1, &c.s1, &c.s1)
	c.mulMod(&c.d1, &c.d1, &c.d1)
	c.mulMod(&r.x, &c.s1, &diff.z)
	c.mulMod(&r.z, &c.d1, &diff.x)
}

// ladder sets r1 = k*p and r2 = (k+1)*p for k >= 1 with the Montgomery
// ladder. p must not alias r1 or r2.
func (c *ecCurve) ladder(r1, r2, p *xzPoint, k uint64) {
	r1.set(p)
	c.double(r2, p)
	top := 63
	for k>>uint(top)&1 == 0 {
		top--
	}
	for bit := top - 1; bit >= 0; bit-- {
		if k>>uint(bit)&1 == 1 {
			c.add(r1, r1, r2, p)
			c.double(r2, r2)
		} else {
			c.add(r2, r1, r2, p)
			c.double(r1, r1)
		}
	}
}

// mul sets r = k*p for k >= 1. r may alias p.
func (c *ecCurve) mul(r, p *xzPoint, k uint64) {
	c.ladder(&c.l1, &c.l2, p, k)
	r.set(&c.l1)
}

// normalise sets x = p.x/p.z. If p.z is not invertible it returns its gcd
// with n.
func (c *ecCurve) normalise(x *big.Int, p *xzPoint) *big.Int {
	var inv big.Int
	if g := c.inverse(&inv, &p.z); g != nil {
		return g
	}
	c.mulMod(x, &p.x, &inv)
	return nil
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

// proper returns g, or nil if g is n and the curve has no factor to give.
func (c *ecCurve) proper(g *big.Int) *big.Int {
	if g == c.n {
		return nil
	}
	return g
}

func gcdInt(a, b int) int {
	for b != 0 {
		a, b = b, a%b
	}
	return a
}

// primeBits returns a bitset with bit p set for each prime p <= max.
func primeBits(max int) []uint64 {
	bits := make([]uint64, max/64+1)
	if max < 2 {
		return bits
	}
	// Mark the odd numbers >= 3, then clear the odd composites.
	for i := range bits {
		bits[i] = 0xAAAAAAAAAAAAAAAA
	}
	bits[0] = bits[0]&^(1<<1) | 1<<2
	bits[len(bits)-1] &= 1<<(max%64+1) - 1
	for i := 3; i*i <= max; i += 2 {
		if bits[i/64]&(1<<(i%64)) == 0 {
			continue
		}
		for j := i * i; j <= max; j += 2 * i {
			bits[j/64] &^= 1 << (j % 64)
		}
	}
	return bits
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

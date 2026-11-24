package qs

import "math/bits"

// mont does arithmetic mod an odd n in Montgomery form, x·2⁶⁴ mod n, which
// needs no division. bits.Div64 is software on arm64.
type mont struct {
	n, nInv uint64 // nInv = −n⁻¹ mod 2⁶⁴
	one, r2 uint64 // 2⁶⁴ and 2¹²⁸ mod n
}

func newMont(n uint64) mont {
	inv := n // correct to 3 bits for odd n, and each step doubles that
	for i := 0; i < 5; i++ {
		inv *= 2 - n*inv
	}
	one := -n % n
	_, r2 := bits.Div64(one, 0, n)
	return mont{n: n, nInv: -inv, one: one, r2: r2}
}

// mul returns a·b·2⁻⁶⁴ mod n, for a and b less than n.
func (m mont) mul(a, b uint64) uint64 {
	hi, lo := bits.Mul64(a, b)
	mh, ml := bits.Mul64(lo*m.nInv, m.n)
	_, carry := bits.Add64(lo, ml, 0)
	r, c := bits.Add64(hi, mh, carry)
	if c != 0 || r >= m.n {
		r -= m.n
	}
	return r
}

func (m mont) to(a uint64) uint64 { return m.mul(a%m.n, m.r2) }

func (m mont) pow(base, exp uint64) uint64 {
	res := m.one
	for ; exp > 0; exp >>= 1 {
		if exp&1 != 0 {
			res = m.mul(res, base)
		}
		base = m.mul(base, base)
	}
	return res
}

func gcd64(a, b uint64) uint64 {
	for b != 0 {
		a, b = b, a%b
	}
	return a
}

// isPrime64 is Miller–Rabin with bases that make it exact for every uint64,
// from https://miller-rabin.appspot.com (Jim Sinclair, 2011).
func isPrime64(n uint64) bool {
	if n < 2 {
		return false
	}
	for _, p := range [...]uint64{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37} {
		if n%p == 0 {
			return n == p
		}
	}
	m := newMont(n)
	minusOne := n - m.one
	d := n - 1
	s := bits.TrailingZeros64(d)
	d >>= uint(s)
	for _, a := range [...]uint64{2, 325, 9375, 28178, 450775, 9780504, 1795265022} {
		if a%n == 0 {
			continue
		}
		x := m.pow(m.to(a), d)
		if x == m.one || x == minusOne {
			continue
		}
		composite := true
		for r := 1; r < s; r++ {
			x = m.mul(x, x)
			if x == minusOne {
				composite = false
				break
			}
		}
		if composite {
			return false
		}
	}
	return true
}

// rho64 returns a non-trivial factor of the odd composite n, found with
// Brent's variant of Pollard's rho, or 0 if it gives up. It works in
// Montgomery form, which doesn't change any gcd since 2⁶⁴ is coprime to n.
func rho64(n uint64) uint64 {
	m := newMont(n)
	for c := uint64(1); c < 16; c++ {
		f := func(v uint64) uint64 {
			v = m.mul(v, v) + c
			if v >= n {
				v -= n
			}
			return v
		}
		y, r, q, g := uint64(2), uint64(1), uint64(1), uint64(1)
		var x, ys uint64
		for g == 1 {
			x = y
			for i := uint64(0); i < r; i++ {
				y = f(y)
			}
			// Multiply up to 64 differences before each gcd.
			for k := uint64(0); k < r && g == 1; k += 64 {
				ys = y
				for i := uint64(0); i < 64 && k+i < r; i++ {
					y = f(y)
					q = m.mul(q, absDiff(x, y))
				}
				g = gcd64(q, n)
			}
			r <<= 1
		}
		if g == n {
			// The batch overshot, so retry its steps one at a time.
			for g = 1; g == 1; {
				ys = f(ys)
				g = gcd64(absDiff(x, ys), n)
			}
		}
		if g > 1 && g < n {
			return g
		}
	}
	return 0
}

func absDiff(a, b uint64) uint64 {
	if a > b {
		return a - b
	}
	return b - a
}

// splitCofactor returns the two primes, each at most lpBound, whose product
// is c. Every prime up to pMax must already be divided out of c.
func splitCofactor(c, pMax, lpBound uint64) (lp1, lp2 uint64, ok bool) {
	if c < pMax*pMax || isPrime64(c) {
		return 0, 0, false // a prime too large to keep
	}
	d := rho64(c)
	if d == 0 {
		return 0, 0, false
	}
	e := c / d
	if d > e {
		d, e = e, d
	}
	// Both factors are above pMax, so each is prime if below pMax².
	if e > lpBound || (d >= pMax*pMax && !isPrime64(d)) || (e >= pMax*pMax && !isPrime64(e)) {
		return 0, 0, false
	}
	return d, e, true
}

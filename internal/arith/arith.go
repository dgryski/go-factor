// Package arith provides 64-bit modular arithmetic, primality testing,
// cofactor splitting, and polynomial arithmetic over finite fields.
package arith

import (
	"math"
	"math/big"
	"math/bits"
)

// Mont does arithmetic mod an odd n in Montgomery form, x·2⁶⁴ mod n, which
// needs no division. bits.Div64 is software on arm64.
type Mont struct {
	N, nInv uint64 // nInv = −n⁻¹ mod 2⁶⁴
	One, r2 uint64 // 2⁶⁴ and 2¹²⁸ mod n
}

// NewMont returns a Montgomery context for the odd modulus n.
func NewMont(n uint64) Mont {
	inv := n // correct to 3 bits for odd n, and each step doubles that
	for i := 0; i < 5; i++ {
		inv *= 2 - n*inv
	}
	one := -n % n
	_, r2 := bits.Div64(one, 0, n)
	return Mont{N: n, nInv: -inv, One: one, r2: r2}
}

// Mul returns a·b·2⁻⁶⁴ mod n, for a and b less than n.
func (m Mont) Mul(a, b uint64) uint64 {
	hi, lo := bits.Mul64(a, b)
	mh, ml := bits.Mul64(lo*m.nInv, m.N)
	_, carry := bits.Add64(lo, ml, 0)
	r, c := bits.Add64(hi, mh, carry)
	if c != 0 || r >= m.N {
		r -= m.N
	}
	return r
}

// To converts a into Montgomery form.
func (m Mont) To(a uint64) uint64 { return m.Mul(a%m.N, m.r2) }

// Pow computes base^exp mod n in Montgomery form.
func (m Mont) Pow(base, exp uint64) uint64 {
	res := m.One
	for ; exp > 0; exp >>= 1 {
		if exp&1 != 0 {
			res = m.Mul(res, base)
		}
		base = m.Mul(base, base)
	}
	return res
}

// Gcd64 returns the greatest common divisor of a and b.
func Gcd64(a, b uint64) uint64 {
	for b != 0 {
		a, b = b, a%b
	}
	return a
}

// IsPrime64 is Miller–Rabin with bases that make it exact for every uint64,
// from https://miller-rabin.appspot.com (Jim Sinclair, 2011).
func IsPrime64(n uint64) bool {
	if n < 2 {
		return false
	}
	for _, p := range [...]uint64{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37} {
		if n%p == 0 {
			return n == p
		}
	}
	m := NewMont(n)
	minusOne := n - m.One
	d := n - 1
	s := bits.TrailingZeros64(d)
	d >>= uint(s)
	for _, a := range [...]uint64{2, 325, 9375, 28178, 450775, 9780504, 1795265022} {
		if a%n == 0 {
			continue
		}
		x := m.Pow(m.To(a), d)
		if x == m.One || x == minusOne {
			continue
		}
		composite := true
		for r := 1; r < s; r++ {
			x = m.Mul(x, x)
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

// Rho64 returns a non-trivial factor of the composite n, found with Brent's
// variant of Pollard's rho, or 0 if it gives up. It works in Montgomery
// form, which doesn't change any gcd since 2⁶⁴ is coprime to n. Each
// constant c runs at most about 2¹⁹ steps, so a prime n can't make it spin.
func Rho64(n uint64) uint64 {
	if n <= 1 {
		return 0
	}
	if n&1 == 0 {
		return 2 // Montgomery form needs n odd
	}
	m := NewMont(n)
	// c < n keeps f(v) below 2n, so one subtraction reduces it.
	for c := uint64(1); c < 16 && c < n; c++ {
		f := func(v uint64) uint64 {
			v = m.Mul(v, v) + c
			if v >= n {
				v -= n
			}
			return v
		}
		y, r, q, g := uint64(2), uint64(1), uint64(1), uint64(1)
		var x, ys uint64
		for g == 1 && r <= 1<<18 {
			x = y
			for i := uint64(0); i < r; i++ {
				y = f(y)
			}
			// Multiply up to 64 differences before each gcd.
			for k := uint64(0); k < r && g == 1; k += 64 {
				ys = y
				for i := uint64(0); i < 64 && k+i < r; i++ {
					y = f(y)
					q = m.Mul(q, absDiff(x, y))
				}
				g = Gcd64(q, n)
			}
			r <<= 1
		}
		if g == n {
			// The batch overshot, so retry its steps one at a time.
			for {
				ys = f(ys)
				if g = Gcd64(absDiff(x, ys), n); g > 1 {
					break
				}
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

// SplitCofactor returns the two primes, each at most lpBound, whose product
// is c. Every prime up to pMax must already be divided out of c.
func SplitCofactor(c, pMax, lpBound uint64) (lp1, lp2 uint64, ok bool) {
	if hi, lo := bits.Mul64(lpBound, lpBound); hi == 0 && c > lo {
		return 0, 0, false // too large for two primes up to lpBound
	}
	if c < pMax*pMax || IsPrime64(c) {
		return 0, 0, false // a prime too large to keep
	}
	d := Rho64(c)
	if d == 0 {
		return 0, 0, false
	}
	e := c / d
	if d > e {
		d, e = e, d
	}
	// Both factors are above pMax, so each is prime if below pMax².
	if e > lpBound || (d >= pMax*pMax && !IsPrime64(d)) || (e >= pMax*pMax && !IsPrime64(e)) {
		return 0, 0, false
	}
	return d, e, true
}

// PrimesUpTo returns all primes <= b.
func PrimesUpTo(b int64) []int64 {
	if b < 2 {
		return nil
	}
	composite := make([]bool, b+1)
	var primes []int64
	for p := int64(2); p <= b; p++ {
		if composite[p] {
			continue
		}
		primes = append(primes, p)
		for i := p * p; i <= b; i += p {
			composite[i] = true
		}
	}
	return primes
}

// MulMod64 returns (a * b) % p without overflowing uint64.
func MulMod64(a, b, p uint64) uint64 {
	if a >= p {
		a %= p
	}
	if b >= p {
		b %= p
	}
	if p <= 1<<32 {
		return (a * b) % p
	}
	hi, lo := bits.Mul64(a, b)
	_, rem := bits.Div64(hi, lo, p)
	return rem
}

// PowMod64 returns base^exp mod p for any uint64 p >= 1.
func PowMod64(base, exp, p uint64) uint64 {
	if p == 1 {
		return 0
	}
	res := uint64(1)
	base %= p
	if p <= 1<<32 {
		for exp > 0 {
			if exp&1 != 0 {
				res = res * base % p
			}
			base = base * base % p
			exp >>= 1
		}
		return res
	}
	for exp > 0 {
		if exp&1 != 0 {
			res = MulMod64(res, base, p)
		}
		base = MulMod64(base, base, p)
		exp >>= 1
	}
	return res
}

// InvMod returns the inverse of a mod the prime p. a must not be 0 mod p.
func InvMod(a, p uint64) uint64 {
	return PowMod64(a, p-2, p)
}

// InvMod64 returns the inverse of a mod m (for any m > 1 with gcd(a, m) == 1),
// or (0, false) if gcd(a, m) != 1.
func InvMod64(a, m uint64) (uint64, bool) {
	if m <= 1 {
		return 0, false
	}
	a %= m
	if a == 0 {
		return 0, false
	}
	var t, newT int64 = 0, 1
	r, newR := m, a
	for newR != 0 {
		q := r / newR
		t, newT = newT, t-int64(q)*newT
		r, newR = newR, r-q*newR
	}
	if r != 1 {
		return 0, false
	}
	if t < 0 {
		t += int64(m)
	}
	return uint64(t), true
}

// SqrtModPrime returns t with t² ≡ n (mod p) using Tonelli–Shanks. n must be
// a quadratic residue mod the prime p.
func SqrtModPrime(n, p uint32) uint32 {
	return uint32(SqrtModPrime64(uint64(n), uint64(p)))
}

// SqrtModPrime64 returns t with t² ≡ n (mod p) using Tonelli–Shanks for any
// 64-bit prime p. n must be a quadratic residue mod p.
func SqrtModPrime64(n, p uint64) uint64 {
	n %= p
	if n == 0 || p == 2 {
		return n
	}
	if p&3 == 3 {
		return PowMod64(n, (p+1)>>2, p)
	}

	// p − 1 = q·2^s with q odd.
	q, s := p-1, uint32(0)
	for q&1 == 0 {
		q >>= 1
		s++
	}

	// Find a quadratic non-residue z.
	z := uint64(2)
	for PowMod64(z, (p-1)>>1, p) != p-1 {
		z++
	}

	m := s
	c := PowMod64(z, q, p)
	t := PowMod64(n, q, p)
	r := PowMod64(n, (q+1)>>1, p)
	for t != 1 {
		// Find the least i with t^(2^i) = 1.
		i := uint32(1)
		for tt := MulMod64(t, t, p); tt != 1; tt = MulMod64(tt, tt, p) {
			i++
			if i == m {
				return 0 // not a quadratic residue
			}
		}
		b := c
		for j := uint32(0); j < m-i-1; j++ {
			b = MulMod64(b, b, p)
		}
		m = i
		c = MulMod64(b, b, p)
		t = MulMod64(t, c, p)
		r = MulMod64(r, b, p)
	}
	return r
}

// Log2Big returns log2(v) for v > 0. It reads v's top 53 bits from its
// words, so it doesn't allocate on 64-bit platforms.
func Log2Big(v *big.Int) float64 {
	if bits.UintSize != 64 {
		shift := v.BitLen() - 53
		if shift <= 0 {
			return math.Log2(float64(v.Uint64()))
		}
		top := new(big.Int).Rsh(v, uint(shift)).Uint64()
		return math.Log2(float64(top)) + float64(shift)
	}
	words := v.Bits()
	if len(words) == 0 {
		return math.Inf(-1)
	}
	top := uint64(words[len(words)-1])
	topBits := bits.Len64(top)
	bitLen := (len(words)-1)*64 + topBits
	if bitLen <= 53 {
		return math.Log2(float64(top))
	}
	shift := topBits - 53
	sig := top >> uint(shift)
	if shift < 0 {
		prev := uint64(words[len(words)-2])
		sig = top<<uint(-shift) | prev>>uint(64+shift)
	}
	return math.Log2(float64(sig)) + float64(bitLen-53)
}

// NthRoot returns floor(n^(1/d)) for n >= 0 and d >= 1.
func NthRoot(n *big.Int, d int) *big.Int {
	if n.Sign() == 0 || d == 1 {
		return new(big.Int).Set(n)
	}
	if d == 2 {
		return new(big.Int).Sqrt(n)
	}
	// Initial overestimate: 2^ceil(bitLen/d).
	shift := (n.BitLen() + d - 1) / d
	x := new(big.Int).Lsh(big.NewInt(1), uint(shift))
	dBig := big.NewInt(int64(d))
	d1Big := big.NewInt(int64(d - 1))
	pow := new(big.Int)
	num := new(big.Int)
	next := new(big.Int)
	for {
		pow.Exp(x, d1Big, nil)
		num.Mul(x, d1Big)
		num.Add(num, new(big.Int).Quo(n, pow))
		next.Quo(num, dBig)
		if next.Cmp(x) >= 0 {
			return x
		}
		x.Set(next)
	}
}

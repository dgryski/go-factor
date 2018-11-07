package factor

import (
	"math/big"
)

const squfofMaxIters = 1000000

var maxSqufof = big.NewInt(squfofMaxIters)

func Squfof(n *big.Int, k int64) (primes, composites []*big.Int) {

	if n.ProbablyPrime(10) {
		return []*big.Int{n}, nil
	}

	// temporaries
	a := big.NewInt(0)
	b := big.NewInt(0)

	// apply small multiplier
	N := big.NewInt(0).Mul(n, big.NewInt(k))

	// s = sqrt(N)
	s := big.NewInt(0).Sqrt(N)

	// research says "50 should be big enough" -- need link to paper
	list := make([]*big.Int, 0, 50)

	Qprev := big.NewInt(1)

	// P = s
	P := big.NewInt(0).Set(s)

	// L = floor ( 2 * sqrt(2 * s)) == sqrt(8 * s)
	// B = 3*L == 3*sqrt(8*s)
	a.Mul(s, big.NewInt(8))
	L := big.NewInt(0).Sqrt(a)
	lOverTwo := big.NewInt(0).Rsh(L, 1)
	B := big.NewInt(0).Mul(L, big.NewInt(3))

	if B.Cmp(maxSqufof) >= 0 {
		B.Set(maxSqufof)
	}
	B64 := B.Uint64()

	// Q = N - P*P
	Q := big.NewInt(0).Set(N)
	a.Mul(P, P)
	Q.Sub(Q, a)

	// needed in loop
	q := big.NewInt(0)
	Pnext := big.NewInt(0)
	t := big.NewInt(0)

	var foundSquare bool
	for i := uint64(0); i < B64; i++ {

		/* q = (s + P) / Q */
		q.Add(s, P)
		q.Div(q, Q)

		/* Pnext = q * Q - P */
		a.Mul(q, Q)
		Pnext.Sub(a, P)

		if Q.Cmp(L) <= 0 {
			e := big.NewInt(0).Set(Q)
			if e.Bit(0) == 0 {
				e.Rsh(e, 1)
				list = append(list, e)
			} else if e.Cmp(lOverTwo) <= 0 {
				list = append(list, e)
			}
		}

		/* t = Qprev + q * (P - Pnext) */
		a.Sub(P, Pnext)
		b.Mul(q, a)
		t.Add(Qprev, b)

		Qprev.Set(Q)
		Q.Set(t)
		P.Set(Pnext)

		if i&1 == 0 {
			if r, ok := isSquare(Q); ok && r.Cmp(one) > 0 && !inList(r, list) {
				// save for outside the loop
				Qprev.Set(r)
				foundSquare = true
				break
			}
		}
	}

	if !foundSquare {
		return primes, []*big.Int{n}
	}

	// alias so the formulas stay correct
	r := Qprev

	// P = P + r * floor((s-P)/r)
	b.Sub(s, P)
	a.Div(b, r)
	b.Mul(a, r)
	P.Add(P, b)

	// Q = (N - P*P) / Qprev
	b.Mul(P, P)
	a.Sub(N, b)
	Q.Div(a, Qprev)

	for {
		/* q = (s + P) / Q */
		a.Add(s, P)
		q.Div(a, Q)

		/* Pnext = q * Q - P */
		a.Mul(q, Q)
		Pnext.Sub(a, P)

		if P.Cmp(Pnext) == 0 {
			break
		}

		/* t = Qprev + q * (P - Pnext) */
		a.Sub(P, Pnext)
		b.Mul(q, a)
		t.Add(Qprev, b)

		Qprev.Set(Q)
		Q.Set(t)
		P.Set(Pnext)
	}

	if newN, newG, ok := checkGCD(n, Q); ok {
		pr, co := Squfof(newN, k)
		primes, composites = append(primes, pr...), append(composites, co...)
		pr, co = Squfof(newG, k)
		primes, composites = append(primes, pr...), append(composites, co...)
		return primes, composites
	}

	return nil, []*big.Int{n}
}

func isSquare(q *big.Int) (*big.Int, bool) {
	// https://www.johndcook.com/blog/2008/11/17/fast-way-to-test-whether-a-number-is-a-square/
	if h := q.Uint64() & 0x0f; h != 0 && h != 1 && h != 4 && h != 9 {
		return nil, false
	}

	r := big.NewInt(0).Sqrt(q)
	return r, big.NewInt(0).Mul(r, r).Cmp(q) == 0
}

func inList(r *big.Int, l []*big.Int) bool {
	for _, v := range l {
		if r.Cmp(v) == 0 {
			return true
		}
	}
	return false
}

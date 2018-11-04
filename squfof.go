package factor

import (
	"math/big"
)

const squfofMaxIters = 100000

func Squfof(n *big.Int, k int64) (primes, composites []*big.Int) {

	if n.ProbablyPrime(10) {
		return []*big.Int{n}, nil
	}

	// temporaries
	a := big.NewInt(0)
	b := big.NewInt(0)

	d := big.NewInt(0).Mul(n, big.NewInt(k))

	isqrtn := big.NewInt(0).Sqrt(d)
	bi := big.NewInt(0)
	pim1 := big.NewInt(0).Set(isqrtn)
	pi := big.NewInt(0)

	qi := big.NewInt(0).Set(d)
	a.Mul(pim1, pim1)
	qi.Sub(qi, a)

	qim1 := big.NewInt(1)
	qip1 := big.NewInt(0)

	var qiSquare bool

	for i := 0; i < squfofMaxIters; i++ {

		if isSquare(qi) {
			qiSquare = true
			break
		}

		/* bi = (isqrtn + pim1) / qi; */
		bi.Add(isqrtn, pim1)
		bi.Div(bi, qi)

		/* pi = bi * qi - pim1; */
		a.Mul(bi, qi)
		pi.Sub(a, pim1)

		/* qip1 = qim1 + bi * (pim1 - pi); */
		a.Sub(pim1, pi)
		b.Mul(bi, a)
		qip1.Add(qim1, b)

		qim1.Set(qi)
		qi.Set(qip1)
		pim1.Set(pi)
	}

	if !qiSquare {
		return primes, []*big.Int{n}
	}

	/* bi = (isqrtn - pim1) / sqrt(qi); */
	b.Sqrt(qi)
	a.Sub(isqrtn, pim1)
	bi.Div(a, b)

	/* pi = bi * sqrt(qi) + pim1; */
	a.Mul(bi, b)
	pi.Add(a, pim1)

	pim1.Set(pi)
	qim1.Set(b)

	/* qi = (n - pi * pi) / qim1; */
	a.Set(d)
	b.Mul(pi, pi)
	a.Sub(a, b)
	qi.Div(a, qim1)

	for i := 0; i < squfofMaxIters; i++ {

		/* bi = (isqrtn + pim1) / qi; */
		bi.Add(isqrtn, pim1)
		bi.Div(bi, qi)

		/* pi = bi * qi - pim1; */
		a.Mul(bi, qi)
		pi.Sub(a, pim1)

		if pim1.Cmp(pi) == 0 {
			// found a factor --
			break
		}

		/* qip1 = qim1 + bi * (pim1 - pi); */
		a.Sub(pim1, pi)
		b.Mul(bi, a)
		qip1.Add(qim1, b)

		qim1.Set(qi)
		qi.Set(qip1)
		pim1.Set(pi)
	}

	if newN, newG, ok := checkGCD(n, pi); ok {
		pr, co := Squfof(newN, k)
		primes, composites = append(primes, pr...), append(composites, co...)
		pr, co = Squfof(newG, k)
		primes, composites = append(primes, pr...), append(composites, co...)
	}

	return primes, composites
}

func isSquare(n *big.Int) bool {
	t := big.NewInt(0)
	t = t.Sqrt(n)
	t.Mul(t, t)
	return t.Cmp(n) == 0
}

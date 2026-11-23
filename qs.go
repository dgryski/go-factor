package factor

import (
	"math/big"
	"math/bits"
)

// qsExtraRelations is how many relations to collect beyond the number of
// matrix columns. Each extra relation guarantees at least one more
// dependency, and each dependency has roughly a 50% chance of giving a
// non-trivial factor.
const qsExtraRelations = 20

// legendre calculates the Legendre symbol (a/p) for an odd prime p.
func legendre(a *big.Int, p int64) int {
	bigP := big.NewInt(p)
	pow := new(big.Int).Sub(bigP, one)
	pow.Div(pow, two)

	res := new(big.Int).Exp(a, pow, bigP)
	if res.Sign() == 0 {
		return 0
	}
	if res.Cmp(one) == 0 {
		return 1
	}
	return -1
}

// primesUpTo returns all primes <= b.
func primesUpTo(b int64) []int64 {
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

// buildFactorBase returns 2 and the odd primes in primes for which n is a
// quadratic residue. n must be odd and have no factors in primes.
func buildFactorBase(n *big.Int, primes []int64) []int64 {
	base := []int64{2}
	for _, p := range primes {
		if p > 2 && legendre(n, p) == 1 {
			base = append(base, p)
		}
	}
	return base
}

// qsRelation stores a smooth relation X^2 ≡ Y (mod N), where Y is recorded
// only as its exponent vector over the factor base.
type qsRelation struct {
	X       *big.Int
	Factors []int // index 0 is the sign of Y, index i+1 is the exponent of factorBase[i]
}

// tryFactorize attempts to fully factorize y over the factor base. y must be
// non-zero.
func tryFactorize(y *big.Int, factorBase []int64) ([]int, bool) {
	temp := new(big.Int).Abs(y)
	exponents := make([]int, len(factorBase)+1)
	if y.Sign() < 0 {
		exponents[0] = 1
	}

	var bigP, q, r big.Int
	for i, p := range factorBase {
		bigP.SetInt64(p)
		for {
			q.QuoRem(temp, &bigP, &r)
			if r.Sign() != 0 {
				break
			}
			exponents[i+1]++
			temp.Set(&q)
		}
	}

	if temp.Cmp(one) == 0 {
		return exponents, true
	}
	return nil, false
}

// findDependencies returns subsets of rows whose exponent vectors sum to zero
// mod 2. Each subset is a bitset over row indices.
func findDependencies(relations []qsRelation, numCols int) [][]uint64 {
	numRows := len(relations)
	colWords := (numCols + 63) / 64
	rowWords := (numRows + 63) / 64

	// Each row holds its exponent vector mod 2 and a record of which
	// original rows were XORed into it, starting with just itself.
	vec := make([][]uint64, numRows)
	combo := make([][]uint64, numRows)
	for i, rel := range relations {
		vec[i] = make([]uint64, colWords)
		for j, e := range rel.Factors {
			if e&1 != 0 {
				vec[i][j/64] |= 1 << (j % 64)
			}
		}
		combo[i] = make([]uint64, rowWords)
		combo[i][i/64] |= 1 << (i % 64)
	}

	isPivot := make([]bool, numRows)
	for c := 0; c < numCols; c++ {
		w, bit := c/64, uint64(1)<<(c%64)
		pivot := -1
		for i := 0; i < numRows; i++ {
			if !isPivot[i] && vec[i][w]&bit != 0 {
				pivot = i
				break
			}
		}
		if pivot == -1 {
			continue
		}
		isPivot[pivot] = true
		for i := 0; i < numRows; i++ {
			if i != pivot && vec[i][w]&bit != 0 {
				for k := range vec[i] {
					vec[i][k] ^= vec[pivot][k]
				}
				for k := range combo[i] {
					combo[i][k] ^= combo[pivot][k]
				}
			}
		}
	}

	// Every column with a set bit in a non-pivot row was cleared when its
	// pivot was chosen, so non-pivot rows are now all zero.
	var deps [][]uint64
	for i := 0; i < numRows; i++ {
		if !isPivot[i] {
			deps = append(deps, combo[i])
		}
	}
	return deps
}

// tryDependency combines the relations selected by dep into a congruence of
// squares a^2 ≡ b^2 (mod n) and returns a non-trivial factor of n, or nil.
func tryDependency(n *big.Int, factorBase []int64, relations []qsRelation, dep []uint64) *big.Int {
	a := big.NewInt(1)
	exps := make([]int, len(factorBase)+1)
	for w, word := range dep {
		for word != 0 {
			i := w*64 + bits.TrailingZeros64(word)
			word &= word - 1
			a.Mul(a, relations[i].X)
			a.Mod(a, n)
			for j, e := range relations[i].Factors {
				exps[j] += e
			}
		}
	}

	// b = sqrt(∏ Y_i). The sign exponent is even, so the product is positive.
	b := big.NewInt(1)
	var p, e, term big.Int
	for j := 1; j < len(exps); j++ {
		if exps[j] == 0 {
			continue
		}
		p.SetInt64(factorBase[j-1])
		e.SetInt64(int64(exps[j] / 2))
		term.Exp(&p, &e, n)
		b.Mul(b, &term)
		b.Mod(b, n)
	}

	var t, g big.Int
	for _, op := range []func(x, y *big.Int) *big.Int{t.Sub, t.Add} {
		op(a, b)
		g.GCD(nil, nil, t.Abs(&t), n)
		if g.Cmp(one) > 0 && g.Cmp(n) < 0 {
			return new(big.Int).Set(&g)
		}
	}
	return nil
}

// QuadraticSieve factors n by finding congruences of squares X^2 ≡ Y^2 (mod n).
// Candidates are trial-divided over the factor base rather than sieved, so it
// is only practical for n up to about 60 bits.
func QuadraticSieve(n *big.Int) (primes, composites []*big.Int) {
	if n.ProbablyPrime(10) {
		return []*big.Int{n}, nil
	}

	p, q := qsSplit(n)
	if p == nil {
		return nil, []*big.Int{n}
	}

	pr, co := QuadraticSieve(p)
	primes, composites = append(primes, pr...), append(composites, co...)
	pr, co = QuadraticSieve(q)
	primes, composites = append(primes, pr...), append(composites, co...)
	return primes, composites
}

// qsSplit attempts to split n into two non-trivial factors. It returns
// nil, nil if n is prime, too small, or no factor was found.
func qsSplit(n *big.Int) (*big.Int, *big.Int) {
	if n.Cmp(big.NewInt(4)) < 0 || n.ProbablyPrime(20) {
		return nil, nil
	}

	// Perfect squares would make some Y = X^2 - N zero.
	sqrtN := new(big.Int).Sqrt(n)
	if new(big.Int).Mul(sqrtN, sqrtN).Cmp(n) == 0 {
		return sqrtN, new(big.Int).Set(sqrtN)
	}
	// Round up so that X = sqrtN gives the smallest positive Y.
	sqrtN.Add(sqrtN, one)

	var b, sieveRange int64
	switch bitLen := n.BitLen(); {
	case bitLen < 32:
		b, sieveRange = 250, 2000
	case bitLen < 48:
		b, sieveRange = 1000, 20000
	default:
		b, sieveRange = 5000, 100000
	}

	// Remove small factors first: primes that divide n have Legendre
	// symbol 0 and would otherwise be silently dropped from the base.
	primes := primesUpTo(b)
	var m big.Int
	for _, p := range primes {
		if m.Mod(n, big.NewInt(p)).Sign() == 0 {
			f := big.NewInt(p)
			return f, new(big.Int).Div(n, f)
		}
	}

	factorBase := buildFactorBase(n, primes)
	numCols := len(factorBase) + 1 // +1 for the sign of Y
	requiredRelations := numCols + qsExtraRelations

	var relations []qsRelation
	var xSquared, y big.Int
	try := func(x *big.Int) {
		xSquared.Mul(x, x)
		y.Sub(&xSquared, n)
		if exponents, smooth := tryFactorize(&y, factorBase); smooth {
			relations = append(relations, qsRelation{X: x, Factors: exponents})
		}
	}

	// Test X = sqrtN + i and X = sqrtN - i. The negative side stops before
	// X reaches zero, since X and -X give the same Y and would only add
	// duplicate rows (and so trivial dependencies).
	for i := int64(0); i < sieveRange && len(relations) < requiredRelations; i++ {
		off := big.NewInt(i)
		try(new(big.Int).Add(sqrtN, off))
		if i > 0 && len(relations) < requiredRelations {
			if x := new(big.Int).Sub(sqrtN, off); x.Sign() > 0 {
				try(x)
			}
		}
	}

	if len(relations) <= numCols {
		return nil, nil // not enough relations to guarantee a dependency
	}

	for _, dep := range findDependencies(relations, numCols) {
		if f := tryDependency(n, factorBase, relations, dep); f != nil {
			return f, new(big.Int).Div(n, f)
		}
	}
	return nil, nil
}

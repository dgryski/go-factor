package factor

import (
	"math"
	"math/big"
	"math/bits"
)

// qsExtraRelations is how many relations to collect beyond the number of
// matrix columns. Each extra relation guarantees at least one more
// dependency, and each dependency has roughly a 50% chance of giving a
// non-trivial factor.
const qsExtraRelations = 20

const (
	qsBlockSize      = 32 * 1024 // sieve block size, sized for L1 cache
	qsMinSievePrime  = 30        // smaller primes are not sieved, only allowed for in the threshold
	qsThresholdSlack = 1.75      // allow this many factors of log2(largest prime) to be missing
)

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

// qsRoot holds the sieve state for one factor base prime p. Q(x) = (x+m)^2 - N
// is divisible by p exactly when x ≡ r1 or r2 (mod p), where r1, r2 = ±t - m
// and t^2 ≡ N (mod p).
type qsRoot struct {
	p          uint32
	logP       byte   // rounded log2(p)
	r1, r2     uint32 // roots of Q(x) mod p; equal for p = 2
	off1, off2 uint32 // offset of the next hit for r1, r2 in the current block
}

// powMod64 returns base^exp mod p. p must be less than 2^32.
func powMod64(base, exp, p uint64) uint64 {
	res := uint64(1)
	base %= p
	for exp > 0 {
		if exp&1 != 0 {
			res = res * base % p
		}
		base = base * base % p
		exp >>= 1
	}
	return res
}

// sqrtModPrime returns t with t^2 ≡ n (mod p) using Tonelli-Shanks. n must be
// a quadratic residue mod the prime p.
func sqrtModPrime(n, p uint32) uint32 {
	p64 := uint64(p)
	n64 := uint64(n) % p64
	if n64 == 0 || p == 2 {
		return uint32(n64)
	}
	if p&3 == 3 {
		return uint32(powMod64(n64, (p64+1)>>2, p64))
	}

	// p - 1 = q * 2^s with q odd.
	q, s := p64-1, uint32(0)
	for q&1 == 0 {
		q >>= 1
		s++
	}

	// Find a quadratic non-residue z.
	z := uint64(2)
	for powMod64(z, (p64-1)>>1, p64) != p64-1 {
		z++
	}

	m := s
	c := powMod64(z, q, p64)
	t := powMod64(n64, q, p64)
	r := powMod64(n64, (q+1)>>1, p64)
	for t != 1 {
		// Find the least i with t^(2^i) = 1.
		i := uint32(1)
		for tt := t * t % p64; tt != 1; tt = tt * tt % p64 {
			i++
			if i == m {
				return 0 // not a quadratic residue
			}
		}
		b := c
		for j := uint32(0); j < m-i-1; j++ {
			b = b * b % p64
		}
		m = i
		c = b * b % p64
		t = t * c % p64
		r = r * b % p64
	}
	return uint32(r)
}

// buildRootTable computes the sieve roots of Q(x) = (x+m)^2 - n for each
// prime in factorBase. No prime in factorBase may divide n.
func buildRootTable(n, m *big.Int, factorBase []int64) []qsRoot {
	roots := make([]qsRoot, len(factorBase))
	var bigP, r big.Int
	for i, p := range factorBase {
		bigP.SetInt64(p)
		t := uint64(sqrtModPrime(uint32(r.Mod(n, &bigP).Uint64()), uint32(p)))
		mModP := r.Mod(m, &bigP).Uint64()
		p64 := uint64(p)
		roots[i] = qsRoot{
			p:    uint32(p),
			logP: byte(math.Round(math.Log2(float64(p)))),
			r1:   uint32((t + p64 - mModP) % p64),
			r2:   uint32((2*p64 - t - mModP) % p64),
		}
	}
	return roots
}

// nextOffset returns the smallest off >= 0 with start+off ≡ root (mod p).
func nextOffset(start int64, root, p uint32) uint32 {
	off := (int64(root) - start) % int64(p)
	if off < 0 {
		off += int64(p)
	}
	return uint32(off)
}

// tryFactorize attempts to fully factor y = Q(x) over the factor base,
// trial dividing only by primes whose sieve roots match x. y must be
// non-zero.
func tryFactorize(x int64, y *big.Int, roots []qsRoot) ([]int, bool) {
	temp := new(big.Int).Abs(y)
	exponents := make([]int, len(roots)+1)
	if y.Sign() < 0 {
		exponents[0] = 1
	}

	var bigP, q, r big.Int
	for i := range roots {
		root := &roots[i]
		xr := x % int64(root.p)
		if xr < 0 {
			xr += int64(root.p)
		}
		if uint32(xr) != root.r1 && uint32(xr) != root.r2 {
			continue
		}
		bigP.SetUint64(uint64(root.p))
		for {
			q.QuoRem(temp, &bigP, &r)
			if r.Sign() != 0 {
				break
			}
			exponents[i+1]++
			temp.Set(&q)
		}
		if temp.Cmp(one) == 0 {
			return exponents, true
		}
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
// It sieves a single polynomial, so it is only practical for n up to about
// 120 bits.
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

	b, sieveRange := qsParams(n)

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
	relations := sieve(n, sqrtN, factorBase, sieveRange, requiredRelations)

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

// qsParams chooses the factor base bound b ≈ 3 L(n)^(1/2), where
// L(n) = exp(sqrt(ln n ln ln n)), and the sieve half-width. Without the
// factor of 3 the base is too small to find enough relations below about
// 80 bits.
func qsParams(n *big.Int) (b, sieveRange int64) {
	lnN := float64(n.BitLen()) * math.Ln2
	b = int64(3 * math.Sqrt(math.Exp(math.Sqrt(lnN*math.Log(lnN)))))
	if b < 500 {
		b = 500 // smaller bases often run out of sieve range before finding enough relations
	}
	mult := b / 3
	if mult < 256 {
		mult = 256
	}
	return b, mult * b
}

// sieve collects up to want relations X^2 ≡ Y (mod n) with Y smooth over
// factorBase, for X = m + x and x in [-sieveRange, sieveRange). Each block
// adds log2(p) at every x where p divides Y, and only positions whose total
// is close to log2|Y| are trial divided.
func sieve(n, m *big.Int, factorBase []int64, sieveRange int64, want int) []qsRelation {
	roots := buildRootTable(n, m, factorBase)

	// Stop before X reaches zero, since X and -X give the same Y and would
	// only add duplicate rows (and so trivial dependencies).
	lo := -sieveRange
	if m.IsInt64() && 1-m.Int64() > lo {
		lo = 1 - m.Int64()
	}

	// Skip the small primes, which hit often but add little, and lower the
	// threshold by their expected contribution of about 2 log2(p) / p.
	start, smallSlack := 0, 0.0
	for start < len(roots) && roots[start].p < qsMinSievePrime {
		r := roots[start]
		if r.p == 2 {
			smallSlack += float64(r.logP) / 2
		} else {
			smallSlack += 2 / float64(r.p) * float64(r.logP)
		}
		start++
	}
	if len(roots)-start < 8 {
		start, smallSlack = 0, 0
	}

	for i := start; i < len(roots); i++ {
		r := &roots[i]
		r.off1 = nextOffset(lo, r.r1, r.p)
		r.off2 = nextOffset(lo, r.r2, r.p)
	}

	mf, _ := new(big.Float).SetInt(m).Float64()
	slack := qsThresholdSlack*math.Log2(float64(roots[len(roots)-1].p)) + smallSlack

	var relations []qsRelation
	var buf [qsBlockSize]byte
	var x, y big.Int
	for blockStart := lo; blockStart < sieveRange && len(relations) < want; blockStart += qsBlockSize {
		blockLen := int64(qsBlockSize)
		if sieveRange-blockStart < blockLen {
			blockLen = sieveRange - blockStart
		}
		block := buf[:blockLen]
		for i := range block {
			block[i] = 0
		}

		for i := start; i < len(roots); i++ {
			r := &roots[i]
			p := int(r.p)
			j := int(r.off1)
			for ; j < len(block); j += p {
				block[j] += r.logP
			}
			r.off1 = uint32(j - len(block))
			if r.r2 != r.r1 {
				j = int(r.off2)
				for ; j < len(block); j += p {
					block[j] += r.logP
				}
				r.off2 = uint32(j - len(block))
			}
		}

		// |Y| = |x| (2m + x) is largest at whichever end of the block is
		// furthest from zero.
		maxX := float64(blockStart + blockLen - 1)
		if a := float64(-blockStart); a > maxX {
			maxX = a
		}
		if maxX < 1 {
			maxX = 1
		}
		threshold := math.Log2(maxX*(2*mf+maxX)) - slack

		for i, score := range block {
			if float64(score) < threshold {
				continue
			}
			off := blockStart + int64(i)
			x.SetInt64(off)
			x.Add(&x, m)
			y.Mul(&x, &x)
			y.Sub(&y, n)
			if exponents, smooth := tryFactorize(off, &y, roots); smooth {
				relations = append(relations, qsRelation{X: new(big.Int).Set(&x), Factors: exponents})
				if len(relations) >= want {
					break
				}
			}
		}
	}
	return relations
}

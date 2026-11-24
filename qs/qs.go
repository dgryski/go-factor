// Package qs implements the quadratic sieve.
//
// Small inputs sieve the single polynomial (x + ⌈√kN⌉)² − kN. Larger inputs
// use the self-initializing quadratic sieve (SIQS), which switches between
// many polynomials (Ax + B)² − kN so that the values being sieved stay small.
// Both use a Knuth–Schroeppel multiplier k and accept relations with one
// large prime.
package qs

import (
	"math"
	"math/big"
	"math/bits"
	"math/rand"
	"sort"
	"time"
)

var one = big.NewInt(1)

// extraRelations is how many relations to collect beyond the number of
// matrix columns. Each extra relation guarantees at least one more
// dependency, and each dependency has roughly a 50% chance of giving a
// non-trivial factor.
const extraRelations = 20

const (
	blockSize     = 32 * 1024 // sieve block size, sized for L1 cache
	minSievePrime = 30        // smaller primes are not sieved, only allowed for in the threshold

	// A candidate may score this many times log2(largest base prime) below
	// log2|Y|.
	singleSlack = 1.75
	siqsSlack   = 1.25

	// SIQS also keeps partial relations with one prime up to this many times
	// the largest base prime, and lowers the threshold by half of that
	// prime's size. Lowering it by all of it adds too many candidates that
	// fail trial division. With a single polynomial, partials don't pay.
	largePrimeMult = 64

	siqsMinBits = 60   // below this, a single polynomial is faster
	siqsMinB    = 1000 // smaller bases leave too few primes to build enough different A values
)

// QuadraticSieve factors n by finding congruences of squares X² ≡ Y² (mod n).
func QuadraticSieve(n *big.Int) (primes, composites []*big.Int) {
	if n.ProbablyPrime(10) {
		return []*big.Int{n}, nil
	}

	p, q := split(n)
	if p == nil {
		return nil, []*big.Int{n}
	}

	pr, co := QuadraticSieve(p)
	primes, composites = append(primes, pr...), append(composites, co...)
	pr, co = QuadraticSieve(q)
	primes, composites = append(primes, pr...), append(composites, co...)
	return primes, composites
}

// split attempts to split n into two non-trivial factors. It returns
// nil, nil if n is prime, too small, or no factor was found.
func split(n *big.Int) (*big.Int, *big.Int) {
	if n.Cmp(big.NewInt(4)) < 0 || n.ProbablyPrime(20) {
		return nil, nil
	}

	// Perfect squares would make some Y zero.
	sqrtN := new(big.Int).Sqrt(n)
	if new(big.Int).Mul(sqrtN, sqrtN).Cmp(n) == 0 {
		return sqrtN, new(big.Int).Set(sqrtN)
	}

	b, siqs, sieveRange := params(n)

	// Remove small factors first: primes that divide n would otherwise be
	// silently dropped from the base, and the multiplier must be coprime to n.
	primes := primesUpTo(b)
	var r, bigP big.Int
	for _, p := range primes {
		bigP.SetInt64(p)
		if r.Mod(n, &bigP).Sign() == 0 {
			return new(big.Int).Set(&bigP), new(big.Int).Div(n, &bigP)
		}
	}

	k := knuthSchroeppel(n)
	kn := new(big.Int).Mul(n, big.NewInt(k))
	factorBase := buildFactorBase(kn, primes)
	s := newSieve(n, kn, factorBase, siqs)

	method := "single polynomial"
	if siqs {
		method = "SIQS"
	}
	logf("qs: %d bits, %s, k=%d, B=%d (%d primes), M=%d, want %d relations",
		n.BitLen(), method, k, b, len(factorBase), sieveRange, s.want)

	var f *big.Int
	if siqs {
		f = s.siqs(sieveRange)
	} else {
		f = s.singlePoly(sieveRange)
	}
	s.logProgress()
	if f != nil {
		logf("qs: found factor %v while sieving", f)
		return f, new(big.Int).Div(n, f)
	}

	if len(s.relations) <= s.numCols {
		logf("qs: not enough relations, giving up")
		return nil, nil // not enough relations to guarantee a dependency
	}
	start := time.Now()
	deps := findDependencies(s.relations, s.numCols)
	logf("qs: matrix %d×%d, %d dependencies in %v",
		len(s.relations), s.numCols, len(deps), time.Since(start).Round(time.Millisecond))
	for i, dep := range deps {
		if f := tryDependency(n, factorBase, s.relations, dep); f != nil {
			logf("qs: dependency %d gave factor %v", i+1, f)
			return f, new(big.Int).Div(n, f)
		}
	}
	logf("qs: all dependencies were trivial")
	return nil, nil
}

// Logf, if set, receives progress messages while the sieve runs.
var Logf func(format string, args ...interface{})

func logf(format string, args ...interface{}) {
	if Logf != nil {
		Logf(format, args...)
	}
}

// logProgress reports how many relations the sieve has collected.
func (s *sieve) logProgress() {
	full := len(s.relations) - s.combined
	logf("qs: %d/%d relations (%d full, %d from partials), %d partials, %d polynomials, %v",
		len(s.relations), s.want, full, s.combined, len(s.partials), s.polys,
		time.Since(s.start).Round(time.Millisecond))
}

// params chooses the factor base bound b and the sieve half-width, and
// whether to use SIQS. L(n) = exp(sqrt(ln n ln ln n)).
func params(n *big.Int) (b int64, siqs bool, sieveRange int64) {
	lnN := float64(n.BitLen()) * math.Ln2
	l := math.Exp(math.Sqrt(lnN * math.Log(lnN)))

	if n.BitLen() < siqsMinBits {
		// A single polynomial's values grow with x, so it needs a larger
		// base. Without the factor of 3 there aren't enough relations below
		// about 80 bits.
		b = int64(3 * math.Sqrt(l))
		if b < 500 {
			b = 500 // smaller bases often run out of sieve range before finding enough relations
		}
		mult := b / 3
		if mult < 256 {
			mult = 256
		}
		return b, false, mult * b
	}

	b = int64(math.Pow(l, 0.46))
	if b < siqsMinB {
		b = siqsMinB
	}
	switch {
	case b < 10000:
		sieveRange = 32 * 1024
	case b < 100000:
		sieveRange = 64 * 1024
	default:
		sieveRange = 128 * 1024
	}
	return b, true, sieveRange
}

// knuthSchroeppel picks a small odd square-free multiplier k so that kN has
// many small primes as quadratic residues. n must have no factors below 100.
func knuthSchroeppel(n *big.Int) int64 {
	smallPrimes := primesUpTo(97)[1:]
	nModP := make([]uint64, len(smallPrimes))
	var r, bigP big.Int
	for i, p := range smallPrimes {
		bigP.SetInt64(p)
		nModP[i] = r.Mod(n, &bigP).Uint64()
	}
	nMod8 := r.And(n, big.NewInt(7)).Uint64()

	bestK, bestScore := int64(1), math.Inf(-1)
	for k := int64(1); k < 100; k += 2 {
		if !squareFree(k) {
			continue
		}

		// Expected contribution of each small prime to log(Y), less the cost
		// of Y being sqrt(k) times larger.
		score := -0.5 * math.Log(float64(k))
		switch uint64(k) * nMod8 & 7 {
		case 1:
			score += 2 * math.Ln2
		case 5:
			score += math.Ln2
		default:
			score += 0.5 * math.Ln2
		}
		for i, p := range smallPrimes {
			knp := uint64(k) * nModP[i] % uint64(p)
			lnP := math.Log(float64(p))
			if knp == 0 {
				score += lnP / float64(p)
			} else if powMod64(knp, uint64(p-1)/2, uint64(p)) == 1 {
				score += 2 * lnP / float64(p-1)
			}
		}

		if score > bestScore {
			bestK, bestScore = k, score
		}
	}
	return bestK
}

func squareFree(k int64) bool {
	for p := int64(2); p*p <= k; p++ {
		if k%(p*p) == 0 {
			return false
		}
	}
	return true
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

// buildFactorBase returns the primes in primes for which kn is a square
// mod p: 2, the primes dividing the multiplier, and the odd primes with
// Legendre symbol 1.
func buildFactorBase(kn *big.Int, primes []int64) []int64 {
	var base []int64
	var r, bigP big.Int
	for _, p := range primes {
		bigP.SetInt64(p)
		knp := r.Mod(kn, &bigP).Uint64()
		if p == 2 || knp == 0 || powMod64(knp, uint64(p-1)/2, uint64(p)) == 1 {
			base = append(base, p)
		}
	}
	return base
}

// powMod64 returns base^exp mod p. p must be less than 2^32.
func powMod64(base, exp, p uint64) uint64 {
	res := uint64(1) % p
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

// invMod returns the inverse of a mod the prime p. a must not be 0 mod p.
func invMod(a, p uint64) uint64 {
	return powMod64(a, p-2, p)
}

// sqrtModPrime returns t with t² ≡ n (mod p) using Tonelli–Shanks. n must be
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

	// p − 1 = q·2^s with q odd.
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

// log2Big returns log2(v) for v > 0.
func log2Big(v *big.Int) float64 {
	shift := v.BitLen() - 53
	if shift <= 0 {
		return math.Log2(float64(v.Uint64()))
	}
	top := new(big.Int).Rsh(v, uint(shift)).Uint64()
	return math.Log2(float64(top)) + float64(shift)
}

// relation stores X² ≡ Y (mod N), with Y recorded only as its factors.
type relation struct {
	X *big.Int // reduced mod N
	// Factors has one entry per prime factor of Y, counted with
	// multiplicity: 0 for −1, and i+1 for factorBase[i].
	Factors []int32
}

// root holds the sieve state for one factor base prime p. For the
// polynomial Y(x) = ((Ax + B)² − kN) / A, p divides Y(x) exactly when
// x ≡ r1 or r2 (mod p), where r1, r2 = (±t − B)·A⁻¹ and t² ≡ kN (mod p).
type root struct {
	p          uint32
	t          uint32 // 0 if p divides the multiplier, giving a single root
	logP       byte   // rounded log2(p)
	inA        bool   // p divides the current A, so it isn't sieved
	r1, r2     uint32 // equal for p = 2 and for p dividing the multiplier
	off1, off2 uint32 // offset of the next hit for r1, r2 in the current block
}

type sieve struct {
	n, kn     *big.Int
	log2KN    float64
	roots     []root
	numCols   int
	want      int
	sieveFrom int     // index of the first sieved prime
	slack     float64 // how far below log2|Y| a candidate may score
	lpBound   uint64

	relations []relation
	partials  map[uint64]relation // partial relations by their large prime
	seenX     map[string]bool     // min(X, n−X) for every relation and partial
	scratch   []int32
	buf       [blockSize]byte

	// For progress logging.
	combined int // relations made from two partials
	polys    int // polynomials sieved
	start    time.Time
	lastLog  time.Time
}

func newSieve(n, kn *big.Int, factorBase []int64, siqs bool) *sieve {
	s := &sieve{
		n:        n,
		kn:       kn,
		log2KN:   log2Big(kn),
		roots:    make([]root, len(factorBase)),
		numCols:  len(factorBase) + 1, // +1 for the sign of Y
		partials: make(map[uint64]relation),
		seenX:    make(map[string]bool),
		start:    time.Now(),
	}
	s.want = s.numCols + extraRelations

	var r, bigP big.Int
	for i, p := range factorBase {
		bigP.SetInt64(p)
		s.roots[i] = root{
			p:    uint32(p),
			t:    sqrtModPrime(uint32(r.Mod(kn, &bigP).Uint64()), uint32(p)),
			logP: byte(math.Round(math.Log2(float64(p)))),
		}
	}

	// Skip the small primes, which hit often but add little, and lower the
	// threshold by their expected contribution of about 2·log2(p)/p.
	smallSlack := 0.0
	for s.sieveFrom < len(s.roots) && s.roots[s.sieveFrom].p < minSievePrime {
		r := s.roots[s.sieveFrom]
		hits := 2.0
		if r.p == 2 || r.t == 0 {
			hits = 1
		}
		smallSlack += hits / float64(r.p) * float64(r.logP)
		s.sieveFrom++
	}
	if len(s.roots)-s.sieveFrom < 8 {
		s.sieveFrom, smallSlack = 0, 0
	}

	pMax := uint64(s.roots[len(s.roots)-1].p)
	s.slack = singleSlack*math.Log2(float64(pMax)) + smallSlack
	if siqs {
		s.lpBound = largePrimeMult * pMax
		s.slack = siqsSlack*math.Log2(float64(pMax)) + math.Log2(float64(s.lpBound))/2 + smallSlack
	}
	return s
}

// singlePoly sieves Y(x) = (x + m)² − kN, with m = ⌈√kN⌉, for x in
// [−sieveRange, sieveRange). It returns a factor of n if it finds one by
// accident.
func (s *sieve) singlePoly(sieveRange int64) *big.Int {
	m := new(big.Int).Sqrt(s.kn)
	m.Add(m, one)
	s.setRoots(one, m)

	// Stop before X reaches zero, since X and −X give the same Y and would
	// only add duplicate rows (and so trivial dependencies).
	lo := -sieveRange
	if m.IsInt64() && 1-m.Int64() > lo {
		lo = 1 - m.Int64()
	}
	s.polys = 1
	return s.sievePoly(one, m, nil, lo, sieveRange)
}

// siqs sieves x in [−sieveRange, sieveRange) for each of a sequence of
// polynomials until it has enough relations or runs out of polynomials. It
// returns a factor of n if it finds one by accident.
func (s *sieve) siqs(sieveRange int64) *big.Int {
	// With A ≈ √(2kN) / M, the largest |Y| at x = 0 and at x = ±M are
	// about the same.
	log2Target := (s.log2KN+1)/2 - math.Log2(float64(sieveRange))

	var pool []int
	for i, r := range s.roots {
		if r.p > minSievePrime && r.t != 0 {
			pool = append(pool, i)
		}
	}

	rng := rand.New(rand.NewSource(1))
	seen := make(map[string]bool)
	var bigP, r, aq, twoBv big.Int
	aInv := make([]uint32, len(s.roots))
	var deltas [][]uint32

	for len(s.relations) < s.want {
		a, aIdx := chooseA(rng, s.roots, pool, log2Target, seen)
		if a == nil {
			logf("qs: ran out of polynomials")
			return nil
		}
		for _, idx := range aIdx {
			s.roots[idx].inA = true
		}

		// B = Σ B_j with B_j ≡ ±t_j (mod q_j) and B_j ≡ 0 (mod A/q_j), so
		// B² ≡ kN (mod A). Flipping the signs of B_0..B_{s−2} gives
		// 2^(s−1) polynomials for this A.
		bj := make([]*big.Int, len(aIdx))
		b := new(big.Int)
		for j, idx := range aIdx {
			q := uint64(s.roots[idx].p)
			bigP.SetUint64(q)
			aq.Quo(a, &bigP)
			gamma := uint64(s.roots[idx].t) * invMod(r.Mod(&aq, &bigP).Uint64(), q) % q
			if gamma > q/2 {
				gamma = q - gamma
			}
			bj[j] = new(big.Int).Mul(new(big.Int).SetUint64(gamma), &aq)
			b.Add(b, bj[j])
		}

		// When the sign of B_v flips, every root moves by 2·B_v·A⁻¹ (mod p).
		for len(deltas) < len(aIdx)-1 {
			deltas = append(deltas, make([]uint32, len(s.roots)))
		}
		for i := range s.roots {
			rt := &s.roots[i]
			if rt.inA {
				continue
			}
			p := uint64(rt.p)
			bigP.SetUint64(p)
			aInv[i] = uint32(invMod(r.Mod(a, &bigP).Uint64(), p))
			for v := 0; v < len(aIdx)-1; v++ {
				deltas[v][i] = uint32(2 * r.Mod(bj[v], &bigP).Uint64() * uint64(aInv[i]) % p)
			}
		}
		s.setRootsInv(b, aInv)

		for poly := 1; ; poly++ {
			if f := s.sievePoly(a, b, aIdx, -sieveRange, sieveRange); f != nil {
				return f
			}
			s.polys++
			if Logf != nil && time.Since(s.lastLog) > 5*time.Second {
				s.logProgress()
				s.lastLog = time.Now()
			}
			if poly == 1<<uint(len(aIdx)-1) || len(s.relations) >= s.want {
				break
			}

			// Gray code order: step poly flips the sign of B_v.
			v := bits.TrailingZeros(uint(poly))
			twoBv.Lsh(bj[v], 1)
			if poly>>uint(v+1)&1 == 0 {
				b.Sub(b, &twoBv)
				s.moveRoots(deltas[v], true)
			} else {
				b.Add(b, &twoBv)
				s.moveRoots(deltas[v], false)
			}
		}

		for _, idx := range aIdx {
			s.roots[idx].inA = false
		}
	}
	return nil
}

// chooseA returns a product A of distinct primes from pool, with log2(A)
// close to log2Target, and the indices of its primes. It returns nil if it
// can't find an A that isn't in seen.
func chooseA(rng *rand.Rand, roots []root, pool []int, log2Target float64, seen map[string]bool) (*big.Int, []int) {
	if len(pool) < 4 {
		return nil, nil
	}
	poolP := func(i int) float64 { return float64(roots[pool[i]].p) }

	// Aim for primes around 2000, as far as the pool and target allow, and
	// at least two of them so that there is some choice.
	n := int(math.Round(log2Target / 11))
	for n < len(pool)/2 && log2Target/float64(n) > math.Log2(poolP(len(pool)-1)) {
		n++
	}
	for n > 2 && log2Target/float64(n) < math.Log2(poolP(0))+1 {
		n--
	}
	if n < 2 {
		n = 2
	}

	// Pick all but the last prime at random from those within a factor of
	// 1.5 of the ideal size, widening the range if it's too small.
	q := math.Exp2(log2Target / float64(n))
	lo := sort.Search(len(pool), func(i int) bool { return poolP(i) >= q/1.5 })
	hi := sort.Search(len(pool), func(i int) bool { return poolP(i) > q*1.5 })
	for hi-lo < 2*n+8 && (lo > 0 || hi < len(pool)) {
		if lo > 0 {
			lo--
		}
		if hi < len(pool) {
			hi++
		}
	}
	if hi-lo < n {
		return nil, nil
	}

	chosen := make([]int, n)
	used := make(map[int]bool, n)
	for attempt := 0; attempt < 100; attempt++ {
		for k := range used {
			delete(used, k)
		}
		log2A := 0.0
		for j := 0; j < n-1; j++ {
			pi := lo + rng.Intn(hi-lo)
			for used[pi] {
				pi = lo + rng.Intn(hi-lo)
			}
			used[pi] = true
			chosen[j] = pi
			log2A += math.Log2(poolP(pi))
		}

		// The last prime is the unused one closest to the remainder.
		want := math.Exp2(log2Target - log2A)
		pos := sort.Search(len(pool), func(i int) bool { return poolP(i) >= want })
		last := -1
		for d := 0; d < len(pool) && last < 0; d++ {
			if i := pos + d; i < len(pool) && !used[i] {
				last = i
			} else if i := pos - 1 - d; i >= 0 && !used[i] {
				last = i
			}
		}
		chosen[n-1] = last

		aIdx := make([]int, n)
		a := big.NewInt(1)
		var bigP big.Int
		for j, pi := range chosen {
			aIdx[j] = pool[pi]
			a.Mul(a, bigP.SetUint64(uint64(roots[pool[pi]].p)))
		}
		if key := a.String(); !seen[key] {
			seen[key] = true
			sort.Ints(aIdx)
			return a, aIdx
		}
	}
	return nil, nil
}

// setRoots computes the roots of every prime not in A for the polynomial
// with coefficients a and b.
func (s *sieve) setRoots(a, b *big.Int) {
	aInv := make([]uint32, len(s.roots))
	var r, bigP big.Int
	for i := range s.roots {
		if !s.roots[i].inA {
			p := uint64(s.roots[i].p)
			aInv[i] = uint32(invMod(r.Mod(a, bigP.SetUint64(p)).Uint64(), p))
		}
	}
	s.setRootsInv(b, aInv)
}

// setRootsInv is setRoots with A⁻¹ mod p already computed for each prime.
func (s *sieve) setRootsInv(b *big.Int, aInv []uint32) {
	var r, bigP big.Int
	for i := range s.roots {
		rt := &s.roots[i]
		if rt.inA {
			continue
		}
		p, t := uint64(rt.p), uint64(rt.t)
		bModP := r.Mod(b, bigP.SetUint64(p)).Uint64()
		rt.r1 = uint32((t + p - bModP) % p * uint64(aInv[i]) % p)
		rt.r2 = uint32((2*p - t - bModP) % p * uint64(aInv[i]) % p)
	}
}

// moveRoots adds delta[i] to both roots of each prime, or subtracts it.
func (s *sieve) moveRoots(delta []uint32, add bool) {
	for i := range s.roots {
		rt := &s.roots[i]
		d := delta[i]
		if rt.inA || d == 0 {
			continue
		}
		if !add {
			d = rt.p - d
		}
		rt.r1 = addMod(rt.r1, d, rt.p)
		rt.r2 = addMod(rt.r2, d, rt.p)
	}
}

func addMod(a, b, p uint32) uint32 {
	if a >= p-b {
		return a - (p - b)
	}
	return a + b
}

// sievePoly sieves Y(x) = ((ax + b)² − kN) / a for x in [lo, hi) and adds
// the relations it finds. aIdx lists the primes dividing a. It returns a
// factor of n if it finds one by accident.
func (s *sieve) sievePoly(a, b *big.Int, aIdx []int, lo, hi int64) *big.Int {
	for i := s.sieveFrom; i < len(s.roots); i++ {
		rt := &s.roots[i]
		rt.off1 = nextOffset(lo, rt.r1, rt.p)
		rt.off2 = nextOffset(lo, rt.r2, rt.p)
	}

	var x, y big.Int
	for blockStart := lo; blockStart < hi && len(s.relations) < s.want; blockStart += blockSize {
		blockLen := int64(blockSize)
		if hi-blockStart < blockLen {
			blockLen = hi - blockStart
		}
		block := s.buf[:blockLen]
		for i := range block {
			block[i] = 0
		}

		for i := s.sieveFrom; i < len(s.roots); i++ {
			rt := &s.roots[i]
			if rt.inA {
				continue
			}
			p := int(rt.p)
			j := int(rt.off1)
			for ; j < len(block); j += p {
				block[j] += rt.logP
			}
			rt.off1 = uint32(j - len(block))
			if rt.r2 != rt.r1 {
				j = int(rt.off2)
				for ; j < len(block); j += p {
					block[j] += rt.logP
				}
				rt.off2 = uint32(j - len(block))
			}
		}

		threshold := s.maxLog2Y(a, b, blockStart, blockStart+blockLen) - s.slack
		for i, score := range block {
			if float64(score) < threshold {
				continue
			}
			off := blockStart + int64(i)
			x.SetInt64(off)
			x.Mul(&x, a)
			x.Add(&x, b)
			y.Mul(&x, &x)
			y.Sub(&y, s.kn)
			y.Quo(&y, a)
			if y.Sign() == 0 {
				continue
			}
			factors, lp, ok := s.tryFactorize(off, &y, aIdx)
			if !ok {
				continue
			}
			if f := s.addRelation(new(big.Int).Mod(&x, s.n), factors, lp); f != nil {
				return f
			}
			if len(s.relations) >= s.want {
				break
			}
		}
	}
	return nil
}

// maxLog2Y returns log2 of the largest |Y(x)| for x in [from, to). |Y| is
// largest either at one end, or where ax + b = 0 if that lies in between.
func (s *sieve) maxLog2Y(a, b *big.Int, from, to int64) float64 {
	var x, y big.Int
	best, signs := 0.0, 0
	for _, off := range [2]int64{from, to - 1} {
		x.SetInt64(off)
		x.Mul(&x, a)
		x.Add(&x, b)
		signs += x.Sign()
		y.Mul(&x, &x)
		y.Sub(&y, s.kn)
		if y.Sign() != 0 {
			best = math.Max(best, log2Big(y.Abs(&y)))
		}
	}
	if signs == 0 {
		best = math.Max(best, s.log2KN)
	}
	return best - log2Big(a)
}

// nextOffset returns the smallest off >= 0 with start+off ≡ root (mod p).
func nextOffset(start int64, root, p uint32) uint32 {
	off := (int64(root) - start) % int64(p)
	if off < 0 {
		off += int64(p)
	}
	return uint32(off)
}

// tryFactorize attempts to factor y = Y(x) over the factor base, trial
// dividing only by primes whose roots match x. A relation's Y is a·Y(x), so
// each prime of a is counted once more. A cofactor that is a prime up to the
// large prime bound is returned as lp.
func (s *sieve) tryFactorize(x int64, y *big.Int, aIdx []int) (factors []int32, lp uint64, ok bool) {
	factors = s.scratch[:0]
	temp := new(big.Int).Abs(y)
	if y.Sign() < 0 {
		factors = append(factors, 0)
	}

	var bigP, q, r big.Int
	divide := func(i int) {
		bigP.SetUint64(uint64(s.roots[i].p))
		for {
			q.QuoRem(temp, &bigP, &r)
			if r.Sign() != 0 {
				break
			}
			factors = append(factors, int32(i+1))
			temp.Set(&q)
		}
	}

	for _, i := range aIdx {
		factors = append(factors, int32(i+1))
		divide(i)
	}
	for i := range s.roots {
		rt := &s.roots[i]
		if rt.inA {
			continue
		}
		xr := x % int64(rt.p)
		if xr < 0 {
			xr += int64(rt.p)
		}
		if uint32(xr) == rt.r1 || uint32(xr) == rt.r2 {
			if divide(i); temp.Cmp(one) == 0 {
				break
			}
		}
	}
	s.scratch = factors[:0]

	// Every prime up to the largest base prime has been divided out, so a
	// cofactor below its square is prime.
	if temp.Cmp(one) != 0 {
		if !temp.IsUint64() || temp.Uint64() > s.lpBound {
			return nil, 0, false
		}
		lp = temp.Uint64()
	}
	return append([]int32(nil), factors...), lp, true
}

// addRelation records a relation, combining it with an earlier partial
// relation if it has the same large prime lp. It returns a factor of n if lp
// divides n.
func (s *sieve) addRelation(x *big.Int, factors []int32, lp uint64) *big.Int {
	// Every polynomial gives Y = X² − kN, so different polynomials can find
	// the same X, and −X gives the same Y. Either way the copy only adds a
	// trivial dependency.
	key := new(big.Int).Lsh(x, 1)
	if key.Cmp(s.n) > 0 {
		key.Sub(s.n, x)
	} else {
		key.Set(x)
	}
	if k := string(key.Bytes()); s.seenX[k] {
		return nil
	} else {
		s.seenX[k] = true
	}

	if lp == 0 {
		s.relations = append(s.relations, relation{X: x, Factors: factors})
		return nil
	}

	prev, ok := s.partials[lp]
	if !ok {
		s.partials[lp] = relation{X: x, Factors: factors}
		return nil
	}
	// X1² ≡ Y1·L and X2² ≡ Y2·L give (X1·X2/L)² ≡ Y1·Y2.
	l := new(big.Int).SetUint64(lp)
	lInv := new(big.Int).ModInverse(l, s.n)
	if lInv == nil {
		return l // a prime that isn't invertible divides n
	}
	x.Mul(x, prev.X)
	x.Mul(x, lInv)
	x.Mod(x, s.n)
	combined := append(append([]int32(nil), prev.Factors...), factors...)
	s.relations = append(s.relations, relation{X: x, Factors: combined})
	s.combined++
	return nil
}

// findDependencies returns subsets of rows whose exponent vectors sum to zero
// mod 2. Each subset is a bitset over row indices.
func findDependencies(relations []relation, numCols int) [][]uint64 {
	numRows := len(relations)
	colWords := (numCols + 63) / 64
	rowWords := (numRows + 63) / 64

	// Each row holds its exponent vector mod 2 and a record of which
	// original rows were XORed into it, starting with just itself. The rows
	// share two flat backing arrays.
	vecBuf := make([]uint64, numRows*colWords)
	comboBuf := make([]uint64, numRows*rowWords)
	vec := make([][]uint64, numRows)
	combo := make([][]uint64, numRows)
	for i, rel := range relations {
		vec[i] = vecBuf[i*colWords : (i+1)*colWords]
		for _, c := range rel.Factors {
			vec[i][c/64] ^= 1 << (uint(c) % 64)
		}
		combo[i] = comboBuf[i*rowWords : (i+1)*rowWords]
		combo[i][i/64] |= 1 << (uint(i) % 64)
	}

	// Reduce to row echelon form, swapping each pivot up so that only the
	// rows below it need clearing. Words before the pivot's are already
	// zero in every row below.
	pivots := 0
	for c := 0; c < numCols && pivots < numRows; c++ {
		w, bit := c/64, uint64(1)<<(uint(c)%64)
		p := pivots
		for p < numRows && vec[p][w]&bit == 0 {
			p++
		}
		if p == numRows {
			continue
		}
		vec[pivots], vec[p] = vec[p], vec[pivots]
		combo[pivots], combo[p] = combo[p], combo[pivots]

		pv, pc := vec[pivots], combo[pivots]
		for i := pivots + 1; i < numRows; i++ {
			if vec[i][w]&bit == 0 {
				continue
			}
			for k := w; k < colWords; k++ {
				vec[i][k] ^= pv[k]
			}
			for k := range pc {
				combo[i][k] ^= pc[k]
			}
		}
		pivots++
	}

	// The rows below the last pivot are now all zero.
	return combo[pivots:]
}

// tryDependency combines the relations selected by dep into a congruence of
// squares a² ≡ b² (mod n) and returns a non-trivial factor of n, or nil.
func tryDependency(n *big.Int, factorBase []int64, relations []relation, dep []uint64) *big.Int {
	a := big.NewInt(1)
	exps := make([]int, len(factorBase)+1)
	for w, word := range dep {
		for word != 0 {
			i := w*64 + bits.TrailingZeros64(word)
			word &= word - 1
			a.Mul(a, relations[i].X)
			a.Mod(a, n)
			for _, c := range relations[i].Factors {
				exps[c]++
			}
		}
	}

	// b = √(∏ Y_i). The sign exponent is even, so the product is positive.
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

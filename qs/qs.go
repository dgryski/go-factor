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
	"runtime"
	"sort"
	"sync"
	"sync/atomic"
	"time"

	"github.com/dgryski/go-factor/internal/progress"
)

var one = big.NewInt(1)

// extraRelations is how many relations to collect beyond the number of
// matrix columns. Each extra relation guarantees at least one more
// dependency, and each dependency has roughly a 50% chance of giving a
// non-trivial factor.
const extraRelations = 20

const (
	blockShift    = 15 // sieve blocks are 32 KB, sized for L1 cache
	blockSize     = 1 << blockShift
	bucketCut     = 16384 // primes from here up are sieved through per-block buckets
	resieveCut    = 2048  // primes from here to bucketCut are found again by resieving
	minSievePrime = 30    // smaller primes are not sieved, only allowed for in the threshold

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
	setup := newSetupLog()

	// Remove small factors first: primes that divide n would otherwise be
	// silently dropped from the base, and the multiplier must be coprime to n.
	primes := primesUpTo(b)
	setup.step("found %d primes up to %d", len(primes), b)
	var r, bigP big.Int
	for _, p := range primes {
		bigP.SetInt64(p)
		if r.Mod(n, &bigP).Sign() == 0 {
			return new(big.Int).Set(&bigP), new(big.Int).Div(n, &bigP)
		}
	}
	setup.step("trial divided by them")

	k := knuthSchroeppel(n)
	kn := new(big.Int).Mul(n, big.NewInt(k))
	factorBase := buildFactorBase(kn, primes)
	setup.step("built a factor base of %d primes", len(factorBase))
	s := newSieve(n, kn, factorBase, siqs)
	setup.step("found the square roots of kN mod each prime")

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

	primes = s.finish(factorBase)
	numCols := len(primes) + 1
	if len(s.relations) <= numCols {
		logf("qs: not enough relations, giving up")
		return nil, nil // not enough relations to guarantee a dependency
	}
	if SaveDir != "" {
		saveRelations(n, primes, s.relations)
	}
	start := time.Now()
	deps := findDependencies(s.relations, numCols)
	logf("qs: matrix %d×%d, %d dependencies in %v",
		len(s.relations), numCols, len(deps), time.Since(start).Round(time.Millisecond))
	for i, dep := range deps {
		if f := tryDependency(n, primes, s.relations, dep); f != nil {
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

// setupLog times the steps before sieving. Small inputs set up in well
// under a second, so it only logs once the setup has run for 5 seconds.
type setupLog struct{ start, last time.Time }

func newSetupLog() *setupLog {
	now := time.Now()
	return &setupLog{now, now}
}

// step logs what the step that just finished did and how long it took.
func (l *setupLog) step(format string, args ...interface{}) {
	now := time.Now()
	if Logf != nil && now.Sub(l.start) > 5*time.Second {
		logf("qs: setup: "+format+" in %v, %v so far", append(args,
			now.Sub(l.last).Round(time.Millisecond), now.Sub(l.start).Round(time.Millisecond))...)
	}
	l.last = now
}

// logProgress reports how many relations the sieve has collected and how
// long the rest should take at the recent rate. The first message comes
// after the first batch, and the time before it is mostly setup, so it
// only sets the rate's starting point.
func (s *sieve) logProgress() {
	now := time.Now()
	have := len(s.relations) + s.partialExcess
	s.rate.Sample(now, have)
	s.lastLog = now

	left := ""
	if have < s.want {
		left = ", " + s.rate.Left(have, s.want)
	}
	logf("qs: %d/%d relations (%d full, %d from partials), %d partials, %d polynomials, %v%s",
		have, s.want, len(s.relations), s.partialExcess, len(s.partials), s.polys,
		now.Sub(s.start).Round(time.Millisecond), left)
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
		sieveRange = largeSieveRange
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

// log2Big returns log2(v) for v > 0. It reads v's top 53 bits from its
// words, so it doesn't allocate on 64-bit platforms.
func log2Big(v *big.Int) float64 {
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
	roots     []root // the roots before any polynomial, copied by each worker
	numCols   int
	want      int
	sieveFrom int     // index of the first sieved prime
	slack     float64 // how far below log2|Y| a candidate may score
	lpBound   uint64
	pMax      uint64
	dlpCap    uint64 // cofactors up to this are split into two large primes

	// Primes from bucketFrom up are sieved by walking each polynomial's
	// interval once and filing the hits by block. A bucket entry packs
	// log2(p) << 48 | prime index << 16 | offset in the block. Primes from
	// resieveFrom to bucketFrom are sieved block by block, and found again
	// for the block's candidates by walking back from their next offsets.
	resieveFrom int
	bucketFrom  int

	// The rest is only used by the goroutine collecting relations.
	relations []relation
	seenX     map[string]bool // min(X, n−X) for every relation and partial
	key       big.Int
	keyBuf    []byte

	// Partial relations have one or two large primes up to lpBound left over.
	// Large primes that occur in more than one partial become matrix columns
	// in finish.
	partials      []partial
	partialExcess int // surviving partials minus large primes at the last count
	nextCount     int // len(partials) at which to count again

	// For progress logging.
	polys   int // polynomials sieved
	start   time.Time
	lastLog time.Time
	rate    progress.Rate // relations per second
}

// worker sieves polynomials. Each worker has its own copy of the roots and
// its own buffers, so several can run at once. It passes the relations it
// finds to emit, and stops early once stop is set.
type worker struct {
	s     *sieve // read only
	roots []root
	emit  func(found)
	stop  *atomic.Bool
	poly  int // index of the polynomial being sieved within its A

	aInv    []uint32
	deltas  [][]uint32
	buf     [blockSize]byte
	buckets [][]uint64
	small   []smallRoot // primes below resieveCut not dividing a

	// The current block's candidates, and the primes from resieveFrom up
	// hitting each.
	cands    []int
	candHits [][]int32
	candBits [blockSize / 64]uint64
	candSlot [blockSize]uint16 // index in cands of each candidate's offset

	scratch []int32
	temp    big.Int
}

func (s *sieve) newWorker(stop *atomic.Bool, emit func(found)) *worker {
	return &worker{
		s:     s,
		roots: append([]root(nil), s.roots...),
		emit:  emit,
		stop:  stop,
		aInv:  make([]uint32, len(s.roots)),
	}
}

// found is a relation from a worker, not yet checked against the others.
type found struct {
	x        *big.Int // reduced mod n
	factors  []int32
	lp1, lp2 uint64
	poly     int // index of its polynomial within its A
}

// Workers is the most goroutines that sieve SIQS polynomials. 0 means
// runtime.GOMAXPROCS(0).
var Workers = 0

// workers returns how many goroutines to sieve with. Small inputs need too
// few polynomials to keep many workers busy, so extra ones mostly sieve past
// the point where sieving stops. The limits are from timings on a 15-core
// laptop.
func (s *sieve) workers() int {
	n := runtime.GOMAXPROCS(0)
	if Workers > 0 {
		n = Workers
	}
	limit := n
	switch bits := s.n.BitLen(); {
	case bits < 80:
		limit = 1
	case bits < 130:
		limit = 4
	case bits < 170:
		limit = 8
	}
	if n > limit {
		n = limit
	}
	return n
}

func newSieve(n, kn *big.Int, factorBase []int64, siqs bool) *sieve {
	s := &sieve{
		n:       n,
		kn:      kn,
		log2KN:  log2Big(kn),
		roots:   make([]root, len(factorBase)),
		numCols: len(factorBase) + 1, // +1 for the sign of Y
		seenX:   make(map[string]bool),
		start:   time.Now(),
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

	s.resieveFrom = sort.Search(len(s.roots), func(i int) bool { return s.roots[i].p >= resieveCut })
	if s.resieveFrom < s.sieveFrom {
		s.resieveFrom = s.sieveFrom
	}
	s.bucketFrom = sort.Search(len(s.roots), func(i int) bool { return s.roots[i].p >= bucketCut })
	if s.bucketFrom < s.resieveFrom {
		s.bucketFrom = s.resieveFrom
	}

	pMax := uint64(s.roots[len(s.roots)-1].p)
	s.pMax = pMax
	s.slack = singleSlack*math.Log2(float64(pMax)) + smallSlack
	if siqs {
		s.lpBound = largePrimeMult * pMax
		s.slack = siqsSlack*math.Log2(float64(pMax)) + math.Log2(float64(s.lpBound))/2 + smallSlack
		if n.BitLen() >= dlpMinBits {
			s.dlpCap = s.lpBound * s.lpBound
			if s.dlpCap > dlpMaxCofactor {
				s.dlpCap = dlpMaxCofactor
			}
			s.slack = siqsSlack*math.Log2(float64(pMax)) + dlpSlackFrac*math.Log2(float64(s.dlpCap)) + smallSlack
		}
	}
	return s
}

// singlePoly sieves Y(x) = (x + m)² − kN, with m = ⌈√kN⌉, for x in
// [−sieveRange, sieveRange). It returns a factor of n if it finds one by
// accident.
func (s *sieve) singlePoly(sieveRange int64) *big.Int {
	m := new(big.Int).Sqrt(s.kn)
	m.Add(m, one)
	var stop atomic.Bool
	var factor *big.Int
	w := s.newWorker(&stop, func(f found) {
		if factor = s.addRelation(f.x, f.factors, f.lp1, f.lp2); factor != nil || s.done() {
			stop.Store(true)
		}
	})
	w.setRoots(one, m)

	// Stop before X reaches zero, since X and −X give the same Y and would
	// only add duplicate rows (and so trivial dependencies).
	lo := -sieveRange
	if m.IsInt64() && 1-m.Int64() > lo {
		lo = 1 - m.Int64()
	}
	s.polys = 1
	w.sievePoly(one, m, nil, lo, sieveRange)
	return factor
}

// siqs sieves x in [−sieveRange, sieveRange) for each of a sequence of
// polynomials until it has enough relations or runs out of polynomials. It
// returns a factor of n if it finds one by accident.
//
// One goroutine chooses the A values, workers sieve all the polynomials for
// one A at a time, and the calling goroutine adds their relations in the
// order of the A values. The relations are the same for any number of
// workers. Workers only do extra work past the point where sieving stops.
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

	type job struct {
		seq  int
		a    *big.Int
		aIdx []int
	}
	type batch struct {
		seq, polys int
		rels       []found
	}
	nw := s.workers()
	jobs := make(chan job, nw)
	results := make(chan batch, nw)
	quit := make(chan struct{})
	var stop atomic.Bool

	go func() {
		defer close(jobs)
		rng := rand.New(rand.NewSource(1))
		seen := make(map[string]bool)
		for seq := 0; ; seq++ {
			a, aIdx := chooseA(rng, s.roots, pool, log2Target, seen)
			if a == nil {
				return
			}
			select {
			case jobs <- job{seq, a, aIdx}:
			case <-quit:
				return
			}
		}
	}()

	var wg sync.WaitGroup
	for i := 0; i < nw; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			var b batch
			w := s.newWorker(&stop, func(f found) { b.rels = append(b.rels, f) })
			for j := range jobs {
				b = batch{seq: j.seq}
				b.polys = w.sieveA(j.a, j.aIdx, sieveRange)
				results <- b
			}
		}()
	}
	go func() {
		wg.Wait()
		close(results)
	}()

	var factor *big.Int
	finished := false
	pending := make(map[int]batch)
	next := 0
	for b := range results {
		if finished {
			continue // drain until the workers exit
		}
		pending[b.seq] = b
		for !finished {
			b, ok := pending[next]
			if !ok {
				break
			}
			delete(pending, next)
			next++
			polys := b.polys
			for _, f := range b.rels {
				if factor = s.addRelation(f.x, f.factors, f.lp1, f.lp2); factor != nil || s.done() {
					polys, finished = f.poly+1, true
					break
				}
			}
			s.polys += polys
			if Logf != nil && time.Since(s.lastLog) > 5*time.Second {
				s.logProgress()
			}
		}
		if finished {
			stop.Store(true)
			close(quit)
		}
	}
	if !finished {
		logf("qs: ran out of polynomials")
	}
	return factor
}

// sieveA sieves the 2^(s−1) polynomials for a, whose primes are aIdx, in
// Gray code order, until they're done or w.stop is set. It returns how many
// it sieved.
func (w *worker) sieveA(a *big.Int, aIdx []int, sieveRange int64) int {
	if w.stop.Load() {
		return 0
	}
	for _, idx := range aIdx {
		w.roots[idx].inA = true
	}
	defer func() {
		for _, idx := range aIdx {
			w.roots[idx].inA = false
		}
	}()

	// B = Σ B_j with B_j ≡ ±t_j (mod q_j) and B_j ≡ 0 (mod A/q_j), so
	// B² ≡ kN (mod A). Flipping the signs of B_0..B_{s−2} gives
	// 2^(s−1) polynomials for this A.
	var bigP, r, aq, twoBv big.Int
	bj := make([]*big.Int, len(aIdx))
	b := new(big.Int)
	for j, idx := range aIdx {
		q := uint64(w.roots[idx].p)
		bigP.SetUint64(q)
		aq.Quo(a, &bigP)
		gamma := uint64(w.roots[idx].t) * invMod(r.Mod(&aq, &bigP).Uint64(), q) % q
		if gamma > q/2 {
			gamma = q - gamma
		}
		bj[j] = new(big.Int).Mul(new(big.Int).SetUint64(gamma), &aq)
		b.Add(b, bj[j])
	}

	// When the sign of B_v flips, every root moves by 2·B_v·A⁻¹ (mod p).
	for len(w.deltas) < len(aIdx)-1 {
		w.deltas = append(w.deltas, make([]uint32, len(w.roots)))
	}
	for i := range w.roots {
		rt := &w.roots[i]
		if rt.inA {
			continue
		}
		p := uint64(rt.p)
		bigP.SetUint64(p)
		w.aInv[i] = uint32(invMod(r.Mod(a, &bigP).Uint64(), p))
		for v := 0; v < len(aIdx)-1; v++ {
			w.deltas[v][i] = uint32(2 * r.Mod(bj[v], &bigP).Uint64() * uint64(w.aInv[i]) % p)
		}
	}
	w.setRootsInv(b, w.aInv)

	for poly := 1; ; poly++ {
		w.poly = poly - 1
		w.sievePoly(a, b, aIdx, -sieveRange, sieveRange)
		if poly == 1<<uint(len(aIdx)-1) || w.stop.Load() {
			return poly
		}

		// Gray code order: step poly flips the sign of B_v.
		v := bits.TrailingZeros(uint(poly))
		twoBv.Lsh(bj[v], 1)
		if poly>>uint(v+1)&1 == 0 {
			b.Sub(b, &twoBv)
			w.moveRoots(w.deltas[v], true)
		} else {
			b.Add(b, &twoBv)
			w.moveRoots(w.deltas[v], false)
		}
	}
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
func (w *worker) setRoots(a, b *big.Int) {
	var r, bigP big.Int
	for i := range w.roots {
		if !w.roots[i].inA {
			p := uint64(w.roots[i].p)
			w.aInv[i] = uint32(invMod(r.Mod(a, bigP.SetUint64(p)).Uint64(), p))
		}
	}
	w.setRootsInv(b, w.aInv)
}

// setRootsInv is setRoots with A⁻¹ mod p already computed for each prime.
func (w *worker) setRootsInv(b *big.Int, aInv []uint32) {
	var r, bigP big.Int
	for i := range w.roots {
		rt := &w.roots[i]
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
func (w *worker) moveRoots(delta []uint32, add bool) {
	for i := range w.roots {
		rt := &w.roots[i]
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

// sievePoly sieves Y(x) = ((ax + b)² − kN) / a for x in [lo, hi) and passes
// the relations it finds to w.emit. aIdx lists the primes dividing a.
func (w *worker) sievePoly(a, b *big.Int, aIdx []int, lo, hi int64) {
	s := w.s
	w.small = w.small[:0]
	for i := 0; i < s.resieveFrom; i++ {
		if rt := &w.roots[i]; !rt.inA {
			w.small = append(w.small, smallRoot{idx: int32(i), p: rt.p, r1: rt.r1, r2: rt.r2})
		}
	}
	for i := s.sieveFrom; i < s.bucketFrom; i++ {
		rt := &w.roots[i]
		rt.off1 = nextOffset(lo, rt.r1, rt.p)
		rt.off2 = nextOffset(lo, rt.r2, rt.p)
	}
	w.fillBuckets(lo, hi)

	var x, y big.Int
	for blockStart, bi := lo, 0; blockStart < hi && !w.stop.Load(); blockStart, bi = blockStart+blockSize, bi+1 {
		blockLen := int64(blockSize)
		if hi-blockStart < blockLen {
			blockLen = hi - blockStart
		}
		block := w.buf[:blockLen]
		for i := range block {
			block[i] = 0
		}

		for i := s.sieveFrom; i < s.bucketFrom; i++ {
			rt := &w.roots[i]
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
		bucket := w.buckets[bi]
		for _, e := range bucket {
			block[uint16(e)] += byte(e >> 48)
		}

		thresh := int(math.Ceil(w.maxLog2Y(a, b, blockStart, blockStart+blockLen) - s.slack))
		if thresh < 0 {
			thresh = 0
		}
		w.cands = w.cands[:0]
		for i, score := range block {
			if int(score) >= thresh {
				w.cands = append(w.cands, i)
			}
		}
		if len(w.cands) == 0 {
			continue
		}
		w.bucketHits(bucket, len(block))

		for k, i := range w.cands {
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
			factors, lp1, lp2, ok := w.tryFactorize(off, &y, aIdx, w.candHits[k])
			if !ok {
				continue
			}
			w.emit(found{new(big.Int).Mod(&x, s.n), factors, lp1, lp2, w.poly})
			if w.stop.Load() {
				break
			}
		}
	}
}

// bucketHits sets w.candHits to the primes from resieveFrom up hitting
// each candidate. It must run after the block of length blockLen is sieved,
// when off1 and off2 hold the offsets of the next block's first hits.
func (w *worker) bucketHits(bucket []uint64, blockLen int) {
	s := w.s
	for len(w.candHits) < len(w.cands) {
		w.candHits = append(w.candHits, nil)
	}
	for k, i := range w.cands {
		w.candBits[i/64] |= 1 << (uint(i) % 64)
		w.candSlot[i] = uint16(k)
		w.candHits[k] = w.candHits[k][:0]
	}
	for i := s.resieveFrom; i < s.bucketFrom; i++ {
		rt := &w.roots[i]
		if rt.inA {
			continue
		}
		p := int(rt.p)
		for j := int(rt.off1) + blockLen - p; j >= 0; j -= p {
			if w.candBits[uint(j)/64]&(1<<(uint(j)%64)) != 0 {
				k := w.candSlot[j]
				w.candHits[k] = append(w.candHits[k], int32(i))
			}
		}
		if rt.r2 != rt.r1 {
			for j := int(rt.off2) + blockLen - p; j >= 0; j -= p {
				if w.candBits[uint(j)/64]&(1<<(uint(j)%64)) != 0 {
					k := w.candSlot[j]
					w.candHits[k] = append(w.candHits[k], int32(i))
				}
			}
		}
	}
	for _, e := range bucket {
		if j := uint16(e); w.candBits[j/64]&(1<<(j%64)) != 0 {
			k := w.candSlot[j]
			w.candHits[k] = append(w.candHits[k], int32(uint32(e>>16)))
		}
	}
	for _, i := range w.cands {
		w.candBits[i/64] = 0
	}
}

// fillBuckets files every hit in [lo, hi) of the primes from bucketFrom up
// into the bucket of its block.
func (w *worker) fillBuckets(lo, hi int64) {
	s := w.s
	n := int((hi - lo + blockSize - 1) >> blockShift)
	for len(w.buckets) < n {
		w.buckets = append(w.buckets, nil)
	}
	for i := range w.buckets[:n] {
		w.buckets[i] = w.buckets[i][:0]
	}
	total := uint64(hi - lo)
	for i := s.bucketFrom; i < len(w.roots); i++ {
		rt := &w.roots[i]
		if rt.inA {
			continue
		}
		p := uint64(rt.p)
		tag := uint64(rt.logP)<<48 | uint64(i)<<16
		for off := uint64(nextOffset(lo, rt.r1, rt.p)); off < total; off += p {
			b := off >> blockShift
			w.buckets[b] = append(w.buckets[b], tag|off&(blockSize-1))
		}
		if rt.r2 != rt.r1 {
			for off := uint64(nextOffset(lo, rt.r2, rt.p)); off < total; off += p {
				b := off >> blockShift
				w.buckets[b] = append(w.buckets[b], tag|off&(blockSize-1))
			}
		}
	}
}

// maxLog2Y returns log2 of the largest |Y(x)| for x in [from, to). |Y| is
// largest either at one end, or where ax + b = 0 if that lies in between.
func (w *worker) maxLog2Y(a, b *big.Int, from, to int64) float64 {
	s := w.s
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

// tryFactorize attempts to factor y = Y(x) over the factor base. It divides
// by the small primes whose roots match x and by the larger primes in hits.
// A relation's Y is a·Y(x), so each prime of a is counted once more. A
// cofactor made of one or two primes up to lpBound is returned as lp1, lp2.
func (w *worker) tryFactorize(x int64, y *big.Int, aIdx []int, hits []int32) (factors []int32, lp1, lp2 uint64, ok bool) {
	s := w.s
	factors = w.scratch[:0]
	temp := w.temp.Abs(y)
	if y.Sign() < 0 {
		factors = append(factors, 0)
	}

	var bigP, q, r big.Int
	divide := func(i int32) {
		bigP.SetUint64(uint64(w.roots[i].p))
		for {
			q.QuoRem(temp, &bigP, &r)
			if r.Sign() != 0 {
				break
			}
			factors = append(factors, i+1)
			temp.Set(&q)
		}
	}

	for _, i := range aIdx {
		factors = append(factors, int32(i+1))
		divide(int32(i))
	}
	for i := range w.small {
		sr := &w.small[i]
		xr := x % int64(sr.p)
		if xr < 0 {
			xr += int64(sr.p)
		}
		if uint32(xr) == sr.r1 || uint32(xr) == sr.r2 {
			divide(sr.idx)
		}
	}
	for _, i := range hits {
		divide(i)
	}
	w.scratch = factors[:0]

	// Every prime up to the largest base prime has been divided out, so a
	// cofactor below its square is prime.
	if temp.Cmp(one) != 0 {
		if !temp.IsUint64() {
			return nil, 0, 0, false
		}
		c := temp.Uint64()
		switch {
		case c <= s.lpBound && (c < s.pMax*s.pMax || isPrime64(c)):
			lp1 = c
		case c <= s.dlpCap:
			if lp1, lp2, ok = splitCofactor(c, s.pMax, s.lpBound); !ok {
				return nil, 0, 0, false
			}
		default:
			return nil, 0, 0, false
		}
	}
	return append([]int32(nil), factors...), lp1, lp2, true
}

// smallRoot packs what tryFactorize needs to test a small prime.
type smallRoot struct {
	idx       int32
	p, r1, r2 uint32
}

// addRelation records a full relation, or a partial one with large primes
// lp1 and lp2 (0 if absent). The relation keeps x, and may change it. It
// returns a factor of n if a large prime divides n.
func (s *sieve) addRelation(x *big.Int, factors []int32, lp1, lp2 uint64) *big.Int {
	// Every polynomial gives Y = X² − kN, so different polynomials can find
	// the same X, and −X gives the same Y. Either way the copy only adds a
	// trivial dependency.
	key := s.key.Lsh(x, 1)
	if key.Cmp(s.n) > 0 {
		key.Sub(s.n, x)
	} else {
		key.Set(x)
	}
	kb := s.keyBuf[:0]
	if n := (key.BitLen() + 7) / 8; cap(kb) >= n {
		kb = key.FillBytes(kb[:n])
	} else {
		kb = key.Bytes()
	}
	s.keyBuf = kb
	if s.seenX[string(kb)] {
		return nil
	}
	s.seenX[string(kb)] = true

	var r, l big.Int
	for _, lp := range [2]uint64{lp1, lp2} {
		if lp != 0 && r.Mod(s.n, l.SetUint64(lp)).Sign() == 0 {
			return new(big.Int).SetUint64(lp)
		}
	}
	if lp1 != 0 && lp1 == lp2 {
		// X² ≡ Y·L² gives (X/L)² ≡ Y.
		x.Mul(x, new(big.Int).ModInverse(l.SetUint64(lp1), s.n))
		x.Mod(x, s.n)
		lp1, lp2 = 0, 0
	}
	if lp1 == 0 {
		s.relations = append(s.relations, relation{X: x, Factors: factors})
		return nil
	}
	s.partials = append(s.partials, partial{relation{X: x, Factors: factors}, lp1, lp2})
	if len(s.partials) >= s.nextCount {
		s.partialExcess, _, _ = s.countPartials()
		s.nextCount = len(s.partials) + len(s.partials)/20 + 1000
	}
	return nil
}

// done reports whether the relations found so far, counting the partials
// that survive pruning, give a matrix with at least s.want surplus rows.
func (s *sieve) done() bool {
	return len(s.relations)+s.partialExcess >= s.want
}

// partial is a relation X² ≡ Y·lp1·lp2 (mod n), with lp2 = 0 if absent.
type partial struct {
	relation
	lp1, lp2 uint64
}

// countPartials repeatedly drops the partials containing a large prime that
// no other surviving partial has, since they can't be in any dependency. It
// returns the surviving partials minus the large primes they use, which
// partials survive, and the column of each surviving large prime counted
// from 0.
func (s *sieve) countPartials() (excess int, alive []bool, cols map[uint64]int32) {
	degree := make(map[uint64]int, len(s.partials))
	byPrime := make(map[uint64][]int32, len(s.partials))
	for i, p := range s.partials {
		for _, lp := range [2]uint64{p.lp1, p.lp2} {
			if lp != 0 {
				degree[lp]++
				byPrime[lp] = append(byPrime[lp], int32(i))
			}
		}
	}
	alive = make([]bool, len(s.partials))
	for i := range alive {
		alive[i] = true
	}
	var queue []uint64
	for lp, d := range degree {
		if d == 1 {
			queue = append(queue, lp)
		}
	}
	for len(queue) > 0 {
		lp := queue[len(queue)-1]
		queue = queue[:len(queue)-1]
		if degree[lp] != 1 {
			continue
		}
		for _, i := range byPrime[lp] {
			if !alive[i] {
				continue
			}
			alive[i] = false
			for _, q := range [2]uint64{s.partials[i].lp1, s.partials[i].lp2} {
				if q != 0 {
					if degree[q]--; degree[q] == 1 {
						queue = append(queue, q)
					}
				}
			}
		}
	}
	cols = make(map[uint64]int32)
	for _, lp := range sortedKeys(degree) {
		if degree[lp] > 0 {
			cols[lp] = int32(len(cols))
		}
	}
	for _, a := range alive {
		if a {
			excess++
		}
	}
	return excess - len(cols), alive, cols
}

func sortedKeys(m map[uint64]int) []uint64 {
	keys := make([]uint64, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}
	sort.Slice(keys, func(i, j int) bool { return keys[i] < keys[j] })
	return keys
}

// finish adds the surviving partials to s.relations, with a column for each
// large prime they use, and returns the prime of each column after the sign:
// the factor base, then the large primes.
func (s *sieve) finish(factorBase []int64) []int64 {
	_, alive, cols := s.countPartials()
	primes := make([]int64, len(factorBase)+len(cols))
	copy(primes, factorBase)
	for lp, c := range cols {
		primes[len(factorBase)+int(c)] = int64(lp)
	}
	base := int32(len(factorBase) + 1)
	for i, p := range s.partials {
		if !alive[i] {
			continue
		}
		f := append(p.Factors, base+cols[p.lp1])
		if p.lp2 != 0 {
			f = append(f, base+cols[p.lp2])
		}
		s.relations = append(s.relations, relation{X: p.X, Factors: f})
	}
	s.partials = nil
	return primes
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

// largeSieveRange is the SIQS sieve half-width for factor bases from 100,000
// primes up. It's a var while it's being tuned.
var largeSieveRange int64 = 128 * 1024

// Double large primes are still being tuned, so they're off by default. Set
// dlpMinBits to turn them on for inputs of at least that many bits.
var (
	dlpMinBits            = 1 << 30
	dlpSlackFrac          = 0.5 // threshold allowance, as a fraction of log2(dlpCap)
	dlpMaxCofactor uint64 = 1 << 50
)

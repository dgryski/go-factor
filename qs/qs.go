// Package qs implements the quadratic sieve.
//
// Small inputs sieve the single polynomial (x + ⌈√kN⌉)² − kN. Larger inputs
// use the self-initializing quadratic sieve (SIQS), which switches between
// many polynomials (Ax + B)² − kN so that the values being sieved stay small.
// Both use a Knuth–Schroeppel multiplier k and accept relations with one
// large prime.
package qs

import (
	"encoding/binary"
	"math"
	"math/big"
	"math/bits"
	"math/rand"
	"runtime"
	"sort"
	"sync"
	"sync/atomic"
	"time"

	"github.com/dgryski/go-factor/internal/arith"
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
	blockMask     = blockSize - 1
	bucketCut     = 16384 // primes from here up are sieved through per-block buckets
	resieveCut     = 2048 // primes from here to bucketCut are found again by resieving
	minSievePrime  = 30   // smaller primes are not sieved, only allowed for in the threshold
	largePrimeMult = 64

	// A candidate may score this many times log2(largest base prime) below
	// log2|Y|.
	singleSlack = 1.75
	siqsSlack   = 1.25

	siqsMinBits = 60   // below this, a single polynomial is faster
	siqsMinB    = 1000 // smaller bases leave too few primes to build enough different A values
)

// QuadraticSieve factors n by finding congruences of squares X² ≡ Y² (mod n).
func QuadraticSieve(n *big.Int) (primes, composites []*big.Int) {
	if n.ProbablyPrime(10) {
		return []*big.Int{n}, nil
	}

	parts := splitAll(n)
	if len(parts) < 2 {
		return nil, []*big.Int{n}
	}
	for _, f := range parts {
		pr, co := QuadraticSieve(f)
		primes, composites = append(primes, pr...), append(composites, co...)
	}
	return primes, composites
}

// split attempts to split n into two non-trivial factors. It returns
// nil, nil if n is prime, too small, or no factor was found.
func split(n *big.Int) (*big.Int, *big.Int) {
	parts := splitAll(n)
	if len(parts) < 2 {
		return nil, nil
	}
	return parts[0], new(big.Int).Div(n, parts[0])
}

// maxSieveRetries is how many times splitAll sieves more relations when
// every dependency is trivial.
const maxSieveRetries = 3

// splitAll returns non-trivial factors of n whose product is n, refined
// across every dependency, or nil if none was found.
func splitAll(n *big.Int) []*big.Int {
	if n.Cmp(big.NewInt(4)) < 0 || n.ProbablyPrime(20) {
		return nil
	}

	// Perfect squares would make some Y zero.
	sqrtN := new(big.Int).Sqrt(n)
	if new(big.Int).Mul(sqrtN, sqrtN).Cmp(n) == 0 {
		return []*big.Int{sqrtN, new(big.Int).Set(sqrtN)}
	}

	b, siqs, sieveRange := params(n)
	setup := newSetupLog()

	// Remove small factors first so the multiplier is coprime to n.
	primes := arith.PrimesUpTo(b)
	setup.step("found %d primes up to %d", len(primes), b)
	var r, bigP big.Int
	for _, p := range primes {
		bigP.SetInt64(p)
		if r.Mod(n, &bigP).Sign() == 0 {
			return []*big.Int{new(big.Int).Set(&bigP), new(big.Int).Div(n, &bigP)}
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

	for try := 0; ; try++ {
		var f *big.Int
		if siqs {
			f = s.siqs(sieveRange)
		} else {
			f = s.singlePoly(sieveRange)
		}
		s.logProgress()
		if f != nil {
			logf("qs: found factor %v while sieving", f)
			return []*big.Int{f, new(big.Int).Div(n, f)}
		}

		full, partials := len(s.relations), s.partials
		cols := s.finish(factorBase)
		numCols := len(cols) + 1
		if len(s.relations) <= numCols {
			logf("qs: not enough relations, giving up")
			return nil // not enough relations to guarantee a dependency
		}
		if SaveDir != "" {
			saveRelations(n, cols, s.relations)
		}
		start := time.Now()
		deps := findDependencies(s.relations, numCols)
		logf("qs: matrix %d×%d, %d dependencies in %v",
			len(s.relations), numCols, len(deps), time.Since(start).Round(time.Millisecond))
		if parts := refineDependencies(n, cols, s.relations, deps); parts != nil {
			return parts
		}
		if !siqs || try >= maxSieveRetries {
			logf("qs: all dependencies were trivial")
			return nil
		}

		// Restore the relations and partials from before finish and sieve more.
		s.relations, s.partials = s.relations[:full], partials
		extra := s.want / 50
		if extra < 100 {
			extra = 100
		}
		s.want += extra
		s.nextCount = len(s.partials)
		logf("qs: all dependencies were trivial, sieving %d more relations", extra)
	}
}

// refineDependencies splits n across all dependencies until every part is prime.
func refineDependencies(n *big.Int, cols []int64, relations []relation, deps [][]uint64) []*big.Int {
	parts := []*big.Int{n}
	for i, dep := range deps {
		g := tryDependency(n, cols, relations, dep)
		if g == nil {
			continue
		}
		refined := false
		var d big.Int
		for j := 0; j < len(parts); j++ {
			f := parts[j]
			d.GCD(nil, nil, f, g)
			if d.Cmp(one) > 0 && d.Cmp(f) < 0 {
				parts[j] = new(big.Int).Set(&d)
				parts = append(parts, new(big.Int).Div(f, &d))
				refined = true
			}
		}
		if !refined {
			continue
		}
		logf("qs: dependency %d gave factor %v", i+1, g)
		allPrime := true
		for _, f := range parts {
			if !f.ProbablyPrime(20) {
				allPrime = false
				break
			}
		}
		if allPrime {
			break
		}
	}
	if len(parts) < 2 {
		return nil
	}
	return parts
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
	smallPrimes := arith.PrimesUpTo(97)[1:]
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
			} else if arith.PowMod64(knp, uint64(p-1)/2, uint64(p)) == 1 {
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

// buildFactorBase returns the primes in primes for which kn is a square
// mod p: 2, the primes dividing the multiplier, and the odd primes with
// Legendre symbol 1.
func buildFactorBase(kn *big.Int, primes []int64) []int64 {
	var base []int64
	var r, bigP big.Int
	for _, p := range primes {
		bigP.SetInt64(p)
		knp := r.Mod(kn, &bigP).Uint64()
		if p == 2 || knp == 0 || arith.PowMod64(knp, uint64(p-1)/2, uint64(p)) == 1 {
			base = append(base, p)
		}
	}
	return base
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
	inv2ModN  *big.Int
	log2KN    float64
	q2        bool   // sieve ((2Ax + B)² − kN) / (4A) with odd B when kN ≡ 1 (mod 8)
	roots     []root // the roots before any polynomial, copied by each worker
	numCols   int
	want      int
	sieveFrom int     // index of the first sieved prime
	slack     float64 // how far below log2|Y| a candidate may score
	lpBound   uint64
	pMax      uint64
	dlpCap    uint64 // cofactors up to this are split into two large primes

	// Primes from bucketFrom up are sieved by walking each polynomial's
	// interval once and filing hits into a flat uint32 arena grouped by logP.
	resieveFrom   int
	bucketFrom    int
	bucketGroups  []logGroup
	expectedBKHits int

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

	// A value generator state, kept so sieving can resume.
	aRng  *rand.Rand
	aSeen map[string]bool
}

type bkRoot struct {
	p, r1, r2 uint32
}

// worker sieves polynomials. Each worker has its own copy of the roots and
// its own buffers, so several can run at once. It passes the relations it
// finds to emit, and stops early once stop is set.
type worker struct {
	s     *sieve // read only
	roots []root
	bk    []bkRoot
	emit  func(found)
	stop  *atomic.Bool
	poly  int // index of the polynomial being sieved within its A

	deltas    [][]uint32
	buf       [blockSize]byte
	arena     *qsArena
	log2Denom float64
	twoA      big.Int
	fourA     big.Int

	// The current block's candidates, and the primes from resieveFrom up
	// hitting each.
	cands    []int
	candHits [][]int32
	candBits [blockSize / 64]uint64
	candSlot [blockSize]uint16 // index in cands of each candidate's offset

	scratch            []int32
	temp, sx, sy       big.Int
	mx, my             big.Int
}

func (s *sieve) newWorker(stop *atomic.Bool, emit func(found)) *worker {
	bk := make([]bkRoot, len(s.roots)-s.bucketFrom)
	for i := s.bucketFrom; i < len(s.roots); i++ {
		bk[i-s.bucketFrom].p = s.roots[i].p
	}
	return &worker{
		s:     s,
		roots: append([]root(nil), s.roots...),
		bk:    bk,
		emit:  emit,
		stop:  stop,
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
		log2KN:  arith.Log2Big(kn),
		roots:   make([]root, len(factorBase)),
		numCols: len(factorBase) + 1, // +1 for the sign of Y
		seenX:   make(map[string]bool),
		start:   time.Now(),
	}
	s.want = s.numCols + extraRelations
	if siqs && len(factorBase) >= 16 && n.Bit(0) == 1 && kn.Bit(0) == 1 && kn.Bit(1) == 0 && kn.Bit(2) == 0 {
		s.q2 = true
		s.inv2ModN = new(big.Int).Add(n, one)
		s.inv2ModN.Rsh(s.inv2ModN, 1)
	}

	var r, bigP big.Int
	for i, p := range factorBase {
		bigP.SetInt64(p)
		s.roots[i] = root{
			p:    uint32(p),
			t:    arith.SqrtModPrime(uint32(r.Mod(kn, &bigP).Uint64()), uint32(p)),
			logP: byte(math.Round(math.Log2(float64(p)))),
		}
	}

	// Skip the small primes, which hit often but add little, and lower the
	// threshold by their expected contribution of about 2·log2(p)/p.
	smallSlack := 0.0
	for s.sieveFrom < len(s.roots) && s.roots[s.sieveFrom].p < minSievePrime {
		r := s.roots[s.sieveFrom]
		if r.p == 2 && s.q2 {
			smallSlack += 1.0
			s.sieveFrom++
			continue
		}
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
	s.bucketGroups = buildLogGroups(s.roots, s.bucketFrom)
	var invSum float64
	for i := s.bucketFrom; i < len(s.roots); i++ {
		invSum += 2.0 / float64(s.roots[i].p)
	}
	s.expectedBKHits = int(float64(blockSize)*invSum) + 256

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
	// Stop before X reaches zero, since X and −X give the same Y and would
	// only add duplicate rows (and so trivial dependencies).
	lo := -sieveRange
	if m.IsInt64() && 1-m.Int64() > lo {
		lo = 1 - m.Int64()
	}
	w.log2Denom = 0
	w.setRoots(one, m, lo)
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
	// about the same, and Q2 polynomials halve the target A.
	log2Target := (s.log2KN+1)/2 - math.Log2(float64(sieveRange))
	if s.q2 {
		log2Target -= 1
	}

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

	if s.aRng == nil {
		s.aRng = rand.New(rand.NewSource(1))
		s.aSeen = make(map[string]bool)
	}
	go func() {
		defer close(jobs)
		rng, seen := s.aRng, s.aSeen
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
	sLen := len(aIdx)
	var bigP, r, aq, twoBv big.Int
	bj := make([]*big.Int, sLen)
	qSlice := make([]uint64, sLen)
	gammaSlice := make([]uint64, sLen)
	b := new(big.Int)
	q2 := w.s.q2
	for j, idx := range aIdx {
		q := uint64(w.roots[idx].p)
		qSlice[j] = q
		bigP.SetUint64(q)
		aq.Quo(a, &bigP)
		gamma := uint64(w.roots[idx].t) * uint64(invMod32(uint32(r.Mod(&aq, &bigP).Uint64()), uint32(q))) % q
		if gamma > q/2 {
			gamma = q - gamma
		}
		if q2 && j == sLen-1 && (b.Bit(0)^uint(gamma&1)) == 0 {
			gamma = q - gamma
		}
		gammaSlice[j] = gamma
		bj[j] = new(big.Int).Mul(new(big.Int).SetUint64(gamma), &aq)
		b.Add(b, bj[j])
	}

	for len(w.deltas) < sLen-1 {
		w.deltas = append(w.deltas, make([]uint32, len(w.roots)))
	}
	var qmBuf, prefBuf [32]uint64
	qm := qmBuf[:0]
	pref := prefBuf[:0]
	if sLen <= len(qmBuf) {
		qm = qmBuf[:sLen]
		pref = prefBuf[:sLen]
	} else {
		qm = make([]uint64, sLen)
		pref = make([]uint64, sLen)
	}

	// Batch invert q_j mod p via prefix products since B_j*A^-1 ≡ gamma_j*q_j^-1 (mod p).
	w.log2Denom = arith.Log2Big(a)
	if q2 {
		w.log2Denom += 2
		w.twoA.Lsh(a, 1)
		w.fourA.Lsh(a, 2)
	}
	bkFrom := w.s.bucketFrom
	for _, idx := range aIdx {
		for v := 0; v < sLen-1; v++ {
			w.deltas[v][idx] = 0
		}
		if idx >= bkFrom {
			w.bk[idx-bkFrom].r1 = math.MaxUint32
			w.bk[idx-bkFrom].r2 = math.MaxUint32
		}
	}
	uRange := uint64(sieveRange)
	for i := range w.roots {
		rt := &w.roots[i]
		if rt.inA {
			continue
		}
		if q2 && rt.p == 2 {
			rt.r1, rt.r2 = 0, 1
			for v := 0; v < sLen-1; v++ {
				w.deltas[v][i] = 0
			}
			continue
		}
		p := uint64(rt.p)
		for j := 0; j < sLen; j++ {
			q := qSlice[j]
			if q >= p {
				q %= p
			}
			qm[j] = q
			if j == 0 {
				pref[0] = q
			} else {
				pref[j] = pref[j-1] * q % p
			}
		}
		inv := uint64(invMod32(uint32(pref[sLen-1]), rt.p))
		aInv := inv
		var sumBeta uint64
		for j := sLen - 1; j >= 0; j-- {
			var qInv uint64
			if j > 0 {
				qInv = inv * pref[j-1] % p
				inv = inv * qm[j] % p
			} else {
				qInv = inv
			}
			g := gammaSlice[j]
			if g >= p {
				g %= p
			}
			beta := g * qInv % p
			sumBeta += beta
			if j < sLen-1 {
				d := beta
				if !q2 {
					d <<= 1
					if d >= p {
						d -= p
					}
				}
				w.deltas[j][i] = uint32(d)
			}
		}
		sumBeta %= p
		shift := uRange
		if shift >= p {
			shift %= p
		}
		if q2 {
			halfP := (p + 1) >> 1
			aInv = aInv * halfP % p
			sumBeta = sumBeta * halfP % p
		}
		tA := uint64(rt.t) * aInv % p
		r1 := uint32((tA + p - sumBeta + shift) % p)
		r2 := uint32((2*p - tA - sumBeta + shift) % p)
		if i < bkFrom {
			rt.r1, rt.r2 = r1, r2
		} else {
			w.bk[i-bkFrom].r1 = r1
			w.bk[i-bkFrom].r2 = r2
		}
	}

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

func invMod32(a, p uint32) uint32 {
	if p >= 1<<31 {
		return uint32(arith.InvMod(uint64(a), uint64(p)))
	}
	t, newT := int32(0), int32(1)
	r, newR := int32(p), int32(a)
	for newR != 0 {
		q := r / newR
		t, newT = newT, t-q*newT
		r, newR = newR, r-q*newR
	}
	if t < 0 {
		t += int32(p)
	}
	return uint32(t)
}

func negMod(lo int64, p uint32) uint32 {
	m := (-lo) % int64(p)
	if m < 0 {
		m += int64(p)
	}
	return uint32(m)
}

// setRoots sets r1 and r2 of every prime not in A to its first offset >= 0 from lo.
func (w *worker) setRoots(a, b *big.Int, lo int64) {
	bkFrom := w.s.bucketFrom
	var r, bigP big.Int
	for i := range w.roots {
		rt := &w.roots[i]
		if rt.inA {
			if i >= bkFrom {
				w.bk[i-bkFrom].r1 = math.MaxUint32
				w.bk[i-bkFrom].r2 = math.MaxUint32
			}
			continue
		}
		p := uint64(rt.p)
		bigP.SetUint64(p)
		aInv := uint64(invMod32(uint32(r.Mod(a, &bigP).Uint64()), rt.p))
		bModP := r.Mod(b, &bigP).Uint64()
		shift := uint64(negMod(lo, rt.p))
		t := uint64(rt.t)
		r1 := uint32(((t+p-bModP)%p*aInv + shift) % p)
		r2 := uint32(((2*p-t-bModP)%p*aInv + shift) % p)
		if i < bkFrom {
			rt.r1, rt.r2 = r1, r2
		} else {
			w.bk[i-bkFrom].r1 = r1
			w.bk[i-bkFrom].r2 = r2
		}
	}
}

// moveRoots adds delta[i] to both roots of each prime, or subtracts it,
// using branchless modular reduction since p < 2^31.
func (w *worker) moveRoots(delta []uint32, add bool) {
	bkFrom := w.s.bucketFrom
	if add {
		for i := 0; i < bkFrom; i++ {
			rt := &w.roots[i]
			if rt.inA {
				continue
			}
			p := rt.p
			d := delta[i]
			d1 := int32(rt.r1 + d - p)
			rt.r1 = uint32(d1) + (p & uint32(d1>>31))
			d2 := int32(rt.r2 + d - p)
			rt.r2 = uint32(d2) + (p & uint32(d2>>31))
		}
		bks := w.bk
		del := delta[bkFrom:]
		for k := range bks {
			p := bks[k].p
			d := del[k]
			d1 := int32(bks[k].r1 + d - p)
			bks[k].r1 = uint32(d1) + (p & uint32(d1>>31))
			d2 := int32(bks[k].r2 + d - p)
			bks[k].r2 = uint32(d2) + (p & uint32(d2>>31))
		}
	} else {
		for i := 0; i < bkFrom; i++ {
			rt := &w.roots[i]
			if rt.inA {
				continue
			}
			p := rt.p
			d := delta[i]
			d1 := int32(rt.r1 - d)
			rt.r1 = uint32(d1) + (p & uint32(d1>>31))
			d2 := int32(rt.r2 - d)
			rt.r2 = uint32(d2) + (p & uint32(d2>>31))
		}
		bks := w.bk
		del := delta[bkFrom:]
		for k := range bks {
			p := bks[k].p
			d := del[k]
			d1 := int32(bks[k].r1 - d)
			bks[k].r1 = uint32(d1) + (p & uint32(d1>>31))
			d2 := int32(bks[k].r2 - d)
			bks[k].r2 = uint32(d2) + (p & uint32(d2>>31))
		}
	}
}

type logGroup struct {
	logP       byte
	start, end int
}

type overHit struct {
	idx  int32
	off  uint16
	logP byte
}

// qsArena stores bucket hits as packed uint32 entries (idx << blockShift | off)
// grouped by logP in a flat preallocated slab per block.
type qsArena struct {
	numBlocks int
	data      []uint32
	counts    []uint32
	groupEnds [][]uint32
	overflow  [][]overHit
	hasOver   bool
	stride    uint32
	shift     uint
}

func buildLogGroups(roots []root, from int) []logGroup {
	if from >= len(roots) {
		return nil
	}
	var out []logGroup
	curLog := roots[from].logP
	start := from
	for i := from + 1; i < len(roots); i++ {
		lg := roots[i].logP
		if lg != curLog {
			out = append(out, logGroup{logP: curLog, start: start, end: i})
			curLog = lg
			start = i
		}
	}
	return append(out, logGroup{logP: curLog, start: start, end: len(roots)})
}

func newQSArena(numBlocks, numGroups, expectedPerBlock int) *qsArena {
	want := expectedPerBlock * 3
	if want < 2048 {
		want = 2048
	}
	shift := uint(11)
	for (1 << shift) < want && shift < 17 {
		shift++
	}
	ge := make([][]uint32, numGroups)
	if numGroups > 0 && numBlocks > 0 {
		geBacking := make([]uint32, numGroups*numBlocks)
		for g := 0; g < numGroups; g++ {
			ge[g] = geBacking[g*numBlocks : (g+1)*numBlocks]
		}
	}
	return &qsArena{
		numBlocks: numBlocks,
		data:      make([]uint32, numBlocks<<shift),
		counts:    make([]uint32, numBlocks),
		groupEnds: ge,
		overflow:  make([][]overHit, numBlocks),
		stride:    uint32(1 << shift),
		shift:     shift,
	}
}

func (a *qsArena) reset(n int) {
	clear(a.counts[:n])
	if a.hasOver {
		for i := 0; i < n; i++ {
			a.overflow[i] = a.overflow[i][:0]
		}
		a.hasOver = false
	}
}

func (a *qsArena) addOverflow(bi int, idx int32, off uint16, logP byte) {
	a.hasOver = true
	a.overflow[bi] = append(a.overflow[bi], overHit{idx: idx, off: off, logP: logP})
}

func (a *qsArena) endGroup(g, n int) {
	copy(a.groupEnds[g][:n], a.counts[:n])
}

func (a *qsArena) drainBlock(bi int, block []byte, groups []logGroup) {
	base := bi << a.shift
	var cur uint32
	for g, grp := range groups {
		end := a.groupEnds[g][bi]
		if end == cur {
			continue
		}
		lg := grp.logP
		s := a.data[base+int(cur) : base+int(end)]
		i := 0
		for ; i+3 < len(s); i += 4 {
			block[s[i]&blockMask] += lg
			block[s[i+1]&blockMask] += lg
			block[s[i+2]&blockMask] += lg
			block[s[i+3]&blockMask] += lg
		}
		for ; i < len(s); i++ {
			block[s[i]&blockMask] += lg
		}
		cur = end
	}
	if a.hasOver {
		for _, ov := range a.overflow[bi] {
			block[ov.off] += ov.logP
		}
	}
}

func (a *qsArena) gatherHits(bi int, candBits *[blockSize / 64]uint64, candSlot *[blockSize]uint16, dst [][]int32) {
	base := bi << a.shift
	s := a.data[base : base+int(a.counts[bi])]
	for _, e := range s {
		j := e & blockMask
		if candBits[j>>6]&(1<<(j&63)) != 0 {
			k := candSlot[j]
			dst[k] = append(dst[k], int32(e>>blockShift))
		}
	}
	if a.hasOver {
		for _, ov := range a.overflow[bi] {
			if candBits[ov.off>>6]&(1<<(ov.off&63)) != 0 {
				k := candSlot[ov.off]
				dst[k] = append(dst[k], ov.idx)
			}
		}
	}
}

func addSieveStride(block []byte, j, p int, logP byte) uint32 {
	n := len(block)
	p2 := p << 1
	p3 := p2 + p
	p4 := p << 2
	for ; j+p3 < n; j += p4 {
		block[j] += logP
		block[j+p] += logP
		block[j+p2] += logP
		block[j+p3] += logP
	}
	for ; j < n; j += p {
		block[j] += logP
	}
	return uint32(j - n)
}

// appendGe appends indices off where buf[off] >= thresh, scanning 8 bytes per word.
func appendGe(dst []int, buf []byte, thresh int) []int {
	n := len(buf)
	if thresh > 255 || n == 0 {
		return dst
	}
	if thresh <= 0 {
		for off := 0; off < n; off++ {
			dst = append(dst, off)
		}
		return dst
	}
	const low7 = 0x7f7f7f7f7f7f7f7f
	const high1 = 0x8080808080808080
	n8 := n &^ 7
	if thresh <= 128 {
		bias := uint64(128-thresh) * 0x0101010101010101
		for off := 0; off < n8; off += 8 {
			w := binary.LittleEndian.Uint64(buf[off:])
			m := (((w & low7) + bias) | w) & high1
			for m != 0 {
				lane := bits.TrailingZeros64(m) >> 3
				dst = append(dst, off+lane)
				m &= m - 1
			}
		}
	} else {
		bias := uint64(256-thresh) * 0x0101010101010101
		for off := 0; off < n8; off += 8 {
			w := binary.LittleEndian.Uint64(buf[off:])
			m := (((w & low7) + bias) & w) & high1
			for m != 0 {
				lane := bits.TrailingZeros64(m) >> 3
				dst = append(dst, off+lane)
				m &= m - 1
			}
		}
	}
	for off := n8; off < n; off++ {
		if int(buf[off]) >= thresh {
			dst = append(dst, off)
		}
	}
	return dst
}

// sievePoly sieves Y(x) = ((ax + b)² − kN) / a for x in [lo, hi) and passes
// the relations it finds to w.emit. aIdx lists the primes dividing a.
func (w *worker) sievePoly(a, b *big.Int, aIdx []int, lo, hi int64) {
	s := w.s
	for i := s.sieveFrom; i < s.bucketFrom; i++ {
		rt := &w.roots[i]
		rt.off1 = rt.r1
		rt.off2 = rt.r2
	}
	w.fillBuckets(lo, hi, aIdx)

	stepA, denomA := a, a
	if s.q2 && aIdx != nil {
		stepA, denomA = &w.twoA, &w.fourA
	}
	x, y := &w.sx, &w.sy
	for blockStart, bi := lo, 0; blockStart < hi && !w.stop.Load(); blockStart, bi = blockStart+blockSize, bi+1 {
		blockLen := int64(blockSize)
		if hi-blockStart < blockLen {
			blockLen = hi - blockStart
		}
		block := w.buf[:blockLen]
		clear(block)

		for i := s.sieveFrom; i < s.bucketFrom; i++ {
			rt := &w.roots[i]
			if rt.inA {
				continue
			}
			p := int(rt.p)
			rt.off1 = addSieveStride(block, int(rt.off1), p, rt.logP)
			if rt.r2 != rt.r1 {
				rt.off2 = addSieveStride(block, int(rt.off2), p, rt.logP)
			}
		}
		if w.arena != nil {
			w.arena.drainBlock(bi, block, s.bucketGroups)
		}

		thresh := int(math.Ceil(w.maxLog2Y(stepA, b, blockStart, blockStart+blockLen) - s.slack))
		if thresh < 0 {
			thresh = 0
		}
		w.cands = appendGe(w.cands[:0], block, thresh)
		if len(w.cands) == 0 {
			continue
		}
		w.bucketHits(bi, len(block))

		for k, i := range w.cands {
			off := blockStart + int64(i)
			x.SetInt64(off)
			x.Mul(x, stepA)
			x.Add(x, b)
			y.Mul(x, x)
			y.Sub(y, s.kn)
			y.Quo(y, denomA)
			if y.Sign() == 0 {
				continue
			}
			factors, lp1, lp2, ok := w.tryFactorize(uint64(off-lo), y, aIdx, w.candHits[k])
			if !ok {
				continue
			}
			xMod := new(big.Int)
			if s.q2 && aIdx != nil {
				xMod.Mul(x, s.inv2ModN)
				xMod.Mod(xMod, s.n)
			} else {
				xMod.Mod(x, s.n)
			}
			w.emit(found{xMod, factors, lp1, lp2, w.poly})
			if w.stop.Load() {
				break
			}
		}
	}
}

// bucketHits sets w.candHits to the primes from resieveFrom up hitting
// each candidate. It must run after the block of length blockLen is sieved,
// when off1 and off2 hold the offsets of the next block's first hits.
func (w *worker) bucketHits(bi, blockLen int) {
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
	if w.arena != nil {
		w.arena.gatherHits(bi, &w.candBits, &w.candSlot, w.candHits)
	}
	for _, i := range w.cands {
		w.candBits[i/64] = 0
	}
}

// fillBuckets files every hit in [lo, hi) of the primes from bucketFrom up
// into the flat uint32 arena of its block.
func (w *worker) fillBuckets(lo, hi int64, aIdx []int) {
	s := w.s
	if len(s.bucketGroups) == 0 {
		return
	}
	n := int((hi - lo + blockSize - 1) >> blockShift)
	if w.arena == nil || w.arena.numBlocks < n {
		w.arena = newQSArena(n, len(s.bucketGroups), s.expectedBKHits)
	} else {
		w.arena.reset(n)
	}
	bkFrom := s.bucketFrom
	bks := w.bk
	for _, idx := range aIdx {
		if idx >= bkFrom {
			bks[idx-bkFrom].r1 = math.MaxUint32
			bks[idx-bkFrom].r2 = math.MaxUint32
		}
	}
	a := w.arena
	counts := a.counts
	data := a.data
	baseStride := a.stride
	shift := a.shift
	total := uint64(hi - lo)

	for g, grp := range s.bucketGroups {
		lg := grp.logP
		mid := grp.start + sort.Search(grp.end-grp.start, func(k int) bool {
			return uint64(bks[grp.start-bkFrom+k].p) >= total
		})
		for i := grp.start; i < mid; i++ {
			bk := &bks[i-bkFrom]
			p := uint64(bk.p)
			r1, r2 := bk.r1, bk.r2
			stride := baseStride
			if i >= 1<<(32-blockShift) {
				stride = 0
			}
			tag := uint32(i) << blockShift
			for off := uint64(r1); off < total; off += p {
				bi := int(off >> blockShift)
				bo := uint32(off & blockMask)
				if c := counts[bi]; c < stride {
					data[bi<<shift|int(c)] = tag | bo
					counts[bi] = c + 1
				} else {
					a.addOverflow(bi, int32(i), uint16(bo), lg)
				}
			}
			if r2 != r1 {
				for off := uint64(r2); off < total; off += p {
					bi := int(off >> blockShift)
					bo := uint32(off & blockMask)
					if c := counts[bi]; c < stride {
						data[bi<<shift|int(c)] = tag | bo
						counts[bi] = c + 1
					} else {
						a.addOverflow(bi, int32(i), uint16(bo), lg)
					}
				}
			}
		}
		for i := mid; i < grp.end; i++ {
			bk := &bks[i-bkFrom]
			r1, r2 := bk.r1, bk.r2
			stride := baseStride
			if i >= 1<<(32-blockShift) {
				stride = 0
			}
			tag := uint32(i) << blockShift
			if off := uint64(r1); off < total {
				bi := int(off >> blockShift)
				bo := uint32(off & blockMask)
				if c := counts[bi]; c < stride {
					data[bi<<shift|int(c)] = tag | bo
					counts[bi] = c + 1
				} else {
					a.addOverflow(bi, int32(i), uint16(bo), lg)
				}
			}
			if r2 != r1 {
				if off := uint64(r2); off < total {
					bi := int(off >> blockShift)
					bo := uint32(off & blockMask)
					if c := counts[bi]; c < stride {
						data[bi<<shift|int(c)] = tag | bo
						counts[bi] = c + 1
					} else {
						a.addOverflow(bi, int32(i), uint16(bo), lg)
					}
				}
			}
		}
		a.endGroup(g, n)
	}
}

// maxLog2Y returns log2 of the largest |Y(x)| for x in [from, to). |Y| is
// largest either at one end, or where stepA*x + b = 0 if that lies in between.
func (w *worker) maxLog2Y(stepA, b *big.Int, from, to int64) float64 {
	s := w.s
	x, y := &w.mx, &w.my
	best, signs := 0.0, 0
	for _, off := range [2]int64{from, to - 1} {
		x.SetInt64(off)
		x.Mul(x, stepA)
		x.Add(x, b)
		signs += x.Sign()
		y.Mul(x, x)
		y.Sub(y, s.kn)
		if y.Sign() != 0 {
			best = math.Max(best, arith.Log2Big(y.Abs(y)))
		}
	}
	if signs == 0 {
		best = math.Max(best, s.log2KN)
	}
	return best - w.log2Denom
}

// tryFactorize attempts to factor y = Y(x) over the factor base, where z = x - lo.
// A relation's Y is a·Y(x), so each prime of a is counted once more.
func (w *worker) tryFactorize(z uint64, y *big.Int, aIdx []int, hits []int32) (factors []int32, lp1, lp2 uint64, ok bool) {
	s := w.s
	factors = w.scratch[:0]
	if y.Sign() < 0 {
		factors = append(factors, 0)
	}

	if y.BitLen() > 384 {
		return w.tryFactorizeBig(z, y, aIdx, hits, factors)
	}

	var w32 [12]uint32
	n32 := 0
	words := y.Bits()
	if bits.UintSize == 64 {
		for idx, wd := range words {
			u := uint64(wd)
			w32[2*idx] = uint32(u)
			w32[2*idx+1] = uint32(u >> 32)
		}
		n32 = 2 * len(words)
		if n32 > 0 && w32[n32-1] == 0 {
			n32--
		}
	} else {
		for idx, wd := range words {
			w32[idx] = uint32(wd)
		}
		n32 = len(words)
	}

	in64 := false
	var v64 uint64
	if n32 <= 2 {
		in64 = true
		if n32 >= 1 {
			v64 = uint64(w32[0])
		}
		if n32 == 2 {
			v64 |= uint64(w32[1]) << 32
		}
	}

	divide := func(i int32) {
		p64 := uint64(w.roots[i].p)
		if in64 {
			for v64 >= p64 && v64%p64 == 0 {
				factors = append(factors, i+1)
				v64 /= p64
			}
			return
		}
		for {
			var rem uint64
			var q32 [12]uint32
			for k := n32 - 1; k >= 0; k-- {
				cur := (rem << 32) | uint64(w32[k])
				q32[k] = uint32(cur / p64)
				rem = cur % p64
			}
			if rem != 0 {
				return
			}
			factors = append(factors, i+1)
			for n32 > 0 && q32[n32-1] == 0 {
				n32--
			}
			if n32 <= 2 {
				in64 = true
				v64 = uint64(q32[0])
				if n32 == 2 {
					v64 |= uint64(q32[1]) << 32
				}
				for v64 >= p64 && v64%p64 == 0 {
					factors = append(factors, i+1)
					v64 /= p64
				}
				return
			}
			for k := 0; k < n32; k++ {
				w32[k] = q32[k]
			}
		}
	}

	for _, i := range aIdx {
		factors = append(factors, int32(i+1))
		divide(int32(i))
	}
	if z <= math.MaxUint32 {
		z32 := uint32(z)
		for i := 0; i < s.resieveFrom && !(in64 && v64 == 1); i++ {
			rt := &w.roots[i]
			if rt.inA {
				continue
			}
			zm := z32 % rt.p
			if zm == rt.r1 || zm == rt.r2 {
				divide(int32(i))
			}
		}
	} else {
		for i := 0; i < s.resieveFrom && !(in64 && v64 == 1); i++ {
			rt := &w.roots[i]
			if rt.inA {
				continue
			}
			zm := uint32(z % uint64(rt.p))
			if zm == rt.r1 || zm == rt.r2 {
				divide(int32(i))
			}
		}
	}
	for _, i := range hits {
		if in64 && v64 == 1 {
			break
		}
		divide(i)
	}
	w.scratch = factors[:0]

	if !in64 {
		return nil, 0, 0, false
	}
	if v64 != 1 {
		c := v64
		switch {
		case c <= s.lpBound && (c < s.pMax*s.pMax || arith.IsPrime64(c)):
			lp1 = c
		case c <= s.dlpCap:
			if lp1, lp2, ok = arith.SplitCofactor(c, s.pMax, s.lpBound); !ok {
				return nil, 0, 0, false
			}
		default:
			return nil, 0, 0, false
		}
	}
	return append([]int32(nil), factors...), lp1, lp2, true
}

func (w *worker) tryFactorizeBig(z uint64, y *big.Int, aIdx []int, hits []int32, factors []int32) ([]int32, uint64, uint64, bool) {
	s := w.s
	temp := w.temp.Abs(y)
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
	for i := 0; i < s.resieveFrom; i++ {
		rt := &w.roots[i]
		if rt.inA {
			continue
		}
		zm := uint32(z % uint64(rt.p))
		if zm == rt.r1 || zm == rt.r2 {
			divide(int32(i))
		}
	}
	for _, i := range hits {
		divide(i)
	}
	w.scratch = factors[:0]
	var lp1, lp2 uint64
	if temp.Cmp(one) != 0 {
		if !temp.IsUint64() {
			return nil, 0, 0, false
		}
		c := temp.Uint64()
		switch {
		case c <= s.lpBound && (c < s.pMax*s.pMax || arith.IsPrime64(c)):
			lp1 = c
		case c <= s.dlpCap:
			var ok bool
			if lp1, lp2, ok = arith.SplitCofactor(c, s.pMax, s.lpBound); !ok {
				return nil, 0, 0, false
			}
		default:
			return nil, 0, 0, false
		}
	}
	return append([]int32(nil), factors...), lp1, lp2, true
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
		have := len(s.relations) + s.partialExcess
		need := s.want - have
		switch {
		case need <= 0:
			s.nextCount = len(s.partials) + 1
		case need <= 2000:
			step := need * 2
			if step < 100 {
				step = 100
			}
			if step > 1000 {
				step = 1000
			}
			s.nextCount = len(s.partials) + step
		case have >= s.want/2:
			s.nextCount = len(s.partials) + len(s.partials)/50 + 300
		default:
			s.nextCount = len(s.partials) + len(s.partials)/15 + 1000
		}
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

package gnfs

import (
	"math/big"
	"os"
	"strconv"
	"time"

	"github.com/dgryski/go-factor/internal/arith"
	"github.com/dgryski/go-factor/internal/linalg"
	"github.com/dgryski/go-factor/internal/progress"
)

const extraRelations = 32

// Logf, if set, receives progress messages while GNFS runs.
var Logf func(format string, args ...interface{})

// Workers is the maximum number of goroutines used by sieving. 0 means GOMAXPROCS.
var Workers = 0

func logf(format string, args ...interface{}) {
	if Logf != nil {
		Logf(format, args...)
	}
}

// Factor factors n using the General Number Field Sieve, returning prime and
// composite factors whose product is n.
func Factor(n *big.Int) (primes, composites []*big.Int) {
	if n.ProbablyPrime(20) {
		return []*big.Int{n}, nil
	}
	if n.IsUint64() {
		return factor64(n.Uint64())
	}

	p, q := split(n)
	if p == nil {
		return nil, []*big.Int{n}
	}

	pr, co := Factor(p)
	primes, composites = append(primes, pr...), append(composites, co...)
	pr, co = Factor(q)
	primes, composites = append(primes, pr...), append(composites, co...)
	return primes, composites
}

// factor64 factors a composite that fits in 64 bits with rho, which is far
// cheaper than a sieve at this size.
func factor64(n uint64) (primes, composites []*big.Int) {
	if n == 1 {
		return nil, nil
	}
	if arith.IsPrime64(n) {
		return []*big.Int{new(big.Int).SetUint64(n)}, nil
	}
	d := arith.Rho64(n)
	if d == 0 {
		return nil, []*big.Int{new(big.Int).SetUint64(n)}
	}
	pr, co := factor64(d)
	primes, composites = append(primes, pr...), append(composites, co...)
	pr, co = factor64(n / d)
	primes, composites = append(primes, pr...), append(composites, co...)
	return primes, composites
}

func defaultParams(n *big.Int) *Params {
	bits := n.BitLen()
	switch {
	case bits < 115:
		return &Params{
			Degree: 3, Br: 1200, Ba: 1800,
			LpBoundRat: 32000, LpBoundAlg: 48000,
			MfbRat: 30, MfbAlg: 32, MaxLPRat: 2, MaxLPAlg: 2,
			UseLattice: false, LineA: 65536, LineB: 1000,
			SlackRat: 34, SlackAlg: 36, NumQCs: 64,
			RelTargetRatio: 0.80,
		}
	case bits < 150:
		return &Params{
			Degree: 4, Br: 8000, Ba: 8000,
			LpBoundRat: 1 << 17, LpBoundAlg: 1 << 17,
			MfbRat: 34, MfbAlg: 34, MaxLPRat: 2, MaxLPAlg: 2,
			UseLattice: true, QMin: 8000, QMax: 60000, IHalf: 512, JMax: 512,
			SlackRat: 38, SlackAlg: 38, NumQCs: 64,
			RelTargetRatio: 0.80,
		}
	case bits < 200:
		return &Params{
			Degree: 4, Br: 30000, Ba: 30000,
			LpBoundRat: 1 << 18, LpBoundAlg: 1 << 18,
			MfbRat: 36, MfbAlg: 36, MaxLPRat: 2, MaxLPAlg: 2,
			UseLattice: true, QMin: 20000, QMax: 120000, IHalf: 1024, JMax: 512,
			SlackRat: 38, SlackAlg: 38, NumQCs: 64,
			RelTargetRatio: 0.85,
		}
	case bits < 250:
		return &Params{
			Degree: 4, Br: 100000, Ba: 100000,
			LpBoundRat: 1 << 20, LpBoundAlg: 1 << 20,
			MfbRat: 40, MfbAlg: 40, MaxLPRat: 2, MaxLPAlg: 2,
			UseLattice: true, QMin: 60000, QMax: 400000, IHalf: 1024, JMax: 1024,
			SlackRat: 41, SlackAlg: 41, NumQCs: 64,
			RelTargetRatio: 0.82,
		}
	case bits < 275:
		return &Params{
			Degree: 4, Br: 250000, Ba: 250000,
			LpBoundRat: 1 << 21, LpBoundAlg: 1 << 21,
			MfbRat: 42, MfbAlg: 42, MaxLPRat: 2, MaxLPAlg: 2,
			UseLattice: true, QMin: 80000, QMax: 800000, IHalf: 2048, JMax: 1024,
			SlackRat: 43, SlackAlg: 43, NumQCs: 64,
			RelTargetRatio: 0.82,
		}
	case bits < 305:
		return &Params{
			Degree: 4, Br: 400000, Ba: 400000,
			LpBoundRat: 1 << 22, LpBoundAlg: 1 << 22,
			MfbRat: 44, MfbAlg: 44, MaxLPRat: 2, MaxLPAlg: 2,
			UseLattice: true, QMin: 150000, QMax: 2000000, IHalf: 2048, JMax: 2048,
			SlackRat: 45, SlackAlg: 45, NumQCs: 64,
			RelTargetRatio: 0.84,
		}
	case bits < 325:
		return &Params{
			Degree: 4, Br: 550000, Ba: 600000,
			LpBoundRat: 1 << 23, LpBoundAlg: 1 << 23,
			MfbRat: 46, MfbAlg: 46, MaxLPRat: 2, MaxLPAlg: 2,
			UseLattice: true, QMin: 220000, QMax: 3500000, IHalf: 2048, JMax: 2048,
			SlackRat: 47, SlackAlg: 47, NumQCs: 64,
			RelTargetRatio: 0.77,
		}
	case bits < 365:
		return &Params{
			Degree: 4, Br: 700000, Ba: 800000,
			LpBoundRat: 12000000, LpBoundAlg: 12000000,
			MfbRat: 47, MfbAlg: 47, MaxLPRat: 2, MaxLPAlg: 2,
			UseLattice: true, QMin: 280000, QMax: 5000000, IHalf: 2048, JMax: 2048,
			SlackRat: 48, SlackAlg: 48, NumQCs: 64,
			RelTargetRatio: 0.80,
		}
	default:
		return &Params{
			Degree: 5, Br: 500000, Ba: 800000,
			LpBoundRat: 1 << 26, LpBoundAlg: 1 << 26,
			MfbRat: 52, MfbAlg: 52, MaxLPRat: 2, MaxLPAlg: 2,
			UseLattice: true, QMin: 400000, QMax: 10000000, IHalf: 4096, JMax: 2048,
			SlackRat: 53, SlackAlg: 53, NumQCs: 64,
			RelTargetRatio: 1.05,
		}
	}
}

func computeQCMin(params *Params) uint64 {
	qcMin := params.LpBoundAlg
	if uint64(params.Ba) > qcMin {
		qcMin = uint64(params.Ba)
	}
	if params.UseLattice {
		// Keep QC primes above extended special-q ranges while staying below 2^31.
		ext := params.QMax * 64
		const maxQ = (1 << 31) - 10000
		if ext > maxQ {
			ext = maxQ
		}
		if ext > qcMin {
			qcMin = ext
		}
	}
	return qcMin
}

func split(n *big.Int) (*big.Int, *big.Int) {
	return splitWithParams(n, nil)
}

func splitWithParams(n *big.Int, params *Params) (*big.Int, *big.Int) {
	if n.Cmp(big.NewInt(4)) < 0 || n.ProbablyPrime(20) {
		return nil, nil
	}

	// Check perfect power / square first.
	sqrtN := new(big.Int).Sqrt(n)
	if new(big.Int).Mul(sqrtN, sqrtN).Cmp(n) == 0 {
		return sqrtN, new(big.Int).Set(sqrtN)
	}

	if params == nil {
		params = defaultParams(n)
	}

	// Trial divide up to max(Br, Ba) first so factor bases have no prime factors of N.
	maxB := params.Br
	if params.Ba > maxB {
		maxB = params.Ba
	}
	var bigP, r big.Int
	for _, p := range arith.PrimesUpTo(maxB) {
		bigP.SetInt64(p)
		if r.Mod(n, &bigP).Sign() == 0 {
			return new(big.Int).Set(&bigP), new(big.Int).Div(n, &bigP)
		}
	}

	qcMin := computeQCMin(params)

	var poly *Poly
	var savedPairs [][2]int64
	var savedFree []bool
	if SaveDir != "" {
		path := relationsPath(n)
		if f, err := os.Open(path); err == nil {
			nSaved, pSaved, pairs, isFree, rerr := readPolyAndPairs(f)
			_ = f.Close()
			if rerr == nil && nSaved.Cmp(n) == 0 {
				poly = pSaved
				savedPairs = pairs
				savedFree = isFree
				logf("gnfs: resuming %d bits (%d digits) with %d saved pairs from %s",
					n.BitLen(), len(n.String()), len(pairs), path)
			} else if rerr != nil {
				logf("gnfs: ignoring %s: %v", path, rerr)
			}
		}
	}

	if poly == nil {
		polyStart := time.Now()
		var fac *big.Int
		poly, fac = selectPolyWithParams(n, params)
		if fac != nil {
			logf("gnfs: found factor %v during polynomial selection", fac)
			return fac, new(big.Int).Div(n, fac)
		}
		if poly == nil {
			logf("gnfs: polynomial selection failed")
			return nil, nil
		}
		var g big.Int
		g.GCD(nil, nil, poly.A1, n)
		if g.Cmp(bigOne) > 0 && g.Cmp(n) < 0 {
			return &g, new(big.Int).Div(n, &g)
		}
		logf("gnfs: %d bits (%d digits), deg=%d, c_d=%v, a_1=%v, skew=%.1f, alpha=%.2f in %v",
			n.BitLen(), len(n.String()), poly.Degree, poly.Coeffs[poly.Degree], poly.A1,
			poly.Skew, poly.Alpha, time.Since(polyStart).Round(time.Millisecond))
	}

	fb := BuildFactorBases(poly, params.Br, params.Ba, qcMin, params.NumQCs)
	freeRels := freeRelations(poly, fb)
	target := newRelTarget(params, extraRelations)
	logf("gnfs: RFB=%d (Br=%d), AFB=%d (Ba=%d, bad=%d), QCB=%d (>%d), free=%d, cols=%d, want excess=%d",
		len(fb.RFB), params.Br, len(fb.AFB), params.Ba, fb.BadCount, len(fb.QCB), qcMin, len(freeRels), fb.NumFixed, extraRelations)

	initCap := target.initial
	if len(savedPairs) > initCap {
		initCap = len(savedPairs)
	}
	seen := make(map[[2]int64]struct{}, initCap)
	rels := make([]Relation, 0, initCap)
	filt := newLPFilter(initCap)
	var fullCount int
	excess := -fb.NumFixed

	addRel := func(rel Relation) bool {
		key := [2]int64{rel.A, rel.B}
		if _, ok := seen[key]; ok {
			return false
		}
		seen[key] = struct{}{}
		rels = append(rels, rel)
		filt.add(&rels[len(rels)-1])
		if rel.RatLPs[0] == 0 && len(rel.AlgLPs) == 0 {
			fullCount++
		}
		return true
	}

	if len(savedPairs) > 0 {
		loadStart := time.Now()
		for _, rel := range refactorPairs(poly, fb, params, qcMin, savedPairs, savedFree) {
			addRel(rel)
		}
		excess = filt.excess(fb.NumFixed)
		logf("gnfs: refactored %d/%d saved relations (excess %d/%d) in %v",
			len(rels), len(savedPairs), excess, extraRelations, time.Since(loadStart).Round(time.Millisecond))
	}
	loadedCount := len(rels)

	for _, rel := range freeRels {
		addRel(rel)
	}
	if len(rels) > loadedCount {
		excess = filt.excess(fb.NumFixed)
	}

	shouldStop := func() bool {
		return excess >= extraRelations
	}

	if !shouldStop() {
		logf("gnfs: expecting about %d relations", target.initial)

		nextCheck := fb.NumFixed / 2
		if nextCheck < 500 {
			nextCheck = 500
		}
		if len(rels) >= nextCheck {
			nextCheck = len(rels) + 500
		}
		sieveStart := time.Now()
		lastLog := sieveStart
		var rate progress.Rate
		if len(rels) > 0 {
			rate.Sample(sieveStart, len(rels))
		}

		onRel := func(rel Relation) {
			if !addRel(rel) {
				return
			}
			if len(rels) >= nextCheck {
				excess = filt.excess(fb.NumFixed)
				// Check roughly every 7% before the 2-core climbs, every 1%
				// once climbing, and tightly once close to the goal.
				switch {
				case excess >= -2000:
					step := (extraRelations - excess) * 2
					if step < 100 {
						step = 100
					}
					if step > 1000 {
						step = 1000
					}
					nextCheck = len(rels) + step
				case excess >= -fb.NumFixed/2:
					nextCheck = len(rels) + len(rels)/100 + 200
				default:
					nextCheck = len(rels) + len(rels)/15 + 500
				}
				want := target.update(len(rels), excess)
				now := time.Now()
				if now.Sub(lastLog) >= 5*time.Second {
					rate.Sample(now, len(rels))
					lastLog = now
					wantStr, left := "?", "time left unknown"
					if want > 0 {
						wantStr, left = strconv.Itoa(want), rate.Left(len(rels), want)
					}
					logf("gnfs: %d/%s relations (excess %d/%d, %d full), %v, %s",
						len(rels), wantStr, excess, extraRelations, fullCount,
						now.Sub(sieveStart).Round(time.Millisecond), left)
				}
			}
		}

		savedSoFar := loadedCount
		if !params.UseLattice {
			curParams := *params
			var prevA, prevB int64
			for round := 0; round < 8 && !shouldStop(); round++ {
				sieveLine(poly, fb, &curParams, prevA, prevB, onRel, shouldStop)
				excess = filt.excess(fb.NumFixed)
				if !shouldStop() {
					prevA, prevB = curParams.LineA, curParams.LineB
					curParams.LineA *= 2
					if curParams.LineB < 8000 {
						curParams.LineB *= 2
					}
				}
			}
		} else {
			curParams := *params
			span := curParams.QMax - curParams.QMin
			if span < 50000 {
				span = 50000
			}
			for !shouldStop() && curParams.QMin < qcMin {
				if curParams.QMax >= qcMin {
					curParams.QMax = qcMin - 1
				}
				sieveLattice(poly, fb, &curParams, onRel, shouldStop)
				excess = filt.excess(fb.NumFixed)
				if SaveDir != "" && !shouldStop() && len(rels) > savedSoFar {
					saveRelations(n, poly, rels)
					savedSoFar = len(rels)
				}
				if !shouldStop() {
					curParams.QMin = curParams.QMax + 1
					if span < 2000000 {
						span *= 2
						if span > 2000000 {
							span = 2000000
						}
					}
					curParams.QMax = curParams.QMin + span
				}
			}
		}

		excess = filt.excess(fb.NumFixed)
		logf("gnfs: sieved %d unique relations (%d full, excess %d) in %v",
			len(rels), fullCount, excess, time.Since(sieveStart).Round(time.Millisecond))
	}

	if excess <= 0 {
		logf("gnfs: not enough relations (excess=%d), giving up", excess)
		return nil, nil
	}

	if SaveDir != "" && len(rels) > loadedCount {
		saveRelations(n, poly, rels)
	}

	return solveRelations(n, poly, fb, rels)
}

func solveRelations(n *big.Int, poly *Poly, fb *FactorBases, rels []Relation) (*big.Int, *big.Int) {
	filtStart := time.Now()
	logf("gnfs: starting filtering on %d relations", len(rels))
	surviving, rows, numCols := buildMatrix(poly, fb, rels)
	logf("gnfs: filtered matrix %d×%d in %v, starting linear algebra", len(rows), numCols, time.Since(filtStart).Round(time.Millisecond))
	if len(rows) <= numCols {
		logf("gnfs: matrix has no surplus rows, giving up")
		return nil, nil
	}

	laStart := time.Now()
	deps := linalg.FindDependencies(rows, numCols, func(f string, a ...interface{}) {
		logf("gnfs: "+f, a...)
	})
	logf("gnfs: %d dependencies in %v, starting square root", len(deps), time.Since(laStart).Round(time.Millisecond))

	sqrtStart := time.Now()
	badSq := 0
	for i, dep := range deps {
		logf("gnfs: trying dependency %d/%d", i+1, len(deps))
		fac, sqOk := tryDependency(n, poly, fb, surviving, dep)
		if !sqOk {
			badSq++
			continue
		}
		if fac != nil {
			logf("gnfs: dependency %d gave factor %v (sqrt time %v, %d bad squares)",
				i+1, fac, time.Since(sqrtStart).Round(time.Millisecond), badSq)
			return fac, new(big.Int).Div(n, fac)
		}
	}
	logf("gnfs: all %d dependencies failed or trivial (%d bad squares)", len(deps), badSq)
	return nil, nil
}

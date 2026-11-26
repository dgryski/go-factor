package gnfs

import (
	"math/big"
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
			Degree: 3, Br: 800, Ba: 1200,
			LpBoundRat: 24000, LpBoundAlg: 32000,
			MfbRat: 28, MfbAlg: 30, MaxLPRat: 2, MaxLPAlg: 2,
			UseLattice: false, LineA: 65536, LineB: 1000,
			SlackRat: 34, SlackAlg: 36, NumQCs: 64,
		}
	case bits < 150:
		return &Params{
			Degree: 4, Br: 8000, Ba: 8000,
			LpBoundRat: 1 << 17, LpBoundAlg: 1 << 17,
			MfbRat: 34, MfbAlg: 34, MaxLPRat: 2, MaxLPAlg: 2,
			UseLattice: true, QMin: 8000, QMax: 60000, IHalf: 512, JMax: 512,
			SlackRat: 38, SlackAlg: 38, NumQCs: 64,
		}
	case bits < 200:
		return &Params{
			Degree: 4, Br: 30000, Ba: 30000,
			LpBoundRat: 1 << 18, LpBoundAlg: 1 << 18,
			MfbRat: 36, MfbAlg: 36, MaxLPRat: 2, MaxLPAlg: 2,
			UseLattice: true, QMin: 30000, QMax: 120000, IHalf: 1024, JMax: 512,
			SlackRat: 40, SlackAlg: 40, NumQCs: 64,
		}
	case bits < 250:
		return &Params{
			Degree: 4, Br: 100000, Ba: 100000,
			LpBoundRat: 1 << 20, LpBoundAlg: 1 << 20,
			MfbRat: 40, MfbAlg: 40, MaxLPRat: 2, MaxLPAlg: 2,
			UseLattice: true, QMin: 100000, QMax: 400000, IHalf: 1024, JMax: 1024,
			SlackRat: 46, SlackAlg: 46, NumQCs: 64,
		}
	case bits < 300:
		return &Params{
			Degree: 4, Br: 40000, Ba: 80000,
			LpBoundRat: 1 << 24, LpBoundAlg: 1 << 24,
			MfbRat: 48, MfbAlg: 48, MaxLPRat: 2, MaxLPAlg: 2,
			UseLattice: true, QMin: 80000, QMax: 800000, IHalf: 2048, JMax: 1024,
			SlackRat: 54, SlackAlg: 54, NumQCs: 64,
		}
	case bits < 365:
		return &Params{
			Degree: 4, Br: 150000, Ba: 300000,
			LpBoundRat: 1 << 26, LpBoundAlg: 1 << 26,
			MfbRat: 52, MfbAlg: 52, MaxLPRat: 2, MaxLPAlg: 2,
			UseLattice: true, QMin: 300000, QMax: 2500000, IHalf: 2048, JMax: 2048,
			SlackRat: 58, SlackAlg: 58, NumQCs: 64,
		}
	default:
		return &Params{
			Degree: 5, Br: 400000, Ba: 800000,
			LpBoundRat: 1 << 28, LpBoundAlg: 1 << 28,
			MfbRat: 56, MfbAlg: 56, MaxLPRat: 2, MaxLPAlg: 2,
			UseLattice: true, QMin: 800000, QMax: 10000000, IHalf: 4096, JMax: 2048,
			SlackRat: 62, SlackAlg: 62, NumQCs: 64,
		}
	}
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

	polyStart := time.Now()
	poly, fac := SelectPoly(n, params.Degree)
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
	logf("gnfs: %d bits (%d digits), deg=%d, c_d=%v, skew=%.1f, alpha=%.2f in %v",
		n.BitLen(), len(n.String()), poly.Degree, poly.Coeffs[poly.Degree],
		poly.Skew, poly.Alpha, time.Since(polyStart).Round(time.Millisecond))

	// Ensure QC primes exceed every algebraic prime factor and extended special-q range.
	qcMin := params.LpBoundAlg
	if uint64(params.Ba) > qcMin {
		qcMin = uint64(params.Ba)
	}
	if params.UseLattice && params.QMax*8 > qcMin {
		qcMin = params.QMax * 8
	}
	fb := BuildFactorBases(poly, params.Br, params.Ba, qcMin, params.NumQCs)
	target := newRelTarget(params, extraRelations)
	logf("gnfs: RFB=%d (Br=%d), AFB=%d (Ba=%d, bad=%d), QCB=%d (>%d), cols=%d, want excess=%d",
		len(fb.RFB), params.Br, len(fb.AFB), params.Ba, fb.BadCount, len(fb.QCB), qcMin, fb.NumFixed, extraRelations)
	logf("gnfs: expecting about %d relations", target.initial)

	seen := make(map[[2]int64]struct{})
	var rels []Relation
	var fullCount int
	excess := 0
	nextCheck := fb.NumFixed / 2
	if nextCheck < 500 {
		nextCheck = 500
	}
	sieveStart := time.Now()
	lastLog := sieveStart
	var rate progress.Rate

	onRel := func(rel Relation) {
		key := [2]int64{rel.A, rel.B}
		if _, ok := seen[key]; ok {
			return
		}
		seen[key] = struct{}{}
		rels = append(rels, rel)
		if rel.RatLPs[0] == 0 && len(rel.AlgLPs) == 0 {
			fullCount++
		}
		if len(rels) >= nextCheck {
			ex, _, _ := pruneSingletons(rels, fb)
			excess = ex
			// Once cycles form, excess climbs by thousands per few percent
			// more relations, so check often to stop near the goal.
			if excess > -fb.NumFixed {
				nextCheck = len(rels) + len(rels)/100 + 500
			} else {
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
	shouldStop := func() bool {
		return excess >= extraRelations
	}

	if !params.UseLattice {
		curParams := *params
		var prevA, prevB int64
		for round := 0; round < 5 && !shouldStop(); round++ {
			sieveLine(poly, fb, &curParams, prevA, prevB, onRel, shouldStop)
			ex, _, _ := pruneSingletons(rels, fb)
			excess = ex
			if !shouldStop() {
				prevA, prevB = curParams.LineA, curParams.LineB
				curParams.LineA *= 2
				if curParams.LineB < 4000 {
					curParams.LineB *= 2
				}
			}
		}
	} else {
		curParams := *params
		for round := 0; round < 4 && !shouldStop() && curParams.QMax < qcMin; round++ {
			sieveLattice(poly, fb, &curParams, onRel, shouldStop)
			ex, _, _ := pruneSingletons(rels, fb)
			excess = ex
			if !shouldStop() {
				curParams.QMin = curParams.QMax + 1
				curParams.QMax *= 2
				if curParams.QMax >= qcMin {
					curParams.QMax = qcMin - 1
				}
			}
		}
	}

	ex, _, _ := pruneSingletons(rels, fb)
	excess = ex
	logf("gnfs: sieved %d unique relations (%d full, excess %d) in %v",
		len(rels), fullCount, excess, time.Since(sieveStart).Round(time.Millisecond))

	if excess <= 0 {
		logf("gnfs: not enough relations (excess=%d), giving up", excess)
		return nil, nil
	}

	if SaveDir != "" {
		saveRelations(n, poly, rels)
	}

	surviving, rows, numCols := buildMatrix(poly, fb, rels)
	logf("gnfs: filtered matrix %d×%d", len(rows), numCols)
	if len(rows) <= numCols {
		logf("gnfs: matrix has no surplus rows, giving up")
		return nil, nil
	}

	laStart := time.Now()
	deps := linalg.FindDependencies(rows, numCols, func(f string, a ...interface{}) {
		logf("gnfs: "+f, a...)
	})
	logf("gnfs: %d dependencies in %v", len(deps), time.Since(laStart).Round(time.Millisecond))

	sqrtStart := time.Now()
	badSq := 0
	for i, dep := range deps {
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

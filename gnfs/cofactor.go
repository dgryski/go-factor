package gnfs

import (
	"math/big"
	"math/bits"

	"github.com/dgryski/go-factor/internal/arith"
)

// LargeIdeal is a large prime on the algebraic side. R = A * B^-1 mod P,
// or R = P for projective (P | c_d).
type LargeIdeal struct {
	P, R uint64
}

// Relation represents a smooth (or partial) pair (A, B) with gcd(|A|, B) == 1, B > 0.
// A free relation has A = p and B = 0, records p in RatCols or RatLPs like any
// other rational prime, and records its d algebraic ideals in AlgCols.
type Relation struct {
	A, B    int64
	RatCols []int32      // RFB indices (with multiplicity)
	AlgCols []int32      // AFB indices (with multiplicity, including special-q if <= Ba)
	RatLPs  [2]uint64    // rational large primes (0 if absent)
	AlgLPs  []LargeIdeal // algebraic large ideals (including special-q if > Ba)
	IsFree  bool
}

// Params holds the tuning parameters for a GNFS factorization.
type Params struct {
	Degree     int
	Br, Ba     int64  // rational and algebraic factor base bounds
	LpBoundRat uint64 // rational large prime bound
	LpBoundAlg uint64 // algebraic large prime bound
	MfbRat     int    // max cofactor bits on rational side
	MfbAlg     int    // max cofactor bits on algebraic side
	MaxLPRat   int    // 1 or 2 large primes on rational side
	MaxLPAlg   int    // 1 or 2 large primes on algebraic side
	QMin, QMax uint64 // special-q range (if lattice sieve is used)
	IHalf      int    // lattice sieve i in [-IHalf, IHalf)
	JMax       int    // lattice sieve j in [1, JMax)
	UseLattice bool   // false = line sieve, true = special-q lattice sieve
	LineA      int64  // line sieve half-width a in [-LineA, LineA)
	LineB      int64  // line sieve b in [1, LineB]
	SlackRat   float64
	SlackAlg   float64
	NumQCs     int
}

// splitCofactor splits c (already divided by all primes <= pMax) into 1 or 2
// primes each <= lpBound, returning them sorted ascending.
func splitCofactor(c, pMax, lpBound uint64, maxLP int, mfb int) ([2]uint64, bool) {
	if c == 1 {
		return [2]uint64{}, true
	}
	if mfb > 0 && bits.Len64(c) > mfb {
		return [2]uint64{}, false
	}
	if c <= lpBound && (c < pMax*pMax || arith.IsPrime64(c)) {
		return [2]uint64{c, 0}, true
	}
	if maxLP < 2 {
		return [2]uint64{}, false
	}
	lp1, lp2, ok := arith.SplitCofactor(c, pMax, lpBound)
	if !ok {
		return [2]uint64{}, false
	}
	return [2]uint64{lp1, lp2}, true
}

// cofactorCandidate attempts to factor G(a, b) and F(a, b) for coprime (a, b) with b > 0.
// ratHits are RFB indices >= smallRatCut hitting (a, b); algHits are AFB indices >= smallAlgCut.
// specialQ (if sq.Q != 0) is known to divide F(a, b).
func cofactorCandidate(
	poly *Poly,
	fb *FactorBases,
	p *Params,
	a, b int64,
	ratSmall []int32, ratHits []int32,
	algSmall []int32, algHits []int32,
	sq *LatticeQ,
	scratch *cofactorScratch,
) (Relation, bool) {
	// Rational side first (much smaller numbers than F(a, b)).
	gVal := evalGInto(&scratch.val, poly, a, b)
	if gVal.Sign() == 0 {
		return Relation{}, false
	}
	temp := scratch.temp.Abs(gVal)
	ratCols := scratch.ratCols[:0]

	var bigP, q, rem big.Int
	divRat := func(idx int32) {
		pr := fb.RFB[idx].P
		bigP.SetUint64(uint64(pr))
		for {
			q.QuoRem(temp, &bigP, &rem)
			if rem.Sign() != 0 {
				break
			}
			ratCols = append(ratCols, idx)
			temp.Set(&q)
		}
	}

	for _, idx := range ratSmall {
		id := fb.RFB[idx]
		pr := uint64(id.P)
		if id.R == id.P {
			if uint64(b)%pr == 0 {
				divRat(idx)
			}
		} else {
			am := modInt64(a, pr)
			bm := uint64(b) % pr
			if am == bm*uint64(id.R)%pr {
				divRat(idx)
			}
		}
	}
	for _, idx := range ratHits {
		divRat(idx)
	}

	var ratLPs [2]uint64
	if temp.Cmp(bigOne) != 0 {
		if !temp.IsUint64() {
			return Relation{}, false
		}
		var ok bool
		ratLPs, ok = splitCofactor(temp.Uint64(), uint64(p.Br), p.LpBoundRat, p.MaxLPRat, p.MfbRat)
		if !ok {
			return Relation{}, false
		}
	}

	// Algebraic side.
	fVal := evalFInto(&scratch.val, poly, a, b, scratch.aPows, scratch.bPows)
	if fVal.Sign() == 0 {
		return Relation{}, false
	}
	temp.Abs(fVal)
	algCols := scratch.algCols[:0]
	var algLPs []LargeIdeal

	divAlg := func(idx int32) {
		pr := fb.AFB[idx].P
		bigP.SetUint64(uint64(pr))
		for {
			q.QuoRem(temp, &bigP, &rem)
			if rem.Sign() != 0 {
				break
			}
			algCols = append(algCols, idx)
			temp.Set(&q)
		}
	}

	// Divide out special-q if present.
	if sq != nil && sq.Q != 0 {
		bigP.SetUint64(uint64(sq.Q))
		if sq.AFBIdx >= 0 {
			for {
				q.QuoRem(temp, &bigP, &rem)
				if rem.Sign() != 0 {
					break
				}
				algCols = append(algCols, sq.AFBIdx)
				temp.Set(&q)
			}
		} else {
			sqIdeal := LargeIdeal{P: uint64(sq.Q), R: uint64(sq.R)}
			for {
				q.QuoRem(temp, &bigP, &rem)
				if rem.Sign() != 0 {
					break
				}
				algLPs = append(algLPs, sqIdeal)
				temp.Set(&q)
			}
		}
	}

	for _, idx := range algSmall {
		id := fb.AFB[idx]
		pr := uint64(id.P)
		if id.R == id.P {
			if uint64(b)%pr == 0 {
				divAlg(idx)
			}
		} else {
			am := modInt64(a, pr)
			bm := uint64(b) % pr
			if am == bm*uint64(id.R)%pr {
				divAlg(idx)
			}
		}
	}
	for _, idx := range algHits {
		divAlg(idx)
	}

	if temp.Cmp(bigOne) != 0 {
		if !temp.IsUint64() {
			return Relation{}, false
		}
		lps, ok := splitCofactor(temp.Uint64(), uint64(p.Ba), p.LpBoundAlg, p.MaxLPAlg, p.MfbAlg)
		if !ok {
			return Relation{}, false
		}
		for _, lp := range lps {
			if lp == 0 {
				continue
			}
			bm := uint64(b) % lp
			if bm == 0 {
				// Projective large prime (lp | c_d); shouldn't happen when c_d is smooth <= Ba.
				algLPs = append(algLPs, LargeIdeal{P: lp, R: lp})
			} else {
				invB, ok := arith.InvMod64(bm, lp)
				if !ok {
					return Relation{}, false
				}
				r := arith.MulMod64(modInt64(a, lp), invB, lp)
				if isMultipleRoot(poly, r, lp) {
					// v_lp(F) need not be the ideal's valuation, and no
					// character covers it.
					return Relation{}, false
				}
				algLPs = append(algLPs, LargeIdeal{P: lp, R: r})
			}
		}
	}

	rel := Relation{
		A:       a,
		B:       b,
		RatCols: append([]int32(nil), ratCols...),
		AlgCols: append([]int32(nil), algCols...),
		RatLPs:  ratLPs,
		AlgLPs:  algLPs,
	}
	return rel, true
}

var bigOne = big.NewInt(1)

// isMultipleRoot reports whether r is a multiple root of f mod p.
func isMultipleRoot(poly *Poly, r, p uint64) bool {
	dfp := arith.PolyDerivModP(poly.CoeffsModP(p), p)
	return arith.PolyEvalModP(dfp, r, p) == 0
}

type cofactorScratch struct {
	val, temp big.Int
	aPows     []big.Int
	bPows     []big.Int
	ratCols   []int32
	algCols   []int32
}

func newCofactorScratch(deg int) *cofactorScratch {
	return &cofactorScratch{
		aPows: make([]big.Int, deg+1),
		bPows: make([]big.Int, deg+1),
	}
}

func evalGInto(dst *big.Int, poly *Poly, a, b int64) *big.Int {
	var t1, t2 big.Int
	t1.SetInt64(a)
	t1.Mul(&t1, poly.A1)
	t2.SetInt64(b)
	t2.Mul(&t2, poly.A0)
	return dst.Sub(&t1, &t2)
}

func evalFInto(dst *big.Int, poly *Poly, a, b int64, aPows, bPows []big.Int) *big.Int {
	d := poly.Degree
	aPows[0].SetInt64(1)
	bPows[0].SetInt64(1)
	var bigA, bigB, term big.Int
	bigA.SetInt64(a)
	bigB.SetInt64(b)
	for i := 1; i <= d; i++ {
		aPows[i].Mul(&aPows[i-1], &bigA)
		bPows[i].Mul(&bPows[i-1], &bigB)
	}
	dst.SetInt64(0)
	for i := 0; i <= d; i++ {
		term.Mul(poly.Coeffs[i], &aPows[i])
		term.Mul(&term, &bPows[d-i])
		dst.Add(dst, &term)
	}
	return dst
}

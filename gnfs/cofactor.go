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
	Degree         int
	Br, Ba         int64  // rational and algebraic factor base bounds
	LpBoundRat     uint64 // rational large prime bound
	LpBoundAlg     uint64 // algebraic large prime bound
	MfbRat         int    // max cofactor bits on rational side
	MfbAlg         int    // max cofactor bits on algebraic side
	MaxLPRat       int    // 1 or 2 large primes on rational side
	MaxLPAlg       int    // 1 or 2 large primes on algebraic side
	QMin, QMax     uint64 // special-q range (if lattice sieve is used)
	IHalf          int    // lattice sieve i in [-IHalf, IHalf)
	JMax           int    // lattice sieve j in [1, JMax)
	UseLattice     bool   // false = line sieve, true = special-q lattice sieve
	LineA          int64  // line sieve half-width a in [-LineA, LineA)
	LineB          int64  // line sieve b in [1, LineB]
	SlackRat       float64
	SlackAlg       float64
	NumQCs         int
	RelTargetRatio float64 // expected relations / LP ideals (0 uses defaultRelTargetRatio)
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

// divRem128 divides (hi, lo) by d (d > 0) and returns quotient (qHi, qLo) and remainder r.
func divRem128(hi, lo, d uint64) (qHi, qLo, r uint64) {
	if hi == 0 {
		return 0, lo / d, lo % d
	}
	if hi < d {
		qLo, r = bits.Div64(hi, lo, d)
		return 0, qLo, r
	}
	qHi = hi / d
	qLo, r = bits.Div64(hi%d, lo, d)
	return qHi, qLo, r
}

func bigToU128(v *big.Int) (hi, lo uint64, ok bool) {
	if bits.UintSize == 64 {
		words := v.Bits()
		switch len(words) {
		case 0:
			return 0, 0, true
		case 1:
			return 0, uint64(words[0]), true
		case 2:
			return uint64(words[1]), uint64(words[0]), true
		default:
			return 0, 0, false
		}
	}
	if v.BitLen() > 128 {
		return 0, 0, false
	}
	lo = v.Uint64()
	var top big.Int
	top.Rsh(v, 64)
	return top.Uint64(), lo, true
}

// normState holds a positive integer being trial-divided, switching from big.Int
// to 128-bit and 64-bit arithmetic as soon as it fits.
type normState struct {
	hi, lo uint64
	use128 bool
}

func initNorm(n *normState, v *big.Int) {
	if hi, lo, ok := bigToU128(v); ok {
		n.hi, n.lo, n.use128 = hi, lo, true
	} else {
		n.use128 = false
	}
}

func (n *normState) divAll(p uint64, temp *big.Int, sc *cofactorScratch) int {
	mult := 0
	if !n.use128 {
		sc.bigP.SetUint64(p)
		for {
			sc.q.QuoRem(temp, &sc.bigP, &sc.rem)
			if sc.rem.Sign() != 0 {
				break
			}
			mult++
			temp.Set(&sc.q)
			if hi, lo, ok := bigToU128(temp); ok {
				n.hi, n.lo, n.use128 = hi, lo, true
				break
			}
		}
		if !n.use128 {
			return mult
		}
	}
	for n.hi != 0 {
		qh, ql, r := divRem128(n.hi, n.lo, p)
		if r != 0 {
			return mult
		}
		n.hi, n.lo = qh, ql
		mult++
	}
	for n.lo >= p && n.lo%p == 0 {
		n.lo /= p
		mult++
	}
	return mult
}

func (n *normState) uint64Val() (uint64, bool) {
	if !n.use128 || n.hi != 0 {
		return 0, false
	}
	return n.lo, true
}

// cofactorRatLattice factors G(a, b) for a lattice-sieve candidate at (i, j)
// with normalized coprime (a, b), b > 0.
func cofactorRatLattice(
	poly *Poly,
	fb *FactorBases,
	p *Params,
	i, j int,
	a, b int64,
	ratSmall []int32,
	ratLines []linePrime,
	ratVerts []vertPrime,
	ratAll []int32,
	ratHits []int32,
	sc *cofactorScratch,
) (ratCols []int32, ratLPs [2]uint64, ok bool) {
	gVal := evalGInto(&sc.val, poly, a, b, sc)
	if gVal.Sign() == 0 {
		return nil, [2]uint64{}, false
	}
	temp := sc.temp.Abs(gVal)
	var ns normState
	initNorm(&ns, temp)
	cols := sc.ratCols[:0]

	addMult := func(idx int32, pr uint64) {
		for m := ns.divAll(pr, temp, sc); m > 0; m-- {
			cols = append(cols, idx)
		}
	}

	for _, idx := range ratAll {
		addMult(idx, uint64(fb.RFB[idx].P))
	}
	for _, idx := range ratSmall {
		id := fb.RFB[idx]
		pr := uint64(id.P)
		if id.R == id.P {
			if uint64(b)%pr == 0 {
				addMult(idx, pr)
			}
		} else if modInt64(a, pr) == (uint64(b)%pr)*uint64(id.R)%pr {
			addMult(idx, pr)
		}
	}
	if i >= -16384 && i <= 16384 && j > 0 && j <= 16384 {
		i32, j32 := uint32(int32(i)), uint32(j)
		for _, vp := range ratVerts {
			if vp.pInv != 0 {
				if j32*vp.pInv <= vp.pLim {
					addMult(vp.idx, uint64(vp.p))
				}
			} else if j%int(vp.p) == 0 {
				addMult(vp.idx, uint64(vp.p))
			}
		}
		for k := range ratLines {
			lp := &ratLines[k]
			u := j32*lp.lambda + lp.bias - i32
			if u*lp.pInv <= lp.pLim {
				addMult(lp.idx, uint64(lp.p))
			}
		}
	} else {
		i64, j64 := int64(i), int64(j)
		for _, vp := range ratVerts {
			if j%int(vp.p) == 0 {
				addMult(vp.idx, uint64(vp.p))
			}
		}
		for k := range ratLines {
			lp := &ratLines[k]
			if (i64-j64*int64(lp.lambda))%int64(lp.p) == 0 {
				addMult(lp.idx, uint64(lp.p))
			}
		}
	}
	for _, idx := range ratHits {
		addMult(idx, uint64(fb.RFB[idx].P))
	}

	rem64, fits := ns.uint64Val()
	if !fits {
		return nil, [2]uint64{}, false
	}
	if rem64 != 1 {
		ratLPs, ok = splitCofactor(rem64, uint64(p.Br), p.LpBoundRat, p.MaxLPRat, p.MfbRat)
		if !ok {
			return nil, [2]uint64{}, false
		}
	}
	return cols, ratLPs, true
}

// cofactorAlgLattice factors F(a, b) for a lattice-sieve candidate that already
// passed rational cofactorization.
func cofactorAlgLattice(
	poly *Poly,
	fb *FactorBases,
	p *Params,
	i, j int,
	a, b int64,
	ratCols []int32,
	ratLPs [2]uint64,
	algSmall []int32,
	algLines []linePrime,
	algVerts []vertPrime,
	algAll []int32,
	algHits []int32,
	sq *LatticeQ,
	sc *cofactorScratch,
) (Relation, bool) {
	fVal := evalFInto(&sc.val, poly, a, b, sc)
	if fVal.Sign() == 0 {
		return Relation{}, false
	}
	temp := sc.temp.Abs(fVal)
	var ns normState
	initNorm(&ns, temp)
	algCols := sc.algCols[:0]
	var algLPs []LargeIdeal

	addMult := func(idx int32, pr uint64) {
		for m := ns.divAll(pr, temp, sc); m > 0; m-- {
			algCols = append(algCols, idx)
		}
	}

	if sq != nil && sq.Q != 0 {
		m := ns.divAll(uint64(sq.Q), temp, sc)
		if sq.AFBIdx >= 0 {
			for ; m > 0; m-- {
				algCols = append(algCols, sq.AFBIdx)
			}
		} else {
			sqIdeal := LargeIdeal{P: uint64(sq.Q), R: uint64(sq.R)}
			for ; m > 0; m-- {
				algLPs = append(algLPs, sqIdeal)
			}
		}
	}

	for _, idx := range algAll {
		addMult(idx, uint64(fb.AFB[idx].P))
	}
	for _, idx := range algSmall {
		id := fb.AFB[idx]
		pr := uint64(id.P)
		if id.R == id.P {
			if uint64(b)%pr == 0 {
				addMult(idx, pr)
			}
		} else if modInt64(a, pr) == (uint64(b)%pr)*uint64(id.R)%pr {
			addMult(idx, pr)
		}
	}
	if i >= -16384 && i <= 16384 && j > 0 && j <= 16384 {
		i32, j32 := uint32(int32(i)), uint32(j)
		for _, vp := range algVerts {
			if vp.pInv != 0 {
				if j32*vp.pInv <= vp.pLim {
					addMult(vp.idx, uint64(vp.p))
				}
			} else if j%int(vp.p) == 0 {
				addMult(vp.idx, uint64(vp.p))
			}
		}
		for k := range algLines {
			lp := &algLines[k]
			u := j32*lp.lambda + lp.bias - i32
			if u*lp.pInv <= lp.pLim {
				addMult(lp.idx, uint64(lp.p))
			}
		}
	} else {
		i64, j64 := int64(i), int64(j)
		for _, vp := range algVerts {
			if j%int(vp.p) == 0 {
				addMult(vp.idx, uint64(vp.p))
			}
		}
		for k := range algLines {
			lp := &algLines[k]
			if (i64-j64*int64(lp.lambda))%int64(lp.p) == 0 {
				addMult(lp.idx, uint64(lp.p))
			}
		}
	}
	for _, idx := range algHits {
		addMult(idx, uint64(fb.AFB[idx].P))
	}

	rem64, fits := ns.uint64Val()
	if !fits {
		return Relation{}, false
	}
	if rem64 != 1 {
		lps, ok := splitCofactor(rem64, uint64(p.Ba), p.LpBoundAlg, p.MaxLPAlg, p.MfbAlg)
		if !ok {
			return Relation{}, false
		}
		for _, lp := range lps {
			if lp == 0 {
				continue
			}
			ideal, ok := makeLargeIdeal(poly, a, b, lp)
			if !ok {
				return Relation{}, false
			}
			algLPs = append(algLPs, ideal)
		}
	}

	return Relation{
		A:       a,
		B:       b,
		RatCols: append([]int32(nil), ratCols...),
		AlgCols: append([]int32(nil), algCols...),
		RatLPs:  ratLPs,
		AlgLPs:  algLPs,
	}, true
}

// cofactorCandidate attempts to factor G(a, b) and F(a, b) for coprime (a, b) with b > 0.
// Used by the line sieve.
func cofactorCandidate(
	poly *Poly,
	fb *FactorBases,
	p *Params,
	a, b int64,
	ratSmall []int32, ratHits []int32,
	algSmall []int32, algHits []int32,
	sq *LatticeQ,
	sc *cofactorScratch,
) (Relation, bool) {
	gVal := evalGInto(&sc.val, poly, a, b, sc)
	if gVal.Sign() == 0 {
		return Relation{}, false
	}
	temp := sc.temp.Abs(gVal)
	var ns normState
	initNorm(&ns, temp)
	ratCols := sc.ratCols[:0]

	divRat := func(idx int32, pr uint64) {
		for m := ns.divAll(pr, temp, sc); m > 0; m-- {
			ratCols = append(ratCols, idx)
		}
	}

	for _, idx := range ratSmall {
		id := fb.RFB[idx]
		pr := uint64(id.P)
		if id.R == id.P {
			if uint64(b)%pr == 0 {
				divRat(idx, pr)
			}
		} else if modInt64(a, pr) == (uint64(b)%pr)*uint64(id.R)%pr {
			divRat(idx, pr)
		}
	}
	for _, idx := range ratHits {
		divRat(idx, uint64(fb.RFB[idx].P))
	}

	remRat, fits := ns.uint64Val()
	if !fits {
		return Relation{}, false
	}
	var ratLPs [2]uint64
	if remRat != 1 {
		var ok bool
		ratLPs, ok = splitCofactor(remRat, uint64(p.Br), p.LpBoundRat, p.MaxLPRat, p.MfbRat)
		if !ok {
			return Relation{}, false
		}
	}

	fVal := evalFInto(&sc.val, poly, a, b, sc)
	if fVal.Sign() == 0 {
		return Relation{}, false
	}
	temp.Abs(fVal)
	initNorm(&ns, temp)
	algCols := sc.algCols[:0]
	var algLPs []LargeIdeal

	divAlg := func(idx int32, pr uint64) {
		for m := ns.divAll(pr, temp, sc); m > 0; m-- {
			algCols = append(algCols, idx)
		}
	}

	if sq != nil && sq.Q != 0 {
		m := ns.divAll(uint64(sq.Q), temp, sc)
		if sq.AFBIdx >= 0 {
			for ; m > 0; m-- {
				algCols = append(algCols, sq.AFBIdx)
			}
		} else {
			sqIdeal := LargeIdeal{P: uint64(sq.Q), R: uint64(sq.R)}
			for ; m > 0; m-- {
				algLPs = append(algLPs, sqIdeal)
			}
		}
	}

	for _, idx := range algSmall {
		id := fb.AFB[idx]
		pr := uint64(id.P)
		if id.R == id.P {
			if uint64(b)%pr == 0 {
				divAlg(idx, pr)
			}
		} else if modInt64(a, pr) == (uint64(b)%pr)*uint64(id.R)%pr {
			divAlg(idx, pr)
		}
	}
	for _, idx := range algHits {
		divAlg(idx, uint64(fb.AFB[idx].P))
	}

	remAlg, fits := ns.uint64Val()
	if !fits {
		return Relation{}, false
	}
	if remAlg != 1 {
		lps, ok := splitCofactor(remAlg, uint64(p.Ba), p.LpBoundAlg, p.MaxLPAlg, p.MfbAlg)
		if !ok {
			return Relation{}, false
		}
		for _, lp := range lps {
			if lp == 0 {
				continue
			}
			ideal, ok := makeLargeIdeal(poly, a, b, lp)
			if !ok {
				return Relation{}, false
			}
			algLPs = append(algLPs, ideal)
		}
	}

	return Relation{
		A:       a,
		B:       b,
		RatCols: append([]int32(nil), ratCols...),
		AlgCols: append([]int32(nil), algCols...),
		RatLPs:  ratLPs,
		AlgLPs:  algLPs,
	}, true
}

func makeLargeIdeal(poly *Poly, a, b int64, lp uint64) (LargeIdeal, bool) {
	bm := uint64(b) % lp
	if bm == 0 {
		return LargeIdeal{P: lp, R: lp}, true
	}
	invB, ok := arith.InvMod64(bm, lp)
	if !ok {
		return LargeIdeal{}, false
	}
	r := arith.MulMod64(modInt64(a, lp), invB, lp)
	if isMultipleRoot(poly, r, lp) {
		return LargeIdeal{}, false
	}
	return LargeIdeal{P: lp, R: r}, true
}

// cofactorPair reconstructs a Relation from a saved (a, b) pair without sieve
// hits, accepting algebraic cofactors that include a special-q below qcMin.
func cofactorPair(
	poly *Poly,
	fb *FactorBases,
	p *Params,
	qcMin uint64,
	a, b int64,
	isFree bool,
	sc *cofactorScratch,
) (Relation, bool) {
	if isFree {
		if a <= 1 || b != 0 || !arith.IsPrime64(uint64(a)) {
			return Relation{}, false
		}
		up := uint64(a)
		var ratCols []int32
		var ratLPs [2]uint64
		if idx, ok := fb.RatLookup[uint32(up)]; ok && up <= uint64(p.Br) {
			ratCols = []int32{idx}
		} else if up <= p.LpBoundRat {
			ratLPs = [2]uint64{up, 0}
		} else {
			return Relation{}, false
		}
		fp := poly.CoeffsModP(up)
		if len(fp)-1 != poly.Degree {
			return Relation{}, false
		}
		roots := arith.PolyRootsModP(fp, up)
		if len(roots) != poly.Degree {
			return Relation{}, false
		}
		var algCols []int32
		var algLPs []LargeIdeal
		for _, rt := range roots {
			if isMultipleRoot(poly, rt, up) {
				return Relation{}, false
			}
			if up <= uint64(p.Ba) {
				idx, ok := fb.AlgLookup[algKey(uint32(up), uint32(rt))]
				if !ok {
					return Relation{}, false
				}
				algCols = append(algCols, idx)
			} else if up < qcMin {
				algLPs = append(algLPs, LargeIdeal{P: up, R: rt})
			} else {
				return Relation{}, false
			}
		}
		return Relation{
			A:       a,
			B:       0,
			RatCols: ratCols,
			AlgCols: algCols,
			RatLPs:  ratLPs,
			AlgLPs:  algLPs,
			IsFree:  true,
		}, true
	}

	if b <= 0 {
		return Relation{}, false
	}
	ua := uint64(a)
	if a < 0 {
		ua = uint64(-a)
	}
	if arith.Gcd64(ua, uint64(b)) != 1 {
		return Relation{}, false
	}

	gVal := evalGInto(&sc.val, poly, a, b, sc)
	if gVal.Sign() == 0 {
		return Relation{}, false
	}
	temp := sc.temp.Abs(gVal)
	var ns normState
	initNorm(&ns, temp)
	ratCols := sc.ratCols[:0]

	for idx, id := range fb.RFB {
		pr := uint64(id.P)
		hit := false
		if id.R == id.P {
			hit = uint64(b)%pr == 0
		} else {
			hit = modInt64(a, pr) == (uint64(b)%pr)*uint64(id.R)%pr
		}
		if hit {
			for m := ns.divAll(pr, temp, sc); m > 0; m-- {
				ratCols = append(ratCols, int32(idx))
			}
		}
	}
	remRat, fits := ns.uint64Val()
	if !fits {
		return Relation{}, false
	}
	var ratLPs [2]uint64
	if remRat != 1 {
		var ok bool
		ratLPs, ok = splitCofactor(remRat, uint64(p.Br), p.LpBoundRat, p.MaxLPRat, p.MfbRat)
		if !ok {
			return Relation{}, false
		}
	}

	fVal := evalFInto(&sc.val, poly, a, b, sc)
	if fVal.Sign() == 0 {
		return Relation{}, false
	}
	temp.Abs(fVal)
	initNorm(&ns, temp)
	algCols := sc.algCols[:0]

	for idx, id := range fb.AFB {
		pr := uint64(id.P)
		hit := false
		if id.R == id.P {
			hit = uint64(b)%pr == 0
		} else {
			hit = modInt64(a, pr) == (uint64(b)%pr)*uint64(id.R)%pr
		}
		if hit {
			for m := ns.divAll(pr, temp, sc); m > 0; m-- {
				algCols = append(algCols, int32(idx))
			}
		}
	}

	lps, ok := factorSavedAlgCofactor(&ns, temp, uint64(p.Ba), p.LpBoundAlg, qcMin, p.MaxLPAlg)
	if !ok {
		return Relation{}, false
	}
	var algLPs []LargeIdeal
	for _, lp := range lps {
		ideal, ok := makeLargeIdeal(poly, a, b, lp)
		if !ok {
			return Relation{}, false
		}
		algLPs = append(algLPs, ideal)
	}

	return Relation{
		A:       a,
		B:       b,
		RatCols: append([]int32(nil), ratCols...),
		AlgCols: append([]int32(nil), algCols...),
		RatLPs:  ratLPs,
		AlgLPs:  algLPs,
	}, true
}

// factorSavedAlgCofactor factors the algebraic cofactor of a saved relation,
// allowing up to maxLP primes <= lpBound plus one optional special-q < qcMin.
func factorSavedAlgCofactor(ns *normState, temp *big.Int, ba, lpBound, qcMin uint64, maxLP int) ([]uint64, bool) {
	if u, ok := ns.uint64Val(); ok {
		if u == 1 {
			return nil, true
		}
		if lps, ok := splitCofactor(u, ba, lpBound, maxLP, 0); ok {
			if lps[1] != 0 {
				return []uint64{lps[0], lps[1]}, true
			}
			return []uint64{lps[0]}, true
		}
		var facs []uint64
		if !factorUint64Into(u, &facs) {
			return nil, false
		}
		return checkSavedAlgPrimes(facs, lpBound, qcMin, maxLP)
	}

	// Cofactor > 64 bits (only possible when special-q * lp1 * lp2 >= 2^64).
	var cur big.Int
	if ns.use128 {
		cur.SetUint64(ns.hi)
		cur.Lsh(&cur, 64)
		cur.Or(&cur, new(big.Int).SetUint64(ns.lo))
	} else {
		cur.Set(temp)
	}
	if cur.BitLen() > 96 {
		return nil, false
	}
	var facs []uint64
	for !cur.IsUint64() {
		d := rhoBig(&cur)
		if d == nil || !d.IsUint64() {
			return nil, false
		}
		du := d.Uint64()
		if !factorUint64Into(du, &facs) {
			return nil, false
		}
		cur.Quo(&cur, d)
	}
	if cur.Uint64() > 1 {
		if !factorUint64Into(cur.Uint64(), &facs) {
			return nil, false
		}
	}
	return checkSavedAlgPrimes(facs, lpBound, qcMin, maxLP)
}

func factorUint64Into(n uint64, dst *[]uint64) bool {
	if n == 1 {
		return true
	}
	if arith.IsPrime64(n) {
		*dst = append(*dst, n)
		return true
	}
	d := arith.Rho64(n)
	if d == 0 {
		return false
	}
	return factorUint64Into(d, dst) && factorUint64Into(n/d, dst)
}

func checkSavedAlgPrimes(facs []uint64, lpBound, qcMin uint64, maxLP int) ([]uint64, bool) {
	if len(facs) > maxLP+1 {
		return nil, false
	}
	aboveLP := 0
	for _, f := range facs {
		if f >= qcMin {
			return nil, false
		}
		if f > lpBound {
			aboveLP++
		}
	}
	if aboveLP > 1 {
		return nil, false
	}
	return facs, true
}

func rhoBig(n *big.Int) *big.Int {
	if n.ProbablyPrime(20) {
		return new(big.Int).Set(n)
	}
	var x, y, c, diff, g big.Int
	for ci := int64(1); ci <= 8; ci++ {
		c.SetInt64(ci)
		x.SetInt64(2)
		y.SetInt64(2)
		g.SetInt64(1)
		steps := 0
		for g.Cmp(bigOne) == 0 && steps < 1<<16 {
			x.Mul(&x, &x).Add(&x, &c).Mod(&x, n)
			y.Mul(&y, &y).Add(&y, &c).Mod(&y, n)
			y.Mul(&y, &y).Add(&y, &c).Mod(&y, n)
			diff.Sub(&x, &y).Abs(&diff)
			g.GCD(nil, nil, &diff, n)
			steps++
		}
		if g.Cmp(bigOne) > 0 && g.Cmp(n) < 0 {
			return new(big.Int).Set(&g)
		}
	}
	return nil
}

var bigOne = big.NewInt(1)

// isMultipleRoot reports whether r is a multiple root of f mod p.
func isMultipleRoot(poly *Poly, r, p uint64) bool {
	dfp := arith.PolyDerivModP(poly.CoeffsModP(p), p)
	return arith.PolyEvalModP(dfp, r, p) == 0
}

type cofactorScratch struct {
	val, temp    big.Int
	t1, t2       big.Int
	bigA, bigB   big.Int
	term         big.Int
	bigP, q, rem big.Int
	aPows        []big.Int
	bPows        []big.Int
	ratCols      []int32
	algCols      []int32
}

func newCofactorScratch(deg int) *cofactorScratch {
	return &cofactorScratch{
		aPows: make([]big.Int, deg+1),
		bPows: make([]big.Int, deg+1),
	}
}

func evalGInto(dst *big.Int, poly *Poly, a, b int64, sc *cofactorScratch) *big.Int {
	sc.t1.SetInt64(a)
	sc.t1.Mul(&sc.t1, poly.A1)
	sc.t2.SetInt64(b)
	sc.t2.Mul(&sc.t2, poly.A0)
	return dst.Sub(&sc.t1, &sc.t2)
}

func evalFInto(dst *big.Int, poly *Poly, a, b int64, sc *cofactorScratch) *big.Int {
	d := poly.Degree
	sc.aPows[0].SetInt64(1)
	sc.bPows[0].SetInt64(1)
	sc.bigA.SetInt64(a)
	sc.bigB.SetInt64(b)
	for i := 1; i <= d; i++ {
		sc.aPows[i].Mul(&sc.aPows[i-1], &sc.bigA)
		sc.bPows[i].Mul(&sc.bPows[i-1], &sc.bigB)
	}
	dst.SetInt64(0)
	for i := 0; i <= d; i++ {
		sc.term.Mul(poly.Coeffs[i], &sc.aPows[i])
		sc.term.Mul(&sc.term, &sc.bPows[d-i])
		dst.Add(dst, &sc.term)
	}
	return dst
}

package gnfs

import (
	"math"
	"math/big"

	"github.com/dgryski/go-factor/internal/arith"
)

const (
	// ColAllOnes is the matrix column forcing |S| even.
	ColAllOnes = 0
	// ColRatSign is the matrix column for G(a, b) < 0.
	ColRatSign = 1
	// ColAlgSign is the matrix column for F(a, b) < 0.
	ColAlgSign = 2
	// NumSignCols is the number of fixed parity/sign columns before RFB.
	NumSignCols = 3
)

// RatIdeal represents a prime ideal on the rational side.
type RatIdeal struct {
	P    uint32
	R    uint32 // A0 * A1^-1 mod P, or P for projective (P | A1)
	LogP byte
}

// AlgIdeal represents a degree-1 prime ideal (p, r) on the algebraic side.
type AlgIdeal struct {
	P    uint32
	R    uint32 // root mod P, or P for the projective root (P | c_d)
	LogP byte
}

// QChar is a quadratic character (q, s) where f(s) ≡ 0 (mod q), f'(s) ≠ 0 (mod q).
type QChar struct {
	Q uint64
	S uint64
}

// FactorBases holds the rational factor base, algebraic factor base,
// and quadratic character base for a polynomial pair.
type FactorBases struct {
	RFB       []RatIdeal
	AFB       []AlgIdeal
	QCB       []QChar
	BadCount  int
	ColRFB    int
	ColAFB    int
	ColQCB    int
	NumFixed  int
	LogScale  float64
	RatLookup map[uint32]int32 // P -> RFB index
	AlgLookup map[uint64]int32 // (uint64(P)<<32)|uint64(R) -> AFB index
}

func algKey(p, r uint32) uint64 {
	return uint64(p)<<32 | uint64(r)
}

// BuildFactorBases builds the rational factor base (primes <= br), algebraic
// factor base (primes <= ba), and quadratic character base (primes > qcMin).
func BuildFactorBases(poly *Poly, br, ba int64, qcMin uint64, numQCs int) *FactorBases {
	fb := &FactorBases{
		LogScale:  1.0,
		RatLookup: make(map[uint32]int32),
		AlgLookup: make(map[uint64]int32),
	}

	maxB := br
	if ba > maxB {
		maxB = ba
	}
	primes := arith.PrimesUpTo(maxB)

	var bigP, r, a0Mod, a1Mod big.Int
	for _, prime := range primes {
		up := uint64(prime)
		p32 := uint32(prime)
		logP := byte(math.Round(math.Log2(float64(prime))))
		if logP == 0 {
			logP = 1
		}
		bigP.SetUint64(up)

		if prime <= br {
			a1m := a1Mod.Mod(poly.A1, &bigP).Uint64()
			var root uint32
			if a1m == 0 {
				root = p32 // projective (P | A1)
			} else {
				a0m := a0Mod.Mod(poly.A0, &bigP).Uint64()
				root = uint32(arith.MulMod64(a0m, arith.InvMod(a1m, up), up))
			}
			fb.RatLookup[p32] = int32(len(fb.RFB))
			fb.RFB = append(fb.RFB, RatIdeal{P: p32, R: root, LogP: logP})
		}

		if prime <= ba {
			fp := poly.CoeffsModP(up)
			roots := arith.PolyRootsModP(fp, up)
			cdDivides := r.Mod(poly.Coeffs[poly.Degree], &bigP).Sign() == 0
			var dfp []uint64
			if !cdDivides {
				dfp = arith.PolyDerivModP(fp, up)
			}
			if cdDivides {
				fb.AlgLookup[algKey(p32, p32)] = int32(len(fb.AFB))
				fb.AFB = append(fb.AFB, AlgIdeal{P: p32, R: p32, LogP: logP})
				fb.BadCount++
			}
			for _, rt := range roots {
				r32 := uint32(rt)
				fb.AlgLookup[algKey(p32, r32)] = int32(len(fb.AFB))
				fb.AFB = append(fb.AFB, AlgIdeal{P: p32, R: r32, LogP: logP})
				if cdDivides || arith.PolyEvalModP(dfp, rt, up) == 0 {
					fb.BadCount++
				}
			}
		}
	}

	wantQCs := numQCs + fb.BadCount
	if wantQCs < 64 {
		wantQCs = 64
	}
	if wantQCs > 128 {
		wantQCs = 128
	}

	q := qcMin + 1
	if q < 3 {
		q = 3
	}
	if q&1 == 0 {
		q++
	}
	for len(fb.QCB) < wantQCs {
		if arith.IsPrime64(q) {
			bigP.SetUint64(q)
			if r.Mod(poly.Coeffs[poly.Degree], &bigP).Sign() != 0 {
				fp := poly.CoeffsModP(q)
				dfp := arith.PolyDerivModP(fp, q)
				for _, s := range arith.PolyRootsModP(fp, q) {
					if arith.PolyEvalModP(dfp, s, q) != 0 {
						fb.QCB = append(fb.QCB, QChar{Q: q, S: s})
						if len(fb.QCB) >= wantQCs {
							break
						}
					}
				}
			}
		}
		q += 2
	}

	fb.ColRFB = NumSignCols
	fb.ColAFB = fb.ColRFB + len(fb.RFB)
	fb.ColQCB = fb.ColAFB + len(fb.AFB)
	fb.NumFixed = fb.ColQCB + len(fb.QCB)
	return fb
}

// QCBits appends the matrix column indices for quadratic characters of (a, b)
// whose Legendre symbol ((a - b*s)/q) == -1.
func (fb *FactorBases) QCBits(a, b int64, dst []int32) []int32 {
	for i, qc := range fb.QCB {
		q := qc.Q
		am := modInt64(a, q)
		bm := modInt64(b, q)
		bs := arith.MulMod64(bm, qc.S, q)
		diff := am + q - bs
		if diff >= q {
			diff -= q
		}
		if diff != 0 && arith.PowMod64(diff, (q-1)>>1, q) == q-1 {
			dst = append(dst, int32(fb.ColQCB+i))
		}
	}
	return dst
}

func modInt64(a int64, m uint64) uint64 {
	if a >= 0 {
		return uint64(a) % m
	}
	r := uint64(-a) % m
	if r == 0 {
		return 0
	}
	return m - r
}

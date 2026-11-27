package gnfs

import (
	"math/big"
	"math/rand"
	"testing"

	"github.com/dgryski/go-factor/internal/arith"
	"github.com/dgryski/go-factor/internal/testutil"
)

func TestBuildFactorBasesAndQCB(t *testing.T) {
	rnd := rand.New(rand.NewSource(11))
	n := new(big.Int).Mul(testutil.RandomPrime(rnd, 60), testutil.RandomPrime(rnd, 60))
	poly, fac := BaseMPoly(n, 3, 6)
	if fac != nil {
		t.Skip("got lucky factor")
	}

	const br, ba, qcMin = 500, 500, 2000
	fb := BuildFactorBases(poly, br, ba, qcMin, 64)
	if len(fb.RFB) == 0 || len(fb.AFB) == 0 || len(fb.QCB) < 64 {
		t.Fatalf("unexpected base sizes: rfb=%d afb=%d qcb=%d", len(fb.RFB), len(fb.AFB), len(fb.QCB))
	}

	// Verify every AFB entry is a genuine root (or projective root) of f mod p.
	for _, id := range fb.AFB {
		p := uint64(id.P)
		fp := poly.CoeffsModP(p)
		if id.R == id.P {
			if len(fp)-1 >= poly.Degree {
				t.Fatalf("projective root recorded for p=%d not dividing c_d", p)
			}
		} else {
			if arith.PolyEvalModP(fp, uint64(id.R), p) != 0 {
				t.Fatalf("AFB (%d, %d) is not a root", id.P, id.R)
			}
		}
	}

	// Verify every QCB entry has q > qcMin, simple root, and q !| c_d.
	for _, qc := range fb.QCB {
		if qc.Q <= qcMin {
			t.Fatalf("QC prime %d <= %d", qc.Q, qcMin)
		}
		fp := poly.CoeffsModP(qc.Q)
		if arith.PolyEvalModP(fp, qc.S, qc.Q) != 0 {
			t.Fatalf("QC (%d, %d) is not a root", qc.Q, qc.S)
		}
		dfp := arith.PolyDerivModP(fp, qc.Q)
		if arith.PolyEvalModP(dfp, qc.S, qc.Q) == 0 {
			t.Fatalf("QC (%d, %d) is not simple", qc.Q, qc.S)
		}
	}
}

func TestExactBadIdealsRound2(t *testing.T) {
	// Cover non maximal orders, mixed inertia/ramification, and non monic projective primes.
	cases := []*Poly{
		{
			Degree: 3,
			Coeffs: []*big.Int{big.NewInt(-4), big.NewInt(0), big.NewInt(0), big.NewInt(1)},
			A1:     big.NewInt(1),
			A0:     big.NewInt(7),
		},
		{
			Degree: 4,
			Coeffs: []*big.Int{big.NewInt(12), big.NewInt(0), big.NewInt(0), big.NewInt(0), big.NewInt(1)},
			A1:     big.NewInt(1),
			A0:     big.NewInt(11),
		},
		{
			Degree: 4,
			Coeffs: []*big.Int{big.NewInt(7), big.NewInt(-2), big.NewInt(3), big.NewInt(5), big.NewInt(6)},
			A1:     big.NewInt(1),
			A0:     big.NewInt(13),
		},
		{
			Degree: 5,
			Coeffs: []*big.Int{big.NewInt(3), big.NewInt(2), big.NewInt(-1), big.NewInt(1), big.NewInt(-2), big.NewInt(1)},
			A1:     big.NewInt(1),
			A0:     big.NewInt(17),
		},
	}

	rnd := rand.New(rand.NewSource(313))
	for ci, poly := range cases {
		if !IsIrreducible(poly) {
			t.Fatalf("case %d polynomial is reducible", ci)
		}
		fb := BuildFactorBases(poly, 100, 100, 1000, 64)
		if len(fb.BadDecomps) == 0 {
			t.Fatalf("case %d expected bad prime decompositions", ci)
		}
		d := poly.Degree
		cd := poly.Coeffs[d]
		cdPow := new(big.Int).Exp(cd, big.NewInt(int64(d-1)), nil)

		for _, bd := range fb.BadDecomps {
			efSum := 0
			for _, id := range bd.Ideals {
				if id.E <= 0 || id.F <= 0 {
					t.Fatalf("case %d p=%d invalid (e,f)=(%d,%d)", ci, bd.P, id.E, id.F)
				}
				efSum += id.E * id.F
			}
			if efSum != d {
				t.Fatalf("case %d p=%d sum(e*f)=%d, want %d", ci, bd.P, efSum, d)
			}

			bigP := big.NewInt(int64(bd.P))
			var vals []int
			var rem big.Int
			for trial := 0; trial < 120; trial++ {
				a := int64(rnd.Intn(401) - 200)
				b := int64(1 + rnd.Intn(200))
				ua := uint64(a)
				if a < 0 {
					ua = uint64(-a)
				}
				if arith.Gcd64(ua, uint64(b)) != 1 {
					continue
				}
				vals = bd.valuations(a, b, vals[:0])
				gotNormVal := 0
				for j, v := range vals {
					gotNormVal += bd.Ideals[j].F * v
				}
				norm := new(big.Int).Mul(cdPow, poly.EvalF(a, b))
				norm.Abs(norm)
				wantNormVal := 0
				for norm.Sign() > 0 {
					norm.QuoRem(norm, bigP, &rem)
					if rem.Sign() != 0 {
						break
					}
					wantNormVal++
				}
				if gotNormVal != wantNormVal {
					t.Fatalf("case %d p=%d (a,b)=(%d,%d): sum(f*v)=%d, v_p(norm)=%d", ci, bd.P, a, b, gotNormVal, wantNormVal)
				}
			}
		}
	}
}

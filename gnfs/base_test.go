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

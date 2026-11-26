package gnfs

import (
	"math/big"
	"math/rand"
	"testing"

	"github.com/dgryski/go-factor/internal/arith"
	"github.com/dgryski/go-factor/internal/testutil"
)

func TestFrankeKleinjungMatchesBruteForce(t *testing.T) {
	const I, J = 32, 24
	W := 2 * I
	total := (J - 1) * W
	numBlocks := (total + blockSize - 1) >> blockShift

	primes := arith.PrimesUpTo(400)
	for _, p64 := range primes {
		p := uint32(p64)
		if int(p) < W {
			continue
		}
		for lam := uint32(0); lam < p; lam += 7 {
			buckets := make([][]uint64, numBlocks)
			addFKHits(p, lam, I, J, uint64(p)<<16, buckets)
			got := make(map[[2]int]bool)
			for bi, bkt := range buckets {
				for _, e := range bkt {
					gPos := bi*blockSize + int(uint16(e))
					j := 1 + gPos/W
					i := (gPos % W) - I
					pos := [2]int{i, j}
					if got[pos] {
						t.Fatalf("duplicate hit (%d,%d) for p=%d lam=%d", i, j, p, lam)
					}
					got[pos] = true
				}
			}

			want := make(map[[2]int]bool)
			for j := 1; j < J; j++ {
				for i := -I; i < I; i++ {
					diff := (int64(i) - int64(lam)*int64(j)) % int64(p)
					if diff == 0 {
						want[[2]int{i, j}] = true
					}
				}
			}
			if len(got) != len(want) {
				t.Fatalf("p=%d lam=%d: got %d hits, want %d", p, lam, len(got), len(want))
			}
			for pos := range want {
				if !got[pos] {
					t.Fatalf("p=%d lam=%d: missing hit %v", p, lam, pos)
				}
			}
		}
	}
}

func TestLatticeRoundTripAndRelations(t *testing.T) {
	rnd := rand.New(rand.NewSource(13))
	n := new(big.Int).Mul(testutil.RandomPrime(rnd, 50), testutil.RandomPrime(rnd, 50))
	poly, fac := BaseMPoly(n, 3, 6)
	if fac != nil {
		t.Skip("got lucky factor")
	}
	fb := BuildFactorBases(poly, 800, 1200, 50000, 64)

	// Pick a special-q > 128 and check reduced lattice basis round trip and hits.
	const I, J = 64, 64
	var lq LatticeQ
	foundQ := false
	for _, id := range fb.AFB {
		if id.P >= 150 && id.R != id.P {
			lq = ReduceLattice(id.P, id.R, poly.Skew, I, J, fb.AlgLookup[algKey(id.P, id.R)])
			foundQ = true
			break
		}
	}
	if !foundQ {
		t.Fatal("no special-q found")
	}
	det := lq.U1*lq.V2 - lq.U2*lq.V1
	if det != int64(lq.Q) && det != -int64(lq.Q) {
		t.Fatalf("bad det=%d for q=%d", det, lq.Q)
	}

	// Check lattice round trip for random (i, j) and root transform accuracy.
	for k := 0; k < 50; k++ {
		i := int64(rnd.Intn(2*I) - I)
		j := int64(1 + rnd.Intn(J-1))
		a := i*lq.U1 + j*lq.V1
		b := i*lq.U2 + j*lq.V2

		// Invert reduced basis: (i, j) = (v2*a - v1*b, -u2*a + u1*b) / det.
		iBack := (lq.V2*a - lq.V1*b) / det
		jBack := (-lq.U2*a + lq.U1*b) / det
		if iBack != i || jBack != j {
			t.Fatalf("round trip failed: (%d,%d) -> (%d,%d) -> (%d,%d)", i, j, a, b, iBack, jBack)
		}

		// Verify q | F(a, b).
		fVal := poly.EvalF(a, b)
		if new(big.Int).Mod(fVal, big.NewInt(int64(lq.Q))).Sign() != 0 {
			t.Fatalf("q=%d does not divide F(%d,%d)", lq.Q, a, b)
		}

		// Verify transformed root matches p | F(a, b) for factor base primes != q.
		for _, id := range fb.AFB[:40] {
			if id.P == lq.Q {
				continue
			}
			lam, vert, _ := transformRoot(id.P, id.R, &lq)
			p64 := int64(id.P)
			var hitTrans bool
			if vert {
				hitTrans = j%p64 == 0
			} else {
				hitTrans = (i-int64(lam)*j)%p64 == 0
			}
			var hitDirect bool
			if id.R == id.P {
				hitDirect = b%p64 == 0
			} else {
				hitDirect = (a-int64(id.R)*b)%p64 == 0
			}
			if hitTrans != hitDirect {
				t.Fatalf("ideal (%d,%d) mismatch at (i,j)=(%d,%d): hitTrans=%v hitDirect=%v", id.P, id.R, i, j, hitTrans, hitDirect)
			}
		}
	}

	// Verify line sieve and lattice sieve produce exact factorizations.
	verifyRel := func(rel Relation) {
		gProd := big.NewInt(1)
		for _, idx := range rel.RatCols {
			gProd.Mul(gProd, big.NewInt(int64(fb.RFB[idx].P)))
		}
		for _, lp := range rel.RatLPs {
			if lp != 0 {
				gProd.Mul(gProd, new(big.Int).SetUint64(lp))
			}
		}
		gWant := new(big.Int).Abs(poly.EvalG(rel.A, rel.B))
		if gProd.Cmp(gWant) != 0 {
			t.Fatalf("G(%d,%d) product=%v, want %v", rel.A, rel.B, gProd, gWant)
		}

		fProd := big.NewInt(1)
		for _, idx := range rel.AlgCols {
			fProd.Mul(fProd, big.NewInt(int64(fb.AFB[idx].P)))
		}
		for _, lp := range rel.AlgLPs {
			fProd.Mul(fProd, new(big.Int).SetUint64(lp.P))
		}
		fWant := new(big.Int).Abs(poly.EvalF(rel.A, rel.B))
		if fProd.Cmp(fWant) != 0 {
			t.Fatalf("F(%d,%d) product=%v, want %v", rel.A, rel.B, fProd, fWant)
		}
	}

	pLine := &Params{
		Degree: 3, Br: 800, Ba: 1200,
		LpBoundRat: 20000, LpBoundAlg: 30000,
		MfbRat: 30, MfbAlg: 30, MaxLPRat: 2, MaxLPAlg: 2,
		LineA: 8192, LineB: 200,
		SlackRat: 22, SlackAlg: 24,
	}
	var lineRels []Relation
	sieveLine(poly, fb, pLine, 0, 0, func(r Relation) { lineRels = append(lineRels, r) }, func() bool { return len(lineRels) >= 50 })
	if len(lineRels) == 0 {
		t.Fatal("line sieve found 0 relations")
	}
	for _, r := range lineRels {
		verifyRel(r)
	}

	pLat := &Params{
		Degree: 3, Br: 800, Ba: 1200,
		LpBoundRat: 20000, LpBoundAlg: 30000,
		MfbRat: 30, MfbAlg: 30, MaxLPRat: 2, MaxLPAlg: 2,
		QMin: 500, QMax: 1500, IHalf: 256, JMax: 256,
		SlackRat: 22, SlackAlg: 24,
	}
	var latRels []Relation
	sieveLattice(poly, fb, pLat, func(r Relation) { latRels = append(latRels, r) }, func() bool { return len(latRels) >= 50 })
	if len(latRels) == 0 {
		t.Fatal("lattice sieve found 0 relations")
	}
	for _, r := range latRels {
		verifyRel(r)
	}
}

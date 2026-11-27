package gnfs

import (
	"math/big"
	"math/rand"
	"testing"

	"github.com/dgryski/go-factor/internal/arith"
	"github.com/dgryski/go-factor/internal/testutil"
)

func TestFrankeKleinjungMatchesBruteForce(t *testing.T) {
	for _, I := range []int{30, 31, 32} {
		const J = 24
		W := 2 * I
		total := (J - 1) * W
		numBlocks := (total + blockSize - 1) >> blockShift

		primes := arith.PrimesUpTo(400)
		for _, forceOverflow := range []bool{false, true} {
			for _, p64 := range primes {
				p := uint32(p64)
				if int(p) < W {
					continue
				}
				for lam := uint32(0); lam < p; lam += 7 {
					arena := newFKArena(numBlocks, 1, 64)
					if forceOverflow {
						arena.stride = 2
					}
					addFKHitsArena(p, lam, I, J, 7, 5, arena)
					arena.endGroup(0)

					got := make(map[[2]int]bool)
					for bi := 0; bi < numBlocks; bi++ {
						base := bi << arena.shift
						for _, e := range arena.data[base : base+int(arena.counts[bi])] {
							if int32(e>>blockShift) != 7 {
								t.Fatalf("bad idx %d", e>>blockShift)
							}
							gPos := bi*blockSize + int(e&blockMask)
							j := 1 + gPos/W
							i := (gPos % W) - I
							if (i|j)&1 == 0 {
								t.Fatalf("emitted both-even (%d,%d) for I=%d p=%d lam=%d", i, j, I, p, lam)
							}
							pos := [2]int{i, j}
							if got[pos] {
								t.Fatalf("duplicate hit (%d,%d) for I=%d p=%d lam=%d", i, j, I, p, lam)
							}
							got[pos] = true
						}
						for _, ov := range arena.overflow[bi] {
							if ov.idx != 7 || ov.logP != 5 {
								t.Fatalf("bad overflow entry %+v", ov)
							}
							gPos := bi*blockSize + int(ov.off)
							j := 1 + gPos/W
							i := (gPos % W) - I
							if (i|j)&1 == 0 {
								t.Fatalf("emitted both-even overflow (%d,%d)", i, j)
							}
							pos := [2]int{i, j}
							if got[pos] {
								t.Fatalf("duplicate overflow hit (%d,%d)", i, j)
							}
							got[pos] = true
						}
					}

					want := make(map[[2]int]bool)
					for j := 1; j < J; j++ {
						for i := -I; i < I; i++ {
							if (i|j)&1 == 0 {
								continue
							}
							diff := (int64(i) - int64(lam)*int64(j)) % int64(p)
							if diff == 0 {
								want[[2]int{i, j}] = true
							}
						}
					}
					if len(got) != len(want) {
						t.Fatalf("I=%d p=%d lam=%d over=%v: got %d hits, want %d", I, p, lam, forceOverflow, len(got), len(want))
					}
					for pos := range want {
						if !got[pos] {
							t.Fatalf("I=%d p=%d lam=%d: missing hit %v", I, p, lam, pos)
						}
					}
				}
			}
		}
	}
}

func TestSWARByteScans(t *testing.T) {
	rng := rand.New(rand.NewSource(99))
	ab := make([]byte, 259)
	rb := make([]byte, 259)
	for i := range ab {
		ab[i] = byte(rng.Intn(256))
		rb[i] = byte(rng.Intn(256))
	}
	ab[0], rb[0] = 0, 255
	ab[7], rb[7] = 128, 128
	ab[8], rb[8] = 127, 129
	ab[258], rb[258] = 255, 255

	var got []int
	for aT := -1; aT <= 257; aT += 3 {
		wantHas := false
		for _, v := range ab {
			if int(v) >= aT {
				wantHas = true
				break
			}
		}
		if hasByteGe(ab, aT) != wantHas {
			t.Fatalf("hasByteGe(aT=%d) = %v, want %v", aT, !wantHas, wantHas)
		}
		for rT := -1; rT <= 257; rT += 5 {
			got = appendBothGe(got[:0], ab, rb, aT, rT)
			var want []int
			for idx := range ab {
				if int(ab[idx]) >= aT && int(rb[idx]) >= rT {
					want = append(want, idx)
				}
			}
			if len(got) != len(want) {
				t.Fatalf("appendBothGe(%d,%d): got %v, want %v", aT, rT, got, want)
			}
			for k := range got {
				if got[k] != want[k] {
					t.Fatalf("appendBothGe(%d,%d)[%d] = %d, want %d", aT, rT, k, got[k], want[k])
				}
			}
		}
	}
}

func TestSieveBlockLinesMatchesBruteForce(t *testing.T) {
	lines := []linePrime{
		{idx: 0, p: 29, lambda: 0, logP: 5},
		{idx: 1, p: 31, lambda: 17, logP: 5},
		{idx: 2, p: 37, lambda: 36, logP: 5},
	}
	verts := []vertPrime{
		{idx: 3, p: 5, logP: 2},
	}
	for _, I := range []int{64, 60} { // 64 hits aligned fast path, 60 hits general path
		W := 2 * I
		const J = 20
		total := (J - 1) * W
		got := make([]byte, total)
		step := 512
		for bStart := 0; bStart < total; bStart += step {
			bEnd := bStart + step
			if bEnd > total {
				bEnd = total
			}
			sieveBlockLines(got[bStart:bEnd], bStart, bEnd, I, W, lines, verts, 3)
		}

		for j := 1; j < J; j++ {
			for i := -I; i < I; i++ {
				pos := (j-1)*W + (i + I)
				if (i|j)&1 != 0 {
					want := byte(3)
					for _, vp := range verts {
						if j%int(vp.p) == 0 {
							want += vp.logP
						}
					}
					for _, lp := range lines {
						if (int64(i)-int64(j)*int64(lp.lambda))%int64(lp.p) == 0 {
							want += lp.logP
						}
					}
					if got[pos] != want {
						t.Fatalf("I=%d (i,j)=(%d,%d): got %d, want %d", I, i, j, got[pos], want)
					}
				} else if got[pos] != 3 {
					t.Fatalf("I=%d both-even (i,j)=(%d,%d) was sieved: got %d, want 3", I, i, j, got[pos])
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

func TestSubSegmentNormInit(t *testing.T) {
	rnd := rand.New(rand.NewSource(41))
	n := new(big.Int).Mul(testutil.RandomPrime(rnd, 130), testutil.RandomPrime(rnd, 130))
	poly, fac := SelectPoly(n, 4)
	if fac != nil || poly == nil {
		t.Skip("lucky factor or nil poly")
	}
	fp := newFastNormPoly(poly)

	for trial := 0; trial < 40; trial++ {
		a := int64(rnd.Intn(200000) - 100000)
		b := int64(1 + rnd.Intn(50000))
		fExact := new(big.Int).Abs(poly.EvalF(a, b))
		if fExact.Sign() > 0 {
			_, gotF := fp.evalLog2F(float64(a), float64(b))
			wantF := arith.Log2Big(fExact)
			if diff := gotF - wantF; diff < -0.02 || diff > 0.02 {
				t.Fatalf("evalLog2F(%d,%d) = %.4f, want %.4f", a, b, gotF, wantF)
			}
		}
		gExact := new(big.Int).Abs(poly.EvalG(a, b))
		if gExact.Sign() > 0 {
			_, gotG := fp.evalLog2G(float64(a), float64(b))
			wantG := arith.Log2Big(gExact)
			if diff := gotG - wantG; diff < -0.02 || diff > 0.02 {
				t.Fatalf("evalLog2G(%d,%d) = %.4f, want %.4f", a, b, gotG, wantG)
			}
		}
	}

	// Verify interior offsets receive positive baseline credit relative to block corners.
	lq := ReduceLattice(100003, 12345, poly.Skew, 2048, 1024, -1)
	var ab [blockSize]byte
	const W = 4096
	aThresh := initLatticeNormsAlg(ab[:], &fp, &lq, 0, blockSize, 2048, W, 16.6, 48.0, 2)
	if aThresh <= 0 {
		t.Fatalf("unexpected non-positive block threshold %d", aThresh)
	}
	center := ab[2048]
	corner := ab[W-1]
	if center <= corner {
		t.Fatalf("expected center baseline (%d) > corner baseline (%d) on j=1", center, corner)
	}
}

func TestOddModInv32AndPowerSieving(t *testing.T) {
	for _, p64 := range arith.PrimesUpTo(5000) {
		if p64 == 2 {
			continue
		}
		p := uint32(p64)
		inv, lim := oddModInv32(p)
		for u := uint32(0); u < 20000; u += 3 {
			got := u*inv <= lim
			want := u%p == 0
			if got != want {
				t.Fatalf("oddModInv32(p=%d) at u=%d: got %v, want %v", p, u, got, want)
			}
		}
	}

	powers := []sievePower{
		{pk: 16, lambda: 3, logP: 1, evenP: true},
		{pk: 16, lambda: 4, logP: 1, evenP: true},
		{pk: 25, lambda: 11, logP: 2, evenP: false},
		{pk: 27, lambda: 0, logP: 2, evenP: false},
	}
	for _, I := range []int{64, 60} {
		W := 2 * I
		const J = 20
		total := (J - 1) * W
		got := make([]byte, total)
		step := 384
		for bStart := 0; bStart < total; bStart += step {
			bEnd := bStart + step
			if bEnd > total {
				bEnd = total
			}
			sieveBlockPowersAdd(got[bStart:bEnd], bStart, bEnd, I, W, powers)
		}
		for j := 1; j < J; j++ {
			for i := -I; i < I; i++ {
				pos := (j-1)*W + (i + I)
				if (i|j)&1 == 0 {
					if got[pos] != 0 {
						t.Fatalf("I=%d both-even (%d,%d) got %d, want 0", I, i, j, got[pos])
					}
					continue
				}
				var want byte
				for _, sp := range powers {
					if (int64(i)-int64(j)*int64(sp.lambda))%int64(sp.pk) == 0 {
						want += sp.logP
					}
				}
				if got[pos] != want {
					t.Fatalf("I=%d (%d,%d): got %d, want %d", I, i, j, got[pos], want)
				}
			}
		}
	}

	lines := []linePrime{
		{idx: 10, p: 29, lambda: 0, logP: 5},
		{idx: 11, p: 31, lambda: 17, logP: 5},
		{idx: 12, p: 37, lambda: 36, logP: 5},
	}
	verts := []vertPrime{{idx: 13, p: 5, logP: 2}}
	for _, I := range []int{64, 60} {
		W := 2 * I
		const J = 16
		total := (J - 1) * W
		var candBits [blockSize / 64]uint64
		var candSlot [blockSize]uint16
		var cands []int
		for pos := 0; pos < total; pos += 13 {
			j := 1 + pos/W
			i := (pos % W) - I
			if (i|j)&1 == 0 {
				continue
			}
			candSlot[pos] = uint16(len(cands))
			candBits[pos>>6] |= 1 << (uint(pos) & 63)
			cands = append(cands, pos)
		}
		hits := make([][]int32, len(cands))
		gatherBlockLineHits(0, total, I, W, lines, verts, &candBits, &candSlot, hits)
		for k, pos := range cands {
			j := 1 + pos/W
			i := (pos % W) - I
			gotSet := make(map[int32]bool)
			for _, id := range hits[k] {
				gotSet[id] = true
			}
			for _, vp := range verts {
				if (j%int(vp.p) == 0) != gotSet[vp.idx] {
					t.Fatalf("I=%d vert (%d,%d) mismatch", I, i, j)
				}
			}
			for _, lp := range lines {
				want := (int64(i)-int64(j)*int64(lp.lambda))%int64(lp.p) == 0
				if want != gotSet[lp.idx] {
					t.Fatalf("I=%d line p=%d (%d,%d) mismatch: got %v, want %v", I, lp.p, i, j, gotSet[lp.idx], want)
				}
			}
		}
	}
}

func TestLatticeMultiStripMatchesSingleStrip(t *testing.T) {
	rnd := rand.New(rand.NewSource(77))
	n := new(big.Int).Mul(testutil.RandomPrime(rnd, 58), testutil.RandomPrime(rnd, 58))
	poly, fac := BaseMPoly(n, 3, 6)
	if fac != nil {
		t.Skip("got lucky factor")
	}
	fb := BuildFactorBases(poly, 1200, 1500, 80000, 64)

	for _, tc := range []struct {
		I, J  int
		strip int
	}{
		{60, 320, 1},
		{256, 256, 2},
	} {
		params := &Params{
			Degree: 3, Br: 1200, Ba: 1500,
			LpBoundRat: 32000, LpBoundAlg: 32000,
			MfbRat: 30, MfbAlg: 30, MaxLPRat: 2, MaxLPAlg: 2,
			UseLattice: true, QMin: 1200, QMax: 1600,
			IHalf: tc.I, JMax: tc.J,
			SlackRat: 28, SlackAlg: 28,
		}
		sqs := buildSpecialQs(poly, fb, params)
		if len(sqs) > 8 {
			sqs = sqs[:8]
		}
		var single, multi []Relation
		sieveLatticeQsStrip(poly, fb, params, sqs, 1, 0, func(r Relation) {
			single = append(single, r)
		}, func() bool { return false })
		sieveLatticeQsStrip(poly, fb, params, sqs, 1, tc.strip, func(r Relation) {
			multi = append(multi, r)
		}, func() bool { return false })

		if len(single) == 0 || len(single) != len(multi) {
			t.Fatalf("I=%d J=%d strip=%d: single=%d relations, multi=%d", tc.I, tc.J, tc.strip, len(single), len(multi))
		}
		for k := range single {
			a, b := single[k], multi[k]
			if a.A != b.A || a.B != b.B || a.RatLPs != b.RatLPs || len(a.RatCols) != len(b.RatCols) || len(a.AlgCols) != len(b.AlgCols) || len(a.AlgLPs) != len(b.AlgLPs) {
				t.Fatalf("I=%d J=%d relation %d mismatch: single=%+v multi=%+v", tc.I, tc.J, k, a, b)
			}
		}

		// Also verify multi worker multi strip produces the identical deduplicated relation set.
		seenMulti4 := make(map[[2]int64]bool)
		sieveLatticeQsStrip(poly, fb, params, sqs, 4, tc.strip, func(r Relation) {
			seenMulti4[[2]int64{r.A, r.B}] = true
		}, func() bool { return false })
		seenSingle := make(map[[2]int64]bool)
		for _, r := range single {
			seenSingle[[2]int64{r.A, r.B}] = true
		}
		if len(seenMulti4) != len(seenSingle) {
			t.Fatalf("I=%d J=%d 4-worker unique=%d, want %d", tc.I, tc.J, len(seenMulti4), len(seenSingle))
		}
	}
}

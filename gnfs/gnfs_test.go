package gnfs

import (
	"bytes"
	"flag"
	"math/big"
	"math/rand"
	"os"
	"strings"
	"testing"
	"time"

	"github.com/dgryski/go-factor/internal/linalg"
	"github.com/dgryski/go-factor/internal/testutil"
)

var relationsFile = flag.String("relations", "", "relations file saved by factor -gnfs-save, for TestSavedRelations")

func TestSavedRelations(t *testing.T) {
	if *relationsFile == "" {
		t.Skip("no -relations file")
	}
	f, err := os.Open(*relationsFile)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()

	start := time.Now()
	n, poly, pairs, isFree, err := readPolyAndPairs(f)
	if err != nil {
		t.Fatal(err)
	}
	params := defaultParams(n)
	qcMin := computeQCMin(params)
	fb := BuildFactorBases(poly, params.Br, params.Ba, qcMin, params.NumQCs)
	rels := refactorPairs(poly, fb, params, qcMin, pairs, isFree)
	t.Logf("read %d/%d relations for %d-bit n in %v", len(rels), len(pairs), n.BitLen(), time.Since(start))

	Logf = t.Logf
	defer func() { Logf = nil }()
	f1, f2 := solveRelations(n, poly, fb, rels)
	if f1 == nil || new(big.Int).Mul(f1, f2).Cmp(n) != 0 {
		t.Fatalf("solveRelations(%v) = %v, %v", n, f1, f2)
	}
	t.Logf("factors: %v * %v", f1, f2)
}

func TestSaveAndResume(t *testing.T) {
	dir := t.TempDir()
	defer func(prev string) { SaveDir = prev }(SaveDir)
	SaveDir = dir

	rnd := rand.New(rand.NewSource(21))
	n := new(big.Int).Mul(testutil.RandomPrime(rnd, 68), testutil.RandomPrime(rnd, 68))

	f1, f2 := split(n)
	if f1 == nil || new(big.Int).Mul(f1, f2).Cmp(n) != 0 {
		t.Fatalf("initial split(%v) failed: %v, %v", n, f1, f2)
	}
	path := relationsPath(n)
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("missing saved relations %s: %v", path, err)
	}

	var logs []string
	Logf = func(format string, args ...interface{}) {
		logs = append(logs, format)
	}
	r1, r2 := split(n)
	Logf = nil
	if r1 == nil || new(big.Int).Mul(r1, r2).Cmp(n) != 0 {
		t.Fatalf("full resume split(%v) failed: %v, %v", n, r1, r2)
	}
	joined := strings.Join(logs, "\n")
	if !strings.Contains(joined, "resuming") || strings.Contains(joined, "sieved") {
		t.Fatalf("expected full resume to skip sieving, logs:\n%s", joined)
	}

	lines := strings.Split(strings.TrimSpace(string(data)), "\n")
	keep := 3 + (len(lines)-3)/2
	if err := os.WriteFile(path, []byte(strings.Join(lines[:keep], "\n")+"\n"), 0644); err != nil {
		t.Fatal(err)
	}
	logs = nil
	Logf = func(format string, args ...interface{}) {
		logs = append(logs, format)
	}
	p1, p2 := split(n)
	Logf = nil
	if p1 == nil || new(big.Int).Mul(p1, p2).Cmp(n) != 0 {
		t.Fatalf("partial resume split(%v) failed: %v, %v", n, p1, p2)
	}
	joined = strings.Join(logs, "\n")
	if !strings.Contains(joined, "resuming") || !strings.Contains(joined, "sieved") {
		t.Fatalf("expected partial resume to continue sieving, logs:\n%s", joined)
	}
}

func TestFactorLineSieve30To42Digits(t *testing.T) {
	Logf = t.Logf
	defer func() { Logf = nil }()
	rnd := rand.New(rand.NewSource(1))
	for _, bits := range []int{50, 60} { // ~30-36 digits
		p1 := testutil.RandomPrime(rnd, bits)
		p2 := testutil.RandomPrime(rnd, bits)
		n := new(big.Int).Mul(p1, p2)
		orig := new(big.Int).Set(n)
		p, c := Factor(n)
		testutil.VerifyFactoring(t, orig, p, c)
		if len(c) != 0 {
			t.Fatalf("Factor(%v) left composites %v", orig, c)
		}
	}
}

func TestFactorLatticeSieve45To55Digits(t *testing.T) {
	Logf = t.Logf
	defer func() { Logf = nil }()
	rnd := rand.New(rand.NewSource(2))
	for _, bits := range []int{72, 85} { // ~43-51 digits
		p1 := testutil.RandomPrime(rnd, bits)
		p2 := testutil.RandomPrime(rnd, bits)
		n := new(big.Int).Mul(p1, p2)
		orig := new(big.Int).Set(n)
		p, c := Factor(n)
		testutil.VerifyFactoring(t, orig, p, c)
		if len(c) != 0 {
			t.Fatalf("Factor(%v) left composites %v", orig, c)
		}
	}
}

func TestFactorDegree4(t *testing.T) {
	rnd := rand.New(rand.NewSource(3))
	p1 := testutil.RandomPrime(rnd, 72)
	p2 := testutil.RandomPrime(rnd, 72)
	n := new(big.Int).Mul(p1, p2)
	params := &Params{
		Degree: 4, Br: 2500, Ba: 4000,
		LpBoundRat: 1 << 18, LpBoundAlg: 1 << 18,
		MfbRat: 36, MfbAlg: 36, MaxLPRat: 2, MaxLPAlg: 2,
		UseLattice: true, QMin: 4000, QMax: 60000, IHalf: 512, JMax: 512,
		SlackRat: 42, SlackAlg: 42, NumQCs: 64,
	}
	f1, f2 := splitWithParams(n, params)
	if f1 == nil || f2 == nil {
		t.Fatalf("degree-4 splitWithParams(%v) failed", n)
	}
	if new(big.Int).Mul(f1, f2).Cmp(n) != 0 {
		t.Fatalf("%v * %v != %v", f1, f2, n)
	}
}

func TestFactorNonMonicG(t *testing.T) {
	// Test non-monic g(x) = A1*x - A0 (A1 != 1) to exercise the A1^(|S|/2) factor in X_alg/X_rat.
	rnd := rand.New(rand.NewSource(7))
	p1 := testutil.RandomPrime(rnd, 45)
	p2 := testutil.RandomPrime(rnd, 45)
	n := new(big.Int).Mul(p1, p2)

	base, fac := BaseMPoly(n, 3, 1)
	if fac != nil {
		t.Skip("got lucky factor")
	}
	// Transform f_new(y) = sum_i c_i * (A1*y - k)^i * A1^(d-i), g(y) = A1*y - (A1*m + k).
	a1 := big.NewInt(2)
	k := big.NewInt(1)
	d := base.Degree
	newCoeffs := make([]*big.Int, d+1)
	for i := range newCoeffs {
		newCoeffs[i] = new(big.Int)
	}
	lin := []*big.Int{new(big.Int).Neg(k), new(big.Int).Set(a1)} // A1*y - k
	for i := 0; i <= d; i++ {
		term := []*big.Int{new(big.Int).Mul(base.Coeffs[i], new(big.Int).Exp(a1, big.NewInt(int64(d-i)), nil))}
		for j := 0; j < i; j++ {
			// Multiply term by (A1*y - k).
			next := make([]*big.Int, len(term)+1)
			for tIdx := range next {
				next[tIdx] = new(big.Int)
			}
			for tIdx, tc := range term {
				next[tIdx].Add(next[tIdx], new(big.Int).Mul(tc, lin[0]))
				next[tIdx+1].Add(next[tIdx+1], new(big.Int).Mul(tc, lin[1]))
			}
			term = next
		}
		for tIdx, tc := range term {
			newCoeffs[tIdx].Add(newCoeffs[tIdx], tc)
		}
	}

	a0 := new(big.Int).Add(new(big.Int).Mul(a1, base.A0), k)
	a1Inv := new(big.Int).ModInverse(a1, n)
	mNew := new(big.Int).Mul(a0, a1Inv)
	mNew.Mod(mNew, n)

	poly := &Poly{
		Degree: d,
		Coeffs: newCoeffs,
		A1:     a1,
		A0:     a0,
		M:      mNew,
	}
	poly.Skew = OptimizeSkewness(poly)
	poly.Alpha = ComputeAlpha(poly, 200)

	fb := BuildFactorBases(poly, 800, 1200, 200000, 64)
	params := &Params{
		Degree: 3, Br: 800, Ba: 1200,
		LpBoundRat: 32000, LpBoundAlg: 32000,
		MfbRat: 30, MfbAlg: 30, MaxLPRat: 2, MaxLPAlg: 2,
		QMin: 1200, QMax: 30000, IHalf: 256, JMax: 256,
		SlackRat: 36, SlackAlg: 36, NumQCs: 64,
	}
	var rels []Relation
	seen := make(map[[2]int64]struct{})
	sieveLattice(poly, fb, params, func(r Relation) {
		key := [2]int64{r.A, r.B}
		if _, ok := seen[key]; !ok {
			seen[key] = struct{}{}
			rels = append(rels, r)
		}
	}, func() bool {
		if len(rels) < fb.NumFixed+40 {
			return false
		}
		ex, _, _ := pruneSingletons(rels, fb)
		return ex >= 20
	})
	surviving, rows, numCols := buildMatrix(poly, fb, rels)
	if len(rows) <= numCols {
		t.Fatalf("not enough rows: %d <= %d", len(rows), numCols)
	}
	// Run linear algebra and sqrt.
	var gotFac *big.Int
	for _, dep := range linalg.FindDependencies(rows, numCols, nil) {
		if fac, ok := tryDependency(n, poly, fb, surviving, dep); ok && fac != nil {
			gotFac = fac
			break
		}
	}
	if gotFac == nil || (gotFac.Cmp(p1) != 0 && gotFac.Cmp(p2) != 0) {
		t.Fatalf("non-monic g failed: got %v, want %v or %v", gotFac, p1, p2)
	}
}

func TestRelationsRoundTrip(t *testing.T) {
	n := big.NewInt(1000036000099)
	poly, _ := BaseMPoly(n, 3, 1)
	rels := []Relation{{A: -17, B: 3}, {A: 42, B: 1}, {A: 13, B: 0, IsFree: true}}
	var buf bytes.Buffer
	if err := writeRelations(&buf, n, poly, rels); err != nil {
		t.Fatal(err)
	}
	n2, poly2, pairs, isFree, err := readPolyAndPairs(&buf)
	if err != nil {
		t.Fatal(err)
	}
	if n2.Cmp(n) != 0 || poly2.Degree != poly.Degree || len(pairs) != 3 || !isFree[2] {
		t.Fatalf("round trip mismatch: n=%v pairs=%v isFree=%v", n2, pairs, isFree)
	}
}

func TestFactor59Digits(t *testing.T) {
	// 195 bits, the top of the 150-200 bit parameter row. The old factor
	// bases (Br=5000, Ba=8000) found no cycles here at all.
	n, _ := new(big.Int).SetString("31282431280413161904649583977615847079904312904104617466157", 10)
	f1, f2 := split(n)
	if f1 == nil || new(big.Int).Mul(f1, f2).Cmp(n) != 0 {
		t.Fatalf("split(%v) = %v, %v", n, f1, f2)
	}
}

func TestFreeRelationsAndLPFilter(t *testing.T) {
	rnd := rand.New(rand.NewSource(37))
	n := new(big.Int).Mul(testutil.RandomPrime(rnd, 55), testutil.RandomPrime(rnd, 55))
	poly, fac := BaseMPoly(n, 3, 6)
	if fac != nil {
		t.Skip("got lucky factor")
	}
	params := &Params{
		Degree: 3, Br: 2000, Ba: 2000,
		LpBoundRat: 32000, LpBoundAlg: 32000,
		MfbRat: 30, MfbAlg: 30, MaxLPRat: 2, MaxLPAlg: 2,
		NumQCs: 64,
	}
	const qcMin = 100000
	fb := BuildFactorBases(poly, params.Br, params.Ba, qcMin, params.NumQCs)
	free := freeRelations(poly, fb)
	if len(free) == 0 {
		t.Fatalf("expected free relations for Br=Ba=2000")
	}

	sc := newCofactorScratch(poly.Degree)
	for _, rel := range free {
		if !rel.IsFree || rel.B != 0 || len(rel.RatCols) != 1 || len(rel.AlgCols) != poly.Degree {
			t.Fatalf("malformed free relation %+v", rel)
		}
		got, ok := cofactorPair(poly, fb, params, qcMin, rel.A, rel.B, true, sc)
		if !ok || len(got.AlgCols) != poly.Degree {
			t.Fatalf("cofactorPair rejected free relation p=%d", rel.A)
		}
		dup := []Relation{rel, rel}
		surv, rows, numCols := buildMatrix(poly, fb, dup)
		if len(rows) != 2 || len(rows[0]) == 0 {
			t.Fatalf("unexpected rows for free relation p=%d", rel.A)
		}
		deps := linalg.FindDependencies(rows, numCols, nil)
		if len(deps) != 1 {
			t.Fatalf("expected 1 dependency for duplicate free relation, got %d", len(deps))
		}
		if _, sqOk := tryDependency(n, poly, fb, surv, deps[0]); !sqOk {
			t.Fatalf("duplicate free relation p=%d failed square check", rel.A)
		}
	}

	// Verify incremental lpFilter recounts match pruneSingletons across steps.
	synth := []Relation{
		{A: 3, B: 1},
		{A: 5, B: 1, RatLPs: [3]uint64{20011}},
		{A: 7, B: 1, RatLPs: [3]uint64{20011, 20021}},
		{A: 9, B: 1, RatLPs: [3]uint64{20021}, AlgLPs: []LargeIdeal{{P: 20029, R: 4}}},
		{A: 11, B: 1, AlgLPs: []LargeIdeal{{P: 20029, R: 4}}},
		{A: 13, B: 1, AlgLPs: []LargeIdeal{{P: 20029, R: 5}}},
	}
	filt := newLPFilter(len(synth))
	for i := range synth {
		filt.add(&synth[i])
		got := filt.excess(fb.NumFixed)
		want, _, _ := pruneSingletons(synth[:i+1], fb)
		if got != want {
			t.Fatalf("step %d: filt.excess=%d, want %d", i, got, want)
		}
	}
}

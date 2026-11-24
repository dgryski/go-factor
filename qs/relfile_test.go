package qs

import (
	"bytes"
	"flag"
	"math/big"
	"math/rand"
	"os"
	"reflect"
	"testing"
	"time"

	"github.com/dgryski/go-factor/internal/testutil"
)

var relationsFile = flag.String("relations", "", "relations file saved by factor -qs-save, for TestSavedRelations")

// TestSavedRelations runs the matrix and square root steps on a relations
// file, so they can be tuned without sieving again.
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
	n, factorBase, relations, err := readRelations(f)
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("read %d relations over %d primes in %v", len(relations), len(factorBase), time.Since(start))

	Logf = t.Logf
	defer func() { Logf = nil }()
	start = time.Now()
	deps := findDependencies(relations, len(factorBase)+1)
	t.Logf("%d dependencies in %v", len(deps), time.Since(start))
	for i, dep := range deps {
		if f := tryDependency(n, factorBase, relations, dep); f != nil {
			t.Logf("dependency %d gave factor %v", i+1, f)
			return
		}
	}
	t.Errorf("no factor of %v found", n)
}

func TestRelationsRoundTrip(t *testing.T) {
	rnd := rand.New(rand.NewSource(3))
	n := new(big.Int).Mul(testutil.RandomPrime(rnd, 35), testutil.RandomPrime(rnd, 35))
	b, siqs, sieveRange := params(n)
	kn := new(big.Int).Mul(n, big.NewInt(knuthSchroeppel(n)))
	factorBase := buildFactorBase(kn, primesUpTo(b))
	s := newSieve(n, kn, factorBase, siqs)
	s.siqs(sieveRange)
	factorBase = s.finish(factorBase)

	var buf bytes.Buffer
	if err := writeRelations(&buf, n, factorBase, s.relations); err != nil {
		t.Fatal(err)
	}
	n2, fb2, rels2, err := readRelations(&buf)
	if err != nil {
		t.Fatal(err)
	}
	if n2.Cmp(n) != 0 || !reflect.DeepEqual(fb2, factorBase) || len(rels2) != len(s.relations) {
		t.Fatalf("read back n=%v, %d primes, %d relations; want %v, %d, %d",
			n2, len(fb2), len(rels2), n, len(factorBase), len(s.relations))
	}
	for i, rel := range rels2 {
		if rel.X.Cmp(s.relations[i].X) != 0 || !reflect.DeepEqual(rel.Factors, s.relations[i].Factors) {
			t.Fatalf("relation %d differs", i)
		}
	}
}

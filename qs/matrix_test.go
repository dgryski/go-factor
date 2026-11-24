package qs

import (
	"fmt"
	"math/bits"
	"math/rand"
	"testing"
)

// randomRelations returns relations shaped like the sieve's: columns for
// small primes are dense, those for large primes sparse.
func randomRelations(rng *rand.Rand, numRows, numCols int) []relation {
	rels := make([]relation, numRows)
	for i := range rels {
		for c := 0; c < numCols; c++ {
			// Column c is hit with probability about 1/(c+2), like a prime.
			if rng.Intn(c+2) == 0 {
				for e := 1 + rng.Intn(3); e > 0; e-- {
					rels[i].Factors = append(rels[i].Factors, int32(c))
				}
			}
		}
	}
	return rels
}

func rowsOf(rels []relation) [][]int32 {
	rows := make([][]int32, len(rels))
	for i, rel := range rels {
		rows[i] = oddColumns(rel.Factors)
	}
	return rows
}

func TestFindDependencies(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	for _, size := range []struct{ rows, cols int }{{10, 5}, {60, 50}, {300, 280}, {2020, 2000}} {
		rels := randomRelations(rng, size.rows, size.cols)
		rows := rowsOf(rels)
		deps := findDependencies(rels, size.cols)
		want := len(denseDependencies(rows, size.cols))
		if len(deps) != want {
			t.Errorf("%d×%d: %d dependencies, want %d", size.rows, size.cols, len(deps), want)
		}

		// Each dependency must be non-empty and sum to zero.
		depRows := make([][]int32, len(deps))
		for d, dep := range deps {
			sum := make([]bool, size.cols)
			for w, word := range dep {
				for word != 0 {
					i := w*64 + bits.TrailingZeros64(word)
					word &= word - 1
					depRows[d] = append(depRows[d], int32(i))
					for _, c := range rows[i] {
						sum[c] = !sum[c]
					}
				}
			}
			if len(depRows[d]) == 0 {
				t.Fatalf("%d×%d: dependency %d is empty", size.rows, size.cols, d)
			}
			for c, odd := range sum {
				if odd {
					t.Fatalf("%d×%d: dependency %d leaves column %d odd", size.rows, size.cols, d, c)
				}
			}
		}

		// And they must be independent.
		if n := len(denseDependencies(depRows, size.rows)); n != 0 {
			t.Errorf("%d×%d: %d of the dependencies are combinations of the others", size.rows, size.cols, n)
		}
	}
}

func BenchmarkFindDependencies(b *testing.B) {
	for _, numCols := range []int{2000, 10000} {
		rng := rand.New(rand.NewSource(42))
		rels := randomRelations(rng, numCols+extraRelations, numCols)
		b.Run(fmt.Sprintf("cols=%d", numCols), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				findDependencies(rels, numCols)
			}
		})
	}
}

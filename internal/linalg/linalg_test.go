package linalg

import (
	"fmt"
	"math/bits"
	"math/rand"
	"sort"
	"testing"
)

func randomSparseRows(rng *rand.Rand, numRows, numCols, weight int) [][]int32 {
	rows := make([][]int32, numRows)
	for r := range rows {
		seen := make(map[int32]bool)
		for len(rows[r]) < weight {
			c := int32(rng.Intn(numCols))
			if !seen[c] {
				seen[c] = true
				rows[r] = append(rows[r], c)
			}
		}
		sort.Slice(rows[r], func(i, j int) bool { return rows[r][i] < rows[r][j] })
	}
	return rows
}

func TestFindDependencies(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	for _, size := range []struct{ rows, cols int }{{10, 5}, {60, 50}, {300, 280}, {2020, 2000}, {6064, 6000}} {
		rows := make([][]int32, size.rows)
		for i := range rows {
			var f []int32
			for c := 0; c < size.cols; c++ {
				if rng.Intn(c+2) == 0 {
					for e := 1 + rng.Intn(3); e > 0; e-- {
						f = append(f, int32(c))
					}
				}
			}
			rows[i] = OddColumns(f)
		}
		origRows := make([][]int32, len(rows))
		for i, r := range rows {
			origRows[i] = append([]int32(nil), r...)
		}
		deps := FindDependencies(rows, size.cols, nil)
		want := len(DenseDependencies(origRows, size.cols))
		if want > KeepSurplus {
			want = KeepSurplus
		}
		if len(deps) < want {
			t.Errorf("%d×%d: %d dependencies, want at least %d", size.rows, size.cols, len(deps), want)
		}
		depRows := make([][]int32, len(deps))
		for d, dep := range deps {
			sum := make([]bool, size.cols)
			for w, word := range dep {
				for word != 0 {
					i := w*64 + bits.TrailingZeros64(word)
					word &= word - 1
					depRows[d] = append(depRows[d], int32(i))
					for _, c := range origRows[i] {
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
		if n := len(DenseDependencies(depRows, size.rows)); n != 0 {
			t.Errorf("%d×%d: %d of the dependencies are combinations of the others", size.rows, size.cols, n)
		}
	}
}

func TestBlockLanczosParallel(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	for _, tc := range []struct{ rows, cols, weight int }{{5064, 5000, 24}, {10064, 10000, 24}} {
		rows := randomSparseRows(rng, tc.rows, tc.cols, tc.weight)
		deps := BlockLanczos(rows, tc.cols, 1)
		for i, dep := range deps {
			if !IsDependency(rows, tc.cols, dep) {
				t.Fatalf("%d×%d: dependency %d doesn't sum to zero", tc.rows, tc.cols, i)
			}
		}
		indep := IndependentCount(deps, tc.rows)
		if indep < LanczosMinDeps {
			t.Errorf("%d×%d: %d independent dependencies, want at least %d", tc.rows, tc.cols, indep, LanczosMinDeps)
		}
	}
}

func BenchmarkBlockLanczos(b *testing.B) {
	for _, n := range []int{10000, 50000, 150000} {
		rows := randomSparseRows(rand.New(rand.NewSource(99)), n+64, n, 40)
		b.Run(fmt.Sprintf("cols=%d", n), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				BlockLanczos(rows, n, 1)
			}
		})
	}
}

package qs

import (
	"fmt"
	"math/rand"
	"sort"
	"testing"
)

// randomSparseRows returns rows with weight distinct random columns each.
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

func TestFindWinv2x2Pivots(t *testing.T) {
	// Every diagonal entry is 0, which forces 32 2×2 pivots.
	var sym mat64
	for i := 0; i < 64; i += 2 {
		sym[i] |= 1 << uint(i+1)
		sym[i+1] |= 1 << uint(i)
	}
	wInv, sMask := findWinv(sym, 0)
	if sMask != ^uint64(0) {
		t.Fatalf("mask %016x, want full rank", sMask)
	}
	prod := mulMat64(wInv, sym)
	for i := range prod {
		if prod[i] != 1<<uint(i) {
			t.Fatalf("row %d of wInv·t is %016x, want %016x", i, prod[i], uint64(1)<<uint(i))
		}
	}
}

func TestBlockLanczos(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	for _, tc := range []struct{ rows, cols, weight int }{{90, 70, 8}, {320, 300, 14}, {2030, 2000, 20}, {10064, 10000, 24}} {
		rows := randomSparseRows(rng, tc.rows, tc.cols, tc.weight)
		deps := blockLanczos(rows, tc.cols, 1)
		for i, dep := range deps {
			if !isDependency(rows, tc.cols, dep) {
				t.Fatalf("%d×%d: dependency %d doesn't sum to zero", tc.rows, tc.cols, i)
			}
		}
		indep := independentCount(deps, tc.rows)
		if indep < lanczosMinDeps {
			t.Errorf("%d×%d: %d independent dependencies, want at least %d", tc.rows, tc.cols, indep, lanczosMinDeps)
		}
		t.Logf("%d×%d: %d dependencies, %d independent", tc.rows, tc.cols, len(deps), indep)
	}
}

func BenchmarkBlockLanczos(b *testing.B) {
	for _, n := range []int{10000, 50000} {
		rows := randomSparseRows(rand.New(rand.NewSource(99)), n+64, n, 24)
		b.Run(fmt.Sprintf("cols=%d", n), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				blockLanczos(rows, n, 1)
			}
		})
	}
}

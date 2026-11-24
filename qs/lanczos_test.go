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

// rankGap returns rank(B) and rank(B·Bᵀ).
func rankGap(rows [][]int32, numCols int) (rankB, rankA int) {
	colRows := make([][]int32, numCols)
	for r, row := range rows {
		for _, c := range row {
			colRows[c] = append(colRows[c], int32(r))
		}
	}
	aRows := make([][]int32, len(rows))
	for r, row := range rows {
		odd := make(map[int32]bool)
		for _, c := range row {
			for _, s := range colRows[c] {
				odd[s] = !odd[s]
			}
		}
		for s, v := range odd {
			if v {
				aRows[r] = append(aRows[r], s)
			}
		}
	}
	return len(rows) - len(denseDependencies(rows, numCols)), len(rows) - len(denseDependencies(aRows, len(rows)))
}

// TestBlockLanczosLargeGap uses matrices where rank(B) − rank(B·Bᵀ) is over
// 64. Repeating a block of columns makes those columns cancel in B·Bᵀ.
func TestBlockLanczosLargeGap(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	dup := func(m [][]int32, n int, rest [][]int32) [][]int32 {
		rows := make([][]int32, len(m))
		for r := range m {
			rows[r] = append(rows[r], m[r]...)
			for _, c := range m[r] {
				rows[r] = append(rows[r], c+int32(n))
			}
			if rest != nil {
				for _, c := range rest[r] {
					rows[r] = append(rows[r], c+int32(2*n))
				}
			}
		}
		return rows
	}
	for _, tc := range []struct {
		name string
		rows [][]int32
		cols int
	}{
		{"partial", dup(randomSparseRows(rng, 360, 90, 12), 90, randomSparseRows(rng, 360, 120, 14)), 300},
		{"B·Bᵀ = 0", dup(randomSparseRows(rng, 300, 100, 14), 100, nil), 200},
	} {
		rankB, rankA := rankGap(tc.rows, tc.cols)
		if rankB-rankA < 64 {
			t.Fatalf("%s: gap %d, want at least 64", tc.name, rankB-rankA)
		}
		deps := blockLanczos(tc.rows, tc.cols, 1)
		for i, dep := range deps {
			if !isDependency(tc.rows, tc.cols, dep) {
				t.Fatalf("%s: dependency %d doesn't sum to zero", tc.name, i)
			}
		}
		indep := independentCount(deps, len(tc.rows))
		t.Logf("%s: gap %d, %d independent dependencies", tc.name, rankB-rankA, indep)
		if indep < 32 {
			t.Errorf("%s: %d independent dependencies, want at least 32", tc.name, indep)
		}
	}
}

// TestFMatchesMontgomery checks the F term against Montgomery's formula at
// every iteration.
func TestFMatchesMontgomery(t *testing.T) {
	rows := randomSparseRows(rand.New(rand.NewSource(555)), 450, 400, 18)
	checks := 0
	deps := blockLanczosCheck(rows, 400, 999, func(iter int, f, fMont mat64) {
		checks++
		if f != fMont {
			t.Fatalf("iteration %d: F differs from Montgomery's", iter)
		}
	})
	if checks == 0 || len(deps) == 0 {
		t.Fatalf("%d F checks and %d dependencies, want some of each", checks, len(deps))
	}
}

func TestPreconditionerSymmetric(t *testing.T) {
	const c = 257
	rng := rand.New(rand.NewSource(7))
	p := newSymPrecond(c, rng)
	x, y := make([]uint64, c), make([]uint64, c)
	for i := range x {
		x[i], y[i] = rng.Uint64(), rng.Uint64()
	}
	px, py := append([]uint64(nil), x...), append([]uint64(nil), y...)
	p.apply(px)
	p.apply(py)
	if innerProd64(x, py) != innerProd64(px, y) {
		t.Fatal("xᵀ·P·y ≠ (P·x)ᵀ·y")
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

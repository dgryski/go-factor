// Package linalg implements GF(2) structured Gaussian elimination,
// dense elimination, and Montgomery's parallel Block Lanczos algorithm.
package linalg

import (
	"math/bits"
	"sort"
	"time"
)

const (
	// KeepSurplus is the number of rows beyond the column count kept for dependencies.
	KeepSurplus = 64

	// LanczosMinCols is the minimum reduced column count for Block Lanczos.
	LanczosMinCols = 1000
	// LanczosMinDeps is the minimum independent dependencies required from a Lanczos run.
	LanczosMinDeps = 16
	// LanczosSeeds is the maximum number of Lanczos seeds tried before falling back to dense elimination.
	LanczosSeeds = 3
)

// pivotLimits returns the heaviest column the structured elimination pivots
// on, and the most entries a pivot may add.
func pivotLimits(numCols int) (weight, fill int) {
	weight, fill = numCols/1000, numCols/16
	if weight < 6 {
		weight = 6
	} else if weight > 50 {
		weight = 50
	}
	if fill < 64 {
		fill = 64
	} else if fill > 4000 {
		fill = 4000
	}
	return weight, fill
}

// FindDependencies returns subsets of rows whose exponent vectors sum to zero
// mod 2, as bitsets over row indices. Each input row must contain sorted
// distinct column indices (as produced by OddColumns).
func FindDependencies(rows [][]int32, numCols int, logf func(string, ...interface{})) [][]uint64 {
	start := time.Now()
	r := NewReducer(rows, numCols)
	r.Reduce()
	dense, denseCols, nodes := r.Compact()
	if logf != nil {
		logf("reduced matrix to %d×%d in %v", len(dense), denseCols, time.Since(start).Round(time.Millisecond))
	}

	var deps [][]uint64
	for _, node := range r.found {
		deps = append(deps, r.Expand([]int32{node}))
	}
	var roots []int32
	// A failed block Lanczos run costs much less than dense elimination, so
	// try a few seeds first.
	for seed := int64(1); denseCols >= LanczosMinCols && seed <= LanczosSeeds; seed++ {
		start := time.Now()
		bl := BlockLanczos(dense, denseCols, seed)
		indep := IndependentCount(bl, len(dense))
		if logf != nil {
			logf("block Lanczos (seed %d) found %d independent dependencies in %v",
				seed, indep, time.Since(start).Round(time.Millisecond))
		}
		if indep >= LanczosMinDeps {
			for _, dep := range bl {
				roots = roots[:0]
				for _, i := range dep {
					roots = append(roots, nodes[i])
				}
				deps = append(deps, r.Expand(roots))
			}
			return deps
		}
	}
	for _, combo := range DenseDependencies(dense, denseCols) {
		roots = roots[:0]
		for w, word := range combo {
			for word != 0 {
				roots = append(roots, nodes[w*64+bits.TrailingZeros64(word)])
				word &= word - 1
			}
		}
		deps = append(deps, r.Expand(roots))
	}
	return deps
}

// IndependentCount returns how many of deps, subsets of numRows rows, are
// linearly independent.
func IndependentCount(deps [][]int, numRows int) int {
	rows := make([][]int32, len(deps))
	for i, dep := range deps {
		rows[i] = make([]int32, len(dep))
		for j, r := range dep {
			rows[i][j] = int32(r)
		}
	}
	return len(deps) - len(DenseDependencies(rows, numRows))
}

// OddColumns returns the sorted columns that appear an odd number of times
// in factors.
func OddColumns(factors []int32) []int32 {
	sorted := append([]int32(nil), factors...)
	sort.Slice(sorted, func(i, j int) bool { return sorted[i] < sorted[j] })
	var cols []int32
	for i := 0; i < len(sorted); {
		j := i
		for j < len(sorted) && sorted[j] == sorted[i] {
			j++
		}
		if (j-i)&1 != 0 {
			cols = append(cols, sorted[i])
		}
		i = j
	}
	return cols
}

// Reducer performs structured Gaussian elimination on the rows of a GF(2)
// matrix. Every row is a node: the original rows are nodes 0..numOrig−1, and
// each XOR of one row into another adds a node, so that dependencies of the
// reduced matrix can be traced back to the original rows.
type Reducer struct {
	rows    [][]int32  // sorted columns of each row
	alive   []bool     // false once a row is deleted or becomes empty
	node    []int32    // current node of each row
	merges  [][2]int32 // node numOrig+i is the XOR of nodes merges[i]
	numOrig int

	maxWeight, maxFill int // pivot limits from pivotLimits

	colRows [][]int32 // rows that have contained each column, possibly stale
	weight  []int     // number of live rows containing each column
	queue   []int32   // columns whose weight dropped to 1 or 2
	found   []int32   // nodes of rows that became empty, each a dependency
}

// NewReducer creates a structured Gaussian elimination reducer for rows over numCols columns.
func NewReducer(rows [][]int32, numCols int) *Reducer {
	r := &Reducer{
		rows:    rows,
		alive:   make([]bool, len(rows)),
		node:    make([]int32, len(rows)),
		numOrig: len(rows),
		colRows: make([][]int32, numCols),
		weight:  make([]int, numCols),
	}
	r.maxWeight, r.maxFill = pivotLimits(numCols)
	for i, row := range rows {
		r.node[i] = int32(i)
		r.alive[i] = len(row) > 0
		if len(row) == 0 {
			r.found = append(r.found, int32(i))
		}
		for _, c := range row {
			r.colRows[c] = append(r.colRows[c], int32(i))
			r.weight[c]++
		}
	}
	for c := range r.weight {
		r.enqueue(int32(c))
	}
	return r
}

func (r *Reducer) enqueue(c int32) {
	if w := r.weight[c]; w == 1 || w == 2 {
		r.queue = append(r.queue, c)
	}
}

// liveRows returns the live rows containing column c, dropping stale entries.
func (r *Reducer) liveRows(c int32) []int32 {
	list := r.colRows[c]
	n := 0
	for _, row := range list {
		if !r.alive[row] || !hasColumn(r.rows[row], c) {
			continue
		}
		dup := false
		for _, kept := range list[:n] {
			dup = dup || kept == row
		}
		if !dup {
			list[n] = row
			n++
		}
	}
	r.colRows[c] = list[:n]
	return list[:n]
}

func hasColumn(cols []int32, c int32) bool {
	i := sort.Search(len(cols), func(i int) bool { return cols[i] >= c })
	return i < len(cols) && cols[i] == c
}

// deleteRow removes a row that can't be part of any dependency.
func (r *Reducer) deleteRow(row int32) {
	r.alive[row] = false
	for _, c := range r.rows[row] {
		r.weight[c]--
		r.enqueue(c)
	}
}

// xorInto replaces row dst with dst XOR src.
func (r *Reducer) xorInto(dst, src int32) {
	a, b := r.rows[dst], r.rows[src]
	out := make([]int32, 0, len(a)+len(b))
	i, j := 0, 0
	for i < len(a) || j < len(b) {
		switch {
		case j == len(b) || (i < len(a) && a[i] < b[j]):
			out = append(out, a[i])
			i++
		case i == len(a) || b[j] < a[i]:
			c := b[j]
			out = append(out, c)
			r.weight[c]++
			r.colRows[c] = append(r.colRows[c], dst)
			j++
		default:
			c := a[i]
			r.weight[c]--
			r.enqueue(c)
			i++
			j++
		}
	}
	r.rows[dst] = out
	r.merges = append(r.merges, [2]int32{r.node[dst], r.node[src]})
	r.node[dst] = int32(r.numOrig + len(r.merges) - 1)
	if len(out) == 0 {
		r.alive[dst] = false
		r.found = append(r.found, r.node[dst])
	}
}

// Reduce removes rows with a column no other row has, merges the two rows
// of each column of weight 2, pivots on light columns, and deletes surplus
// rows.
func (r *Reducer) Reduce() {
	for {
		r.pivot()
		if !r.pruneSurplus() {
			return
		}
	}
}

// pruneSurplus deletes the heaviest rows while there are more than
// KeepSurplus rows beyond the column count, half the surplus at a time.
func (r *Reducer) pruneSurplus() bool {
	var live []int32
	for i, alive := range r.alive {
		if alive {
			live = append(live, int32(i))
		}
	}
	cols := 0
	for _, w := range r.weight {
		if w > 0 {
			cols++
		}
	}
	surplus := len(live) - cols - KeepSurplus
	if surplus <= 0 {
		return false
	}
	sort.Slice(live, func(i, j int) bool { return len(r.rows[live[i]]) > len(r.rows[live[j]]) })
	for _, row := range live[:(surplus+1)/2] {
		r.deleteRow(row)
	}
	return true
}

func (r *Reducer) pivot() {
	type pivCand struct {
		col  int32
		cost int
	}
	var cands []pivCand
	var others []int32
	for {
		r.drain()
		cands = cands[:0]
		for c, w := range r.weight {
			if w < 3 || w > r.maxWeight {
				continue
			}
			rows := r.liveRows(int32(c))
			if len(rows) < 3 {
				continue
			}
			minLen := len(r.rows[rows[0]])
			for _, rw := range rows[1:] {
				if l := len(r.rows[rw]); l < minLen {
					minLen = l
				}
			}
			cost := (len(rows) - 1) * (minLen - 1)
			if cost <= r.maxFill {
				cands = append(cands, pivCand{col: int32(c), cost: cost})
			}
		}
		if len(cands) == 0 {
			return
		}
		// Pivot in ascending Markowitz fill order (Pomerance and Smith, Exper. Math. 1, 1992).
		sort.Slice(cands, func(i, j int) bool {
			if cands[i].cost != cands[j].cost {
				return cands[i].cost < cands[j].cost
			}
			return cands[i].col < cands[j].col
		})

		progress := false
		for _, cand := range cands {
			w := r.weight[cand.col]
			if w < 3 || w > r.maxWeight {
				continue
			}
			rows := r.liveRows(cand.col)
			if len(rows) < 3 {
				continue
			}
			piv := 0
			for i := range rows {
				if len(r.rows[rows[i]]) < len(r.rows[rows[piv]]) {
					piv = i
				}
			}
			cost := (len(rows) - 1) * (len(r.rows[rows[piv]]) - 1)
			if cost > r.maxFill || cost > cand.cost*2+16 {
				continue
			}
			p := rows[piv]
			others = append(others[:0], rows[:piv]...)
			others = append(others, rows[piv+1:]...)
			for _, row := range others {
				r.xorInto(row, p)
			}
			r.deleteRow(p)
			r.drain()
			progress = true
		}
		if !progress {
			return
		}
	}
}

func (r *Reducer) drain() {
	for len(r.queue) > 0 {
		c := r.queue[len(r.queue)-1]
		r.queue = r.queue[:len(r.queue)-1]
		switch r.weight[c] {
		case 1:
			r.deleteRow(r.liveRows(c)[0])
		case 2:
			rows := r.liveRows(c)
			light, heavy := rows[0], rows[1]
			if len(r.rows[light]) > len(r.rows[heavy]) {
				light, heavy = heavy, light
			}
			r.xorInto(heavy, light)
			r.deleteRow(light)
		}
	}
}

// Compact returns the live rows with their columns renumbered from 0, the
// number of columns left, and the node of each row.
func (r *Reducer) Compact() (rows [][]int32, numCols int, nodes []int32) {
	colMap := make([]int32, len(r.weight))
	for c, w := range r.weight {
		colMap[c] = int32(numCols)
		if w > 0 {
			numCols++
		}
	}
	for i, row := range r.rows {
		if !r.alive[i] {
			continue
		}
		mapped := make([]int32, len(row))
		for j, c := range row {
			mapped[j] = colMap[c]
		}
		rows = append(rows, mapped)
		nodes = append(nodes, r.node[i])
	}
	return rows, numCols, nodes
}

// Expand returns the original rows whose XOR is the XOR of the given nodes,
// as a bitset.
func (r *Reducer) Expand(roots []int32) []uint64 {
	used := make([]bool, r.numOrig+len(r.merges))
	for _, n := range roots {
		used[n] = !used[n]
	}
	for n := len(used) - 1; n >= r.numOrig; n-- {
		if used[n] {
			m := r.merges[n-r.numOrig]
			used[m[0]] = !used[m[0]]
			used[m[1]] = !used[m[1]]
		}
	}
	dep := make([]uint64, (r.numOrig+63)/64)
	for i := 0; i < r.numOrig; i++ {
		if used[i] {
			dep[i/64] |= 1 << (uint(i) % 64)
		}
	}
	return dep
}

// DenseDependencies returns subsets of rows that sum to zero mod 2, as
// bitsets over row indices. Each row lists its columns, each at most once.
func DenseDependencies(rows [][]int32, numCols int) [][]uint64 {
	numRows := len(rows)
	colWords := (numCols + 63) / 64
	rowWords := (numRows + 63) / 64

	vecBuf := make([]uint64, numRows*colWords)
	comboBuf := make([]uint64, numRows*rowWords)
	vec := make([][]uint64, numRows)
	combo := make([][]uint64, numRows)
	for i, row := range rows {
		vec[i] = vecBuf[i*colWords : (i+1)*colWords]
		for _, c := range row {
			vec[i][c/64] |= 1 << (uint(c) % 64)
		}
		combo[i] = comboBuf[i*rowWords : (i+1)*rowWords]
		combo[i][i/64] |= 1 << (uint(i) % 64)
	}

	pivots := 0
	for c := 0; c < numCols && pivots < numRows; c++ {
		w, bit := c/64, uint64(1)<<(uint(c)%64)
		p := pivots
		for p < numRows && vec[p][w]&bit == 0 {
			p++
		}
		if p == numRows {
			continue
		}
		vec[pivots], vec[p] = vec[p], vec[pivots]
		combo[pivots], combo[p] = combo[p], combo[pivots]

		pv, pc := vec[pivots], combo[pivots]
		for i := pivots + 1; i < numRows; i++ {
			if vec[i][w]&bit == 0 {
				continue
			}
			for k := w; k < colWords; k++ {
				vec[i][k] ^= pv[k]
			}
			for k := range pc {
				combo[i][k] ^= pc[k]
			}
		}
		pivots++
	}

	return combo[pivots:]
}

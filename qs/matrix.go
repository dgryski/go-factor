package qs

import (
	"math/bits"
	"sort"
	"time"
)

const (
	maxPivotWeight = 6  // the structured elimination pivots on columns with at most this many rows
	maxPivotFill   = 64 // and skips a pivot that would add more entries than this
)

// findDependencies returns subsets of relations whose exponent vectors sum
// to zero mod 2, as bitsets over relation indices. Structured Gaussian
// elimination shrinks the matrix and dense elimination solves what is left.
func findDependencies(relations []relation, numCols int) [][]uint64 {
	start := time.Now()
	rows := make([][]int32, len(relations))
	for i, rel := range relations {
		rows[i] = oddColumns(rel.Factors)
	}
	r := newReducer(rows, numCols)
	r.reduce()
	dense, denseCols, nodes := r.compact()
	logf("qs: reduced matrix to %d×%d in %v", len(dense), denseCols, time.Since(start).Round(time.Millisecond))

	var deps [][]uint64
	for _, node := range r.found {
		deps = append(deps, r.expand([]int32{node}))
	}
	var roots []int32
	for _, combo := range denseDependencies(dense, denseCols) {
		roots = roots[:0]
		for w, word := range combo {
			for word != 0 {
				roots = append(roots, nodes[w*64+bits.TrailingZeros64(word)])
				word &= word - 1
			}
		}
		deps = append(deps, r.expand(roots))
	}
	return deps
}

// oddColumns returns the sorted columns that appear an odd number of times
// in factors.
func oddColumns(factors []int32) []int32 {
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

// reducer performs structured Gaussian elimination on the rows of a GF(2)
// matrix. Every row is a node: the original rows are nodes 0..numOrig−1, and
// each XOR of one row into another adds a node, so that dependencies of the
// reduced matrix can be traced back to the original rows.
type reducer struct {
	rows    [][]int32  // sorted columns of each row
	alive   []bool     // false once a row is deleted or becomes empty
	node    []int32    // current node of each row
	merges  [][2]int32 // node numOrig+i is the XOR of nodes merges[i]
	numOrig int

	colRows [][]int32 // rows that have contained each column, possibly stale
	weight  []int     // number of live rows containing each column
	queue   []int32   // columns whose weight dropped to 1 or 2
	found   []int32   // nodes of rows that became empty, each a dependency
}

func newReducer(rows [][]int32, numCols int) *reducer {
	r := &reducer{
		rows:    rows,
		alive:   make([]bool, len(rows)),
		node:    make([]int32, len(rows)),
		numOrig: len(rows),
		colRows: make([][]int32, numCols),
		weight:  make([]int, numCols),
	}
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

func (r *reducer) enqueue(c int32) {
	if w := r.weight[c]; w == 1 || w == 2 {
		r.queue = append(r.queue, c)
	}
}

// liveRows returns the live rows containing column c, dropping stale entries.
func (r *reducer) liveRows(c int32) []int32 {
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
func (r *reducer) deleteRow(row int32) {
	r.alive[row] = false
	for _, c := range r.rows[row] {
		r.weight[c]--
		r.enqueue(c)
	}
}

// xorInto replaces row dst with dst XOR src.
func (r *reducer) xorInto(dst, src int32) {
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

// reduce removes rows with a column no other row has, merges the two rows
// of each column of weight 2, and pivots on light columns while the fill-in
// stays small.
func (r *reducer) reduce() {
	for {
		r.drain()
		progress := false
		for c := range r.weight {
			if w := r.weight[c]; w < 3 || w > maxPivotWeight {
				continue
			}
			rows := r.liveRows(int32(c))
			piv := 0
			for i := range rows {
				if len(r.rows[rows[i]]) < len(r.rows[rows[piv]]) {
					piv = i
				}
			}
			if (len(rows)-1)*(len(r.rows[rows[piv]])-1) > maxPivotFill {
				continue
			}
			p := rows[piv]
			others := append([]int32(nil), rows[:piv]...)
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

// drain handles every queued column of weight 1 or 2.
func (r *reducer) drain() {
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

// compact returns the live rows with their columns renumbered from 0, the
// number of columns left, and the node of each row.
func (r *reducer) compact() (rows [][]int32, numCols int, nodes []int32) {
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

// expand returns the original rows whose XOR is the XOR of the given nodes,
// as a bitset.
func (r *reducer) expand(roots []int32) []uint64 {
	used := make([]bool, r.numOrig+len(r.merges))
	for _, n := range roots {
		used[n] = !used[n]
	}
	// A merge node only refers to earlier nodes, so one pass from the top
	// pushes every selection down to the original rows.
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

// denseDependencies returns subsets of rows that sum to zero mod 2, as
// bitsets over row indices. Each row lists its columns, each at most once.
func denseDependencies(rows [][]int32, numCols int) [][]uint64 {
	numRows := len(rows)
	colWords := (numCols + 63) / 64
	rowWords := (numRows + 63) / 64

	// Each row holds its vector and a record of which original rows were
	// XORed into it, starting with just itself.
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

	// Reduce to row echelon form, swapping each pivot up so that only the
	// rows below it need clearing. Words before the pivot's are already
	// zero in every row below.
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

	// The rows below the last pivot are now all zero.
	return combo[pivots:]
}

package gnfs

import (
	"sort"

	"github.com/dgryski/go-factor/internal/linalg"
)

type lpKey struct {
	p uint64
	r uint64 // 0 for rational LP, R+1 for algebraic LP (to distinguish from rational)
}

func ratLPKey(p uint64) lpKey {
	return lpKey{p: p, r: 0}
}

func algLPKey(id LargeIdeal) lpKey {
	return lpKey{p: id.P, r: id.R + 1}
}

func oddLPKeysInto(dst []lpKey, rel *Relation) []lpKey {
	for _, lp := range rel.RatLPs {
		if lp != 0 {
			dst = append(dst, ratLPKey(lp))
		}
	}
	for _, id := range rel.AlgLPs {
		dst = append(dst, algLPKey(id))
	}
	if len(dst) <= 1 {
		return dst
	}
	for i := 1; i < len(dst); i++ {
		v := dst[i]
		j := i - 1
		for j >= 0 && (dst[j].p > v.p || (dst[j].p == v.p && dst[j].r > v.r)) {
			dst[j+1] = dst[j]
			j--
		}
		dst[j+1] = v
	}
	out := 0
	for i := 0; i < len(dst); {
		j := i + 1
		for j < len(dst) && dst[j] == dst[i] {
			j++
		}
		if (j-i)&1 != 0 {
			dst[out] = dst[i]
			out++
		}
		i = j
	}
	return dst[:out]
}

// oddLPKeys returns the distinct large prime keys appearing an odd number of
// times in rel.
func oddLPKeys(rel *Relation) []lpKey {
	var buf [12]lpKey
	keys := oddLPKeysInto(buf[:0], rel)
	if len(keys) == 0 {
		return nil
	}
	return append([]lpKey(nil), keys...)
}

// lpFilter interns large prime ideals to dense IDs so singleton recounts during
// sieving run in flat reusable buffers without heap allocations.
type lpFilter struct {
	idMap    map[lpKey]int32
	keys     []lpKey
	totalDeg []int32
	relOff   []int32
	relIDs   []int32

	deg   []int32
	head  []int32
	cur   []int32
	adj   []int32
	dead  []uint64
	queue []int32
}

func newLPFilter(capHint int) *lpFilter {
	if capHint < 64 {
		capHint = 64
	}
	relOff := make([]int32, 1, capHint+1)
	return &lpFilter{
		idMap:  make(map[lpKey]int32, capHint),
		relOff: relOff,
	}
}

func (f *lpFilter) add(rel *Relation) {
	var buf [12]lpKey
	keys := oddLPKeysInto(buf[:0], rel)
	for _, k := range keys {
		id, ok := f.idMap[k]
		if !ok {
			id = int32(len(f.keys))
			f.idMap[k] = id
			f.keys = append(f.keys, k)
			f.totalDeg = append(f.totalDeg, 0)
		}
		f.totalDeg[id]++
		f.relIDs = append(f.relIDs, id)
	}
	f.relOff = append(f.relOff, int32(len(f.relIDs)))
}

func ensureInt32Len(s []int32, n int) []int32 {
	if cap(s) < n {
		return make([]int32, n)
	}
	return s[:n]
}

func (f *lpFilter) excess(numFixed int) int {
	numRels := len(f.relOff) - 1
	numLP := len(f.keys)
	words := (numRels + 63) >> 6
	if cap(f.dead) < words {
		f.dead = make([]uint64, words)
	} else {
		f.dead = f.dead[:words]
		clear(f.dead)
	}
	if numLP == 0 {
		return numRels - numFixed
	}

	f.deg = ensureInt32Len(f.deg, numLP)
	copy(f.deg, f.totalDeg)
	f.head = ensureInt32Len(f.head, numLP+1)
	f.cur = ensureInt32Len(f.cur, numLP)
	f.adj = ensureInt32Len(f.adj, len(f.relIDs))
	f.queue = f.queue[:0]

	var sum int32
	activeLPs := 0
	for id := 0; id < numLP; id++ {
		f.head[id] = sum
		f.cur[id] = sum
		d := f.deg[id]
		sum += d
		if d > 0 {
			activeLPs++
			if d == 1 {
				f.queue = append(f.queue, int32(id))
			}
		}
	}
	f.head[numLP] = sum

	for r := 0; r < numRels; r++ {
		start, end := f.relOff[r], f.relOff[r+1]
		for p := start; p < end; p++ {
			id := f.relIDs[p]
			pos := f.cur[id]
			f.adj[pos] = int32(r)
			f.cur[id] = pos + 1
		}
	}

	numDead := 0
	for len(f.queue) > 0 {
		id := f.queue[len(f.queue)-1]
		f.queue = f.queue[:len(f.queue)-1]
		if f.deg[id] != 1 {
			continue
		}
		hStart, hEnd := f.head[id], f.head[id+1]
		for pos := hStart; pos < hEnd; pos++ {
			r := f.adj[pos]
			w, bit := r>>6, uint64(1)<<(uint(r)&63)
			if f.dead[w]&bit != 0 {
				continue
			}
			f.dead[w] |= bit
			numDead++
			rStart, rEnd := f.relOff[r], f.relOff[r+1]
			for p := rStart; p < rEnd; p++ {
				q := f.relIDs[p]
				f.deg[q]--
				if f.deg[q] == 0 {
					activeLPs--
				} else if f.deg[q] == 1 {
					f.queue = append(f.queue, q)
				}
			}
		}
	}
	return (numRels - numDead) - activeLPs - numFixed
}

// pruneSingletons iteratively removes singleton ideals and, when surplus exceeds
// KeepSurplus + 16, prunes heavy weight 2 LP cliques (Cavallar, ANTS IV 2000, Section 3.2).
func pruneSingletons(rels []Relation, fb *FactorBases) (excess int, alive []bool, lpCols map[lpKey]int32) {
	f := newLPFilter(len(rels))
	for i := range rels {
		f.add(&rels[i])
	}
	excess = f.excess(fb.NumFixed)
	const targetExcess = linalg.KeepSurplus + 16
	if excess > targetExcess {
		excess = f.pruneCliques(rels, fb.NumFixed, targetExcess)
	}

	alive = make([]bool, len(rels))
	for i := range rels {
		alive[i] = f.dead[i>>6]&(1<<(uint(i)&63)) == 0
	}

	var active []lpKey
	for id, k := range f.keys {
		if f.deg[id] > 0 {
			active = append(active, k)
		}
	}
	sort.Slice(active, func(i, j int) bool {
		if active[i].p != active[j].p {
			return active[i].p < active[j].p
		}
		return active[i].r < active[j].r
	})
	lpCols = make(map[lpKey]int32, len(active))
	for i, k := range active {
		lpCols[k] = int32(fb.NumFixed + i)
	}
	return excess, alive, lpCols
}

func (f *lpFilter) isDead(r int32) bool {
	return f.dead[r>>6]&(1<<(uint(r)&63)) != 0
}

func (f *lpFilter) pruneCliques(rels []Relation, numFixed, targetExcess int) int {
	numRels := len(rels)
	numLP := len(f.keys)
	numDead := 0
	for r := 0; r < numRels; r++ {
		if f.isDead(int32(r)) {
			numDead++
		}
	}
	activeLPs := 0
	for id := 0; id < numLP; id++ {
		if f.deg[id] > 0 {
			activeLPs++
		}
	}
	curExcess := (numRels - numDead) - activeLPs - numFixed
	if curExcess <= targetExcess {
		return curExcess
	}

	// Build adjacency across live relations linked by weight 2 LP ideals.
	relDeg2Count := make([]int32, numRels)
	type edge2 struct{ u, v int32 }
	var edges []edge2
	for id := 0; id < numLP; id++ {
		if f.deg[id] != 2 {
			continue
		}
		var pair [2]int32
		n := 0
		for pos := f.head[id]; pos < f.head[id+1]; pos++ {
			r := f.adj[pos]
			if !f.isDead(r) {
				if n < 2 {
					pair[n] = r
				}
				n++
			}
		}
		if n == 2 && pair[0] != pair[1] {
			edges = append(edges, edge2{pair[0], pair[1]})
			relDeg2Count[pair[0]]++
			relDeg2Count[pair[1]]++
		}
	}
	if len(edges) == 0 {
		return curExcess
	}

	rHead := make([]int32, numRels+1)
	rCur := make([]int32, numRels)
	var totalAdj int32
	for r := 0; r < numRels; r++ {
		rHead[r] = totalAdj
		rCur[r] = totalAdj
		totalAdj += relDeg2Count[r]
	}
	rHead[numRels] = totalAdj
	rAdj := make([]int32, totalAdj)
	for _, e := range edges {
		rAdj[rCur[e.u]] = e.v
		rCur[e.u]++
		rAdj[rCur[e.v]] = e.u
		rCur[e.v]++
	}

	type cliqueCand struct {
		rels   []int32
		weight int32
	}
	visited := make([]bool, numRels)
	var cliques []cliqueCand
	var comp []int32
	for r := 0; r < numRels; r++ {
		if relDeg2Count[r] == 0 || visited[r] {
			continue
		}
		comp = comp[:0]
		visited[r] = true
		comp = append(comp, int32(r))
		var edgeTwice int
		for head := 0; head < len(comp); head++ {
			u := comp[head]
			for pos := rHead[u]; pos < rHead[u+1]; pos++ {
				v := rAdj[pos]
				edgeTwice++
				if !visited[v] {
					visited[v] = true
					comp = append(comp, v)
				}
			}
		}
		if len(comp) >= 2 && len(comp) == edgeTwice/2+1 {
			var w int32
			for _, u := range comp {
				w += int32(len(rels[u].RatCols)+len(rels[u].AlgCols)) + (f.relOff[u+1] - f.relOff[u])
			}
			cliques = append(cliques, cliqueCand{
				rels:   append([]int32(nil), comp...),
				weight: w,
			})
		}
	}
	if len(cliques) == 0 {
		return curExcess
	}

	sort.Slice(cliques, func(i, j int) bool {
		if cliques[i].weight != cliques[j].weight {
			return cliques[i].weight > cliques[j].weight
		}
		return cliques[i].rels[0] < cliques[j].rels[0]
	})

	lpTouch := make([]int32, numLP)
	var touched []int32
	for _, cl := range cliques {
		if curExcess <= targetExcess {
			break
		}
		allAlive := true
		touched = touched[:0]
		for _, r := range cl.rels {
			if f.isDead(r) {
				allAlive = false
				break
			}
			for p := f.relOff[r]; p < f.relOff[r+1]; p++ {
				q := f.relIDs[p]
				if lpTouch[q] == 0 {
					touched = append(touched, q)
				}
				lpTouch[q]++
			}
		}
		safe := allAlive
		if safe {
			for _, q := range touched {
				rem := f.deg[q] - lpTouch[q]
				if rem == 1 || (rem == 0 && f.deg[q] != 2) {
					safe = false
					break
				}
			}
		}
		for _, q := range touched {
			lpTouch[q] = 0
		}
		if !safe {
			continue
		}

		// Delete cl.rels[0] and cascade singletons across the tree.
		r0 := cl.rels[0]
		f.dead[r0>>6] |= 1 << (uint(r0) & 63)
		numDead++
		f.queue = f.queue[:0]
		for p := f.relOff[r0]; p < f.relOff[r0+1]; p++ {
			q := f.relIDs[p]
			f.deg[q]--
			if f.deg[q] == 0 {
				activeLPs--
			} else if f.deg[q] == 1 {
				f.queue = append(f.queue, q)
			}
		}
		for len(f.queue) > 0 {
			id := f.queue[len(f.queue)-1]
			f.queue = f.queue[:len(f.queue)-1]
			if f.deg[id] != 1 {
				continue
			}
			for pos := f.head[id]; pos < f.head[id+1]; pos++ {
				r := f.adj[pos]
				if f.isDead(r) {
					continue
				}
				f.dead[r>>6] |= 1 << (uint(r) & 63)
				numDead++
				for p := f.relOff[r]; p < f.relOff[r+1]; p++ {
					q := f.relIDs[p]
					f.deg[q]--
					if f.deg[q] == 0 {
						activeLPs--
					} else if f.deg[q] == 1 {
						f.queue = append(f.queue, q)
					}
				}
			}
		}
		curExcess = (numRels - numDead) - activeLPs - numFixed
	}
	return curExcess
}

// buildMatrix prunes singletons and constructs the sparse GF(2) rows for SGE + Block Lanczos.
func buildMatrix(poly *Poly, fb *FactorBases, rels []Relation) (surviving []Relation, rows [][]int32, numCols int) {
	_, alive, lpCols := pruneSingletons(rels, fb)
	numCols = fb.NumFixed + len(lpCols)

	var scratchVal = newCofactorScratch(poly.Degree)
	var factors []int32
	var lpBuf [12]lpKey
	badEpoch := make([]int, len(fb.BadDecomps))
	hasAFBCol := len(fb.AFBCol) == len(fb.AFB)
	for i := range rels {
		if !alive[i] {
			continue
		}
		rel := rels[i]
		surviving = append(surviving, rel)

		factors = factors[:0]
		if !rel.IsFree {
			factors = append(factors, ColAllOnes)
			if evalGInto(&scratchVal.val, poly, rel.A, rel.B, scratchVal).Sign() < 0 {
				factors = append(factors, ColRatSign)
			}
			if evalFInto(&scratchVal.val, poly, rel.A, rel.B, scratchVal).Sign() < 0 {
				factors = append(factors, ColAlgSign)
			}
			for bi := range fb.BadDecomps {
				if fb.BadDecomps[bi].DividesCD {
					badEpoch[bi] = i + 1
					factors = fb.BadDecomps[bi].appendOddCols(rel.A, rel.B, factors)
				}
			}
		}
		for _, idx := range rel.RatCols {
			factors = append(factors, int32(fb.ColRFB)+idx)
		}
		for _, idx := range rel.AlgCols {
			if !hasAFBCol {
				factors = append(factors, int32(fb.ColAFB)+idx)
				continue
			}
			if col := fb.AFBCol[idx]; col >= 0 {
				factors = append(factors, col)
			} else if !rel.IsFree {
				bi := fb.BadByP[fb.AFB[idx].P]
				if badEpoch[bi] != i+1 {
					badEpoch[bi] = i + 1
					factors = fb.BadDecomps[bi].appendOddCols(rel.A, rel.B, factors)
				}
			}
		}
		factors = fb.QCBits(rel.A, rel.B, factors)
		for _, k := range oddLPKeysInto(lpBuf[:0], &rel) {
			factors = append(factors, lpCols[k])
		}
		rows = append(rows, linalg.OddColumns(factors))
	}
	return surviving, rows, numCols
}

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

// oddLPKeys returns the distinct large-prime keys appearing an odd number of
// times in rel.
func oddLPKeys(rel *Relation) []lpKey {
	var keys [5]lpKey
	n := 0
	for _, lp := range rel.RatLPs {
		if lp != 0 {
			keys[n] = ratLPKey(lp)
			n++
		}
	}
	for _, id := range rel.AlgLPs {
		if n < len(keys) {
			keys[n] = algLPKey(id)
			n++
		}
	}
	if n == 0 {
		return nil
	}
	if n == 1 {
		return []lpKey{keys[0]}
	}
	if len(rel.AlgLPs)+2 > len(keys) {
		all := make([]lpKey, 0, len(rel.AlgLPs)+2)
		for _, lp := range rel.RatLPs {
			if lp != 0 {
				all = append(all, ratLPKey(lp))
			}
		}
		for _, id := range rel.AlgLPs {
			all = append(all, algLPKey(id))
		}
		return oddLPKeySlice(all)
	}
	return oddLPKeySlice(keys[:n])
}

func oddLPKeySlice(s []lpKey) []lpKey {
	sort.Slice(s, func(i, j int) bool {
		if s[i].p != s[j].p {
			return s[i].p < s[j].p
		}
		return s[i].r < s[j].r
	})
	var out []lpKey
	for i := 0; i < len(s); {
		j := i
		for j < len(s) && s[j] == s[i] {
			j++
		}
		if (j-i)&1 != 0 {
			out = append(out, s[i])
		}
		i = j
	}
	return out
}

// pruneSingletons iteratively removes relations containing a large-prime ideal
// of weight 1 until a 2-core remains. It returns the excess (surviving relations
// minus active LP ideals minus fb.NumFixed), alive flags, and column map for
// surviving LP ideals starting at fb.NumFixed.
func pruneSingletons(rels []Relation, fb *FactorBases) (excess int, alive []bool, lpCols map[lpKey]int32) {
	degree := make(map[lpKey]int, len(rels))
	byKey := make(map[lpKey][]int32, len(rels))
	relKeys := make([][]lpKey, len(rels))

	for i := range rels {
		keys := oddLPKeys(&rels[i])
		relKeys[i] = keys
		for _, k := range keys {
			degree[k]++
			byKey[k] = append(byKey[k], int32(i))
		}
	}

	alive = make([]bool, len(rels))
	for i := range alive {
		alive[i] = true
	}

	var queue []lpKey
	for k, d := range degree {
		if d == 1 {
			queue = append(queue, k)
		}
	}

	for len(queue) > 0 {
		k := queue[len(queue)-1]
		queue = queue[:len(queue)-1]
		if degree[k] != 1 {
			continue
		}
		for _, i := range byKey[k] {
			if !alive[i] {
				continue
			}
			alive[i] = false
			for _, q := range relKeys[i] {
				if degree[q]--; degree[q] == 1 {
					queue = append(queue, q)
				}
			}
		}
	}

	var active []lpKey
	for k, d := range degree {
		if d > 0 {
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

	numAlive := 0
	for _, a := range alive {
		if a {
			numAlive++
		}
	}
	excess = numAlive - len(lpCols) - fb.NumFixed
	return excess, alive, lpCols
}

// buildMatrix prunes singletons and constructs the sparse GF(2) rows for SGE + Block Lanczos.
func buildMatrix(poly *Poly, fb *FactorBases, rels []Relation) (surviving []Relation, rows [][]int32, numCols int) {
	_, alive, lpCols := pruneSingletons(rels, fb)
	numCols = fb.NumFixed + len(lpCols)

	var scratchVal = newCofactorScratch(poly.Degree)
	var factors []int32
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
		}
		for _, idx := range rel.RatCols {
			factors = append(factors, int32(fb.ColRFB)+idx)
		}
		for _, idx := range rel.AlgCols {
			factors = append(factors, int32(fb.ColAFB)+idx)
		}
		factors = fb.QCBits(rel.A, rel.B, factors)
		for _, k := range oddLPKeys(&rel) {
			factors = append(factors, lpCols[k])
		}
		rows = append(rows, linalg.OddColumns(factors))
	}
	return surviving, rows, numCols
}

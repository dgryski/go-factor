package gnfs

import "github.com/dgryski/go-factor/internal/arith"

const defaultRelTargetRatio = 0.95

// relTarget estimates the relation count needed from LP bounds first,
// then extrapolates excess slope once the 2-core starts climbing.
type relTarget struct {
	initial    int
	goal       int
	prevRels   int
	prevExcess int
	havePrev   bool
	want       int
}

func newRelTarget(p *Params, goal int) *relTarget {
	ratio := p.RelTargetRatio
	if ratio <= 0 {
		ratio = defaultRelTargetRatio
	}
	ideals := len(arith.PrimesUpTo(int64(p.LpBoundRat))) + len(arith.PrimesUpTo(int64(p.LpBoundAlg)))
	initial := int(ratio * float64(ideals))
	return &relTarget{initial: initial, goal: goal, want: initial}
}

func (t *relTarget) update(rels, excess int) int {
	rising := t.havePrev && excess > t.prevExcess && rels > t.prevRels
	if rising {
		slope := float64(excess-t.prevExcess) / float64(rels-t.prevRels)
		est := rels + int(float64(t.goal-excess)/slope+0.5)
		// Ignore tiny pre-percolation slopes before reaching the initial target.
		if rels < t.initial && est > t.initial && slope < 0.2 {
			t.want = t.initial
		} else {
			t.want = est
		}
	} else if t.want <= rels {
		t.want = 0
	}
	t.prevRels, t.prevExcess, t.havePrev = rels, excess, true
	return t.want
}

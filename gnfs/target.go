package gnfs

import "github.com/dgryski/go-factor/internal/arith"

// relTargetRatio is the number of relations needed as a fraction of the
// large-prime ideals, π(LpBoundRat) + π(LpBoundAlg), measured on test runs.
const relTargetRatio = 0.9

// relTarget estimates how many relations the sieve needs. Excess makes a
// poor progress count: it stays at -NumFixed until large-prime cycles
// form, dips, then climbs steeply. The number of relations grows steadily,
// so progress is measured in relations against this target instead.
//
// The first estimate comes from the large-prime bounds. Once excess is
// rising, the target extrapolates its slope against the relation count.
// Excess grows faster than linearly, so this overestimates and then
// converges as the sieve nears the goal.
type relTarget struct {
	initial    int
	goal       int // excess wanted
	prevRels   int
	prevExcess int
	havePrev   bool
	want       int // current estimate, 0 if unknown
}

func newRelTarget(p *Params, goal int) *relTarget {
	ideals := len(arith.PrimesUpTo(int64(p.LpBoundRat))) + len(arith.PrimesUpTo(int64(p.LpBoundAlg)))
	initial := int(relTargetRatio * float64(ideals))
	return &relTarget{initial: initial, goal: goal, want: initial}
}

// update records the excess after rels relations and returns the target,
// or 0 if there is no usable estimate.
func (t *relTarget) update(rels, excess int) int {
	rising := t.havePrev && excess > t.prevExcess && rels > t.prevRels
	if rising {
		slope := float64(excess-t.prevExcess) / float64(rels-t.prevRels)
		t.want = rels + int(float64(t.goal-excess)/slope+0.5)
	} else if t.want <= rels {
		// The initial estimate was too low and excess isn't climbing yet.
		t.want = 0
	}
	t.prevRels, t.prevExcess, t.havePrev = rels, excess, true
	return t.want
}

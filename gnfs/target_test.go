package gnfs

import "testing"

func TestRelTarget(t *testing.T) {
	// Excess after each checkpoint of a 230-bit run: flat at -NumFixed,
	// a dip as cycles form, then a steep climb.
	steps := []struct {
		rels, excess int
		want         int // 0: unknown, -1: the initial estimate
	}{
		{33471, -3585, -1},
		{198146, -3585, -1},
		{211855, -14891, -1},
		{226478, -10610, 262829},
		{242076, -5074, 256462},
	}
	tg := &relTarget{initial: 200000, goal: 32, want: 200000}
	for _, s := range steps {
		want := s.want
		if want == -1 {
			want = tg.initial
		}
		if s.rels >= tg.initial && s.want == -1 {
			want = 0 // the initial estimate has been passed
		}
		if got := tg.update(s.rels, s.excess); got != want {
			t.Errorf("update(%d, %d) = %d, want %d", s.rels, s.excess, got, want)
		}
	}
}

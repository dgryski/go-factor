package gnfs

import "testing"

func TestRelTarget(t *testing.T) {
	steps := []struct {
		rels, excess int
		want         int
	}{
		{33471, -19000, -1},
		{43276, -18841, -1},
		{83276, -6841, 106186},
		{150000, -5000, -1},
		{170000, 0, 170128},
		{175000, 0, 0},
		{226478, 10, 339730},
	}
	tg := &relTarget{initial: 200000, goal: 32, want: 200000}
	for _, s := range steps {
		want := s.want
		if want == -1 {
			want = tg.initial
		}
		if got := tg.update(s.rels, s.excess); got != want {
			t.Errorf("update(%d, %d) = %d, want %d", s.rels, s.excess, got, want)
		}
	}
}

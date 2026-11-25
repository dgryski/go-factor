package progress

import (
	"testing"
	"time"
)

// TestRateBursty feeds a count that stalls for 30s and then catches up in
// 5s, at 100 a second on average, starting 12 minutes in. The estimate
// should stay near the average rate instead of swinging with each burst.
func TestRateBursty(t *testing.T) {
	var r Rate
	now := time.Unix(0, 0)
	n := 0
	r.Sample(now, n)
	now = now.Add(12 * time.Minute)
	n = 72000
	r.Sample(now, n)
	const want = 500000
	for i := 0; i < 20; i++ {
		now = now.Add(30 * time.Second)
		n += 500
		r.Sample(now, n)
		now = now.Add(5 * time.Second)
		n += 3000
		r.Sample(now, n)
		if i < 2 {
			continue
		}
		eta := time.Duration(float64(want-n) / r.rate * float64(time.Second))
		avg := time.Duration(want-n) * time.Second / 100
		if eta < avg*8/10 || eta > avg*12/10 {
			t.Errorf("step %d: estimate %v, want about %v", i, eta.Round(time.Second), avg.Round(time.Second))
		}
	}
}

func TestRateFirstSampleOnlyStarts(t *testing.T) {
	var r Rate
	now := time.Unix(0, 0)
	r.Sample(now, 0)
	if got := r.Left(0, 10); got != "time left unknown" {
		t.Errorf("Left after one sample = %q", got)
	}
	r.Sample(now.Add(5*time.Second), 5)
	if got := r.Left(5, 10); got != "about 5s left" {
		t.Errorf("Left = %q, want about 5s left", got)
	}
	if got := r.Left(10, 10); got != "none left" {
		t.Errorf("Left at want = %q", got)
	}
}

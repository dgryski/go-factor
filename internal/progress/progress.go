// Package progress estimates the time left in a long computation from the
// rate at which a count grows between progress messages.
package progress

import (
	"fmt"
	"math"
	"time"
)

// A sample covers at least a tenth of the time since the first, kept
// between these bounds. The count can stall and then jump, as when the
// quadratic sieve waits for a slow batch, and a window that spans a stall
// and the burst after it measures the rate in between.
const (
	minWindow = 5 * time.Second
	maxWindow = 30 * time.Second
)

// Rate is an EWMA of how fast a count grows, in units per second. Each
// sample is weighted by its length, with a time constant of twice the
// window.
type Rate struct {
	start time.Time // the first sample
	last  time.Time // the end of the last window
	lastN int
	rate  float64
}

// Sample records the count n at now. The first call only sets the starting
// point. Later calls fold the growth since the last window into the rate
// once the window is long enough.
func (r *Rate) Sample(now time.Time, n int) {
	if r.start.IsZero() {
		r.start, r.last, r.lastN = now, now, n
		return
	}
	window := min(max(now.Sub(r.start)/10, minWindow), maxWindow)
	dt := now.Sub(r.last)
	if dt < window {
		return
	}
	s := max(0, float64(n-r.lastN)/dt.Seconds())
	if r.rate == 0 {
		r.rate = s
	} else {
		a := 1 - math.Exp(-dt.Seconds()/(2*window).Seconds())
		r.rate = a*s + (1-a)*r.rate
	}
	r.last, r.lastN = now, n
}

// Left describes the time for the count to grow from n to want at the
// current rate.
func (r *Rate) Left(n, want int) string {
	if n >= want {
		return "none left"
	}
	if r.rate == 0 {
		return "time left unknown"
	}
	eta := time.Duration(float64(want-n) / r.rate * float64(time.Second))
	return fmt.Sprintf("about %v left", eta.Round(time.Second))
}

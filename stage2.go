package factor

import "sync"

// stage2D is the giant step of the baby-step giant-step stage 2 used by ECM
// and p+1: 2*3*5*7*11.
const stage2D = 2310

// stage2Plan writes each prime q in (b1, b2] as c±m, with c a multiple of
// stage2D and m an odd baby step at most stage2D/2 and coprime to stage2D.
// Both methods have a value whose difference for c and m is zero mod p
// exactly when c-m or c+m is, so one product term covers both.
type stage2Plan struct {
	b1, b2  int
	isPrime []uint64 // bitset of primes <= b2 + stage2D, from sieve
	babies  []int
	once    sync.Once
}

func newStage2Plan(b1, b2 int) *stage2Plan {
	if b2 < b1 {
		b2 = b1
	}
	s := &stage2Plan{b1: b1, b2: b2}
	for m := 1; m <= stage2D/2; m += 2 {
		if gcdInt(m, stage2D) == 1 {
			s.babies = append(s.babies, m)
		}
	}
	return s
}

// sieve builds the prime bitset. Stage 2 calls it on entry, so a deadline
// that passes in stage 1 doesn't pay for it.
func (s *stage2Plan) sieve() {
	s.once.Do(func() { s.isPrime = primeBits(s.b2 + stage2D) })
}

// prime reports whether q is a prime in (b1, b2]. It needs sieve.
func (s *stage2Plan) prime(q int) bool {
	return q > s.b1 && q <= s.b2 && s.isPrime[q/64]&(1<<(q%64)) != 0
}

// first returns the first giant step: the multiple of stage2D nearest b1,
// at least 1.
func (s *stage2Plan) first() int {
	return max(1, (s.b1+stage2D/2)/stage2D)
}

// last reports whether giant step iv covers b2.
func (s *stage2Plan) last(iv int) bool {
	return iv*stage2D+stage2D/2 >= s.b2
}

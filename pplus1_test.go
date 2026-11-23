package factor

import (
	"math/big"
	"testing"
)

func TestLucasV(t *testing.T) {
	// V_k(3) for Q = 1 is 2, 3, 7, 18, 47, 123, 322, ...
	want := []int64{3, 7, 18, 47, 123, 322, 843, 2207}
	n := big.NewInt(1000003)
	var r big.Int
	for i, w := range want {
		lucasV(&r, big.NewInt(3), uint64(i+1), n)
		if r.Int64() != w {
			t.Errorf("V_%d(3) = %v, want %v", i+1, &r, w)
		}
	}
	// V_mk(P) = V_m(V_k(P))
	var a, b big.Int
	lucasV(&a, big.NewInt(5), 35, n)
	lucasV(&b, big.NewInt(5), 5, n)
	lucasV(&b, &b, 7, n)
	if a.Cmp(&b) != 0 {
		t.Errorf("V_35(5) = %v, V_7(V_5(5)) = %v", &a, &b)
	}
}

// smoothPlus1 returns the n'th prime p, counting from 0, with p+1 = m*s*r for a
// fixed 2000-smooth s and a prime r below 2000.
func smoothPlus1(m int64, n int) *big.Int {
	s := big.NewInt(m * 1024 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23)
	p := new(big.Int)
	for _, r := range smallPrimes[9:] {
		p.Mul(s, big.NewInt(r))
		p.Sub(p, one)
		if p.ProbablyPrime(20) {
			if n == 0 {
				return p
			}
			n--
		}
	}
	panic("no prime found")
}

func TestPPlus1(t *testing.T) {
	// q has q-1 and q+1 with large prime factors.
	q := nextPrime("1000000000000000003")
	tests := []struct {
		name string
		p    *big.Int
	}{
		{"stage 1", smoothPlus1(1, 0)},
		{"stage 2", smoothPlus1(100003, 0)},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			n := new(big.Int).Mul(tt.p, q)
			pr, c := PPlus1(n, 2000, 10)
			verifyFactoring(t, n, pr, c)
			if len(c) != 0 {
				t.Errorf("PPlus1(%v) left composites %v", n, c)
			}
		})
	}
}

func TestPPlus1Backtrack(t *testing.T) {
	// p-1 and p+1 are smooth over the first batch of primes for both, so
	// both are found in that batch whatever the starting value.
	n := big.NewInt(2944999 * 2740219)
	pr, c := PPlus1(n, 2000, 1)
	verifyFactoring(t, big.NewInt(2944999*2740219), pr, c)
	if len(c) != 0 {
		t.Errorf("PPlus1(%v) left composites %v", n, c)
	}
}

func TestPPlus1Small(t *testing.T) {
	for _, v := range []int64{1, 2, 4, 7, 2 * 1000003} {
		n := big.NewInt(v)
		p, c := PPlus1(n, 2000, 4)
		verifyProduct(t, big.NewInt(v), p, c)
		if len(c) != 0 {
			t.Errorf("PPlus1(%v) left composites %v", v, c)
		}
	}
}

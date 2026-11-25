package factor

import (
	"math/big"
	"testing"

	"github.com/dgryski/go-factor/internal/testutil"
)

func nextPrime(s string) *big.Int {
	n, _ := new(big.Int).SetString(s, 10)
	for !n.ProbablyPrime(20) {
		n.Add(n, one)
	}
	return n
}

func product(ps ...*big.Int) *big.Int {
	n := big.NewInt(1)
	for _, p := range ps {
		n.Mul(n, p)
	}
	return n
}

func TestECM(t *testing.T) {
	tests := []struct {
		name    string
		factors []*big.Int
	}{
		{"13+27 digits", []*big.Int{nextPrime("1000000000039"), nextPrime("314159265358979323846264338")}},
		{"three 12 digit", []*big.Int{nextPrime("100000000000"), nextPrime("271828182845"), nextPrime("577215664901")}},
		{"repeated", []*big.Int{nextPrime("1000003"), nextPrime("1000003"), nextPrime("99999999977")}},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			n := product(tt.factors...)
			orig := new(big.Int).Set(n)
			p, c := ECM(n, 2000, 50)
			testutil.VerifyFactoring(t, orig, p, c)
			if len(c) != 0 {
				t.Errorf("composites left: %v", c)
			}
			if len(p) != len(tt.factors) {
				t.Errorf("got %d primes %v, want %d", len(p), p, len(tt.factors))
			}
		})
	}
}

func TestECMPrime(t *testing.T) {
	n := nextPrime("1000000000039")
	p, c := ECM(n, 2000, 10)
	if len(p) != 1 || p[0].Cmp(n) != 0 || len(c) != 0 {
		t.Errorf("ECM(prime) = %v, %v", p, c)
	}
	p, c = ECM(big.NewInt(1), 2000, 10)
	if len(p) != 0 || len(c) != 0 {
		t.Errorf("ECM(1) = %v, %v", p, c)
	}
}

func TestECMTorsion12(t *testing.T) {
	for _, p := range []int64{1009, 1013, 2003, 4001} {
		n := big.NewInt(p)
		for k := int64(2); k < 10; k++ {
			c := &ecCurve{n: n}
			var pt xzPoint
			if g := c.torsion12(&pt, k); g != nil {
				continue
			}

			// Count points on By^2 = x^3 + Ax^2 + x with B chosen so pt is on it.
			var A, f0, f, x big.Int
			A.Lsh(&c.a24, 2)
			A.Sub(&A, two)
			cubic := func(r, x *big.Int) {
				var t big.Int
				t.Add(x, &A)
				t.Mul(&t, x)
				t.Add(&t, one)
				t.Mul(&t, x)
				r.Mod(&t, n)
			}
			cubic(&f0, &pt.x)
			order := int64(1)
			for i := int64(0); i < p; i++ {
				x.SetInt64(i)
				cubic(&f, &x)
				f.Mul(&f, &f0)
				order += int64(1 + big.Jacobi(f.Mod(&f, n), n))
			}

			if order%12 != 0 {
				t.Errorf("p=%d k=%d: order %d not divisible by 12", p, k, order)
			}
			var q xzPoint
			c.mul(&q, &pt, uint64(order))
			if q.z.Sign() != 0 {
				t.Errorf("p=%d k=%d: order*P is not O", p, k)
			}
		}
	}
}

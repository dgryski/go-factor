package factor

import (
	"context"
	"math/big"
	"testing"
	"time"

	"github.com/dgryski/go-factor/internal/testutil"
)

func TestPMinus1(t *testing.T) {
	// 2500069784293916113216294633942411
	// 923874293742934712938741239847123984712394871
	n := big.NewInt(3928471)
	orig := (&big.Int{}).Set(n)
	p, c := PMinus1(n)
	t.Log(p, c)
	testutil.VerifyFactoring(t, orig, p, c)
}

func TestPMinus1SmallFactor(t *testing.T) {
	for _, v := range []int64{7 * 1000003, 11 * 1000003, 2 * 3928471, 4 * 3928471, 2, 4} {
		n := big.NewInt(v)
		orig := new(big.Int).Set(n)
		p, c := PMinus1(n)
		testutil.VerifyFactoring(t, orig, p, c)
		if len(c) != 0 {
			t.Errorf("PMinus1(%v) left composites %v", orig, c)
		}
	}
}

func TestPMinus1LastBatch(t *testing.T) {
	// p-1 has the largest small prime as a factor.
	last := smallPrimes[len(smallPrimes)-1]
	p := new(big.Int)
	for k := int64(1); ; k++ {
		p.SetInt64(2*k*last + 1)
		if p.ProbablyPrime(20) {
			break
		}
	}
	n := new(big.Int).Mul(p, big.NewInt(1000003))
	orig := new(big.Int).Set(n)
	pr, c := PMinus1(n)
	testutil.VerifyFactoring(t, orig, pr, c)
}

func TestPMinus1Backtrack(t *testing.T) {
	// p-1 and q-1 are both 3-smooth, so both are found in the first batch.
	p := new(big.Int).Lsh(big.NewInt(3), 30)
	p.Add(p, one)
	q := new(big.Int).Lsh(big.NewInt(81), 32)
	q.Add(q, one)
	if !p.ProbablyPrime(20) || !q.ProbablyPrime(20) {
		t.Fatal("test factors are not prime")
	}
	n := new(big.Int).Mul(p, q)
	orig := new(big.Int).Set(n)
	pr, c := PMinus1(n)
	testutil.VerifyFactoring(t, orig, pr, c)
	if len(c) != 0 {
		t.Errorf("PMinus1(%v) left composites %v", orig, c)
	}
}

func TestPMinus1Context(t *testing.T) {
	// p-1 = 3*2^30 and q-1 = 3^4*2^32 need only B1=3, with powers of 2 up to n.
	p := new(big.Int).Lsh(big.NewInt(3), 30)
	p.Add(p, one)
	q := new(big.Int).Lsh(big.NewInt(81), 32)
	q.Add(q, one)
	n := new(big.Int).Mul(p, q)
	orig := new(big.Int).Set(n)
	pr, c := PMinus1Context(context.Background(), n, 3)
	testutil.VerifyFactoring(t, orig, pr, c)
	if len(c) != 0 {
		t.Errorf("PMinus1Context(%v, 3) left composites %v", orig, c)
	}

	// 1000003-1 = 2*3*166667, so B1 below 166667 must not find it. The
	// other factor is a safe prime 2r+1 with r > 10^6, out of reach of both.
	r := nextPrime("1000000")
	for {
		q.Lsh(r, 1)
		q.Add(q, one)
		if q.ProbablyPrime(20) {
			break
		}
		r = nextPrime(r.Add(r, one).String())
	}
	n = new(big.Int).Mul(big.NewInt(1000003), q)
	if pr, _ := PMinus1Context(context.Background(), n, 100000); len(pr) != 0 {
		t.Errorf("PMinus1Context(%v, 100000) = %v, want nothing below B1", n, pr)
	}
	if pr, _ := PMinus1Context(context.Background(), n, 166667); len(pr) == 0 {
		t.Errorf("PMinus1Context(%v, 166667) found nothing", n)
	}
}

// TestContextCancel checks that each method returns soon after its deadline
// on a 60 digit semiprime it cannot split.
func TestContextCancel(t *testing.T) {
	n := product(nextPrime("314159265358979323846264338327"), nextPrime("271828182845904523536028747135"))
	methods := []struct {
		name string
		f    func(ctx context.Context) ([]*big.Int, []*big.Int)
	}{
		{"brent", func(ctx context.Context) ([]*big.Int, []*big.Int) { return BrentContext(ctx, n, 2, 1) }},
		{"p-1", func(ctx context.Context) ([]*big.Int, []*big.Int) { return PMinus1Context(ctx, n, 10000000) }},
		{"p+1", func(ctx context.Context) ([]*big.Int, []*big.Int) { return PPlus1Context(ctx, n, 1000000, 10) }},
	}
	for _, m := range methods {
		ctx, cancel := context.WithTimeout(context.Background(), 20*time.Millisecond)
		start := time.Now()
		p, c := m.f(ctx)
		cancel()
		if elapsed := time.Since(start); elapsed > time.Second {
			t.Errorf("%s took %v after a 20ms deadline", m.name, elapsed)
		}
		if len(p) != 0 || len(c) != 1 || c[0].Cmp(n) != 0 {
			t.Errorf("%s = %v, %v, want n back", m.name, p, c)
		}
	}
}

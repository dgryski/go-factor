package factor

import (
	"math/big"
	"testing"
	"time"

	"github.com/dgryski/go-factor/internal/testutil"
)

func TestTrial(t *testing.T) {
	n := big.NewInt(275)
	orig := (&big.Int{}).Set(n)
	p, c := Trial(n)
	t.Log(p, c)
	testutil.VerifyFactoring(t, orig, p, c)
}

func TestInputUnchanged(t *testing.T) {
	for _, tt := range []struct {
		name string
		n    string
		f    func(*big.Int) ([]*big.Int, []*big.Int)
	}{
		{"Trial", "275", Trial},
		{"PMinus1", "3928471", PMinus1},
		{"Brent", "3928471", func(n *big.Int) ([]*big.Int, []*big.Int) { return Brent(n, 19, 1) }},
		{"Squfof", "29384129384129384712938471", func(n *big.Int) ([]*big.Int, []*big.Int) { return Squfof(n, 5) }},
	} {
		n, _ := new(big.Int).SetString(tt.n, 10)
		tt.f(n)
		if n.String() != tt.n {
			t.Errorf("%s changed its input %v to %v", tt.name, tt.n, n)
		}
	}
}

func TestTrialNonPositive(t *testing.T) {
	for _, v := range []int64{0, -1, -12} {
		n := big.NewInt(v)
		var p, c []*big.Int
		testutil.WithTimeout(t, 5*time.Second, "Trial("+n.String()+")", func() {
			p, c = Trial(n)
		})
		testutil.VerifyProduct(t, big.NewInt(v), p, c)
	}
}

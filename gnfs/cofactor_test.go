package gnfs

import (
	"math/big"
	"math/bits"
	"math/rand"
	"testing"

	"github.com/dgryski/go-factor/internal/testutil"
)

func bigFromU128(v u128) *big.Int {
	z := new(big.Int).SetUint64(v.hi)
	z.Lsh(z, 64)
	z.Or(z, new(big.Int).SetUint64(v.lo))
	return z
}

func u128FromBig(z *big.Int) u128 {
	hi, lo, ok := bigToU128(z)
	if !ok {
		panic("big.Int > 128 bits")
	}
	return u128{hi: hi, lo: lo}
}

func TestMont128AndPrimality(t *testing.T) {
	rnd := rand.New(rand.NewSource(101))
	for _, nBits := range []int{65, 80, 100, 127, 128} {
		for trial := 0; trial < 15; trial++ {
			lim := new(big.Int).Lsh(bigOne, uint(nBits))
			nBig := new(big.Int).Rand(rnd, lim)
			nBig.SetBit(nBig, nBits-1, 1)
			nBig.SetBit(nBig, 0, 1)
			n := u128FromBig(nBig)
			m := newMont128(n)

			aBig := new(big.Int).Rand(rnd, nBig)
			bBig := new(big.Int).Rand(rnd, nBig)
			ma := m.to(u128FromBig(aBig))
			mb := m.to(u128FromBig(bBig))
			gotMul := m.mul(ma, mb)
			got := m.mul(gotMul, u128{lo: 1})
			want := new(big.Int).Mul(aBig, bBig)
			want.Mod(want, nBig)
			if bigFromU128(got).Cmp(want) != 0 {
				t.Fatalf("mont128 mul(%d bits): got %v, want %v", nBits, bigFromU128(got), want)
			}

			expBig := new(big.Int).Rand(rnd, lim)
			gotPow := m.mul(m.pow128(ma, u128FromBig(expBig)), u128{lo: 1})
			wantPow := new(big.Int).Exp(aBig, expBig, nBig)
			if bigFromU128(gotPow).Cmp(wantPow) != 0 {
				t.Fatalf("mont128 pow128(%d bits): got %v, want %v", nBits, bigFromU128(gotPow), wantPow)
			}
		}
	}

	for _, pBits := range []int{66, 85, 110} {
		p := testutil.RandomPrime(rnd, pBits)
		if !isPrime128(u128FromBig(p)) {
			t.Fatalf("isPrime128(%v) = false on %d-bit prime", p, pBits)
		}
		comp := new(big.Int).Mul(testutil.RandomPrime(rnd, pBits/2), testutil.RandomPrime(rnd, pBits-pBits/2))
		if isPrime128(u128FromBig(comp)) {
			t.Fatalf("isPrime128(%v) = true on composite", comp)
		}
	}
}

func TestSplitNormCofactor128And3LP(t *testing.T) {
	rnd := rand.New(rand.NewSource(202))

	// Test 2LP cofactors above 64 bits from two 35 bit primes.
	const lpBound2 = uint64(1) << 36
	for i := 0; i < 20; i++ {
		p1 := testutil.RandomPrime(rnd, 35).Uint64()
		p2 := testutil.RandomPrime(rnd, 35).Uint64()
		if p1 > p2 {
			p1, p2 = p2, p1
		}
		hi, lo := bits.Mul64(p1, p2)
		ns := normState{hi: hi, lo: lo, use128: true}
		got, ok := splitNormCofactor(&ns, 500000, lpBound2, 2, 72)
		if !ok || got[0] != p1 || got[1] != p2 || got[2] != 0 {
			t.Fatalf("2LP 70-bit split(%d * %d): got %v ok=%v", p1, p2, got, ok)
		}
	}

	// Test 3LP cofactors across 72 to 87 bits from three 24 to 29 bit primes.
	const lpBound3 = uint64(1) << 30
	for i := 0; i < 30; i++ {
		b1 := 24 + rnd.Intn(6)
		b2 := 24 + rnd.Intn(6)
		b3 := 24 + rnd.Intn(6)
		p1 := testutil.RandomPrime(rnd, b1)
		p2 := testutil.RandomPrime(rnd, b2)
		p3 := testutil.RandomPrime(rnd, b3)
		prod := new(big.Int).Mul(p1, p2)
		prod.Mul(prod, p3)
		var ns normState
		initNorm(&ns, prod)
		got, ok := splitNormCofactor(&ns, 500000, lpBound3, 3, 92)
		if !ok || got[2] == 0 {
			t.Fatalf("3LP %d-bit split(%v * %v * %v): got %v ok=%v", prod.BitLen(), p1, p2, p3, got, ok)
		}
		chk := new(big.Int).SetUint64(got[0])
		chk.Mul(chk, new(big.Int).SetUint64(got[1]))
		chk.Mul(chk, new(big.Int).SetUint64(got[2]))
		if chk.Cmp(prod) != 0 {
			t.Fatalf("3LP product mismatch: got %v, want %v", chk, prod)
		}
	}

	// Reject an 80 bit prime and a 4LP composite when maxLP == 3.
	bigPrime80 := testutil.RandomPrime(rnd, 80)
	var nsP normState
	initNorm(&nsP, bigPrime80)
	if _, ok := splitNormCofactor(&nsP, 500000, lpBound3, 3, 92); ok {
		t.Fatalf("accepted 80-bit inert prime %v", bigPrime80)
	}
	fourLP := big.NewInt(1)
	for k := 0; k < 4; k++ {
		fourLP.Mul(fourLP, testutil.RandomPrime(rnd, 21))
	}
	var ns4 normState
	initNorm(&ns4, fourLP)
	if _, ok := splitNormCofactor(&ns4, 500000, lpBound3, 3, 92); ok {
		t.Fatalf("accepted 4LP composite %v with maxLP=3", fourLP)
	}

	// Verify end to end 3LP GNFS factorization.
	n := new(big.Int).Mul(testutil.RandomPrime(rnd, 72), testutil.RandomPrime(rnd, 72))
	params := &Params{
		Degree: 4, Br: 2500, Ba: 3500,
		LpBoundRat: 1 << 18, LpBoundAlg: 1 << 18,
		MfbRat: 50, MfbAlg: 52, MaxLPRat: 3, MaxLPAlg: 3,
		UseLattice: true, QMin: 3500, QMax: 50000, IHalf: 512, JMax: 512,
		SlackRat: 50, SlackAlg: 52, NumQCs: 64,
		RelTargetRatio: 0.85,
	}
	f1, f2 := splitWithParams(n, params)
	if f1 == nil || new(big.Int).Mul(f1, f2).Cmp(n) != 0 {
		t.Fatalf("3LP end-to-end splitWithParams(%v) = %v, %v", n, f1, f2)
	}
}

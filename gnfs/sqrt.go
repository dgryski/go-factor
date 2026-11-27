package gnfs

import (
	"math"
	"math/big"
	"math/bits"
	"math/cmplx"
	"math/rand"
	"runtime"
	"sync"

	"github.com/dgryski/go-factor/internal/arith"
)

// tryDependency attempts to extract a non-trivial factor of n from the
// dependency bitset dep over surviving relations. It returns (factor, sqOk)
// where sqOk reports whether beta^2 == gamma succeeded.
func tryDependency(n *big.Int, poly *Poly, fb *FactorBases, rels []Relation, dep []uint64) (factor *big.Int, sqOk bool) {
	var selected []Relation
	sievedCount := 0
	for w, word := range dep {
		for word != 0 {
			i := w*64 + bits.TrailingZeros64(word)
			word &= word - 1
			rel := rels[i]
			selected = append(selected, rel)
			if !rel.IsFree {
				sievedCount++
			}
		}
	}
	if len(selected) == 0 || sievedCount%2 != 0 {
		return nil, false
	}

	fHat, fHatDeriv := poly.Monic()
	d := poly.Degree
	cd := poly.Coeffs[d]
	mHat := new(big.Int).Mul(cd, poly.M)
	mHat.Mod(mHat, n)

	// Check rational sign and build rational square root R = prod p^(e_p/2) mod N.
	ratExps := make([]int, len(fb.RFB))
	ratLPExps := make(map[uint64]int)
	negRat := 0
	gSc := newCofactorScratch(poly.Degree)
	for _, rel := range selected {
		if !rel.IsFree && evalGInto(&gSc.val, poly, rel.A, rel.B, gSc).Sign() < 0 {
			negRat++
		}
		for _, idx := range rel.RatCols {
			ratExps[idx]++
		}
		for _, lp := range rel.RatLPs {
			if lp != 0 {
				ratLPExps[lp]++
			}
		}
	}
	if negRat%2 != 0 {
		return nil, false
	}

	R := big.NewInt(1)
	var pBig, eBig, term big.Int
	for idx, exp := range ratExps {
		if exp == 0 {
			continue
		}
		if exp%2 != 0 {
			return nil, false
		}
		pBig.SetUint64(uint64(fb.RFB[idx].P))
		eBig.SetInt64(int64(exp / 2))
		term.Exp(&pBig, &eBig, n)
		R.Mul(R, &term)
		R.Mod(R, n)
	}
	for lp, exp := range ratLPExps {
		if exp%2 != 0 {
			return nil, false
		}
		pBig.SetUint64(lp)
		eBig.SetInt64(int64(exp / 2))
		term.Exp(&pBig, &eBig, n)
		R.Mul(R, &term)
		R.Mod(R, n)
	}

	// Evaluate fHat'(mHat) mod N.
	fDerivAtM := big.NewInt(0)
	for i := len(fHatDeriv) - 1; i >= 0; i-- {
		fDerivAtM.Mul(fDerivAtM, mHat)
		fDerivAtM.Add(fDerivAtM, fHatDeriv[i])
		fDerivAtM.Mod(fDerivAtM, n)
	}

	halfS := big.NewInt(int64(sievedCount / 2))
	cdPow := new(big.Int).Exp(cd, halfS, n)
	xRat := new(big.Int).Mul(fDerivAtM, cdPow)
	xRat.Mod(xRat, n)
	xRat.Mul(xRat, R)
	xRat.Mod(xRat, n)

	// Algebraic side: compute gamma = fHatDeriv^2 * prod(c_d*A - B*alphaHat) in Z[y]/(fHat).
	gamma := computeGamma(selected, cd, fHat, fHatDeriv)

	// Bound |beta_i| via complex embeddings and inverse Vandermonde.
	boundBits := boundBetaBits(selected, cd, fHat, fHatDeriv)

	// Take square root beta in Z[alphaHat] via p-adic Hensel lifting.
	beta, ok := algebraicSqrt(gamma, fHat, cd, boundBits)
	if !ok {
		return nil, false
	}

	// Evaluate X_alg = beta(mHat) * A1^(sievedCount/2) mod N.
	betaAtM := big.NewInt(0)
	for i := len(beta) - 1; i >= 0; i-- {
		betaAtM.Mul(betaAtM, mHat)
		betaAtM.Add(betaAtM, beta[i])
		betaAtM.Mod(betaAtM, n)
	}
	a1Pow := new(big.Int).Exp(poly.A1, halfS, n)
	xAlg := new(big.Int).Mul(betaAtM, a1Pow)
	xAlg.Mod(xAlg, n)

	// Verify X_alg^2 ≡ X_rat^2 (mod N).
	lhs := new(big.Int).Mul(xAlg, xAlg)
	lhs.Mod(lhs, n)
	rhs := new(big.Int).Mul(xRat, xRat)
	rhs.Mod(rhs, n)
	if lhs.Cmp(rhs) != 0 {
		logf("gnfs: WARNING: X_alg^2 != X_rat^2 (mod N) despite beta^2 == gamma")
		return nil, true
	}

	var diff, g big.Int
	one := big.NewInt(1)
	for _, op := range []func(x, y *big.Int) *big.Int{diff.Sub, diff.Add} {
		op(xAlg, xRat)
		g.GCD(nil, nil, diff.Abs(&diff), n)
		if g.Cmp(one) > 0 && g.Cmp(n) < 0 {
			return new(big.Int).Set(&g), true
		}
	}
	return nil, true
}

// computeGamma returns fHat'(alpha_hat)^2 * prod_S (c_d*A - B*alpha_hat) mod fHat(alpha_hat).
func computeGamma(selected []Relation, cd *big.Int, fHat, fHatDeriv []*big.Int) []*big.Int {
	d := len(fHat) - 1
	factors := make([][]*big.Int, len(selected))
	for i, rel := range selected {
		p := make([]*big.Int, d)
		for k := range p {
			p[k] = new(big.Int)
		}
		if rel.IsFree {
			p[0].SetInt64(rel.A)
		} else {
			p[0].Mul(cd, big.NewInt(rel.A))
			p[1].SetInt64(-rel.B)
		}
		factors[i] = p
	}

	var treeMul func(items [][]*big.Int, depth int) []*big.Int
	treeMul = func(items [][]*big.Int, depth int) []*big.Int {
		if len(items) == 1 {
			return items[0]
		}
		if len(items) == 2 {
			return polyAlgMulModPar(items[0], items[1], fHat, depth == 0)
		}
		mid := len(items) / 2
		var left, right []*big.Int
		if depth < 4 && len(items) >= 256 && runtime.GOMAXPROCS(0) > 1 {
			var wg sync.WaitGroup
			wg.Add(1)
			go func() {
				defer wg.Done()
				left = treeMul(items[:mid], depth+1)
			}()
			right = treeMul(items[mid:], depth+1)
			wg.Wait()
		} else {
			left = treeMul(items[:mid], depth+1)
			right = treeMul(items[mid:], depth+1)
		}
		return polyAlgMulModPar(left, right, fHat, depth == 0)
	}

	prod := treeMul(factors, 0)
	fDeriv := make([]*big.Int, d)
	for i := 0; i < d; i++ {
		if i < len(fHatDeriv) {
			fDeriv[i] = new(big.Int).Set(fHatDeriv[i])
		} else {
			fDeriv[i] = new(big.Int)
		}
	}
	fDerivSq := polyAlgMulMod(fDeriv, fDeriv, fHat)
	return polyAlgMulModPar(fDerivSq, prod, fHat, true)
}

func polyAlgMulMod(a, b, fHat []*big.Int) []*big.Int {
	return polyAlgMulModPar(a, b, fHat, false)
}

func polyAlgMulModPar(a, b, fHat []*big.Int, allowPar bool) []*big.Int {
	d := len(fHat) - 1
	prod := polyConvolveBig(a, b, d, allowPar)
	var tmp big.Int
	// Reduce mod monic fHat (fHat[d] == 1).
	for k := len(prod) - 1; k >= d; k-- {
		if prod[k].Sign() == 0 {
			continue
		}
		ck := prod[k]
		shift := k - d
		for j := 0; j < d; j++ {
			if fHat[j].Sign() != 0 {
				tmp.Mul(ck, fHat[j])
				prod[shift+j].Sub(prod[shift+j], &tmp)
			}
		}
	}
	return prod[:d]
}

// polyConvolveBig multiplies degree < d polynomials via 2d - 1 symmetric integer
// point evaluations and Newton forward differences when coefficients >= 4096 bits.
func polyConvolveBig(a, b []*big.Int, d int, allowPar bool) []*big.Int {
	maxBits := 0
	for i := 0; i < d; i++ {
		if bl := a[i].BitLen(); bl > maxBits {
			maxBits = bl
		}
		if bl := b[i].BitLen(); bl > maxBits {
			maxBits = bl
		}
	}
	if d < 3 || maxBits < 4096 {
		return polyConvolveSchoolbook(a, b, d)
	}

	mDeg := 2*d - 2
	nPts := mDeg + 1
	same := len(a) > 0 && len(b) > 0 && a[0] == b[0]
	if same {
		for i := 1; i < d; i++ {
			if a[i] != b[i] {
				same = false
				break
			}
		}
	}

	evalAt := func(poly []*big.Int, x int64) *big.Int {
		if x == 0 {
			return new(big.Int).Set(poly[0])
		}
		res := new(big.Int).Set(poly[d-1])
		var xBig big.Int
		xBig.SetInt64(x)
		for i := d - 2; i >= 0; i-- {
			res.Mul(res, &xBig)
			res.Add(res, poly[i])
		}
		return res
	}

	y := make([]*big.Int, nPts)
	nw := Workers
	if nw <= 0 {
		nw = runtime.GOMAXPROCS(0)
	}
	if allowPar && maxBits >= 32768 && nw > 1 {
		var wg sync.WaitGroup
		for k := 0; k < nPts; k++ {
			wg.Add(1)
			go func(idx int) {
				defer wg.Done()
				xk := int64(idx - (d - 1))
				va := evalAt(a, xk)
				if same {
					y[idx] = va.Mul(va, va)
				} else {
					vb := evalAt(b, xk)
					y[idx] = va.Mul(va, vb)
				}
			}(k)
		}
		wg.Wait()
	} else {
		for k := 0; k < nPts; k++ {
			xk := int64(k - (d - 1))
			va := evalAt(a, xk)
			if same {
				y[k] = va.Mul(va, va)
			} else {
				vb := evalAt(b, xk)
				y[k] = va.Mul(va, vb)
			}
		}
	}

	var stepBig, rem big.Int
	for step := 1; step <= mDeg; step++ {
		stepBig.SetInt64(int64(step))
		for k := mDeg; k >= step; k-- {
			y[k].Sub(y[k], y[k-1])
			y[k].QuoRem(y[k], &stepBig, &rem)
		}
	}

	prod := make([]*big.Int, nPts)
	for i := range prod {
		prod[i] = new(big.Int)
	}
	prod[0].Set(y[mDeg])
	var xjBig, term big.Int
	for j := mDeg - 1; j >= 0; j-- {
		xjBig.SetInt64(int64(j - (d - 1)))
		for k := mDeg - j; k >= 1; k-- {
			term.Mul(prod[k], &xjBig)
			prod[k].Sub(prod[k-1], &term)
		}
		term.Mul(prod[0], &xjBig)
		prod[0].Sub(y[j], &term)
	}
	return prod
}

func polyConvolveSchoolbook(a, b []*big.Int, d int) []*big.Int {
	prod := make([]*big.Int, 2*d-1)
	for i := range prod {
		prod[i] = new(big.Int)
	}
	var tmp big.Int
	for i := 0; i < d; i++ {
		if a[i].Sign() == 0 {
			continue
		}
		for j := 0; j < d; j++ {
			if b[j].Sign() == 0 {
				continue
			}
			tmp.Mul(a[i], b[j])
			prod[i+j].Add(prod[i+j], &tmp)
		}
	}
	return prod
}

// boundBetaBits computes an upper bound in bits for max_i |beta_i|.
func boundBetaBits(selected []Relation, cd *big.Int, fHat, fHatDeriv []*big.Int) int {
	d := len(fHat) - 1
	// Scale fHat by S = 2^scaleShift so roots z_k = alpha_k / S have |z_k| <= 2.
	scaleShift := 0
	for i := 0; i < d; i++ {
		if fHat[i].Sign() == 0 {
			continue
		}
		bl := fHat[i].BitLen()
		s := (bl + (d - i) - 1) / (d - i)
		if s > scaleShift {
			scaleShift = s
		}
	}

	gCoeffs := make([]complex128, d+1)
	for i := 0; i < d; i++ {
		gCoeffs[i] = complex(scaledDownFloat(fHat[i], (d-i)*scaleShift), 0)
	}
	gCoeffs[d] = 1

	zRoots := aberthEhrlich(gCoeffs)
	vInvZ := invertVandermonde(zRoots)
	S := math.Exp2(float64(scaleShift))
	lnS := float64(scaleShift) * math.Ln2

	// Evaluate log |sigma_k(beta)| for each complex root alpha_k = S * z_k.
	cdF, _ := new(big.Float).SetInt(cd).Float64()
	logSigma := make([]float64, d)
	for k, zk := range zRoots {
		// fHat'(S * z_k) = S^(d-1) * g'(z_k)
		gp := complex(0, 0)
		for i := d; i >= 1; i-- {
			gp = gp*zk + complex(float64(i), 0)*gCoeffs[i]
		}
		lnDeriv := float64(d-1)*lnS + math.Log(cmplx.Abs(gp))

		sumLog := 0.0
		for _, rel := range selected {
			if rel.IsFree {
				sumLog += math.Log(math.Abs(float64(rel.A)))
			} else {
				// c_d*A - B*S*z_k = S * ((c_d*A)/S - B*z_k)
				cdaOverS := complex(cdF*float64(rel.A)/S, 0)
				val := cdaOverS - complex(float64(rel.B), 0)*zk
				sumLog += lnS + math.Log(cmplx.Abs(val))
			}
		}
		logSigma[k] = lnDeriv + 0.5*sumLog
	}

	maxLn := math.Inf(-1)
	for i := 0; i < d; i++ {
		for k := 0; k < d; k++ {
			vAbs := cmplx.Abs(vInvZ[i][k])
			if vAbs == 0 {
				continue
			}
			// |(V^-1)_{ik}| = S^(-i) * |(V_z^-1)_{ik}|
			cand := math.Log(vAbs) - float64(i)*lnS + logSigma[k]
			if cand > maxLn {
				maxLn = cand
			}
		}
	}
	if math.IsInf(maxLn, -1) {
		return 64
	}
	bits := int(math.Ceil((maxLn + math.Log(float64(d))) / math.Ln2))
	if bits < 64 {
		bits = 64
	}
	return bits + 64 // 64-bit safety margin
}

func scaledDownFloat(x *big.Int, shift int) float64 {
	if x.Sign() == 0 {
		return 0
	}
	neg := x.Sign() < 0
	ax := new(big.Int).Abs(x)
	// Keep top ~60 bits for float64 conversion.
	bl := ax.BitLen()
	extra := bl - 60
	if extra < 0 {
		extra = 0
	}
	top := new(big.Int).Rsh(ax, uint(extra))
	f := float64(top.Uint64()) * math.Exp2(float64(extra-shift))
	if neg {
		return -f
	}
	return f
}

// aberthEhrlich finds all d complex roots of monic polynomial coeffs[0..d].
func aberthEhrlich(coeffs []complex128) []complex128 {
	d := len(coeffs) - 1
	roots := make([]complex128, d)
	// Initial guesses on a circle of radius 1.5 with asymmetric phase.
	for k := 0; k < d; k++ {
		angle := 2*math.Pi*float64(k)/float64(d) + 0.1234567
		roots[k] = complex(1.5*math.Cos(angle), 1.5*math.Sin(angle))
	}

	eval := func(z complex128) (p, dp complex128) {
		p = coeffs[d]
		dp = 0
		for i := d - 1; i >= 0; i-- {
			dp = dp*z + p
			p = p*z + coeffs[i]
		}
		return p, dp
	}

	for iter := 0; iter < 60; iter++ {
		maxStep := 0.0
		for k := 0; k < d; k++ {
			p, dp := eval(roots[k])
			if cmplx.Abs(p) == 0 {
				continue
			}
			if cmplx.Abs(dp) == 0 {
				dp = 1e-18
			}
			w := p / dp
			sum := complex(0, 0)
			for j := 0; j < d; j++ {
				if j != k {
					diff := roots[k] - roots[j]
					if cmplx.Abs(diff) > 0 {
						sum += 1 / diff
					}
				}
			}
			denom := 1 - w*sum
			if cmplx.Abs(denom) == 0 {
				continue
			}
			step := w / denom
			roots[k] -= step
			if a := cmplx.Abs(step); a > maxStep {
				maxStep = a
			}
		}
		if maxStep < 1e-13 {
			break
		}
	}
	return roots
}

func invertVandermonde(roots []complex128) [][]complex128 {
	d := len(roots)
	mat := make([][]complex128, d)
	inv := make([][]complex128, d)
	for i := 0; i < d; i++ {
		mat[i] = make([]complex128, d)
		inv[i] = make([]complex128, d)
		pow := complex(1, 0)
		for j := 0; j < d; j++ {
			mat[i][j] = pow
			pow *= roots[i]
		}
		inv[i][i] = 1
	}

	for col := 0; col < d; col++ {
		piv := col
		for r := col + 1; r < d; r++ {
			if cmplx.Abs(mat[r][col]) > cmplx.Abs(mat[piv][col]) {
				piv = r
			}
		}
		mat[col], mat[piv] = mat[piv], mat[col]
		inv[col], inv[piv] = inv[piv], inv[col]

		div := mat[col][col]
		for j := 0; j < d; j++ {
			mat[col][j] /= div
			inv[col][j] /= div
		}
		for r := 0; r < d; r++ {
			if r == col {
				continue
			}
			factor := mat[r][col]
			for j := 0; j < d; j++ {
				mat[r][j] -= factor * mat[col][j]
				inv[r][j] -= factor * inv[col][j]
			}
		}
	}
	return inv
}

// algebraicSqrt finds beta in Z[y]/(fHat) with beta^2 == gamma via p-adic Hensel lifting.
func algebraicSqrt(gamma, fHat []*big.Int, cd *big.Int, boundBits int) ([]*big.Int, bool) {
	d := len(fHat) - 1
	fHatP := make([]uint64, d+1)
	gammaP := make([]uint64, d)

	var bigP, rem big.Int
	// Look for an inert prime first (around 2^30); remember a fallback split prime with fewest factors.
	type splitCand struct {
		p    uint64
		facs [][]uint64
	}
	var bestSplit *splitCand

	p := uint64((1 << 30) + 3)
	checked := 0
	var chosenP uint64
	var chosenFacs [][]uint64

	for checked < 80 {
		if arith.IsPrime64(p) {
			bigP.SetUint64(p)
			if rem.Mod(cd, &bigP).Sign() != 0 {
				for i := 0; i <= d; i++ {
					fHatP[i] = rem.Mod(fHat[i], &bigP).Uint64()
				}
				df := arith.PolyDerivModP(fHatP, p)
				if len(arith.PolyGcdModP(fHatP, df, p)) == 1 {
					for i := 0; i < d; i++ {
						gammaP[i] = rem.Mod(gamma[i], &bigP).Uint64()
					}
					if len(arith.PolyGcdModP(fHatP, arith.PolyTrim(gammaP), p)) == 1 {
						checked++
						if arith.PolyIsIrreducibleModP(fHatP, p) {
							chosenP = p
							chosenFacs = [][]uint64{append([]uint64(nil), fHatP...)}
							break
						}
						if bestSplit == nil || len(bestSplit.facs) > 2 {
							facs := arith.PolyFactorModP(fHatP, p)
							if bestSplit == nil || len(facs) < len(bestSplit.facs) {
								bestSplit = &splitCand{p: p, facs: facs}
							}
						}
					}
				}
			}
		}
		p += 2
	}

	if chosenP == 0 {
		if bestSplit == nil {
			return nil, false
		}
		chosenP = bestSplit.p
		chosenFacs = bestSplit.facs
	}

	p = chosenP
	bigP.SetUint64(p)
	for i := 0; i <= d; i++ {
		fHatP[i] = rem.Mod(fHat[i], &bigP).Uint64()
	}
	for i := 0; i < d; i++ {
		gammaP[i] = rem.Mod(gamma[i], &bigP).Uint64()
	}

	// Compute inverse square root s_j ≡ gamma^(-1/2) mod (h_j, p) for each irreducible factor h_j.
	k := len(chosenFacs)
	invSqrts := make([][]uint64, k)
	for j, hj := range chosenFacs {
		gj := arith.PolyModModP(gammaP, hj, p)
		sj, ok := polyInvSqrtModIrred(gj, hj, p)
		if !ok {
			return nil, false // gamma is not a square mod (hj, p)
		}
		invSqrts[j] = sj
	}

	// Precompute CRT idempotents E_j ≡ 1 (mod h_j), ≡ 0 (mod h_i for i != j) mod fHatP.
	crtBasis := make([][]uint64, k)
	if k == 1 {
		crtBasis[0] = []uint64{1}
	} else {
		for j, hj := range chosenFacs {
			other, _ := arith.PolyDivRemModP(fHatP, hj, p)
			invOther, ok := arith.PolyInvModP(other, hj, p)
			if !ok {
				return nil, false
			}
			crtBasis[j] = arith.PolyModModP(arith.PolyMulModP(other, invOther, p), fHatP, p)
		}
	}

	// Build doubling chain of exponents 1 = e_0 < e_1 < ... < e_m = K such that p^K > 2*M_beta.
	log2P := math.Log2(float64(p))
	K := int(math.Ceil(float64(boundBits+2) / log2P))
	if K < 2 {
		K = 2
	}
	var exps []int
	for cur := K; cur > 1; cur = (cur + 1) / 2 {
		exps = append(exps, cur)
	}
	// Reverse exps so it runs from small to K.
	for i, j := 0, len(exps)-1; i < j; i, j = i+1, j-1 {
		exps[i], exps[j] = exps[j], exps[i]
	}

	// Precompute moduli p^{e_i} and reduced gamma mod p^{e_i}.
	moduli := make([]*big.Int, len(exps))
	gammaMods := make([][]*big.Int, len(exps))
	invTwos := make([]*big.Int, len(exps))
	for i, e := range exps {
		m := new(big.Int).Exp(&bigP, big.NewInt(int64(e)), nil)
		moduli[i] = m
		invTwos[i] = new(big.Int).Add(m, bigOne)
		invTwos[i].Rsh(invTwos[i], 1)
		gm := make([]*big.Int, d)
		for j := 0; j < d; j++ {
			gm[j] = new(big.Int).Mod(gamma[j], m)
		}
		gammaMods[i] = gm
	}
	pK := moduli[len(moduli)-1]
	halfPK := new(big.Int).Rsh(pK, 1)

	// Try the 2^(k-1) sign choices for invSqrts[1..k-1].
	for mask := 0; mask < (1 << uint(k-1)); mask++ {
		var x0P []uint64
		if k == 1 {
			x0P = invSqrts[0]
		} else {
			x0P = nil
			for j := 0; j < k; j++ {
				sj := invSqrts[j]
				if j > 0 && mask&(1<<uint(j-1)) != 0 {
					sj = arith.PolySubModP(nil, sj, p)
				}
				term := arith.PolyModModP(arith.PolyMulModP(sj, crtBasis[j], p), fHatP, p)
				x0P = arith.PolyAddModP(x0P, term, p)
			}
		}

		x := make([]*big.Int, d)
		for i := 0; i < d; i++ {
			if i < len(x0P) {
				x[i] = new(big.Int).SetUint64(x0P[i])
			} else {
				x[i] = new(big.Int)
			}
		}

		var beta []*big.Int
		if len(moduli) >= 2 {
			// Lift x = gamma^(-1/2) to half precision, then finish with Karp and Markstein's
			// coupled step (ACM TOMS 23, 1997, Section 3) using only half precision multiplies.
			lastHalf := len(moduli) - 2
			for step := 0; step <= lastHalf; step++ {
				m := moduli[step]
				gm := gammaMods[step]
				x2 := polyAlgMulModModBigPar(x, x, fHat, m, true)
				gx2 := polyAlgMulModModBigPar(gm, x2, fHat, m, true)
				threeMinus := make([]*big.Int, d)
				threeMinus[0] = new(big.Int).Sub(big.NewInt(3), gx2[0])
				threeMinus[0].Mod(threeMinus[0], m)
				for i := 1; i < d; i++ {
					threeMinus[i] = new(big.Int).Neg(gx2[i])
					threeMinus[i].Mod(threeMinus[i], m)
				}
				prod := polyAlgMulModModBigPar(x, threeMinus, fHat, m, true)
				inv2 := invTwos[step]
				for i := 0; i < d; i++ {
					prod[i].Mul(prod[i], inv2)
					prod[i].Mod(prod[i], m)
				}
				x = prod
			}
			m0 := moduli[lastHalf]
			eRem := exps[len(exps)-1] - exps[lastHalf]
			m1 := new(big.Int).Exp(&bigP, big.NewInt(int64(eRem)), nil)
			inv2M1 := new(big.Int).Add(m1, bigOne)
			inv2M1.Rsh(inv2M1, 1)

			b0 := polyAlgMulModModBigPar(gammaMods[lastHalf], x, fHat, m0, true)
			b0Sq := polyAlgMulModModBigPar(b0, b0, fHat, pK, true)
			qPoly := make([]*big.Int, d)
			xLow := make([]*big.Int, d)
			for i := 0; i < d; i++ {
				diff := new(big.Int).Sub(gammaMods[len(moduli)-1][i], b0Sq[i])
				diff.Mod(diff, pK)
				qPoly[i] = new(big.Int).Quo(diff, m0)
				xLow[i] = new(big.Int).Mod(x[i], m1)
			}
			hq := polyAlgMulModModBigPar(xLow, qPoly, fHat, m1, true)
			beta = make([]*big.Int, d)
			for i := 0; i < d; i++ {
				hq[i].Mul(hq[i], inv2M1)
				hq[i].Mod(hq[i], m1)
				beta[i] = new(big.Int).Mul(m0, hq[i])
				beta[i].Add(beta[i], b0[i])
			}
		} else {
			m := moduli[0]
			gm := gammaMods[0]
			x2 := polyAlgMulModModBigPar(x, x, fHat, m, true)
			gx2 := polyAlgMulModModBigPar(gm, x2, fHat, m, true)
			threeMinus := make([]*big.Int, d)
			threeMinus[0] = new(big.Int).Sub(big.NewInt(3), gx2[0])
			threeMinus[0].Mod(threeMinus[0], m)
			for i := 1; i < d; i++ {
				threeMinus[i] = new(big.Int).Neg(gx2[i])
				threeMinus[i].Mod(threeMinus[i], m)
			}
			prod := polyAlgMulModModBigPar(x, threeMinus, fHat, m, true)
			inv2 := invTwos[0]
			for i := 0; i < d; i++ {
				prod[i].Mul(prod[i], inv2)
				prod[i].Mod(prod[i], m)
			}
			beta = polyAlgMulModModBigPar(gm, prod, fHat, pK, true)
		}

		tooBig := false
		for i := 0; i < d; i++ {
			if beta[i].Cmp(halfPK) > 0 {
				beta[i].Sub(beta[i], pK)
			}
			if beta[i].BitLen() > boundBits+1 {
				tooBig = true
			}
		}
		if tooBig && k > 1 {
			continue
		}

		// Verify beta^2 == gamma in Z[y]/(fHat).
		sq := polyAlgMulModPar(beta, beta, fHat, true)
		equal := true
		for i := 0; i < d; i++ {
			if sq[i].Cmp(gamma[i]) != 0 {
				equal = false
				break
			}
		}
		if equal {
			return beta, true
		}
	}
	return nil, false
}

func polyAlgMulModModBig(a, b, fHat []*big.Int, mod *big.Int) []*big.Int {
	return polyAlgMulModModBigPar(a, b, fHat, mod, false)
}

func polyAlgMulModModBigPar(a, b, fHat []*big.Int, mod *big.Int, allowPar bool) []*big.Int {
	d := len(fHat) - 1
	prod := polyConvolveBig(a, b, d, allowPar)
	var tmp big.Int
	for k := len(prod) - 1; k >= d; k-- {
		prod[k].Mod(prod[k], mod)
		if prod[k].Sign() == 0 {
			continue
		}
		ck := prod[k]
		shift := k - d
		for j := 0; j < d; j++ {
			if fHat[j].Sign() != 0 {
				tmp.Mul(ck, fHat[j])
				prod[shift+j].Sub(prod[shift+j], &tmp)
			}
		}
	}
	for i := 0; i < d; i++ {
		prod[i].Mod(prod[i], mod)
	}
	return prod[:d]
}

// polyInvSqrtModIrred computes g^(-1/2) mod (h, p) in F_{p^deg(h)},
// returning false if g is not a square in F_{p^deg(h)}.
func polyInvSqrtModIrred(g, h []uint64, p uint64) ([]uint64, bool) {
	deg := len(h) - 1
	g = arith.PolyModModP(g, h, p)
	if len(g) == 0 {
		return nil, false
	}
	if deg == 1 {
		// Scalar in F_p.
		if arith.PowMod64(g[0], (p-1)>>1, p) != 1 {
			return nil, false
		}
		r := arith.SqrtModPrime64(g[0], p)
		return []uint64{arith.InvMod(r, p)}, true
	}

	// q = p^deg
	q := new(big.Int).Exp(new(big.Int).SetUint64(p), big.NewInt(int64(deg)), nil)
	qMinus1Half := new(big.Int).Sub(q, bigOne)
	qMinus1Half.Rsh(qMinus1Half, 1)

	// Euler criterion.
	if !arith.PolyEqual(arith.PolyPowBigMod(g, qMinus1Half, h, p), []uint64{1}) {
		return nil, false
	}

	// Fast path when q ≡ 3 (mod 4): r = g^((q+1)/4), r^-1 = g^((q-3)/4).
	if new(big.Int).And(q, big.NewInt(3)).Uint64() == 3 {
		exp := new(big.Int).Sub(q, big.NewInt(3))
		exp.Rsh(exp, 2)
		return arith.PolyPowBigMod(g, exp, h, p), true
	}

	// General Tonelli-Shanks in F_q.
	qOdd := new(big.Int).Sub(q, bigOne)
	s := 0
	for qOdd.Bit(0) == 0 {
		qOdd.Rsh(qOdd, 1)
		s++
	}
	rng := rand.New(rand.NewSource(1))
	var z []uint64
	for {
		z = make([]uint64, deg)
		for i := range z {
			z[i] = uint64(rng.Int63n(int64(p)))
		}
		z = arith.PolyTrim(z)
		if len(z) == 0 {
			continue
		}
		e := arith.PolyPowBigMod(z, qMinus1Half, h, p)
		if arith.PolyEqual(e, []uint64{p - 1}) {
			break
		}
	}

	m := s
	c := arith.PolyPowBigMod(z, qOdd, h, p)
	t := arith.PolyPowBigMod(g, qOdd, h, p)
	qOddPlus1Half := new(big.Int).Add(qOdd, bigOne)
	qOddPlus1Half.Rsh(qOddPlus1Half, 1)
	r := arith.PolyPowBigMod(g, qOddPlus1Half, h, p)

	for !arith.PolyEqual(t, []uint64{1}) {
		i := 1
		tt := arith.PolyModModP(arith.PolyMulModP(t, t, p), h, p)
		for !arith.PolyEqual(tt, []uint64{1}) {
			i++
			if i == m {
				return nil, false
			}
			tt = arith.PolyModModP(arith.PolyMulModP(tt, tt, p), h, p)
		}
		b := c
		for j := 0; j < m-i-1; j++ {
			b = arith.PolyModModP(arith.PolyMulModP(b, b, p), h, p)
		}
		m = i
		c = arith.PolyModModP(arith.PolyMulModP(b, b, p), h, p)
		t = arith.PolyModModP(arith.PolyMulModP(t, c, p), h, p)
		r = arith.PolyModModP(arith.PolyMulModP(r, b, p), h, p)
	}
	return arith.PolyInvModP(r, h, p)
}

package gnfs

import (
	"math/big"

	"github.com/dgryski/go-factor/internal/arith"
)

// BadPrimeIdeal represents an exact prime ideal P_j dividing a bad prime p in O_p.
type BadPrimeIdeal struct {
	Col    int32
	E, F   int
	mulMat [][]*big.Int
	mulMod [][]uint64
}

// BadPrimeDecomp holds the p maximal overorder O_p basis and prime ideals above p
// computed via Zassenhaus Round 2 (Cohen, GTM 138, 1993, Algorithms 6.1.8, 6.2.9, 4.8.17).
type BadPrimeDecomp struct {
	P         uint32
	DividesCD bool
	cdA0      []*big.Int
	a1Row     []*big.Int
	Ideals    []BadPrimeIdeal
}

type pOrder struct {
	d      int
	denom  *big.Int
	wMat   [][]*big.Int // omega_j = (1/denom) * sum_i wMat[j][i] * alphaHat^i
	aMat   [][]*big.Int // alphaHat^i = sum_k aMat[i][k] * omega_k
	mulTab [][][]*big.Int
}

func decomposeBadPrime(poly *Poly, fHat []*big.Int, p32 uint32, dividesCD bool, nextCol *int32) BadPrimeDecomp {
	d := poly.Degree
	up := uint64(p32)
	ord := computePMaximalOrder(fHat, up)
	ideals := decomposeOrderModP(ord, up, nextCol)

	cd := poly.Coeffs[d]
	cdA0 := make([]*big.Int, d)
	a1Row := make([]*big.Int, d)
	for k := 0; k < d; k++ {
		cdA0[k] = new(big.Int).Mul(cd, ord.aMat[0][k])
		a1Row[k] = new(big.Int).Set(ord.aMat[1][k])
	}
	return BadPrimeDecomp{
		P:         p32,
		DividesCD: dividesCD,
		cdA0:      cdA0,
		a1Row:     a1Row,
		Ideals:    ideals,
	}
}

// valuations returns exact v_{P_j} valuations of c_d*a - b*alphaHat above P.
func (bd *BadPrimeDecomp) valuations(a, b int64, dst []int) []int {
	d := len(bd.cdA0)
	p := uint64(bd.P)
	var bigP big.Int
	bigP.SetUint64(p)

	g := make([]*big.Int, d)
	var aBig, bBig, term big.Int
	aBig.SetInt64(a)
	bBig.SetInt64(b)
	allZero := true
	for k := 0; k < d; k++ {
		v := new(big.Int).Mul(bd.cdA0[k], &aBig)
		term.Mul(bd.a1Row[k], &bBig)
		v.Sub(v, &term)
		if v.Sign() != 0 {
			allZero = false
		}
		g[k] = v
	}
	dst = dst[:0]
	if allZero {
		for range bd.Ideals {
			dst = append(dst, 0)
		}
		return dst
	}

	var pPow int
	var qTmp, rem big.Int
	for {
		divAll := true
		for k := 0; k < d; k++ {
			if rem.Mod(g[k], &bigP).Sign() != 0 {
				divAll = false
				break
			}
		}
		if !divAll {
			break
		}
		for k := 0; k < d; k++ {
			g[k].Quo(g[k], &bigP)
		}
		pPow++
	}

	w := make([]*big.Int, d)
	nextW := make([]*big.Int, d)
	wMod := make([]uint64, d)
	for i := 0; i < d; i++ {
		w[i] = new(big.Int)
		nextW[i] = new(big.Int)
	}

	for idx := range bd.Ideals {
		id := &bd.Ideals[idx]
		val := pPow * id.E
		for k := 0; k < d; k++ {
			w[k].Set(g[k])
			wMod[k] = rem.Mod(w[k], &bigP).Uint64()
		}
		for {
			zeroModP := true
			for row := 0; row < d; row++ {
				var sum uint64
				mRow := id.mulMod[row]
				for col := 0; col < d; col++ {
					sum += arith.MulMod64(mRow[col], wMod[col], p)
					if sum >= p {
						sum -= p
					}
				}
				if sum != 0 {
					zeroModP = false
					break
				}
			}
			if !zeroModP {
				break
			}
			for row := 0; row < d; row++ {
				nextW[row].SetInt64(0)
				for col := 0; col < d; col++ {
					if id.mulMat[row][col].Sign() == 0 || w[col].Sign() == 0 {
						continue
					}
					qTmp.Mul(id.mulMat[row][col], w[col])
					nextW[row].Add(nextW[row], &qTmp)
				}
				nextW[row].Quo(nextW[row], &bigP)
			}
			for k := 0; k < d; k++ {
				w[k].Set(nextW[k])
				wMod[k] = rem.Mod(w[k], &bigP).Uint64()
			}
			val++
		}
		dst = append(dst, val)
	}
	return dst
}

func (bd *BadPrimeDecomp) appendOddCols(a, b int64, dst []int32) []int32 {
	var vBuf [8]int
	vals := bd.valuations(a, b, vBuf[:0])
	for j, v := range vals {
		if v&1 != 0 {
			dst = append(dst, bd.Ideals[j].Col)
		}
	}
	return dst
}

func computePMaximalOrder(fHat []*big.Int, p uint64) *pOrder {
	d := len(fHat) - 1
	wMat := make([][]*big.Int, d)
	for i := 0; i < d; i++ {
		wMat[i] = make([]*big.Int, d)
		for j := 0; j < d; j++ {
			if i == j {
				wMat[i][j] = big.NewInt(1)
			} else {
				wMat[i][j] = new(big.Int)
			}
		}
	}
	ord := buildPOrder(fHat, big.NewInt(1), wMat)
	bigP := new(big.Int).SetUint64(p)

	for {
		mMod := mulTableModP(ord, p)
		radBasis := nilradicalModP(mMod, d, p)
		r := len(radBasis)
		if r == 0 {
			return ord
		}
		perp := rightNullspaceModP(radBasis, d, p)

		// Enlarge O to the p radical multiplier ring (Cohen, GTM 138, 1993, Algorithm 6.1.8).
		var eqs1 [][]uint64
		for a := 0; a < r; a++ {
			va := radBasis[a]
			prodByB := make([][]uint64, r)
			for b := 0; b < r; b++ {
				prodByB[b] = algMulVecModP(radBasis[b], va, mMod, d, p)
			}
			for coord := 0; coord < d; coord++ {
				row := make([]uint64, r)
				for b := 0; b < r; b++ {
					row[b] = prodByB[b][coord]
				}
				eqs1 = append(eqs1, row)
			}
		}
		lam1 := rightNullspaceModP(eqs1, r, p)
		if len(lam1) == 0 {
			return ord
		}

		s1 := len(lam1)
		zLifts := make([][]uint64, s1)
		for t := 0; t < s1; t++ {
			zv := make([]uint64, d)
			for b := 0; b < r; b++ {
				if lam1[t][b] == 0 {
					continue
				}
				for k := 0; k < d; k++ {
					zv[k] = (zv[k] + arith.MulMod64(lam1[t][b], radBasis[b][k], p)) % p
				}
			}
			zLifts[t] = zv
		}

		var eqs2 [][]uint64
		var rem big.Int
		for a := 0; a < r; a++ {
			va := radBasis[a]
			qMod := make([][]uint64, s1)
			for t := 0; t < s1; t++ {
				exact := algMulExactInt(zLifts[t], va, ord)
				qm := make([]uint64, d)
				for k := 0; k < d; k++ {
					exact[k].Quo(exact[k], bigP)
					qm[k] = rem.Mod(exact[k], bigP).Uint64()
				}
				qMod[t] = qm
			}
			for _, pRow := range perp {
				row := make([]uint64, s1)
				for t := 0; t < s1; t++ {
					var dot uint64
					for k := 0; k < d; k++ {
						dot = (dot + arith.MulMod64(pRow[k], qMod[t][k], p)) % p
					}
					row[t] = dot
				}
				eqs2 = append(eqs2, row)
			}
		}
		lam2 := rightNullspaceModP(eqs2, s1, p)
		if len(lam2) == 0 {
			return ord
		}

		rows := make([][]*big.Int, 0, d+len(lam2))
		for i := 0; i < d; i++ {
			rw := make([]*big.Int, d)
			for j := 0; j < d; j++ {
				rw[j] = new(big.Int).Mul(ord.wMat[i][j], bigP)
			}
			rows = append(rows, rw)
		}
		for _, coeff := range lam2 {
			zv := make([]uint64, d)
			for t := 0; t < s1; t++ {
				if coeff[t] == 0 {
					continue
				}
				for k := 0; k < d; k++ {
					zv[k] = (zv[k] + arith.MulMod64(coeff[t], zLifts[t][k], p)) % p
				}
			}
			rw := make([]*big.Int, d)
			for j := 0; j < d; j++ {
				rw[j] = new(big.Int)
			}
			var zBig, term big.Int
			for k := 0; k < d; k++ {
				if zv[k] == 0 {
					continue
				}
				zBig.SetUint64(zv[k])
				for j := 0; j < d; j++ {
					if ord.wMat[k][j].Sign() != 0 {
						term.Mul(&zBig, ord.wMat[k][j])
						rw[j].Add(rw[j], &term)
					}
				}
			}
			rows = append(rows, rw)
		}

		newDenom := new(big.Int).Mul(ord.denom, bigP)
		newW := hnfUpperSquare(rows, d)
		for new(big.Int).Mod(newDenom, bigP).Sign() == 0 && matDivisibleBy(newW, bigP) {
			newDenom.Quo(newDenom, bigP)
			for i := 0; i < d; i++ {
				for j := 0; j < d; j++ {
					newW[i][j].Quo(newW[i][j], bigP)
				}
			}
		}
		ord = buildPOrder(fHat, newDenom, newW)
	}
}

func buildPOrder(fHat []*big.Int, denom *big.Int, wMat [][]*big.Int) *pOrder {
	d := len(fHat) - 1
	aMat := make([][]*big.Int, d)
	for i := 0; i < d; i++ {
		aMat[i] = make([]*big.Int, d)
		for j := 0; j < d; j++ {
			aMat[i][j] = new(big.Int)
		}
	}
	var prod big.Int
	for col := 0; col < d; col++ {
		for row := d - 1; row >= 0; row-- {
			rhs := new(big.Int)
			if row == col {
				rhs.Set(denom)
			}
			for k := row + 1; k < d; k++ {
				if wMat[row][k].Sign() != 0 && aMat[k][col].Sign() != 0 {
					prod.Mul(wMat[row][k], aMat[k][col])
					rhs.Sub(rhs, &prod)
				}
			}
			aMat[row][col].Quo(rhs, wMat[row][row])
		}
	}

	denomSq := new(big.Int).Mul(denom, denom)
	mulTab := make([][][]*big.Int, d)
	for j := 0; j < d; j++ {
		mulTab[j] = make([][]*big.Int, d)
		for k := 0; k < d; k++ {
			if k < j {
				mulTab[j][k] = mulTab[k][j]
				continue
			}
			pPoly := polyAlgMulMod(wMat[j], wMat[k], fHat)
			coords := make([]*big.Int, d)
			for l := 0; l < d; l++ {
				sum := new(big.Int)
				for i := 0; i < d; i++ {
					if pPoly[i].Sign() != 0 && aMat[i][l].Sign() != 0 {
						prod.Mul(pPoly[i], aMat[i][l])
						sum.Add(sum, &prod)
					}
				}
				sum.Quo(sum, denomSq)
				coords[l] = sum
			}
			mulTab[j][k] = coords
		}
	}

	return &pOrder{
		d:      d,
		denom:  denom,
		wMat:   wMat,
		aMat:   aMat,
		mulTab: mulTab,
	}
}

func decomposeOrderModP(ord *pOrder, p uint64, nextCol *int32) []BadPrimeIdeal {
	d := ord.d
	mMod := mulTableModP(ord, p)
	radBasis := nilradicalModP(mMod, d, p)

	radRREF, pivCols := rrefModP(radBasis, d, p)
	isPiv := make([]bool, d)
	for _, c := range pivCols {
		isPiv[c] = true
	}
	var compCols []int
	colToB := make([]int, d)
	for c := 0; c < d; c++ {
		if isPiv[c] {
			colToB[c] = -1
		} else {
			colToB[c] = len(compCols)
			compCols = append(compCols, c)
		}
	}
	sDim := len(compCols)

	projectToB := func(v []uint64) []uint64 {
		cp := append([]uint64(nil), v...)
		for row, pc := range pivCols {
			f := cp[pc]
			if f == 0 {
				continue
			}
			for c := pc; c < d; c++ {
				sub := arith.MulMod64(f, radRREF[row][c], p)
				if cp[c] >= sub {
					cp[c] -= sub
				} else {
					cp[c] = cp[c] + p - sub
				}
			}
		}
		out := make([]uint64, sDim)
		for bi, c := range compCols {
			out[bi] = cp[c]
		}
		return out
	}

	bMod := make([][][]uint64, sDim)
	for i := 0; i < sDim; i++ {
		bMod[i] = make([][]uint64, sDim)
		vi := make([]uint64, d)
		vi[compCols[i]] = 1
		for j := 0; j < sDim; j++ {
			vj := make([]uint64, d)
			vj[compCols[j]] = 1
			bMod[i][j] = projectToB(algMulVecModP(vi, vj, mMod, d, p))
		}
	}

	bigP := new(big.Int).SetUint64(p)
	var rem big.Int
	oneR := make([]uint64, d)
	for k := 0; k < d; k++ {
		oneR[k] = rem.Mod(ord.aMat[0][k], bigP).Uint64()
	}
	oneB := projectToB(oneR)

	fullBasisB := make([][]uint64, sDim)
	for i := 0; i < sDim; i++ {
		e := make([]uint64, sDim)
		e[i] = 1
		fullBasisB[i] = e
	}
	idempotents := splitEtaleIdempotents(fullBasisB, oneB, bMod, sDim, p)

	var out []BadPrimeIdeal
	for _, idem := range idempotents {
		var eqs [][]uint64
		prodCols := make([][]uint64, d)
		for col := 0; col < d; col++ {
			ej := make([]uint64, d)
			ej[col] = 1
			prodCols[col] = algMulVecModP(projectToB(ej), idem, bMod, sDim, p)
		}
		for coord := 0; coord < sDim; coord++ {
			row := make([]uint64, d)
			for col := 0; col < d; col++ {
				row[col] = prodCols[col][coord]
			}
			eqs = append(eqs, row)
		}
		pIdealBasis := rightNullspaceModP(eqs, d, p)
		fDeg := d - len(pIdealBasis)

		// Compute anti uniformizer beta_j annihilating P_j / p*O_p (Cohen, GTM 138, 1993, Alg. 4.8.17).
		var annEqs [][]uint64
		for _, y := range pIdealBasis {
			yCols := make([][]uint64, d)
			for col := 0; col < d; col++ {
				ec := make([]uint64, d)
				ec[col] = 1
				yCols[col] = algMulVecModP(ec, y, mMod, d, p)
			}
			for coord := 0; coord < d; coord++ {
				row := make([]uint64, d)
				for col := 0; col < d; col++ {
					row[col] = yCols[col][coord]
				}
				annEqs = append(annEqs, row)
			}
		}
		ann := rightNullspaceModP(annEqs, d, p)
		beta := ann[0]

		mulMat := make([][]*big.Int, d)
		mulMod := make([][]uint64, d)
		for row := 0; row < d; row++ {
			mulMat[row] = make([]*big.Int, d)
			mulMod[row] = make([]uint64, d)
			for col := 0; col < d; col++ {
				mulMat[row][col] = new(big.Int)
			}
		}
		var bBig, prod big.Int
		for rIdx := 0; rIdx < d; rIdx++ {
			if beta[rIdx] == 0 {
				continue
			}
			bBig.SetUint64(beta[rIdx])
			for col := 0; col < d; col++ {
				coords := ord.mulTab[rIdx][col]
				for row := 0; row < d; row++ {
					if coords[row].Sign() != 0 {
						prod.Mul(&bBig, coords[row])
						mulMat[row][col].Add(mulMat[row][col], &prod)
					}
				}
			}
		}
		for row := 0; row < d; row++ {
			for col := 0; col < d; col++ {
				mulMod[row][col] = rem.Mod(mulMat[row][col], bigP).Uint64()
			}
		}

		w := make([]*big.Int, d)
		nextW := make([]*big.Int, d)
		for k := 0; k < d; k++ {
			w[k] = new(big.Int).Mul(ord.aMat[0][k], bigP)
			nextW[k] = new(big.Int)
		}
		eIdx := 0
		for {
			divAll := true
			for row := 0; row < d; row++ {
				nextW[row].SetInt64(0)
				for col := 0; col < d; col++ {
					if mulMat[row][col].Sign() != 0 && w[col].Sign() != 0 {
						prod.Mul(mulMat[row][col], w[col])
						nextW[row].Add(nextW[row], &prod)
					}
				}
				if rem.Mod(nextW[row], bigP).Sign() != 0 {
					divAll = false
					break
				}
			}
			if !divAll {
				break
			}
			for row := 0; row < d; row++ {
				w[row].Quo(nextW[row], bigP)
			}
			eIdx++
		}

		col := *nextCol
		*nextCol++
		out = append(out, BadPrimeIdeal{
			Col:    col,
			E:      eIdx,
			F:      fDeg,
			mulMat: mulMat,
			mulMod: mulMod,
		})
	}
	return out
}

// splitEtaleIdempotents splits an etale F_p algebra via Frobenius fixed spaces (Cohen, GTM 138, 1993, Alg. 6.2.9).
func splitEtaleIdempotents(subBasis [][]uint64, oneS []uint64, bMod [][][]uint64, sDim int, p uint64) [][]uint64 {
	kDim := len(subBasis)
	if kDim == 1 {
		return [][]uint64{oneS}
	}
	frobCols := make([][]uint64, kDim)
	for j := 0; j < kDim; j++ {
		xp := algPowVecModP(subBasis[j], oneS, p, bMod, sDim, p)
		diff := make([]uint64, sDim)
		for c := 0; c < sDim; c++ {
			if xp[c] >= subBasis[j][c] {
				diff[c] = xp[c] - subBasis[j][c]
			} else {
				diff[c] = xp[c] + p - subBasis[j][c]
			}
		}
		frobCols[j] = diff
	}
	var eqs [][]uint64
	for c := 0; c < sDim; c++ {
		row := make([]uint64, kDim)
		for j := 0; j < kDim; j++ {
			row[j] = frobCols[j][c]
		}
		eqs = append(eqs, row)
	}
	fixed := rightNullspaceModP(eqs, kDim, p)
	if len(fixed) <= 1 {
		return [][]uint64{oneS}
	}

	var z []uint64
	for _, comb := range fixed {
		cand := make([]uint64, sDim)
		for j := 0; j < kDim; j++ {
			if comb[j] == 0 {
				continue
			}
			for c := 0; c < sDim; c++ {
				cand[c] = (cand[c] + arith.MulMod64(comb[j], subBasis[j][c], p)) % p
			}
		}
		if !vecProportionalModP(cand, oneS, p) {
			z = cand
			break
		}
	}

	minPoly := minPolyModP(z, oneS, bMod, sDim, kDim, p)
	roots := arith.PolyRootsModP(minPoly, p)
	var out [][]uint64
	for idx, r := range roots {
		idem := append([]uint64(nil), oneS...)
		for t, rt := range roots {
			if t == idx {
				continue
			}
			invDiff := arith.InvMod((r+p-rt)%p, p)
			term := make([]uint64, sDim)
			for c := 0; c < sDim; c++ {
				sub := arith.MulMod64(rt, oneS[c], p)
				diff := (z[c] + p - sub) % p
				term[c] = arith.MulMod64(diff, invDiff, p)
			}
			idem = algMulVecModP(idem, term, bMod, sDim, p)
		}
		var childSpan [][]uint64
		for _, bv := range subBasis {
			childSpan = append(childSpan, algMulVecModP(bv, idem, bMod, sDim, p))
		}
		childBasis, _ := rrefModP(childSpan, sDim, p)
		out = append(out, splitEtaleIdempotents(childBasis, idem, bMod, sDim, p)...)
	}
	return out
}

func minPolyModP(z, oneS []uint64, bMod [][][]uint64, sDim, maxDeg int, p uint64) []uint64 {
	pows := make([][]uint64, 0, maxDeg+2)
	pows = append(pows, append([]uint64(nil), oneS...))
	cur := append([]uint64(nil), oneS...)
	for deg := 1; deg <= maxDeg+1; deg++ {
		cur = algMulVecModP(cur, z, bMod, sDim, p)
		pows = append(pows, cur)
		nCols := deg + 1
		var eqs [][]uint64
		for c := 0; c < sDim; c++ {
			row := make([]uint64, nCols)
			for j := 0; j < nCols; j++ {
				row[j] = pows[j][c]
			}
			eqs = append(eqs, row)
		}
		ns := rightNullspaceModP(eqs, nCols, p)
		for _, v := range ns {
			if v[deg] != 0 {
				inv := arith.InvMod(v[deg], p)
				poly := make([]uint64, nCols)
				for j := 0; j < nCols; j++ {
					poly[j] = arith.MulMod64(v[j], inv, p)
				}
				return poly
			}
		}
	}
	return []uint64{0, 1}
}

func vecProportionalModP(a, b []uint64, p uint64) bool {
	piv := -1
	for i := range b {
		if b[i] != 0 {
			piv = i
			break
		}
	}
	if piv < 0 {
		return true
	}
	ratio := arith.MulMod64(a[piv], arith.InvMod(b[piv], p), p)
	for i := range a {
		if a[i] != arith.MulMod64(ratio, b[i], p) {
			return false
		}
	}
	return true
}

func mulTableModP(ord *pOrder, p uint64) [][][]uint64 {
	d := ord.d
	var bigP, rem big.Int
	bigP.SetUint64(p)
	out := make([][][]uint64, d)
	for i := 0; i < d; i++ {
		out[i] = make([][]uint64, d)
		for j := 0; j < d; j++ {
			out[i][j] = make([]uint64, d)
			for k := 0; k < d; k++ {
				out[i][j][k] = rem.Mod(ord.mulTab[i][j][k], &bigP).Uint64()
			}
		}
	}
	return out
}

func algMulVecModP(u, v []uint64, tab [][][]uint64, dim int, p uint64) []uint64 {
	out := make([]uint64, dim)
	for i := 0; i < dim; i++ {
		if u[i] == 0 {
			continue
		}
		for j := 0; j < dim; j++ {
			if v[j] == 0 {
				continue
			}
			uv := arith.MulMod64(u[i], v[j], p)
			coords := tab[i][j]
			for k := 0; k < dim; k++ {
				if coords[k] != 0 {
					out[k] = (out[k] + arith.MulMod64(uv, coords[k], p)) % p
				}
			}
		}
	}
	return out
}

func algPowVecModP(base, one []uint64, exp uint64, tab [][][]uint64, dim int, p uint64) []uint64 {
	res := append([]uint64(nil), one...)
	b := append([]uint64(nil), base...)
	for e := exp; e > 0; e >>= 1 {
		if e&1 != 0 {
			res = algMulVecModP(res, b, tab, dim, p)
		}
		if e > 1 {
			b = algMulVecModP(b, b, tab, dim, p)
		}
	}
	return res
}

func algMulExactInt(u, v []uint64, ord *pOrder) []*big.Int {
	d := ord.d
	out := make([]*big.Int, d)
	for k := 0; k < d; k++ {
		out[k] = new(big.Int)
	}
	var uvBig, term big.Int
	for i := 0; i < d; i++ {
		if u[i] == 0 {
			continue
		}
		for j := 0; j < d; j++ {
			if v[j] == 0 {
				continue
			}
			uvBig.Mul(new(big.Int).SetUint64(u[i]), new(big.Int).SetUint64(v[j]))
			coords := ord.mulTab[i][j]
			for k := 0; k < d; k++ {
				if coords[k].Sign() != 0 {
					term.Mul(&uvBig, coords[k])
					out[k].Add(out[k], &term)
				}
			}
		}
	}
	return out
}

func nilradicalModP(mMod [][][]uint64, d int, p uint64) [][]uint64 {
	q := p
	for q < uint64(d) {
		q *= p
	}
	powQ := func(x []uint64) []uint64 {
		var res []uint64
		b := append([]uint64(nil), x...)
		for e := q; e > 0; e >>= 1 {
			if e&1 != 0 {
				if res == nil {
					res = append([]uint64(nil), b...)
				} else {
					res = algMulVecModP(res, b, mMod, d, p)
				}
			}
			if e > 1 {
				b = algMulVecModP(b, b, mMod, d, p)
			}
		}
		return res
	}
	cols := make([][]uint64, d)
	for j := 0; j < d; j++ {
		ej := make([]uint64, d)
		ej[j] = 1
		cols[j] = powQ(ej)
	}
	eqs := make([][]uint64, d)
	for r := 0; r < d; r++ {
		eqs[r] = make([]uint64, d)
		for c := 0; c < d; c++ {
			eqs[r][c] = cols[c][r]
		}
	}
	return rightNullspaceModP(eqs, d, p)
}

func rrefModP(rows [][]uint64, nCols int, p uint64) (rref [][]uint64, pivCols []int) {
	mat := make([][]uint64, len(rows))
	for i, r := range rows {
		cp := make([]uint64, nCols)
		copy(cp, r)
		mat[i] = cp
	}
	r := 0
	for col := 0; col < nCols && r < len(mat); col++ {
		piv := -1
		for i := r; i < len(mat); i++ {
			if mat[i][col] != 0 {
				piv = i
				break
			}
		}
		if piv < 0 {
			continue
		}
		mat[r], mat[piv] = mat[piv], mat[r]
		inv := arith.InvMod(mat[r][col], p)
		for c := col; c < nCols; c++ {
			mat[r][c] = arith.MulMod64(mat[r][c], inv, p)
		}
		for i := 0; i < len(mat); i++ {
			if i == r || mat[i][col] == 0 {
				continue
			}
			f := mat[i][col]
			for c := col; c < nCols; c++ {
				sub := arith.MulMod64(f, mat[r][c], p)
				if mat[i][c] >= sub {
					mat[i][c] -= sub
				} else {
					mat[i][c] = mat[i][c] + p - sub
				}
			}
		}
		pivCols = append(pivCols, col)
		r++
	}
	return mat[:r], pivCols
}

func rightNullspaceModP(rows [][]uint64, nCols int, p uint64) [][]uint64 {
	rref, pivCols := rrefModP(rows, nCols, p)
	isPiv := make([]bool, nCols)
	for _, c := range pivCols {
		isPiv[c] = true
	}
	var basis [][]uint64
	for freeCol := 0; freeCol < nCols; freeCol++ {
		if isPiv[freeCol] {
			continue
		}
		vec := make([]uint64, nCols)
		vec[freeCol] = 1
		for r, pc := range pivCols {
			if rref[r][freeCol] != 0 {
				vec[pc] = (p - rref[r][freeCol]) % p
			}
		}
		basis = append(basis, vec)
	}
	return basis
}

func hnfUpperSquare(rows [][]*big.Int, d int) [][]*big.Int {
	mat := make([][]*big.Int, len(rows))
	for i, r := range rows {
		cp := make([]*big.Int, d)
		for j := 0; j < d; j++ {
			cp[j] = new(big.Int).Set(r[j])
		}
		mat[i] = cp
	}
	var q, rem, prod big.Int
	for col := 0; col < d; col++ {
		for {
			bestRow := -1
			var bestAbs *big.Int
			for r := col; r < len(mat); r++ {
				if mat[r][col].Sign() == 0 {
					continue
				}
				ab := new(big.Int).Abs(mat[r][col])
				if bestRow < 0 || ab.Cmp(bestAbs) < 0 {
					bestRow = r
					bestAbs = ab
				}
			}
			if bestRow < 0 {
				break
			}
			mat[col], mat[bestRow] = mat[bestRow], mat[col]
			if mat[col][col].Sign() < 0 {
				for c := col; c < d; c++ {
					mat[col][c].Neg(mat[col][c])
				}
			}
			anyNonZero := false
			for r := col + 1; r < len(mat); r++ {
				if mat[r][col].Sign() == 0 {
					continue
				}
				anyNonZero = true
				q.DivMod(mat[r][col], mat[col][col], &rem)
				for c := col; c < d; c++ {
					prod.Mul(&q, mat[col][c])
					mat[r][c].Sub(mat[r][c], &prod)
				}
			}
			if !anyNonZero {
				break
			}
		}
		for r := 0; r < col; r++ {
			if mat[r][col].Sign() == 0 {
				continue
			}
			q.DivMod(mat[r][col], mat[col][col], &rem)
			for c := col; c < d; c++ {
				prod.Mul(&q, mat[col][c])
				mat[r][c].Sub(mat[r][c], &prod)
			}
		}
	}
	return mat[:d]
}

func matDivisibleBy(mat [][]*big.Int, div *big.Int) bool {
	var rem big.Int
	for i := range mat {
		for j := range mat[i] {
			if rem.Mod(mat[i][j], div).Sign() != 0 {
				return false
			}
		}
	}
	return true
}

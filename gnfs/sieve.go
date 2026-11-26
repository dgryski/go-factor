package gnfs

import (
	"encoding/binary"
	"math"
	"math/big"
	"math/bits"
	"runtime"
	"sync"
	"sync/atomic"

	"github.com/dgryski/go-factor/internal/arith"
)

const (
	blockShift    = 15 // 32 KB sieve blocks
	blockSize     = 1 << blockShift
	minSievePrime = 29
)

// LatticeQ represents a special-q algebraic ideal (Q, R) and its Gauss-reduced
// lattice basis u = (U1, U2), v = (V1, V2) with det(u, v) = ±Q.
type LatticeQ struct {
	Q, R   uint32
	AFBIdx int32 // index in fb.AFB if Q <= Ba, else -1
	U1, U2 int64
	V1, V2 int64
	IHalf  int // sieve i in [-IHalf, IHalf)
	JMax   int // sieve j in [1, JMax)
}

type linePrime struct {
	idx    int32
	p      uint32
	lambda uint32
	logP   byte
}

type vertPrime struct {
	idx  int32
	p    uint32
	logP byte
}

// ReduceLattice constructs a Gauss-reduced basis for L = {(a, b) : a ≡ b*r (mod q)}
// under the skewed norm x^2 + skew^2 * y^2.
func ReduceLattice(q, r uint32, skew float64, iHalf, jMax int, afbIdx int32) LatticeQ {
	u1, u2 := float64(q), 0.0
	v1, v2 := float64(r), 1.0
	s2 := skew * skew

	norm2 := func(x, y float64) float64 {
		return x*x + s2*y*y
	}

	for {
		if norm2(u1, u2) > norm2(v1, v2) {
			u1, v1 = v1, u1
			u2, v2 = v2, u2
		}
		den := norm2(u1, u2)
		if den == 0 {
			break
		}
		mu := math.Round((u1*v1 + s2*u2*v2) / den)
		if mu == 0 {
			break
		}
		v1 -= mu * u1
		v2 -= mu * u2
	}

	return LatticeQ{
		Q:      q,
		R:      r,
		AFBIdx: afbIdx,
		U1:     int64(math.Round(u1)),
		U2:     int64(math.Round(u2)),
		V1:     int64(math.Round(v1)),
		V2:     int64(math.Round(v2)),
		IHalf:  iHalf,
		JMax:   jMax,
	}
}

func invMod32(a, p uint32) uint32 {
	if p >= 1<<31 {
		return uint32(arith.InvMod(uint64(a), uint64(p)))
	}
	t, newT := int32(0), int32(1)
	r, newR := int32(p), int32(a)
	for newR != 0 {
		q := r / newR
		t, newT = newT, t-q*newT
		r, newR = newR, r-q*newR
	}
	if t < 0 {
		t += int32(p)
	}
	return uint32(t)
}

// transformRoot computes the lattice-transformed hit condition for ideal (p, r):
//
//	i ≡ lambda * j (mod p), or vertical (j ≡ 0 mod p) if vert is true, or all (everywhere) if all is true.
func transformRoot(p, r uint32, lq *LatticeQ) (lambda uint32, vert, all bool) {
	u1, u2 := lq.U1, lq.U2
	v1, v2 := lq.V1, lq.V2
	p64 := uint64(p)
	var ap, bp uint64
	if r == p {
		// Projective: b ≡ 0 (mod p) => i*u2 + j*v2 ≡ 0 (mod p)
		ap = modInt64(u2, p64)
		bp = modInt64(-v2, p64)
	} else {
		// Affine: a ≡ r*b (mod p) => i*(u1 - r*u2) ≡ j*(r*v2 - v1) (mod p)
		r64 := uint64(r)
		u1m := modInt64(u1, p64)
		u2m := modInt64(u2, p64)
		v1m := modInt64(v1, p64)
		v2m := modInt64(v2, p64)
		ru2 := (r64 * u2m) % p64
		if u1m >= ru2 {
			ap = u1m - ru2
		} else {
			ap = p64 - (ru2 - u1m)
		}
		rv2 := (r64 * v2m) % p64
		if rv2 >= v1m {
			bp = rv2 - v1m
		} else {
			bp = p64 - (v1m - rv2)
		}
	}

	if ap == 0 {
		if bp == 0 {
			return 0, false, true
		}
		return 0, true, false
	}
	invAp := invMod32(uint32(ap), p)
	return uint32((bp * uint64(invAp)) % p64), false, false
}

// fkReduce reduces the p-lattice basis (-p, 0), (lambda, 1) to Franke-Kleinjung
// vectors (a0, b0), (a1, b1) satisfying -W < a0 < 0 < a1 < W and a1 - a0 >= W.
// Requires 0 < lambda < p and p >= W >= 2.
func fkReduce(p, lambda, W int64) (a0, b0, a1, b1 int64) {
	a0, b0 = -p, 0
	a1, b1 = lambda, 1
	for -a0 >= W || a1 >= W {
		if -a0 >= W && a1 < W {
			k := (-a0-W)/a1 + 1
			a0 += k * a1
			b0 += k * b1
			break
		}
		if a1 >= W && -a0 < W {
			k := (a1-W)/(-a0) + 1
			a1 += k * a0
			b1 += k * b0
			break
		}
		if -a0 >= a1 {
			k := (-a0) / a1
			a0 += k * a1
			b0 += k * b1
		} else {
			k := a1 / (-a0)
			a1 += k * a0
			b1 += k * b0
		}
	}
	return a0, b0, a1, b1
}

// addFKHits appends Franke-Kleinjung bucket hits for p >= 2*I while skipping both-even (i, j).
func addFKHits(p, lambda uint32, I, J int, tag uint64, buckets [][]uint64) {
	W := 2 * I
	total := (J - 1) * W
	if total <= 0 {
		return
	}
	if lambda == 0 {
		// Skip even j since i = 0 is even.
		for adjPos := I; adjPos < total; adjPos += 2 * W {
			bi := adjPos >> blockShift
			buckets[bi] = append(buckets[bi], tag|uint64(adjPos&(blockSize-1)))
		}
		return
	}

	a0, b0, a1, b1 := fkReduce(int64(p), int64(lambda), int64(W))
	x := I
	pos := I // starts at (i=0, j=0); subtract W to get adjPos
	thresh0 := int(-a0)
	thresh1 := W - int(a1)
	step0 := int(b0)*W + int(a0)
	step1 := int(b1)*W + int(a1)
	step01 := step0 + step1
	dx0 := int(a0)
	dx1 := int(a1)
	dx01 := int(a0 + a1)
	limit := total + W

	// When W = 2*I is a power of 2 and I is even, bit 0 of pos is (i & 1) and bit W is (j & 1).
	if W&(W-1) == 0 && I&1 == 0 {
		evenMask := W | 1
		for {
			if x >= thresh0 {
				pos += step0
				x += dx0
			} else if x < thresh1 {
				pos += step1
				x += dx1
			} else {
				pos += step01
				x += dx01
			}
			if pos >= limit {
				break
			}
			if pos&evenMask != 0 {
				adjPos := pos - W
				bi := adjPos >> blockShift
				buckets[bi] = append(buckets[bi], tag|uint64(adjPos&(blockSize-1)))
			}
		}
		return
	}

	for {
		if x >= thresh0 {
			pos += step0
			x += dx0
		} else if x < thresh1 {
			pos += step1
			x += dx1
		} else {
			pos += step01
			x += dx01
		}
		if pos >= limit {
			break
		}
		if ((x-I)|(pos/W))&1 != 0 {
			adjPos := pos - W
			bi := adjPos >> blockShift
			buckets[bi] = append(buckets[bi], tag|uint64(adjPos&(blockSize-1)))
		}
	}
}

// hasByteGe reports whether any byte in buf is >= thresh, scanning 8 bytes per word.
func hasByteGe(buf []byte, thresh int) bool {
	if thresh <= 0 {
		return len(buf) > 0
	}
	if thresh > 255 {
		return false
	}
	const low7 = 0x7f7f7f7f7f7f7f7f
	const high1 = 0x8080808080808080
	n8 := len(buf) &^ 7
	if thresh <= 128 {
		bias := uint64(128-thresh) * 0x0101010101010101
		for off := 0; off < n8; off += 8 {
			w := binary.LittleEndian.Uint64(buf[off:])
			if (((w&low7)+bias)|w)&high1 != 0 {
				return true
			}
		}
	} else {
		bias := uint64(256-thresh) * 0x0101010101010101
		for off := 0; off < n8; off += 8 {
			w := binary.LittleEndian.Uint64(buf[off:])
			if (((w&low7)+bias)&w)&high1 != 0 {
				return true
			}
		}
	}
	for off := n8; off < len(buf); off++ {
		if int(buf[off]) >= thresh {
			return true
		}
	}
	return false
}

// appendBothGe appends indices off where ab[off] >= aThresh and rb[off] >= rThresh,
// scanning 8 bytes per word.
func appendBothGe(dst []int, ab, rb []byte, aThresh, rThresh int) []int {
	n := len(ab)
	if aThresh > 255 || rThresh > 255 || n == 0 {
		return dst
	}
	if aThresh <= 0 && rThresh <= 0 {
		for off := 0; off < n; off++ {
			dst = append(dst, off)
		}
		return dst
	}
	const low7 = 0x7f7f7f7f7f7f7f7f
	const high1 = 0x8080808080808080
	n8 := n &^ 7

	if aThresh >= 1 && aThresh <= 128 && rThresh >= 1 && rThresh <= 128 {
		aBias := uint64(128-aThresh) * 0x0101010101010101
		rBias := uint64(128-rThresh) * 0x0101010101010101
		for off := 0; off < n8; off += 8 {
			wa := binary.LittleEndian.Uint64(ab[off:])
			ma := (((wa & low7) + aBias) | wa) & high1
			if ma == 0 {
				continue
			}
			wr := binary.LittleEndian.Uint64(rb[off:])
			m := ma & (((wr & low7) + rBias) | wr)
			for m != 0 {
				lane := bits.TrailingZeros64(m) >> 3
				dst = append(dst, off+lane)
				m &= m - 1
			}
		}
		for off := n8; off < n; off++ {
			if int(ab[off]) >= aThresh && int(rb[off]) >= rThresh {
				dst = append(dst, off)
			}
		}
		return dst
	}

	maskWord := func(w uint64, thresh int) uint64 {
		if thresh <= 0 {
			return high1
		}
		if thresh <= 128 {
			bias := uint64(128-thresh) * 0x0101010101010101
			return (((w & low7) + bias) | w) & high1
		}
		bias := uint64(256-thresh) * 0x0101010101010101
		return (((w & low7) + bias) & w) & high1
	}
	for off := 0; off < n8; off += 8 {
		ma := maskWord(binary.LittleEndian.Uint64(ab[off:]), aThresh)
		if ma == 0 {
			continue
		}
		m := ma & maskWord(binary.LittleEndian.Uint64(rb[off:]), rThresh)
		for m != 0 {
			lane := bits.TrailingZeros64(m) >> 3
			dst = append(dst, off+lane)
			m &= m - 1
		}
	}
	for off := n8; off < n; off++ {
		if int(ab[off]) >= aThresh && int(rb[off]) >= rThresh {
			dst = append(dst, off)
		}
	}
	return dst
}

// sieveLine runs the baseline line sieve over b in [1, LineB] and a in [-LineA, LineA),
// skipping the part with b <= prevB and a in [-prevA, prevA) that an earlier
// round already sieved.
func sieveLine(poly *Poly, fb *FactorBases, params *Params, prevA, prevB int64, onRel func(Relation), shouldStop func() bool) {
	var ratSmall, ratSieve []int32
	for i, id := range fb.RFB {
		if id.P < minSievePrime {
			ratSmall = append(ratSmall, int32(i))
		} else {
			ratSieve = append(ratSieve, int32(i))
		}
	}
	var algSmall, algSieve []int32
	for i, id := range fb.AFB {
		if id.P < minSievePrime {
			algSmall = append(algSmall, int32(i))
		} else {
			algSieve = append(algSieve, int32(i))
		}
	}

	smallSlackRat := expectedSmallLogSumRat(fb, minSievePrime)
	smallSlackAlg := expectedSmallLogSumAlg(fb, minSievePrime)

	nw := Workers
	if nw <= 0 {
		nw = runtime.GOMAXPROCS(0)
	}
	if nw > int(params.LineB) {
		nw = int(params.LineB)
	}
	if nw < 1 {
		nw = 1
	}

	var nextB atomic.Int64
	nextB.Store(1)
	var stop atomic.Bool
	var mu sync.Mutex
	var wg sync.WaitGroup

	for w := 0; w < nw; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			var ratBuf, algBuf [blockSize]byte
			scratch := newCofactorScratch(poly.Degree)
			var cands []int
			var rHits, aHits []int32

			for !stop.Load() {
				b := nextB.Add(1) - 1
				if b > params.LineB {
					return
				}
				evenB := b&1 == 0

				ranges := [2][2]int64{{-params.LineA, params.LineA}}
				if b <= prevB {
					ranges = [2][2]int64{{-params.LineA, -prevA}, {prevA, params.LineA}}
				}
				for _, rg := range ranges {
					hiA := rg[1]
					for bStart := rg[0]; bStart < hiA && !stop.Load(); bStart += blockSize {
						bLen := int(blockSize)
						if hiA-bStart < int64(bLen) {
							bLen = int(hiA - bStart)
						}
						rb := ratBuf[:bLen]
						ab := algBuf[:bLen]
						clear(rb)
						clear(ab)

						for _, idx := range ratSieve {
							id := fb.RFB[idx]
							p := int(id.P)
							if id.R == id.P {
								if b%int64(p) == 0 {
									if evenB {
										j := int(bStart&1) ^ 1
										for ; j < bLen; j += 2 {
											rb[j] += id.LogP
										}
									} else {
										for j := 0; j < bLen; j++ {
											rb[j] += id.LogP
										}
									}
								}
								continue
							}
							target := (b % int64(p)) * int64(id.R) % int64(p)
							j := int(nextOffsetInt64(bStart, target, int64(p)))
							if evenB {
								if (bStart+int64(j))&1 == 0 {
									j += p
								}
								step := p << 1
								for ; j < bLen; j += step {
									rb[j] += id.LogP
								}
							} else {
								for ; j < bLen; j += p {
									rb[j] += id.LogP
								}
							}
						}

						for _, idx := range algSieve {
							id := fb.AFB[idx]
							p := int(id.P)
							if id.R == id.P {
								if b%int64(p) == 0 {
									if evenB {
										j := int(bStart&1) ^ 1
										for ; j < bLen; j += 2 {
											ab[j] += id.LogP
										}
									} else {
										for j := 0; j < bLen; j++ {
											ab[j] += id.LogP
										}
									}
								}
								continue
							}
							target := (b % int64(p)) * int64(id.R) % int64(p)
							j := int(nextOffsetInt64(bStart, target, int64(p)))
							if evenB {
								if (bStart+int64(j))&1 == 0 {
									j += p
								}
								step := p << 1
								for ; j < bLen; j += step {
									ab[j] += id.LogP
								}
							} else {
								for ; j < bLen; j += p {
									ab[j] += id.LogP
								}
							}
						}

						rThresh := blockThresholdRat(poly, bStart, bStart+int64(bLen)-1, b, b, params.SlackRat+smallSlackRat, scratch)
						aThresh := blockThresholdAlg(poly, bStart, bStart+int64(bLen)-1, b, b, 0, params.SlackAlg+smallSlackAlg, scratch)

						cands = appendBothGe(cands[:0], ab, rb, aThresh, rThresh)
						if len(cands) == 0 {
							continue
						}

						for _, j := range cands {
							a := bStart + int64(j)
							if evenB && a&1 == 0 {
								continue
							}
							ua := uint64(a)
							if a < 0 {
								ua = uint64(-a)
							}
							if arith.Gcd64(ua, uint64(b)) != 1 {
								continue
							}
							rHits = rHits[:0]
							for _, idx := range ratSieve {
								id := fb.RFB[idx]
								p := uint64(id.P)
								if id.R == id.P {
									if uint64(b)%p == 0 {
										rHits = append(rHits, idx)
									}
								} else if modInt64(a, p) == (uint64(b)%p)*uint64(id.R)%p {
									rHits = append(rHits, idx)
								}
							}
							aHits = aHits[:0]
							for _, idx := range algSieve {
								id := fb.AFB[idx]
								p := uint64(id.P)
								if id.R == id.P {
									if uint64(b)%p == 0 {
										aHits = append(aHits, idx)
									}
								} else if modInt64(a, p) == (uint64(b)%p)*uint64(id.R)%p {
									aHits = append(aHits, idx)
								}
							}
							rel, ok := cofactorCandidate(poly, fb, params, a, b, ratSmall, rHits, algSmall, aHits, nil, scratch)
							if ok {
								mu.Lock()
								onRel(rel)
								if shouldStop() {
									stop.Store(true)
								}
								mu.Unlock()
								if stop.Load() {
									return
								}
							}
						}
					}
				}
			}
		}()
	}
	wg.Wait()
}

// buildSpecialQs returns the affine algebraic ideals (q, r) with QMin <= q <= QMax,
// leaving out multiple roots of f mod q.
func buildSpecialQs(poly *Poly, fb *FactorBases, params *Params) []LatticeQ {
	var qs []LatticeQ
	var bigP, rem big.Int
	qStart := params.QMin
	if qStart < 3 {
		qStart = 3
	}
	if qStart&1 == 0 {
		qStart++
	}
	for q := qStart; q <= params.QMax; q += 2 {
		if !arith.IsPrime64(q) {
			continue
		}
		bigP.SetUint64(q)
		if rem.Mod(poly.Coeffs[poly.Degree], &bigP).Sign() == 0 {
			continue
		}
		fp := poly.CoeffsModP(q)
		dfp := arith.PolyDerivModP(fp, q)
		for _, r := range arith.PolyRootsModP(fp, q) {
			if arith.PolyEvalModP(dfp, r, q) == 0 {
				continue // multiple root: v_q(F) need not be the ideal's valuation
			}
			afbIdx := int32(-1)
			if q <= uint64(params.Ba) {
				if idx, ok := fb.AlgLookup[algKey(uint32(q), uint32(r))]; ok {
					afbIdx = idx
				}
			}
			qs = append(qs, ReduceLattice(uint32(q), uint32(r), poly.Skew, params.IHalf, params.JMax, afbIdx))
		}
	}
	return qs
}

// sieveLattice runs the special-q lattice sieve with Franke-Kleinjung buckets.
func sieveLattice(poly *Poly, fb *FactorBases, params *Params, onRel func(Relation), shouldStop func() bool) int {
	sqs := buildSpecialQs(poly, fb, params)
	if len(sqs) == 0 {
		return 0
	}

	I := params.IHalf
	J := params.JMax
	W := 2 * I
	total := (J - 1) * W
	numBlocks := (total + blockSize - 1) >> blockShift

	var ratSmall, ratLineIdx, ratFK []int32
	for i, id := range fb.RFB {
		switch {
		case id.P < minSievePrime:
			ratSmall = append(ratSmall, int32(i))
		case int(id.P) < W:
			ratLineIdx = append(ratLineIdx, int32(i))
		default:
			ratFK = append(ratFK, int32(i))
		}
	}
	var algSmall, algLineIdx, algFK []int32
	for i, id := range fb.AFB {
		switch {
		case id.P < minSievePrime:
			algSmall = append(algSmall, int32(i))
		case int(id.P) < W:
			algLineIdx = append(algLineIdx, int32(i))
		default:
			algFK = append(algFK, int32(i))
		}
	}

	smallSlackRat := expectedSmallLogSumRat(fb, minSievePrime)
	smallSlackAlg := expectedSmallLogSumAlg(fb, minSievePrime)

	nw := Workers
	if nw <= 0 {
		nw = runtime.GOMAXPROCS(0)
	}
	if nw > len(sqs) {
		nw = len(sqs)
	}
	if nw < 1 {
		nw = 1
	}

	var nextQ atomic.Int64
	var qDone atomic.Int64
	var stop atomic.Bool
	var mu sync.Mutex
	var wg sync.WaitGroup

	type ratSurv struct {
		off     int
		i, j    int
		a, b    int64
		ratCols []int32
		ratLPs  [2]uint64
	}

	for w := 0; w < nw; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			var ratBuf, algBuf [blockSize]byte
			ratBuckets := make([][]uint64, numBlocks)
			algBuckets := make([][]uint64, numBlocks)
			scratch := newCofactorScratch(poly.Degree)

			ratLines := make([]linePrime, 0, len(ratLineIdx))
			var ratVerts []vertPrime
			var ratAll []int32
			algLines := make([]linePrime, 0, len(algLineIdx))
			var algVerts []vertPrime
			var algAll []int32

			var rawOffs []int
			var candPos []int
			var candI, candJ []int
			var candA, candB []int64
			var candBits [blockSize / 64]uint64
			var candSlot [blockSize]uint16
			var ratCandHits, algCandHits [][]int32
			var survivors []ratSurv

			for !stop.Load() {
				qi := int(nextQ.Add(1) - 1)
				if qi >= len(sqs) {
					return
				}
				lq := &sqs[qi]
				logQ := math.Log2(float64(lq.Q))

				for b := 0; b < numBlocks; b++ {
					ratBuckets[b] = ratBuckets[b][:0]
					algBuckets[b] = algBuckets[b][:0]
				}

				// Transform roots for line primes and fill FK buckets for p >= W.
				ratLines = ratLines[:0]
				ratVerts = ratVerts[:0]
				ratAll = ratAll[:0]
				var ratAllLog byte
				for _, idx := range ratLineIdx {
					id := fb.RFB[idx]
					lam, vert, all := transformRoot(id.P, id.R, lq)
					if all {
						ratAll = append(ratAll, idx)
						ratAllLog += id.LogP
					} else if vert {
						ratVerts = append(ratVerts, vertPrime{idx: idx, p: id.P, logP: id.LogP})
					} else {
						ratLines = append(ratLines, linePrime{idx: idx, p: id.P, lambda: lam, logP: id.LogP})
					}
				}

				algLines = algLines[:0]
				algVerts = algVerts[:0]
				algAll = algAll[:0]
				var algAllLog byte
				for _, idx := range algLineIdx {
					id := fb.AFB[idx]
					if id.P == lq.Q {
						continue
					}
					lam, vert, all := transformRoot(id.P, id.R, lq)
					if all {
						algAll = append(algAll, idx)
						algAllLog += id.LogP
					} else if vert {
						algVerts = append(algVerts, vertPrime{idx: idx, p: id.P, logP: id.LogP})
					} else {
						algLines = append(algLines, linePrime{idx: idx, p: id.P, lambda: lam, logP: id.LogP})
					}
				}

				for _, idx := range ratFK {
					id := fb.RFB[idx]
					lam, vert, all := transformRoot(id.P, id.R, lq)
					tag := uint64(id.LogP)<<48 | uint64(idx)<<16
					if all {
						ratAll = append(ratAll, idx)
						ratAllLog += id.LogP
					} else if vert {
						emitVertBuckets(int(id.P), I, J, W, tag, ratBuckets)
					} else {
						addFKHits(id.P, lam, I, J, tag, ratBuckets)
					}
				}

				for _, idx := range algFK {
					id := fb.AFB[idx]
					if id.P == lq.Q {
						continue
					}
					lam, vert, all := transformRoot(id.P, id.R, lq)
					tag := uint64(id.LogP)<<48 | uint64(idx)<<16
					if all {
						algAll = append(algAll, idx)
						algAllLog += id.LogP
					} else if vert {
						emitVertBuckets(int(id.P), I, J, W, tag, algBuckets)
					} else {
						addFKHits(id.P, lam, I, J, tag, algBuckets)
					}
				}

				// Sieve each block.
				for bi := 0; bi < numBlocks && !stop.Load(); bi++ {
					bStart := bi * blockSize
					bEnd := bStart + blockSize
					if bEnd > total {
						bEnd = total
					}
					bLen := bEnd - bStart

					ab := algBuf[:bLen]
					sieveBlockLines(ab, bStart, bEnd, I, W, algLines, algVerts, algAllLog)
					for _, e := range algBuckets[bi] {
						ab[uint16(e)] += byte(e >> 48)
					}

					aThresh := latticeBlockThresholdAlg(poly, lq, bStart, bEnd, I, W, logQ, params.SlackAlg+smallSlackAlg, scratch)
					if !hasByteGe(ab, aThresh) {
						continue
					}

					rb := ratBuf[:bLen]
					sieveBlockLines(rb, bStart, bEnd, I, W, ratLines, ratVerts, ratAllLog)
					for _, e := range ratBuckets[bi] {
						rb[uint16(e)] += byte(e >> 48)
					}
					rThresh := latticeBlockThresholdRat(poly, lq, bStart, bEnd, I, W, params.SlackRat+smallSlackRat, scratch)

					rawOffs = appendBothGe(rawOffs[:0], ab, rb, aThresh, rThresh)
					if len(rawOffs) == 0 {
						continue
					}

					candPos = candPos[:0]
					candI = candI[:0]
					candJ = candJ[:0]
					candA = candA[:0]
					candB = candB[:0]
					for _, off := range rawOffs {
						gPos := bStart + off
						j := 1 + gPos/W
						i := (gPos % W) - I
						if (i|j)&1 == 0 {
							continue
						}
						ui := uint64(i)
						if i < 0 {
							ui = uint64(-i)
						}
						if arith.Gcd64(ui, uint64(j)) != 1 {
							continue
						}
						a := int64(i)*lq.U1 + int64(j)*lq.V1
						b := int64(i)*lq.U2 + int64(j)*lq.V2
						if b == 0 {
							continue
						}
						if b < 0 {
							a, b = -a, -b
						}
						ua := uint64(a)
						if a < 0 {
							ua = uint64(-a)
						}
						if arith.Gcd64(ua, uint64(b)) != 1 {
							continue
						}
						candPos = append(candPos, off)
						candI = append(candI, i)
						candJ = append(candJ, j)
						candA = append(candA, a)
						candB = append(candB, b)
					}
					if len(candPos) == 0 {
						continue
					}

					// Stage 1 filters G(a, b) before scanning algebraic buckets.
					for len(ratCandHits) < len(candPos) {
						ratCandHits = append(ratCandHits, nil)
					}
					for k, off := range candPos {
						candBits[off/64] |= 1 << (uint(off) % 64)
						candSlot[off] = uint16(k)
						ratCandHits[k] = ratCandHits[k][:0]
					}
					for _, e := range ratBuckets[bi] {
						if off := uint16(e); candBits[off/64]&(1<<(off%64)) != 0 {
							k := candSlot[off]
							ratCandHits[k] = append(ratCandHits[k], int32(uint32(e>>16)))
						}
					}
					for _, off := range candPos {
						candBits[off/64] = 0
					}

					survivors = survivors[:0]
					for k := range candPos {
						ratCols, ratLPs, ok := cofactorRatLattice(
							poly, fb, params,
							candI[k], candJ[k], candA[k], candB[k],
							ratSmall, ratLines, ratVerts, ratAll, ratCandHits[k],
							scratch,
						)
						if ok {
							survivors = append(survivors, ratSurv{
								off:     candPos[k],
								i:       candI[k],
								j:       candJ[k],
								a:       candA[k],
								b:       candB[k],
								ratCols: ratCols,
								ratLPs:  ratLPs,
							})
						}
					}
					if len(survivors) == 0 {
						continue
					}

					// Stage 2 scans algBuckets[bi] only for Stage 1 survivors.
					for len(algCandHits) < len(survivors) {
						algCandHits = append(algCandHits, nil)
					}
					for k, sv := range survivors {
						off := sv.off
						candBits[off/64] |= 1 << (uint(off) % 64)
						candSlot[off] = uint16(k)
						algCandHits[k] = algCandHits[k][:0]
					}
					for _, e := range algBuckets[bi] {
						if off := uint16(e); candBits[off/64]&(1<<(off%64)) != 0 {
							k := candSlot[off]
							algCandHits[k] = append(algCandHits[k], int32(uint32(e>>16)))
						}
					}
					for _, sv := range survivors {
						candBits[sv.off/64] = 0
					}

					for k := range survivors {
						sv := &survivors[k]
						rel, ok := cofactorAlgLattice(
							poly, fb, params,
							sv.i, sv.j, sv.a, sv.b,
							sv.ratCols, sv.ratLPs,
							algSmall, algLines, algVerts, algAll, algCandHits[k],
							lq, scratch,
						)
						if ok {
							mu.Lock()
							onRel(rel)
							if shouldStop() {
								stop.Store(true)
							}
							mu.Unlock()
							if stop.Load() {
								return
							}
						}
					}
				}
				qDone.Add(1)
			}
		}()
	}
	wg.Wait()
	return int(qDone.Load())
}

func emitVertBuckets(p, I, J, W int, tag uint64, buckets [][]uint64) {
	for j := p; j < J; j += p {
		base := (j - 1) * W
		if j&1 == 0 {
			// Skip even i on even j lines where (x ^ I)&1 == 1 selects odd i.
			for x := (I & 1) ^ 1; x < W; x += 2 {
				pos := base + x
				bi := pos >> blockShift
				buckets[bi] = append(buckets[bi], tag|uint64(pos&(blockSize-1)))
			}
		} else {
			for x := 0; x < W; x++ {
				pos := base + x
				bi := pos >> blockShift
				buckets[bi] = append(buckets[bi], tag|uint64(pos&(blockSize-1)))
			}
		}
	}
}

// sieveBlockLines sieves all line-sieve primes (< W) across [bStart, bEnd),
// looping prime-outer over the L1-resident block and skipping both-even (i, j) on even j.
func sieveBlockLines(buf []byte, bStart, bEnd, I, W int, lines []linePrime, verts []vertPrime, allLog byte) {
	if allLog > 0 {
		for i := range buf {
			buf[i] = allLog
		}
	} else {
		clear(buf)
	}

	// Fast path when block boundaries align with full lines of width W (always true for power-of-2 I <= 16384).
	if bStart%W == 0 && bEnd%W == 0 {
		j0 := 1 + bStart/W
		nLines := (bEnd - bStart) / W
		iOddFirst := (I & 1) ^ 1
		for _, vp := range verts {
			vpP := int(vp.p)
			logP := vp.logP
			rem := j0 % vpP
			firstL := 0
			if rem != 0 {
				firstL = vpP - rem
			}
			for l := firstL; l < nLines; l += vpP {
				seg := buf[l*W : (l+1)*W]
				if (j0+l)&1 == 0 {
					for x := iOddFirst; x < W; x += 2 {
						seg[x] += logP
					}
				} else {
					for x := 0; x < W; x++ {
						seg[x] += logP
					}
				}
			}
		}

		i64 := uint64(I)
		j064 := uint64(j0)
		for k := range lines {
			lp := lines[k]
			p := int(lp.p)
			p2 := p << 1
			lam := int(lp.lambda)
			logP := lp.logP
			off := int((j064*uint64(lp.lambda) + i64) % uint64(lp.p))
			for l := 0; l < nLines; l++ {
				seg := buf[l*W : (l+1)*W]
				x := off
				if (j0+l)&1 == 0 {
					if (x-I)&1 == 0 {
						x += p
					}
					for ; x < W; x += p2 {
						seg[x] += logP
					}
				} else {
					for ; x < W; x += p {
						seg[x] += logP
					}
				}
				off += lam
				if off >= p {
					off -= p
				}
			}
		}
		return
	}

	// General path for arbitrary non-aligned W.
	j0 := 1 + bStart/W
	for _, vp := range verts {
		vpP := int(vp.p)
		logP := vp.logP
		for pos := bStart; pos < bEnd; {
			j := 1 + pos/W
			xFrom := pos % W
			segLen := W - xFrom
			if bEnd-pos < segLen {
				segLen = bEnd - pos
			}
			if j%vpP == 0 {
				seg := buf[pos-bStart : pos-bStart+segLen]
				if j&1 == 0 {
					x := ((xFrom - I) & 1) ^ 1
					for ; x < segLen; x += 2 {
						seg[x] += logP
					}
				} else {
					for x := 0; x < segLen; x++ {
						seg[x] += logP
					}
				}
			}
			pos += segLen
		}
	}

	i64 := uint64(I)
	for k := range lines {
		lp := lines[k]
		p := int(lp.p)
		p2 := p << 1
		lam := int(lp.lambda)
		logP := lp.logP
		off := int((uint64(j0)*uint64(lp.lambda) + i64) % uint64(lp.p))
		for pos := bStart; pos < bEnd; {
			j := 1 + pos/W
			xFrom := pos % W
			segLen := W - xFrom
			if bEnd-pos < segLen {
				segLen = bEnd - pos
			}
			seg := buf[pos-bStart : pos-bStart+segLen]
			xGlobal := off
			if xGlobal < xFrom {
				xGlobal += ((xFrom - xGlobal + p - 1) / p) * p
			}
			x := xGlobal - xFrom
			if j&1 == 0 {
				if (xGlobal-I)&1 == 0 {
					x += p
				}
				for ; x < segLen; x += p2 {
					seg[x] += logP
				}
			} else {
				for ; x < segLen; x += p {
					seg[x] += logP
				}
			}
			if xFrom+segLen == W {
				off += lam
				if off >= p {
					off -= p
				}
			}
			pos += segLen
		}
	}
}

func nextOffsetInt64(start, root, p int64) int64 {
	off := (root - start) % p
	if off < 0 {
		off += p
	}
	return off
}

func expectedSmallLogSumRat(fb *FactorBases, cut uint32) float64 {
	sum := 0.0
	for _, id := range fb.RFB {
		if id.P >= cut {
			break
		}
		sum += float64(id.LogP) / float64(id.P+1)
	}
	return sum
}

func expectedSmallLogSumAlg(fb *FactorBases, cut uint32) float64 {
	sum := 0.0
	for _, id := range fb.AFB {
		if id.P >= cut {
			break
		}
		sum += float64(id.LogP) / float64(id.P+1)
	}
	return sum
}

func blockThresholdRat(poly *Poly, a0, a1, b0, b1 int64, slack float64, sc *cofactorScratch) int {
	best := 0.0
	for _, a := range [2]int64{a0, a1} {
		for _, b := range [2]int64{b0, b1} {
			evalGInto(&sc.val, poly, a, b, sc)
			sc.val.Abs(&sc.val)
			if sc.val.Sign() > 0 {
				if l := arith.Log2Big(&sc.val); l > best {
					best = l
				}
			}
		}
	}
	t := int(math.Floor(best - slack))
	if t < 0 {
		t = 0
	}
	return t
}

func blockThresholdAlg(poly *Poly, a0, a1, b0, b1 int64, logQ, slack float64, sc *cofactorScratch) int {
	best := 0.0
	for _, a := range [2]int64{a0, a1} {
		for _, b := range [2]int64{b0, b1} {
			evalFInto(&sc.val, poly, a, b, sc)
			sc.val.Abs(&sc.val)
			if sc.val.Sign() > 0 {
				if l := arith.Log2Big(&sc.val); l > best {
					best = l
				}
			}
		}
	}
	t := int(math.Floor(best - logQ - slack))
	if t < 0 {
		t = 0
	}
	return t
}

func latticeBlockThresholdRat(poly *Poly, lq *LatticeQ, bStart, bEnd, I, W int, slack float64, sc *cofactorScratch) int {
	j0 := int64(1 + bStart/W)
	j1 := int64(1 + (bEnd-1)/W)
	i0 := int64((bStart % W) - I)
	i1 := int64(((bEnd - 1) % W) - I)
	if j1 > j0 {
		i0, i1 = int64(-I), int64(I-1)
	}
	best := 0.0
	for _, ii := range [2]int64{i0, i1} {
		for _, jj := range [2]int64{j0, j1} {
			a := ii*lq.U1 + jj*lq.V1
			b := ii*lq.U2 + jj*lq.V2
			evalGInto(&sc.val, poly, a, b, sc)
			sc.val.Abs(&sc.val)
			if sc.val.Sign() > 0 {
				if l := arith.Log2Big(&sc.val); l > best {
					best = l
				}
			}
		}
	}
	t := int(math.Floor(best - slack))
	if t < 0 {
		t = 0
	}
	return t
}

func latticeBlockThresholdAlg(poly *Poly, lq *LatticeQ, bStart, bEnd, I, W int, logQ, slack float64, sc *cofactorScratch) int {
	j0 := int64(1 + bStart/W)
	j1 := int64(1 + (bEnd-1)/W)
	i0 := int64((bStart % W) - I)
	i1 := int64(((bEnd - 1) % W) - I)
	if j1 > j0 {
		i0, i1 = int64(-I), int64(I-1)
	}
	best := 0.0
	for _, ii := range [2]int64{i0, i1} {
		for _, jj := range [2]int64{j0, j1} {
			a := ii*lq.U1 + jj*lq.V1
			b := ii*lq.U2 + jj*lq.V2
			evalFInto(&sc.val, poly, a, b, sc)
			sc.val.Abs(&sc.val)
			if sc.val.Sign() > 0 {
				if l := arith.Log2Big(&sc.val); l > best {
					best = l
				}
			}
		}
	}
	t := int(math.Floor(best - logQ - slack))
	if t < 0 {
		t = 0
	}
	return t
}

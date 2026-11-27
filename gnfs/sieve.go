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
	blockShift    = 14 // 16 KB sieve blocks keep both sides L1-resident and fit offsets in 14 bits
	blockSize     = 1 << blockShift
	blockMask     = blockSize - 1
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

type logGroup struct {
	logP       byte
	start, end int
}

type overHit struct {
	idx  int32
	off  uint16
	logP byte
}

// fkArena stores Franke-Kleinjung bucket hits as packed uint32 entries
// (idx << blockShift | off) grouped by LogP in a flat slab per block.
type fkArena struct {
	data      []uint32
	counts    []uint32
	groupEnds [][]uint32
	overflow  [][]overHit
	hasOver   bool
	stride    uint32
	shift     uint
}

func buildLogGroupsRat(rfb []RatIdeal, idxs []int32) []logGroup {
	if len(idxs) == 0 {
		return nil
	}
	var out []logGroup
	curLog := rfb[idxs[0]].LogP
	start := 0
	for i := 1; i < len(idxs); i++ {
		lg := rfb[idxs[i]].LogP
		if lg != curLog {
			out = append(out, logGroup{logP: curLog, start: start, end: i})
			curLog = lg
			start = i
		}
	}
	return append(out, logGroup{logP: curLog, start: start, end: len(idxs)})
}

func buildLogGroupsAlg(afb []AlgIdeal, idxs []int32) []logGroup {
	if len(idxs) == 0 {
		return nil
	}
	var out []logGroup
	curLog := afb[idxs[0]].LogP
	start := 0
	for i := 1; i < len(idxs); i++ {
		lg := afb[idxs[i]].LogP
		if lg != curLog {
			out = append(out, logGroup{logP: curLog, start: start, end: i})
			curLog = lg
			start = i
		}
	}
	return append(out, logGroup{logP: curLog, start: start, end: len(idxs)})
}

func newFKArena(numBlocks, numGroups, expectedPerBlock int) *fkArena {
	want := expectedPerBlock * 3
	if want < 2048 {
		want = 2048
	}
	shift := uint(11)
	for (1 << shift) < want && shift < 17 {
		shift++
	}
	stride := uint32(1 << shift)
	ge := make([][]uint32, numGroups)
	geBacking := make([]uint32, numGroups*numBlocks)
	for g := 0; g < numGroups; g++ {
		ge[g] = geBacking[g*numBlocks : (g+1)*numBlocks]
	}
	return &fkArena{
		data:      make([]uint32, numBlocks<<shift),
		counts:    make([]uint32, numBlocks),
		groupEnds: ge,
		overflow:  make([][]overHit, numBlocks),
		stride:    stride,
		shift:     shift,
	}
}

func (a *fkArena) reset() {
	clear(a.counts)
	if a.hasOver {
		for i := range a.overflow {
			a.overflow[i] = a.overflow[i][:0]
		}
		a.hasOver = false
	}
}

func (a *fkArena) addOverflow(bi int, idx int32, off uint16, logP byte) {
	a.hasOver = true
	a.overflow[bi] = append(a.overflow[bi], overHit{idx: idx, off: off, logP: logP})
}

func (a *fkArena) endGroup(g int) {
	copy(a.groupEnds[g], a.counts)
}

func (a *fkArena) drainBlock(bi int, buf []byte, groups []logGroup) {
	base := bi << a.shift
	var cur uint32
	for g, grp := range groups {
		end := a.groupEnds[g][bi]
		if end == cur {
			continue
		}
		lg := grp.logP
		s := a.data[base+int(cur) : base+int(end)]
		i := 0
		for ; i+3 < len(s); i += 4 {
			buf[s[i]&blockMask] += lg
			buf[s[i+1]&blockMask] += lg
			buf[s[i+2]&blockMask] += lg
			buf[s[i+3]&blockMask] += lg
		}
		for ; i < len(s); i++ {
			buf[s[i]&blockMask] += lg
		}
		cur = end
	}
	if a.hasOver {
		for _, ov := range a.overflow[bi] {
			buf[ov.off] += ov.logP
		}
	}
}

func (a *fkArena) gatherHits(bi int, candBits *[blockSize / 64]uint64, candSlot *[blockSize]uint16, dst [][]int32) {
	base := bi << a.shift
	s := a.data[base : base+int(a.counts[bi])]
	for _, e := range s {
		off := e & blockMask
		if candBits[off>>6]&(1<<(off&63)) != 0 {
			k := candSlot[off]
			dst[k] = append(dst[k], int32(e>>blockShift))
		}
	}
	if a.hasOver {
		for _, ov := range a.overflow[bi] {
			if candBits[ov.off>>6]&(1<<(ov.off&63)) != 0 {
				k := candSlot[ov.off]
				dst[k] = append(dst[k], ov.idx)
			}
		}
	}
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

// transformRoot computes the lattice-transformed hit condition for ideal (p, r).
func transformRoot(p, r uint32, lq *LatticeQ) (lambda uint32, vert, all bool) {
	u1, u2 := lq.U1, lq.U2
	v1, v2 := lq.V1, lq.V2
	p64 := uint64(p)
	var ap, bp uint64
	if r == p {
		ap = modInt64(u2, p64)
		bp = modInt64(-v2, p64)
	} else {
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

// addFKHitsArena appends Franke-Kleinjung hits for p >= 2*I into the flat uint32 arena,
// skipping both-even (i, j) pairs.
func addFKHitsArena(p, lambda uint32, I, J int, idx int32, logP byte, a *fkArena) {
	W := 2 * I
	total := (J - 1) * W
	if total <= 0 {
		return
	}
	tag := uint32(idx) << blockShift
	counts := a.counts
	data := a.data
	stride := a.stride
	if idx >= 1<<(32-blockShift) {
		stride = 0
	}
	shift := a.shift

	if lambda == 0 {
		// Skip even j since i = 0 is even.
		for adjPos := I; adjPos < total; adjPos += 2 * W {
			bi := adjPos >> blockShift
			off := uint32(adjPos & blockMask)
			if n := counts[bi]; n < stride {
				data[bi<<shift|int(n)] = tag | off
				counts[bi] = n + 1
			} else {
				a.addOverflow(bi, idx, uint16(off), logP)
			}
		}
		return
	}

	a0, b0, a1, b1 := fkReduce(int64(p), int64(lambda), int64(W))
	x := I
	pos := I
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
				off := uint32(adjPos & blockMask)
				if n := counts[bi]; n < stride {
					data[bi<<shift|int(n)] = tag | off
					counts[bi] = n + 1
				} else {
					a.addOverflow(bi, idx, uint16(off), logP)
				}
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
			off := uint32(adjPos & blockMask)
			if n := counts[bi]; n < stride {
				data[bi<<shift|int(n)] = tag | off
				counts[bi] = n + 1
			} else {
				a.addOverflow(bi, idx, uint16(off), logP)
			}
		}
	}
}

func emitVertArena(p, I, J, W int, idx int32, logP byte, a *fkArena) {
	tag := uint32(idx) << blockShift
	counts := a.counts
	data := a.data
	stride := a.stride
	if idx >= 1<<(32-blockShift) {
		stride = 0
	}
	shift := a.shift
	for j := p; j < J; j += p {
		base := (j - 1) * W
		if j&1 == 0 {
			// Skip even i on even j lines where (x ^ I)&1 == 1 selects odd i.
			for x := (I & 1) ^ 1; x < W; x += 2 {
				pos := base + x
				bi := pos >> blockShift
				off := uint32(pos & blockMask)
				if n := counts[bi]; n < stride {
					data[bi<<shift|int(n)] = tag | off
					counts[bi] = n + 1
				} else {
					a.addOverflow(bi, idx, uint16(off), logP)
				}
			}
		} else {
			for x := 0; x < W; x++ {
				pos := base + x
				bi := pos >> blockShift
				off := uint32(pos & blockMask)
				if n := counts[bi]; n < stride {
					data[bi<<shift|int(n)] = tag | off
					counts[bi] = n + 1
				} else {
					a.addOverflow(bi, idx, uint16(off), logP)
				}
			}
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
	fastPoly := newFastNormPoly(poly)

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
						rThresh := lineBlockNormsRat(rb, &fastPoly, bStart, bLen, b, params.SlackRat+smallSlackRat)
						aThresh := lineBlockNormsAlg(ab, &fastPoly, bStart, bLen, b, params.SlackAlg+smallSlackAlg)

						sieveLineBlock(rb, ab, fb, ratSieve, algSieve, bStart, bLen, b, evenB)

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
				continue
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

func estimateExpectedFKHitsPerBlock(rfb []RatIdeal, afb []AlgIdeal, fkIdxs []int32, isRat bool) int {
	var invSum float64
	if isRat {
		for _, idx := range fkIdxs {
			invSum += 1.0 / float64(rfb[idx].P)
		}
	} else {
		for _, idx := range fkIdxs {
			invSum += 1.0 / float64(afb[idx].P)
		}
	}
	return int(0.8*float64(blockSize)*invSum) + 256
}

// sieveLattice runs the special-q lattice sieve with flat Franke-Kleinjung bucket arenas.
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

	ratGroups := buildLogGroupsRat(fb.RFB, ratFK)
	algGroups := buildLogGroupsAlg(fb.AFB, algFK)
	ratExpHits := estimateExpectedFKHitsPerBlock(fb.RFB, nil, ratFK, true)
	algExpHits := estimateExpectedFKHitsPerBlock(nil, fb.AFB, algFK, false)

	smallSlackRat := expectedSmallLogSumRat(fb, minSievePrime)
	smallSlackAlg := expectedSmallLogSumAlg(fb, minSievePrime)
	fastPoly := newFastNormPoly(poly)

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
			ratArena := newFKArena(numBlocks, len(ratGroups), ratExpHits)
			algArena := newFKArena(numBlocks, len(algGroups), algExpHits)
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

				ratArena.reset()
				algArena.reset()

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

				for g, grp := range ratGroups {
					for _, idx := range ratFK[grp.start:grp.end] {
						id := fb.RFB[idx]
						lam, vert, all := transformRoot(id.P, id.R, lq)
						if all {
							ratAll = append(ratAll, idx)
							ratAllLog += id.LogP
						} else if vert {
							emitVertArena(int(id.P), I, J, W, idx, grp.logP, ratArena)
						} else {
							addFKHitsArena(id.P, lam, I, J, idx, grp.logP, ratArena)
						}
					}
					ratArena.endGroup(g)
				}

				for g, grp := range algGroups {
					for _, idx := range algFK[grp.start:grp.end] {
						id := fb.AFB[idx]
						if id.P == lq.Q {
							continue
						}
						lam, vert, all := transformRoot(id.P, id.R, lq)
						if all {
							algAll = append(algAll, idx)
							algAllLog += id.LogP
						} else if vert {
							emitVertArena(int(id.P), I, J, W, idx, grp.logP, algArena)
						} else {
							addFKHitsArena(id.P, lam, I, J, idx, grp.logP, algArena)
						}
					}
					algArena.endGroup(g)
				}

				for bi := 0; bi < numBlocks && !stop.Load(); bi++ {
					bStart := bi * blockSize
					bEnd := bStart + blockSize
					if bEnd > total {
						bEnd = total
					}
					bLen := bEnd - bStart

					ab := algBuf[:bLen]
					aThresh := initLatticeNormsAlg(ab, &fastPoly, lq, bStart, bEnd, I, W, logQ, params.SlackAlg+smallSlackAlg, algAllLog)
					sieveBlockLinesAdd(ab, bStart, bEnd, I, W, algLines, algVerts)
					algArena.drainBlock(bi, ab, algGroups)

					if !hasByteGe(ab, aThresh) {
						continue
					}

					rb := ratBuf[:bLen]
					rThresh := initLatticeNormsRat(rb, &fastPoly, lq, bStart, bEnd, I, W, params.SlackRat+smallSlackRat, ratAllLog)
					sieveBlockLinesAdd(rb, bStart, bEnd, I, W, ratLines, ratVerts)
					ratArena.drainBlock(bi, rb, ratGroups)

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
					ratArena.gatherHits(bi, &candBits, &candSlot, ratCandHits)
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

					// Stage 2 scans algArena only for Stage 1 survivors.
					for len(algCandHits) < len(survivors) {
						algCandHits = append(algCandHits, nil)
					}
					for k, sv := range survivors {
						off := sv.off
						candBits[off/64] |= 1 << (uint(off) % 64)
						candSlot[off] = uint16(k)
						algCandHits[k] = algCandHits[k][:0]
					}
					algArena.gatherHits(bi, &candBits, &candSlot, algCandHits)
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

func addSieveStride(seg []byte, x, step, limit int, logP byte) {
	step2 := step << 1
	step3 := step2 + step
	step4 := step << 2
	for ; x+step3 < limit; x += step4 {
		seg[x] += logP
		seg[x+step] += logP
		seg[x+step2] += logP
		seg[x+step3] += logP
	}
	for ; x < limit; x += step {
		seg[x] += logP
	}
}

// sieveBlockLines initializes buf to allLog and sieves all line primes (< W) across [bStart, bEnd).
func sieveBlockLines(buf []byte, bStart, bEnd, I, W int, lines []linePrime, verts []vertPrime, allLog byte) {
	if allLog > 0 {
		for i := range buf {
			buf[i] = allLog
		}
	} else {
		clear(buf)
	}
	sieveBlockLinesAdd(buf, bStart, bEnd, I, W, lines, verts)
}

// sieveBlockLinesAdd adds line and vertical prime logP hits onto an initialized block.
func sieveBlockLinesAdd(buf []byte, bStart, bEnd, I, W int, lines []linePrime, verts []vertPrime) {
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
					addSieveStride(seg, iOddFirst, 2, W, logP)
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
					addSieveStride(seg, x, p2, W, logP)
				} else {
					addSieveStride(seg, x, p, W, logP)
				}
				off += lam
				if off >= p {
					off -= p
				}
			}
		}
		return
	}

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
					addSieveStride(seg, ((xFrom-I)&1)^1, 2, segLen, logP)
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
				addSieveStride(seg, x, p2, segLen, logP)
			} else {
				addSieveStride(seg, x, p, segLen, logP)
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

func sieveLineBlock(rb, ab []byte, fb *FactorBases, ratSieve, algSieve []int32, bStart int64, bLen int, b int64, evenB bool) {
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
}

// fastNormPoly stores scaled float64 coefficients of F and G for segment log2 norm evaluation.
type fastNormPoly struct {
	deg    int
	cf     []float64
	shiftF float64
	a1, a0 float64
	shiftG float64
}

func newFastNormPoly(poly *Poly) fastNormPoly {
	maxF := 0
	for _, c := range poly.Coeffs {
		if b := c.BitLen(); b > maxF {
			maxF = b
		}
	}
	shF := maxF - 200
	if shF < 0 {
		shF = 0
	}
	cf := make([]float64, poly.Degree+1)
	for i, c := range poly.Coeffs {
		cf[i] = scaledFloat(c, shF)
	}

	maxG := max(poly.A1.BitLen(), poly.A0.BitLen())
	shG := maxG - 200
	if shG < 0 {
		shG = 0
	}
	return fastNormPoly{
		deg:    poly.Degree,
		cf:     cf,
		shiftF: float64(shF),
		a1:     scaledFloat(poly.A1, shG),
		a0:     scaledFloat(poly.A0, shG),
		shiftG: float64(shG),
	}
}

func (fp *fastNormPoly) evalLog2F(a, b float64) (val, log2Abs float64) {
	d := fp.deg
	val = fp.cf[d]
	bPow := 1.0
	for k := d - 1; k >= 0; k-- {
		bPow *= b
		val = val*a + fp.cf[k]*bPow
	}
	if val == 0 {
		return 0, 0
	}
	return val, math.Log2(math.Abs(val)) + fp.shiftF
}

func (fp *fastNormPoly) evalLog2G(a, b float64) (val, log2Abs float64) {
	val = fp.a1*a - fp.a0*b
	if val == 0 {
		return 0, 0
	}
	return val, math.Log2(math.Abs(val)) + fp.shiftG
}

func segPairThreshold(l0, l1, deduct float64) uint8 {
	lLo, lHi := l0, l1
	if lLo > lHi {
		lLo, lHi = lHi, lLo
	}
	if lLo < lHi-5.0 {
		lLo = lHi - 5.0
	}
	t := int(math.Floor(lLo - deduct))
	if t <= 0 {
		return 0
	}
	if t > 240 {
		return 240
	}
	return uint8(t)
}

func initLatticeNormsAlg(buf []byte, fp *fastNormPoly, lq *LatticeQ, bStart, bEnd, I, W int, logQ, slack float64, allLog byte) int {
	return initLatticeNormsSide(buf, fp, lq, bStart, bEnd, I, W, logQ+slack, slack, allLog, true)
}

func initLatticeNormsRat(buf []byte, fp *fastNormPoly, lq *LatticeQ, bStart, bEnd, I, W int, slack float64, allLog byte) int {
	return initLatticeNormsSide(buf, fp, lq, bStart, bEnd, I, W, slack, slack, allLog, false)
}

func initLatticeNormsSide(buf []byte, fp *fastNormPoly, lq *LatticeQ, bStart, bEnd, I, W int, deduct, slack float64, allLog byte, isAlg bool) int {
	const coarseStep = 64
	const fineStep = 16
	var segEnd [1088]uint16
	var segT [1088]uint8
	nSeg := 0
	var tMax uint8

	u1, u2 := float64(lq.U1), float64(lq.U2)
	v1, v2 := float64(lq.V1), float64(lq.V2)
	eval := fp.evalLog2G
	if isAlg {
		eval = fp.evalLog2F
	}

	appendSeg := func(endOff int, t uint8) {
		if t > tMax {
			tMax = t
		}
		if nSeg > 0 && segT[nSeg-1] == t {
			segEnd[nSeg-1] = uint16(endOff)
			return
		}
		segEnd[nSeg] = uint16(endOff)
		segT[nSeg] = t
		nSeg++
	}

	for pos := bStart; pos < bEnd; {
		j := float64(1 + pos/W)
		xFrom := pos % W
		segLen := W - xFrom
		if bEnd-pos < segLen {
			segLen = bEnd - pos
		}
		xEnd := xFrom + segLen
		baseOff := pos - bStart
		aLine := j*v1 - float64(I)*u1
		bLine := j*v2 - float64(I)*u2

		v0, l0 := eval(aLine+float64(xFrom)*u1, bLine+float64(xFrom)*u2)
		for x0 := xFrom; x0 < xEnd; {
			x1 := x0 + coarseStep
			if x1 > xEnd {
				x1 = xEnd
			}
			xEval := x1
			if x1 == xEnd && xEnd > x0+1 {
				xEval = xEnd - 1
			}
			v1, l1 := eval(aLine+float64(xEval)*u1, bLine+float64(xEval)*u2)
			if (v0*v1 <= 0 || math.Abs(l1-l0) > 4.0) && x1-x0 > fineStep {
				va, la := v0, l0
				for xs := x0; xs < x1; {
					xe := xs + fineStep
					if xe >= x1 {
						appendSeg(baseOff+(x1-xFrom), segPairThreshold(la, l1, deduct))
						break
					}
					vb, lb := eval(aLine+float64(xe)*u1, bLine+float64(xe)*u2)
					appendSeg(baseOff+(xe-xFrom), segPairThreshold(la, lb, deduct))
					va, la = vb, lb
					_ = va
					xs = xe
				}
			} else {
				appendSeg(baseOff+(x1-xFrom), segPairThreshold(l0, l1, deduct))
			}
			v0, l0 = v1, l1
			x0 = x1
		}
		pos += segLen
	}

	return fillNormSegments(buf, segEnd[:nSeg], segT[:nSeg], int(tMax), slack, allLog)
}

func lineBlockNormsAlg(buf []byte, fp *fastNormPoly, bStart int64, bLen int, b int64, slack float64) int {
	return lineBlockNormsSide(buf, fp, bStart, bLen, b, slack, true)
}

func lineBlockNormsRat(buf []byte, fp *fastNormPoly, bStart int64, bLen int, b int64, slack float64) int {
	return lineBlockNormsSide(buf, fp, bStart, bLen, b, slack, false)
}

func lineBlockNormsSide(buf []byte, fp *fastNormPoly, bStart int64, bLen int, b int64, slack float64, isAlg bool) int {
	const coarseStep = 64
	const fineStep = 16
	var segEnd [1088]uint16
	var segT [1088]uint8
	nSeg := 0
	var tMax uint8

	bf := float64(b)
	a0f := float64(bStart)
	eval := fp.evalLog2G
	if isAlg {
		eval = fp.evalLog2F
	}

	appendSeg := func(endOff int, t uint8) {
		if t > tMax {
			tMax = t
		}
		if nSeg > 0 && segT[nSeg-1] == t {
			segEnd[nSeg-1] = uint16(endOff)
			return
		}
		segEnd[nSeg] = uint16(endOff)
		segT[nSeg] = t
		nSeg++
	}

	v0, l0 := eval(a0f, bf)
	for x0 := 0; x0 < bLen; {
		x1 := x0 + coarseStep
		if x1 > bLen {
			x1 = bLen
		}
		xEval := x1
		if x1 == bLen && bLen > x0+1 {
			xEval = bLen - 1
		}
		v1, l1 := eval(a0f+float64(xEval), bf)
		if (v0*v1 <= 0 || math.Abs(l1-l0) > 4.0) && x1-x0 > fineStep {
			la := l0
			for xs := x0; xs < x1; {
				xe := xs + fineStep
				if xe >= x1 {
					appendSeg(x1, segPairThreshold(la, l1, slack))
					break
				}
				_, lb := eval(a0f+float64(xe), bf)
				appendSeg(xe, segPairThreshold(la, lb, slack))
				la = lb
				xs = xe
			}
		} else {
			appendSeg(x1, segPairThreshold(l0, l1, slack))
		}
		v0, l0 = v1, l1
		x0 = x1
	}

	return fillNormSegments(buf, segEnd[:nSeg], segT[:nSeg], int(tMax), slack, 0)
}

func fillNormSegments(buf []byte, segEnd []uint16, segT []uint8, tMax int, slack float64, allLog byte) int {
	maxT := 245 - int(math.Ceil(slack))
	if maxT < int(allLog) {
		maxT = int(allLog)
	}
	effTMax := tMax
	if effTMax > maxT {
		effTMax = maxT
	}
	prev := 0
	for s, end16 := range segEnd {
		end := int(end16)
		initVal := int(allLog) + effTMax - int(segT[s])
		if initVal < int(allLog) {
			initVal = int(allLog)
		} else if initVal > 255 {
			initVal = 255
		}
		bVal := byte(initVal)
		seg := buf[prev:end]
		for i := range seg {
			seg[i] = bVal
		}
		prev = end
	}
	return effTMax
}

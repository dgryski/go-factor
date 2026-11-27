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
	blockShift     = 14 // 16 KB sieve blocks keep both sides L1-resident and fit offsets in 14 bits
	blockSize      = 1 << blockShift
	blockMask      = blockSize - 1
	maxStripBlocks = 256
	minSievePrime  = 29
	minSievePower  = 32
	maxSievePower  = 256
)

// oddModInv32 returns p^-1 mod 2^32 and floor((2^32-1)/p) for exact uint32
// divisibility (Granlund and Montgomery, PLDI 1994, Section 9).
func oddModInv32(p uint32) (inv, lim uint32) {
	inv = p
	for i := 0; i < 4; i++ {
		inv *= 2 - p*inv
	}
	return inv, ^uint32(0) / p
}

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
	bias   uint32
	pInv   uint32
	pLim   uint32
	logP   byte
}

type vertPrime struct {
	idx  int32
	p    uint32
	pInv uint32
	pLim uint32
	logP byte
}

type primePowerIdeal struct {
	pk   uint32
	p    uint32
	a, b uint32
	logP byte
}

type sievePower struct {
	pk     uint32
	lambda uint32
	logP   byte
	evenP  bool
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
	for (1<<shift) < want && shift < 17 {
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

// fkWalker holds per ideal Franke Kleinjung walk state resumed across L2 cache strips.
type fkWalker struct {
	pos     int
	step0   int
	step1   int
	idx     int32
	x       int32
	thresh0 int32
	thresh1 int32
	dx0     int32
	dx1     int32
	vertP   int32
}

func makeFKWalker(p, lam uint32, I, W int, idx int32, vert bool) fkWalker {
	if vert {
		return fkWalker{idx: idx, vertP: int32(p)}
	}
	if lam == 0 {
		return fkWalker{idx: idx, pos: I, step0: 2 * W}
	}
	a0, b0, a1, b1 := fkReduce(int64(p), int64(lam), int64(W))
	t0 := int32(-a0)
	t1 := int32(W - int(a1))
	dx0 := int32(a0)
	dx1 := int32(a1)
	s0 := int(b0)*W + int(a0)
	s1 := int(b1)*W + int(a1)
	pos := I
	x := int32(I)
	if x >= t0 {
		pos += s0
		x += dx0
	} else if x < t1 {
		pos += s1
		x += dx1
	} else {
		pos += s0 + s1
		x += dx0 + dx1
	}
	return fkWalker{
		pos:     pos,
		step0:   s0,
		step1:   s1,
		idx:     idx,
		x:       x,
		thresh0: t0,
		thresh1: t1,
		dx0:     dx0,
		dx1:     dx1,
		vertP:   -1,
	}
}

func emitVertStripArena(p, I, J, W int, idx int32, logP byte, stripBlock0, stripStartPos, stripEndPos int, a *fkArena) {
	if p >= J || stripStartPos >= stripEndPos {
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
	jMin := stripStartPos/W + 1
	jFirst := ((jMin + p - 1) / p) * p
	if jFirst < p {
		jFirst = p
	}
	wantOddX := (I & 1) ^ 1
	for j := jFirst; j < J; j += p {
		base := (j - 1) * W
		if base >= stripEndPos {
			break
		}
		xLo := 0
		if stripStartPos > base {
			xLo = stripStartPos - base
		}
		xHi := W
		if stripEndPos-base < W {
			xHi = stripEndPos - base
		}
		if j&1 == 0 {
			if (xLo & 1) != wantOddX {
				xLo++
			}
			for x := xLo; x < xHi; x += 2 {
				pos := base + x
				bi := (pos >> blockShift) - stripBlock0
				off := uint32(pos & blockMask)
				if n := counts[bi]; n < stride {
					data[bi<<shift|int(n)] = tag | off
					counts[bi] = n + 1
				} else {
					a.addOverflow(bi, idx, uint16(off), logP)
				}
			}
		} else {
			for x := xLo; x < xHi; x++ {
				pos := base + x
				bi := (pos >> blockShift) - stripBlock0
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

func advanceFKWalkers(walkers []fkWalker, I, J, W int, logP byte, stripBlock0, stripStartPos, stripEndPos int, a *fkArena) {
	counts := a.counts
	data := a.data
	baseStride := a.stride
	shift := a.shift
	stripLimit := stripEndPos + W
	pow2W := W&(W-1) == 0 && I&1 == 0
	evenMask := W | 1

	for wi := range walkers {
		w := &walkers[wi]
		if w.vertP > 0 {
			emitVertStripArena(int(w.vertP), I, J, W, w.idx, logP, stripBlock0, stripStartPos, stripEndPos, a)
			continue
		}
		tag := uint32(w.idx) << blockShift
		stride := baseStride
		if w.idx >= 1<<(32-blockShift) {
			stride = 0
		}
		if w.vertP == 0 {
			pos := w.pos
			step := w.step0
			for pos < stripEndPos {
				bi := (pos >> blockShift) - stripBlock0
				off := uint32(pos & blockMask)
				if n := counts[bi]; n < stride {
					data[bi<<shift|int(n)] = tag | off
					counts[bi] = n + 1
				} else {
					a.addOverflow(bi, w.idx, uint16(off), logP)
				}
				pos += step
			}
			w.pos = pos
			continue
		}

		pos := w.pos
		if pos >= stripLimit {
			continue
		}
		x := w.x
		t0, t1 := w.thresh0, w.thresh1
		dx0, dx1 := w.dx0, w.dx1
		dx01 := dx0 + dx1
		s0, s1 := w.step0, w.step1
		s01 := s0 + s1

		if pow2W {
			for pos < stripLimit {
				if pos&evenMask != 0 {
					adjPos := pos - W
					bi := (adjPos >> blockShift) - stripBlock0
					off := uint32(adjPos & blockMask)
					if n := counts[bi]; n < stride {
						data[bi<<shift|int(n)] = tag | off
						counts[bi] = n + 1
					} else {
						a.addOverflow(bi, w.idx, uint16(off), logP)
					}
				}
				if x >= t0 {
					pos += s0
					x += dx0
				} else if x < t1 {
					pos += s1
					x += dx1
				} else {
					pos += s01
					x += dx01
				}
			}
		} else {
			for pos < stripLimit {
				if ((int(x)-I)|(pos/W))&1 != 0 {
					adjPos := pos - W
					bi := (adjPos >> blockShift) - stripBlock0
					off := uint32(adjPos & blockMask)
					if n := counts[bi]; n < stride {
						data[bi<<shift|int(n)] = tag | off
						counts[bi] = n + 1
					} else {
						a.addOverflow(bi, w.idx, uint16(off), logP)
					}
				}
				if x >= t0 {
					pos += s0
					x += dx0
				} else if x < t1 {
					pos += s1
					x += dx1
				} else {
					pos += s01
					x += dx01
				}
			}
		}
		w.pos = pos
		w.x = x
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
	return sieveLatticeQs(poly, fb, params, sqs, onRel, shouldStop)
}

// probeLatticeYield runs a deterministic 30 lattice special q sample across three points
// in [QMin, QMax] and returns unique relations yielded per special q ideal.
func probeLatticeYield(poly *Poly, params *Params) float64 {
	qcMin := computeQCMin(params)
	fb := BuildFactorBases(poly, params.Br, params.Ba, qcMin, -1)
	return probeLatticeYieldWithFB(poly, fb, params)
}

func probeLatticeYieldWithFB(poly *Poly, fb *FactorBases, params *Params) float64 {
	qcMin := computeQCMin(params)
	qMid := uint64(math.Sqrt(float64(params.QMin) * float64(params.QMax)))
	qHi := params.QMin + (params.QMax-params.QMin)/2
	anchors := [3]uint64{params.QMin + 1000, qMid, qHi}

	var sample []LatticeQ
	for _, a := range anchors {
		if a < 3 || a >= qcMin {
			continue
		}
		sub := *params
		sub.QMin = a
		sub.QMax = a + 1500
		if sub.QMax >= qcMin {
			sub.QMax = qcMin - 1
		}
		sqs := buildSpecialQs(poly, fb, &sub)
		if len(sqs) > 10 {
			sqs = sqs[:10]
		}
		sample = append(sample, sqs...)
	}
	if len(sample) == 0 {
		return 0
	}
	seen := make(map[[2]int64]struct{}, len(sample)*24)
	sieveLatticeQsWithWorkers(poly, fb, params, sample, 6, func(r Relation) {
		seen[[2]int64{r.A, r.B}] = struct{}{}
	}, func() bool { return false })
	return float64(len(seen)) / float64(len(sample))
}

func sieveLatticeQs(poly *Poly, fb *FactorBases, params *Params, sqs []LatticeQ, onRel func(Relation), shouldStop func() bool) int {
	return sieveLatticeQsWithWorkers(poly, fb, params, sqs, 0, onRel, shouldStop)
}

func sieveLatticeQsWithWorkers(poly *Poly, fb *FactorBases, params *Params, sqs []LatticeQ, maxWorkers int, onRel func(Relation), shouldStop func() bool) int {
	return sieveLatticeQsStrip(poly, fb, params, sqs, maxWorkers, maxStripBlocks, onRel, shouldStop)
}

func sieveLatticeQsStrip(poly *Poly, fb *FactorBases, params *Params, sqs []LatticeQ, maxWorkers, maxBlocksPerStrip int, onRel func(Relation), shouldStop func() bool) int {
	if len(sqs) == 0 {
		return 0
	}

	I := params.IHalf
	J := params.JMax
	W := 2 * I
	total := (J - 1) * W
	numBlocks := (total + blockSize - 1) >> blockShift
	arenaBlocks := numBlocks
	if maxBlocksPerStrip > 0 && arenaBlocks > maxBlocksPerStrip {
		arenaBlocks = maxBlocksPerStrip
	}
	multiStrip := numBlocks > arenaBlocks

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
	ratPowIdeals := buildRatPowerIdeals(poly, fb, W)
	algPowIdeals := buildAlgPowerIdeals(poly, fb, W)

	smallSlackRat := expectedUnsievedLogSumRat(fb, minSievePrime, ratPowIdeals)
	smallSlackAlg := expectedUnsievedLogSumAlg(fb, minSievePrime, algPowIdeals)
	fastPoly := newFastNormPoly(poly)

	nw := Workers
	if nw <= 0 {
		nw = runtime.GOMAXPROCS(0)
	}
	if maxWorkers > 0 && nw > maxWorkers {
		nw = maxWorkers
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
		off              int
		i, j             int
		a, b             int64
		colStart, colEnd int
		ratLPs           [3]uint64
	}

	for w := 0; w < nw; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			var ratBuf, algBuf [blockSize]byte
			ratArena := newFKArena(arenaBlocks, len(ratGroups), ratExpHits)
			algArena := newFKArena(arenaBlocks, len(algGroups), algExpHits)
			scratch := newCofactorScratch(poly.Degree)

			ratLines := make([]linePrime, 0, len(ratLineIdx))
			var ratVerts []vertPrime
			var ratAll []int32
			ratPowers := make([]sievePower, 0, len(ratPowIdeals))
			algLines := make([]linePrime, 0, len(algLineIdx))
			var algVerts []vertPrime
			var algAll []int32
			algPowers := make([]sievePower, 0, len(algPowIdeals))

			var ratWalkers, algWalkers []fkWalker
			var ratWalkerEnds, algWalkerEnds []int
			if multiStrip {
				ratWalkers = make([]fkWalker, 0, len(ratFK))
				algWalkers = make([]fkWalker, 0, len(algFK))
				ratWalkerEnds = make([]int, len(ratGroups))
				algWalkerEnds = make([]int, len(algGroups))
			}

			var rawOffs []int
			var candPos []int
			var candI, candJ []int
			var candA, candB []int64
			var candBits [blockSize / 64]uint64
			var candSlot [blockSize]uint16
			var ratCandHits, algCandHits [][]int32
			var survivors []ratSurv
			var survColsBuf []int32

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
						inv, lim := uint32(0), uint32(0)
						if id.P&1 != 0 {
							inv, lim = oddModInv32(id.P)
						}
						ratVerts = append(ratVerts, vertPrime{idx: idx, p: id.P, pInv: inv, pLim: lim, logP: id.LogP})
					} else {
						inv, lim := oddModInv32(id.P)
						bias := uint32((I/int(id.P) + 1) * int(id.P))
						ratLines = append(ratLines, linePrime{idx: idx, p: id.P, lambda: lam, bias: bias, pInv: inv, pLim: lim, logP: id.LogP})
					}
				}
				ratPowers = transformPowerIdeals(ratPowers[:0], ratPowIdeals, lq)

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
						inv, lim := uint32(0), uint32(0)
						if id.P&1 != 0 {
							inv, lim = oddModInv32(id.P)
						}
						algVerts = append(algVerts, vertPrime{idx: idx, p: id.P, pInv: inv, pLim: lim, logP: id.LogP})
					} else {
						inv, lim := oddModInv32(id.P)
						bias := uint32((I/int(id.P) + 1) * int(id.P))
						algLines = append(algLines, linePrime{idx: idx, p: id.P, lambda: lam, bias: bias, pInv: inv, pLim: lim, logP: id.LogP})
					}
				}
				algPowers = transformPowerIdeals(algPowers[:0], algPowIdeals, lq)

				if multiStrip {
					ratWalkers = ratWalkers[:0]
					for g, grp := range ratGroups {
						for _, idx := range ratFK[grp.start:grp.end] {
							id := fb.RFB[idx]
							lam, vert, all := transformRoot(id.P, id.R, lq)
							if all {
								ratAll = append(ratAll, idx)
								ratAllLog += id.LogP
							} else {
								ratWalkers = append(ratWalkers, makeFKWalker(id.P, lam, I, W, idx, vert))
							}
						}
						ratWalkerEnds[g] = len(ratWalkers)
					}
					algWalkers = algWalkers[:0]
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
							} else {
								algWalkers = append(algWalkers, makeFKWalker(id.P, lam, I, W, idx, vert))
							}
						}
						algWalkerEnds[g] = len(algWalkers)
					}
				} else {
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
				}

				for stripBlock0 := 0; stripBlock0 < numBlocks && !stop.Load(); stripBlock0 += arenaBlocks {
					stripBlocks := min(arenaBlocks, numBlocks-stripBlock0)
					if multiStrip {
						stripStartPos := stripBlock0 * blockSize
						stripEndPos := min(total, (stripBlock0+stripBlocks)*blockSize)
						ratArena.reset()
						for g, grp := range ratGroups {
							wStart := 0
							if g > 0 {
								wStart = ratWalkerEnds[g-1]
							}
							advanceFKWalkers(ratWalkers[wStart:ratWalkerEnds[g]], I, J, W, grp.logP, stripBlock0, stripStartPos, stripEndPos, ratArena)
							ratArena.endGroup(g)
						}
						algArena.reset()
						for g, grp := range algGroups {
							wStart := 0
							if g > 0 {
								wStart = algWalkerEnds[g-1]
							}
							advanceFKWalkers(algWalkers[wStart:algWalkerEnds[g]], I, J, W, grp.logP, stripBlock0, stripStartPos, stripEndPos, algArena)
							algArena.endGroup(g)
						}
					}

					for localBI := 0; localBI < stripBlocks && !stop.Load(); localBI++ {
						bi := stripBlock0 + localBI
						bStart := bi * blockSize
						bEnd := bStart + blockSize
						if bEnd > total {
							bEnd = total
						}
						bLen := bEnd - bStart

						ab := algBuf[:bLen]
						aThresh := initLatticeNormsAlg(ab, &fastPoly, lq, bStart, bEnd, I, W, logQ, params.SlackAlg+smallSlackAlg, algAllLog)
						sieveBlockLinesAdd(ab, bStart, bEnd, I, W, algLines, algVerts)
						sieveBlockPowersAdd(ab, bStart, bEnd, I, W, algPowers)
						algArena.drainBlock(localBI, ab, algGroups)

						if !hasByteGe(ab, aThresh) {
							continue
						}

						rb := ratBuf[:bLen]
						rThresh := initLatticeNormsRat(rb, &fastPoly, lq, bStart, bEnd, I, W, params.SlackRat+smallSlackRat, ratAllLog)
						sieveBlockLinesAdd(rb, bStart, bEnd, I, W, ratLines, ratVerts)
						sieveBlockPowersAdd(rb, bStart, bEnd, I, W, ratPowers)
						ratArena.drainBlock(localBI, rb, ratGroups)

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
						ratArena.gatherHits(localBI, &candBits, &candSlot, ratCandHits)
						batchRatLines := len(candPos) >= 10
						if batchRatLines {
							gatherBlockLineHits(bStart, bEnd, I, W, ratLines, ratVerts, &candBits, &candSlot, ratCandHits)
						}
						for _, off := range candPos {
							candBits[off/64] = 0
						}

						var ratLinesForCof []linePrime
						var ratVertsForCof []vertPrime
						if !batchRatLines {
							ratLinesForCof = ratLines
							ratVertsForCof = ratVerts
						}

						survivors = survivors[:0]
						survColsBuf = survColsBuf[:0]
						for k := range candPos {
							ratCols, ratLPs, ok := cofactorRatLattice(
								poly, fb, params,
								candI[k], candJ[k], candA[k], candB[k],
								ratSmall, ratLinesForCof, ratVertsForCof, ratAll, ratCandHits[k],
								scratch,
							)
							if ok {
								cStart := len(survColsBuf)
								survColsBuf = append(survColsBuf, ratCols...)
								survivors = append(survivors, ratSurv{
									off:      candPos[k],
									i:        candI[k],
									j:        candJ[k],
									a:        candA[k],
									b:        candB[k],
									colStart: cStart,
									colEnd:   len(survColsBuf),
									ratLPs:   ratLPs,
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
						algArena.gatherHits(localBI, &candBits, &candSlot, algCandHits)
						batchAlgLines := len(survivors) >= 10
						if batchAlgLines {
							gatherBlockLineHits(bStart, bEnd, I, W, algLines, algVerts, &candBits, &candSlot, algCandHits)
						}
						for _, sv := range survivors {
							candBits[sv.off/64] = 0
						}

						var algLinesForCof []linePrime
						var algVertsForCof []vertPrime
						if !batchAlgLines {
							algLinesForCof = algLines
							algVertsForCof = algVerts
						}

						for k := range survivors {
							sv := &survivors[k]
							rel, ok := cofactorAlgLattice(
								poly, fb, params,
								sv.i, sv.j, sv.a, sv.b,
								survColsBuf[sv.colStart:sv.colEnd], sv.ratLPs,
								algSmall, algLinesForCof, algVertsForCof, algAll, algCandHits[k],
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
	if bStart%W == 0 && bEnd%W == 0 && W&63 == 0 {
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

// gatherBlockLineHits scatters line and vertical prime hits in [bStart, bEnd) onto candidates in candBits.
func gatherBlockLineHits(bStart, bEnd, I, W int, lines []linePrime, verts []vertPrime, candBits *[blockSize / 64]uint64, candSlot *[blockSize]uint16, dst [][]int32) {
	j0 := 1 + bStart/W
	if bStart%W == 0 && bEnd%W == 0 && W&63 == 0 {
		nLines := (bEnd - bStart) / W
		for _, vp := range verts {
			vpP := int(vp.p)
			rem := j0 % vpP
			firstL := 0
			if rem != 0 {
				firstL = vpP - rem
			}
			for l := firstL; l < nLines; l += vpP {
				base := l * W
				wStart, wEnd := base>>6, (base+W)>>6
				for w := wStart; w < wEnd; w++ {
					m := candBits[w]
					for m != 0 {
						off := w<<6 | bits.TrailingZeros64(m)
						dst[candSlot[off]] = append(dst[candSlot[off]], vp.idx)
						m &= m - 1
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
			idx := lp.idx
			off := int((j064*uint64(lp.lambda) + i64) % uint64(lp.p))
			for l := 0; l < nLines; l++ {
				base := l * W
				x := off
				step := p
				if (j0+l)&1 == 0 {
					if (x-I)&1 == 0 {
						x += p
					}
					step = p2
				}
				for ; x < W; x += step {
					bOff := base + x
					if candBits[bOff>>6]&(1<<(uint(bOff)&63)) != 0 {
						dst[candSlot[bOff]] = append(dst[candSlot[bOff]], idx)
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

	for _, vp := range verts {
		vpP := int(vp.p)
		for pos := bStart; pos < bEnd; {
			j := 1 + pos/W
			xFrom := pos % W
			segLen := W - xFrom
			if bEnd-pos < segLen {
				segLen = bEnd - pos
			}
			if j%vpP == 0 {
				for x := 0; x < segLen; x++ {
					bOff := pos - bStart + x
					if candBits[bOff>>6]&(1<<(uint(bOff)&63)) != 0 {
						dst[candSlot[bOff]] = append(dst[candSlot[bOff]], vp.idx)
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
		lam := int(lp.lambda)
		idx := lp.idx
		off := int((uint64(j0)*uint64(lp.lambda) + i64) % uint64(lp.p))
		for pos := bStart; pos < bEnd; {
			xFrom := pos % W
			segLen := W - xFrom
			if bEnd-pos < segLen {
				segLen = bEnd - pos
			}
			xGlobal := off
			if xGlobal < xFrom {
				xGlobal += ((xFrom - xGlobal + p - 1) / p) * p
			}
			for x := xGlobal - xFrom; x < segLen; x += p {
				bOff := pos - bStart + x
				if candBits[bOff>>6]&(1<<(uint(bOff)&63)) != 0 {
					dst[candSlot[bOff]] = append(dst[candSlot[bOff]], idx)
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

// buildRatPowerIdeals builds prime power sieve entries p^k in [minSievePower, min(W-1, maxSievePower)]
// for rational primes p.
func buildRatPowerIdeals(poly *Poly, fb *FactorBases, W int) []primePowerIdeal {
	limit := maxSievePower
	if W-1 < limit {
		limit = W - 1
	}
	if limit < minSievePower {
		return nil
	}
	var out []primePowerIdeal
	var bigM, rTmp big.Int
	for _, id := range fb.RFB {
		p := uint64(id.P)
		if p*p > uint64(limit) {
			break
		}
		for pk := p * p; pk <= uint64(limit); pk *= p {
			if int(pk) < minSievePower {
				continue
			}
			bigM.SetUint64(pk)
			a1m := rTmp.Mod(poly.A1, &bigM).Uint64()
			a0m := rTmp.Mod(poly.A0, &bigM).Uint64()
			g := arith.Gcd64(arith.Gcd64(a1m, a0m), pk)
			if g > 1 {
				continue
			}
			out = append(out, primePowerIdeal{
				pk:   uint32(pk),
				p:    id.P,
				a:    uint32(a1m),
				b:    uint32(a0m),
				logP: id.LogP,
			})
		}
	}
	return out
}

// buildAlgPowerIdeals lifts simple roots of f mod p to p^k in [minSievePower, min(W-1, maxSievePower)].
func buildAlgPowerIdeals(poly *Poly, fb *FactorBases, W int) []primePowerIdeal {
	limit := maxSievePower
	if W-1 < limit {
		limit = W - 1
	}
	if limit < minSievePower {
		return nil
	}
	var out []primePowerIdeal
	for _, id := range fb.AFB {
		p := uint64(id.P)
		if p*p > uint64(limit) {
			break
		}
		if id.R == id.P {
			continue
		}
		fp := poly.CoeffsModP(p)
		dfp := arith.PolyDerivModP(fp, p)
		deriv := arith.PolyEvalModP(dfp, uint64(id.R), p)
		if deriv == 0 {
			continue
		}
		invDeriv := arith.InvMod(deriv, p)
		curR := uint64(id.R)
		curPK := p
		for curPK*p <= uint64(limit) {
			nextPK := curPK * p
			fMod := poly.CoeffsModP(nextPK)
			val := arith.PolyEvalModP(fMod, curR, nextPK)
			step := arith.MulMod64((p-(val/curPK)%p)%p, invDeriv, p)
			curR += step * curPK
			curPK = nextPK
			if int(curPK) >= minSievePower {
				out = append(out, primePowerIdeal{
					pk:   uint32(curPK),
					p:    id.P,
					a:    1,
					b:    uint32(curR),
					logP: id.LogP,
				})
			}
		}
	}
	return out
}

func transformPowerIdeals(dst []sievePower, ideals []primePowerIdeal, lq *LatticeQ) []sievePower {
	u1, u2 := lq.U1, lq.U2
	v1, v2 := lq.V1, lq.V2
	for _, id := range ideals {
		if id.p == lq.Q {
			continue
		}
		pk := uint64(id.pk)
		a := uint64(id.a)
		b := uint64(id.b)
		u1m, u2m := modInt64(u1, pk), modInt64(u2, pk)
		v1m, v2m := modInt64(v1, pk), modInt64(v2, pk)
		au1 := (a * u1m) % pk
		bu2 := (b * u2m) % pk
		var ap uint64
		if au1 >= bu2 {
			ap = au1 - bu2
		} else {
			ap = pk - (bu2 - au1)
		}
		if ap%uint64(id.p) == 0 {
			continue
		}
		invAp, ok := arith.InvMod64(ap, pk)
		if !ok {
			continue
		}
		bv2 := (b * v2m) % pk
		av1 := (a * v1m) % pk
		var bp uint64
		if bv2 >= av1 {
			bp = bv2 - av1
		} else {
			bp = pk - (av1 - bv2)
		}
		lam := uint32((bp * invAp) % pk)
		dst = append(dst, sievePower{
			pk:     id.pk,
			lambda: lam,
			logP:   id.logP,
			evenP:  id.p == 2,
		})
	}
	return dst
}

// sieveBlockPowersAdd accumulates logP hits for prime powers pk < W, respecting even j parity for 2^k.
func sieveBlockPowersAdd(buf []byte, bStart, bEnd, I, W int, powers []sievePower) {
	if len(powers) == 0 {
		return
	}
	if bStart%W == 0 && bEnd%W == 0 && W&63 == 0 {
		j0 := 1 + bStart/W
		nLines := (bEnd - bStart) / W
		i64 := uint64(I)
		j064 := uint64(j0)
		for k := range powers {
			sp := powers[k]
			pk := int(sp.pk)
			pk2 := pk << 1
			lam := int(sp.lambda)
			logP := sp.logP
			off := int((j064*uint64(sp.lambda) + i64) % uint64(sp.pk))
			if sp.evenP {
				for l := 0; l < nLines; l++ {
					if (j0+l)&1 != 0 || (off-I)&1 != 0 {
						addSieveStride(buf[l*W:(l+1)*W], off, pk, W, logP)
					}
					off += lam
					if off >= pk {
						off -= pk
					}
				}
			} else {
				for l := 0; l < nLines; l++ {
					seg := buf[l*W : (l+1)*W]
					x := off
					if (j0+l)&1 == 0 {
						if (x-I)&1 == 0 {
							x += pk
						}
						addSieveStride(seg, x, pk2, W, logP)
					} else {
						addSieveStride(seg, x, pk, W, logP)
					}
					off += lam
					if off >= pk {
						off -= pk
					}
				}
			}
		}
		return
	}

	j0 := 1 + bStart/W
	i64 := uint64(I)
	for k := range powers {
		sp := powers[k]
		pk := int(sp.pk)
		pk2 := pk << 1
		lam := int(sp.lambda)
		logP := sp.logP
		off := int((uint64(j0)*uint64(sp.lambda) + i64) % uint64(sp.pk))
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
				xGlobal += ((xFrom - xGlobal + pk - 1) / pk) * pk
			}
			x := xGlobal - xFrom
			if sp.evenP {
				if j&1 != 0 || (xGlobal-I)&1 != 0 {
					addSieveStride(seg, x, pk, segLen, logP)
				}
			} else if j&1 == 0 {
				if (xGlobal-I)&1 == 0 {
					x += pk
				}
				addSieveStride(seg, x, pk2, segLen, logP)
			} else {
				addSieveStride(seg, x, pk, segLen, logP)
			}
			if xFrom+segLen == W {
				off += lam
				if off >= pk {
					off -= pk
				}
			}
			pos += segLen
		}
	}
}

func expectedUnsievedLogSumRat(fb *FactorBases, cut uint32, sievedPows []primePowerIdeal) float64 {
	sieved := make(map[uint32]bool, len(sievedPows))
	for _, sp := range sievedPows {
		sieved[sp.pk] = true
	}
	sum := 0.0
	for _, id := range fb.RFB {
		p := id.P
		if p >= cut && p*p > maxSievePower {
			break
		}
		if p < cut {
			sum += float64(id.LogP) / float64(p+1)
		}
		for pk := uint64(p) * uint64(p); pk <= maxSievePower; pk *= uint64(p) {
			if !sieved[uint32(pk)] {
				sum += float64(id.LogP) * float64(p) / (float64(p+1) * float64(pk))
			}
		}
	}
	return sum
}

func expectedUnsievedLogSumAlg(fb *FactorBases, cut uint32, sievedPows []primePowerIdeal) float64 {
	sievedCount := make(map[uint32]int, len(sievedPows))
	for _, sp := range sievedPows {
		sievedCount[sp.pk]++
	}
	sum := 0.0
	for _, id := range fb.AFB {
		p := id.P
		if p >= cut {
			break
		}
		sum += float64(id.LogP) / float64(p+1)
	}
	// Add unsieved prime powers below minSievePower (such as 4, 8, 16, 9, 27, 25).
	for i := 0; i < len(fb.AFB); {
		p := fb.AFB[i].P
		if uint64(p)*uint64(p) >= minSievePower {
			break
		}
		logP := fb.AFB[i].LogP
		nRoots := 0
		j := i
		for j < len(fb.AFB) && fb.AFB[j].P == p {
			nRoots++
			j++
		}
		for pk := uint64(p) * uint64(p); pk < minSievePower; pk *= uint64(p) {
			unsieved := nRoots - sievedCount[uint32(pk)]
			if unsieved > 0 {
				sum += float64(unsieved) * float64(logP) * float64(p) / (float64(p+1) * float64(pk))
			}
		}
		i = j
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
	const coarseStep = 128
	const fineStep = 16
	var segEnd [1088]uint16
	var segT [1088]uint8
	nSeg := 0
	var tMax uint8

	u1, u2 := float64(lq.U1), float64(lq.U2)
	v1, v2 := float64(lq.V1), float64(lq.V2)

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

		var v0, l0 float64
		if isAlg {
			v0, l0 = fp.evalLog2F(aLine+float64(xFrom)*u1, bLine+float64(xFrom)*u2)
		} else {
			v0, l0 = fp.evalLog2G(aLine+float64(xFrom)*u1, bLine+float64(xFrom)*u2)
		}
		for x0 := xFrom; x0 < xEnd; {
			x1 := x0 + coarseStep
			if x1 > xEnd {
				x1 = xEnd
			}
			xEval := x1
			if x1 == xEnd && xEnd > x0+1 {
				xEval = xEnd - 1
			}
			var v1, l1 float64
			if isAlg {
				v1, l1 = fp.evalLog2F(aLine+float64(xEval)*u1, bLine+float64(xEval)*u2)
			} else {
				v1, l1 = fp.evalLog2G(aLine+float64(xEval)*u1, bLine+float64(xEval)*u2)
			}
			if (v0*v1 <= 0 || math.Abs(l1-l0) > 3.0) && x1-x0 > fineStep {
				la := l0
				for xs := x0; xs < x1; {
					xe := xs + fineStep
					if xe >= x1 {
						t := segPairThreshold(la, l1, deduct)
						if t > tMax {
							tMax = t
						}
						endOff := uint16(baseOff + (x1 - xFrom))
						if nSeg > 0 && segT[nSeg-1] == t {
							segEnd[nSeg-1] = endOff
						} else {
							segEnd[nSeg], segT[nSeg] = endOff, t
							nSeg++
						}
						break
					}
					var lb float64
					if isAlg {
						_, lb = fp.evalLog2F(aLine+float64(xe)*u1, bLine+float64(xe)*u2)
					} else {
						_, lb = fp.evalLog2G(aLine+float64(xe)*u1, bLine+float64(xe)*u2)
					}
					t := segPairThreshold(la, lb, deduct)
					if t > tMax {
						tMax = t
					}
					endOff := uint16(baseOff + (xe - xFrom))
					if nSeg > 0 && segT[nSeg-1] == t {
						segEnd[nSeg-1] = endOff
					} else {
						segEnd[nSeg], segT[nSeg] = endOff, t
						nSeg++
					}
					la = lb
					xs = xe
				}
			} else {
				t := segPairThreshold(l0, l1, deduct)
				if t > tMax {
					tMax = t
				}
				endOff := uint16(baseOff + (x1 - xFrom))
				if nSeg > 0 && segT[nSeg-1] == t {
					segEnd[nSeg-1] = endOff
				} else {
					segEnd[nSeg], segT[nSeg] = endOff, t
					nSeg++
				}
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
		seg := buf[prev:end]
		if initVal == 0 {
			clear(seg)
		} else {
			bVal := byte(initVal)
			wVal := uint64(bVal) * 0x0101010101010101
			n8 := len(seg) &^ 7
			for i := 0; i < n8; i += 8 {
				binary.LittleEndian.PutUint64(seg[i:], wVal)
			}
			for i := n8; i < len(seg); i++ {
				seg[i] = bVal
			}
		}
		prev = end
	}
	return effTMax
}

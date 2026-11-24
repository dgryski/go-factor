# Plan: the quadratic sieve

## Current state (2026-11-24)

Stages 1–3 of the original plan are done, except for a sparse matrix solver.
The code is in package `qs` (`qs/qs.go`). The AI drafts that SIQS, the
multiplier, large primes and the faster elimination started from are in
`notes/_qs-drafts/`. The leading underscore keeps `go build ./...` from
compiling them. See "Fixed from the drafts" below for what was wrong in them.

- `qs.QuadraticSieve(n) (primes, composites []*big.Int)` recurses on the two
  halves found by `split`. `qs.Logf`, if set, gets progress messages.
- `split`:
  1. trial divides up to B
  2. picks a Knuth–Schroeppel multiplier k (odd and square-free, < 100)
  3. builds the factor base over kN: 2, the primes dividing k, and the primes
     with Legendre symbol 1
  4. sieves, using SIQS for n ≥ 60 bits and a single polynomial below that
  5. runs `findDependencies` and `tryDependency`
- **Single polynomial** (`singlePoly`): Y = (x + ⌈√kN⌉)² − kN, with
  B = 3·L(N)^½ (minimum 500) and M = max(B/3, 256)·B. It uses no large primes,
  and the threshold slack is 1.75·log₂(p_max).
- **SIQS** (`siqs`):
  - Y = ((Ax + B)² − kN)/A, x ∈ [−M, M).
  - B = L(N)^0.46 (minimum 1000). M is 32K, 64K or 128K depending on B.
  - A is the product of n primes near 2000 (`chooseA`): n−1 random ones within
    a factor of 1.5 of the ideal size, plus the one that gets closest to the
    target √(2kN)/M.
  - B_j are built by CRT, the 2^(n−1) polynomials are taken in Gray code order,
    and the roots are moved by precomputed 2·B_v·A⁻¹ deltas.
  - Large primes go up to 64·p_max, with slack 1.25·log₂(p_max) plus half of
    log₂(64·p_max).
- **Shared by both:**
  - 32 KB blocks.
  - Primes below 30 aren't sieved, and the threshold allows for their
    contribution.
  - Primes from 16384 up (`bucketCut`) are bucket sieved. `fillBuckets` walks
    each such prime's hits across the polynomial's whole interval once, and
    files each as log₂p << 48 | prime index << 16 | offset in its block's
    bucket. Each block then adds its bucket's entries instead of looping over
    every root.
  - The threshold is per block, from the exact |Y| at the block's ends, and at
    the vertex if it lies inside the block.
  - `tryFactorize` checks primes below 16384 from a packed (p, r1, r2) table
    rebuilt for each polynomial. For larger primes it divides only by the
    bucket entries at the candidate's offset.
- **Relations:**
  - Each relation stores its factors sparsely, as column indices repeated by
    multiplicity, rather than a dense exponent vector.
  - Relations whose ±X was already seen are dropped. Different polynomials do
    find the same X, and before this check those copies made about 1% of
    70-bit SIQS runs fail with only trivial dependencies.
- **`findDependencies`** (`qs/matrix.go`):
  - Structured Gaussian elimination first. It deletes rows with a column no
    other row has, merges the two rows of each weight-2 column, and pivots on
    light columns.
  - The pivot limits scale with the column count (`pivotLimits`): weight
    numCols/1000 between 6 and 50, and fill numCols/16 between 64 and 4000.
    Fixed high limits made 60–120-bit inputs up to 2× slower.
  - Then it deletes the heaviest rows while there are more than 64 surplus
    rows (`keepSurplus`), half the surplus at a time, and repeats. Many primes
    never occur in a relation, so without this the 80-digit matrix kept about
    2,450 surplus rows.
  - Each XOR is recorded as a node, so dependencies of the reduced matrix
    trace back to the original relations.
  - Block Lanczos (`qs/lanczos.go`, from `notes/_qs-drafts4/`) then solves
    the reduced matrix if it has at least 1000 columns. It uses
    A = B·P·Bᵀ with a random symmetric invertible P (see "Block Lanczos"
    below). Dense GF(2) elimination solves smaller matrices, and is the
    fallback if block Lanczos finds fewer than 16 independent dependencies.
  - `TestFindDependencies` checks the dependencies against full dense
    elimination: the count, that each sums to zero, and that they're
    independent.
  - The structured elimination and resieving started from the drafts in
    `notes/_qs-drafts2/`. See "Block Lanczos" below for the draft solver.
  - `factor -qs-save=DIR` saves relations after sieving, and
    `go test ./qs -run SavedRelations -relations=FILE` reruns the matrix and
    square root steps on them. Relations for the 80-digit number and its
    209-bit cofactor are in `~/qs-relations/`.
- **`factor` command:**
  - `-methods=trial,qs` runs only the listed methods, in the fixed order.
  - `-v` turns on QS logging.
  - QS is skipped above 220 bits unless `-methods` names it.
- **Tests:**
  - Tonelli–Shanks against brute force.
  - `TestRelations`, which checks X² ≡ Y (mod n) for every relation from both
    sieves.
  - Random semiprimes from 20 to 150 bits.
  - Benchmarks for `split` and `findDependencies`.
  - The shared test helpers are in `internal/testutil`.

Time for `split` on random balanced semiprimes, from one sample per size:

| bits | first sieve (single poly) | SIQS | SIQS with buckets |
|------|------|-------|------|
| 60   | 0.82 ms | 0.55 ms | 0.49 ms |
| 80   | 5 ms    | 0.93 ms | 0.78 ms |
| 100  | 39 ms   | 2.6 ms  | 2.2 ms |
| 120  | 217 ms  | 7.0 ms  | 6.0 ms |
| 160  | –       | 154 ms  | 143 ms |
| 200  | –       | 2.86 s  | 1.89 s |
| 219  | –       | 11.0 s  | 6.18 s |
| 240  | –       | 40.4 s  | 21.6 s |

The 160–240-bit rows are fixed semiprimes (seed 7), and both columns collected
the same relations.

### Tuning results

These were measured, not guessed. Redo them if the candidate check gets faster.

- **Large primes:**
  - With the threshold lowered by the full log₂ of the large prime bound, they
    made SIQS 2–3× slower up to 220 bits. The extra candidates cost more than
    the partials gave.
  - Lowered by half, they're about 25% faster at 200 bits. Results are flat
    for multipliers 64–128 and slack 1.0–1.25.
  - With a single polynomial, large primes are always slower.
- **SIQS below about 60 bits:** the pool of A primes is too small, so it runs
  out of distinct A values before it has enough relations.
- **`tryFactorize`:** it called `temp.Cmp(one)` for every root, which took
  half the run time. It now checks only after a division.
- **`bucketCut`:** 16384 was fastest among 2048–65536 at 219 and 240 bits.
  2048 was 15% slower and 65536 10% slower.
- **Resieving** (`notes/_qa-drafts5/`) didn't pay and isn't used.
  - It walks each prime from 2048 up back through every block with a
    candidate. At 200 bits that's about 28 µs per block, against about 4.7 µs
    per candidate for the old full scan, so it breaks even at about 7
    candidates per block. Blocks had 3.4 at 200 bits and 1.5 at 240.
  - Looser thresholds gave 40–85 candidates per block but made runs slower,
    because the extra candidates nearly all leave cofactors too big for one
    large prime.
  - Bucket sieving gets most of the benefit instead, because a candidate's
    large primes come from its block's bucket.
  - It may pay again with double large primes, where those candidates become
    partials.

### Fixed from the drafts

- Exponent vectors were dense: about 400 KB allocated per candidate, and about
  20 GB of relations at 80 digits.
- The factor base left out the primes dividing k, and primes dividing N were
  never caught.
- Even multipliers were allowed, which breaks the p = 2 roots.
- The block threshold always included the vertex value kN/A, even when the
  vertex was outside the block.
- `tryDependency` only tried gcd(x − y), not x + y.
- The A = 1 fallback for small n used too short a range.

## Reference: an 80-digit number

73129071533331481181630024590293633952610077040929361500801632604572796624590389
(266 bits) = p18 · p23 · p41:

- 104194448528493551
- 11494857165741861450481
- 61057905163708884083219136437717059663819

GMP-ECM with B1 = 11000 finds p18 within a minute. `factor/msieve.log` has
msieve 1.53 factoring it by SIQS in 1 min 36 s:

- multiplier 5
- 12 blocks of 32 KB per polynomial
- factor base bound 1,301,057 (49,858 primes)
- large prime bound 1.3·10⁸
- 49,954 relations: 26,246 full and 23,708 combined from 267,452 partials
- a sparse matrix solved by block Lanczos

`factor -v -methods=qs` factors it completely in 2 min 42 s, about 1.7×
msieve's time:

- SIQS chose k = 5 and B = 1,570,053 (59,370 primes, M = 128K).
- It collected 59,391 relations in 2 min 34 s, against 5 min before bucket
  sieving.
- Structured elimination reduces the matrix from 59,391×59,371 to
  15,649×15,585 in 0.7 s, and block Lanczos solves that in 1.5 s. The whole
  matrix step takes 2.3 s.
- Earlier versions: 13.6 s with dense elimination on the same reduced matrix,
  and 2 min with the first version of the elimination, which reached
  37,301×34,845.
- The first split found p18. QS then split the 209-bit cofactor into p23 and
  p41 in 14 s.

Without structured elimination, the dense matrix needed about 900 MB and was
still running after 4 minutes.

Pivot limits tried on the saved relations, with the reduced size each gave:

| weight ≤ | fill ≤ | reduced to |
|---|---|---|
| 6 | 64 | 30,939×30,875 |
| 10 | 256 | 20,455×20,391 |
| 20 | 1000 | 16,887×16,823 |
| 50 | 4000 | 15,503×15,439 |
| 1000 | 16000 | 16,909×16,845 (23 s) |

### Block Lanczos

The first two drafts (`notes/_qs-drafts2/`, `notes/_qs-drafts3/`) used
A = B·Bᵀ. Run unchanged apart from inputs and outputs, they gave:

| Matrix | Block Lanczos | Dense |
|---|---|---|
| reduced 80-digit, 15,649×15,585 | 64 valid, independent dependencies in 1.36 s | 68 in 12.4 s |
| reduced 209-bit, 3,707×3,643 | 64 independent of 110 valid, 32 ms | 66 in 153 ms |
| unreduced 80-digit, 59,391×59,371 | none, 5.4 s | – |

Its last step checks every result, so a failure means fewer dependencies,
never wrong ones.

**Why it fails on unreduced matrices.** It solves A·X = A·Y with A = B·Bᵀ,
then looks for dependencies among the 64 columns of X − Y (and those of the
last V, which is zero here). X − Y lies in ker(A). Any part of ker(A) outside
ker(Bᵀ) costs one of those 64 columns. That part has dimension
rank(B) − rank(B·Bᵀ), from row combinations orthogonal to every row.

On the 209-bit relations:

| Matrix | rank(B) | rank(B·Bᵀ) | gap | dependencies found |
|---|---|---|---|---|
| reduced, 3,707×3,643 | 3,641 | 3,641 | 0 | 65 |
| unreduced, 10,328×10,308 | 9,999 | 9,977 | 22 | 42 |

So 64 − gap = 42. The unreduced 80-digit matrix gives none, so its gap is
probably at least 64, but that's an inference. Computing rank(B·Bᵀ) at 59k is
too slow.

Ruled out:
- **Surplus rows:** trimming the unreduced 80-digit matrix to 64 surplus rows
  (57,089×57,025) still gives none.
- **Breakdown:** every failing run ends with V = 0 after about rows/64
  iterations, the expected length.
- **Size:** reduced matrices up to 30,939×30,875, needing 489 iterations,
  work.

The structured elimination's matrices had a gap of 0 in every case tried.
`notes/_qs-drafts3/` also doesn't compile as written:
`int64(0x9e3779b97f4a7c15)` overflows.

**The fix, in `notes/_qs-drafts4/`:** solve with A = B·P·Bᵀ, where P = Lᵀ·L.
L = U1·Perm·L1, with L1 and U1 unit triangular with up to 3 random
off-diagonal entries per row, and Perm a random permutation. B·P·Bᵀ equals
(B·Lᵀ)(B·Lᵀ)ᵀ, so this mixes the columns of B by an invertible map.
ker(Bᵀ) is unchanged, and the row combinations orthogonal to every row
mostly stop being orthogonal. The final step still checks against the
original Bᵀ.

| Matrix | B·P·Bᵀ | B·Bᵀ |
|---|---|---|
| 209-bit unreduced, 10,328×10,308 (gap 22) | 64 independent, 169 ms | 42, 150 ms |
| 80-digit unreduced, 59,391×59,371 | 63 independent, 5.3 s | 0, 4.8 s |
| 80-digit reduced, 15,649×15,585 | 63 independent, 1.6 s | 64, 1.6 s |
| test matrices with gap 90 and 99 | 63 and 64 independent | – |

Applying P costs about 7 passes over a vector per iteration, about 5–10%.
Reducing first is still faster at 80 digits: 2.3 s in all, against 5.3 s for
block Lanczos alone on the unreduced matrix.

`TestFMatchesMontgomery` compares the F term,
`W_{i−2}⁻¹ · (I − S_{i−1}) · T_i · S_i`, with Montgomery's formula at every
iteration on a random 450×400 matrix, and they agree. That supports the
simpler form but doesn't prove it in general.

## Next steps

1. **Double large primes** (`notes/_qs-drafts6/`).
   - `classifyCofactor` splits cofactors up to about 2⁵⁰ with `uint64`
     Miller–Rabin and Brent's rho.
   - Large primes become extra matrix columns after pruning those that appear
     only once, so the structured elimination merges their chains.
   - Still needed: a periodic recount on the large-prime graph to know when
     to stop sieving, and `tryDependency` support for large-prime columns.
   - Then retune slack and the large-prime multiplier upward, and try
     resieving again.
2. **Retune B and M** for SIQS at 70–80 digits. Ours is larger than msieve's.
3. **Allocations.** `runtime.madvise` was 15% of a 200-bit profile.
4. **`factor` defaults.** Raise `qsMaxBits`, since 80 digits now takes under
   3 minutes.

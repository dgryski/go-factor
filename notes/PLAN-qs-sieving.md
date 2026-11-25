# Plan: the quadratic sieve

## Current state (2026-11-24)

Stages 1–3 of the original plan are done, except for a sparse matrix solver.
The code is in package `qs` (`qs/qs.go`). SIQS, the multiplier, large primes
and the faster elimination started from AI drafts. See "Fixed from the
drafts" below for what was wrong in them.

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
    the vertex if it lies inside the block. It's rounded up to an integer
    once per block, so the scan compares bytes to an int.
  - `tryFactorize` checks primes below 2048 (`resieveCut`) from a packed
    (p, r1, r2) table rebuilt for each polynomial. `bucketHits` finds the
    rest for all of a block's candidates at once: primes from 2048 to 16384
    by walking back from their next-block offsets (resieving), and larger
    primes from the block's bucket.
  - `log2Big` reads the top words of a big.Int instead of shifting a copy,
    and `tryFactorize` and `addRelation` use scratch big.Ints on the sieve.
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
  - Block Lanczos (`qs/lanczos.go`) then solves
    the reduced matrix if it has at least 1000 columns. It uses
    A = B·P·Bᵀ with a random symmetric invertible P (see "Block Lanczos"
    below). Dense GF(2) elimination solves smaller matrices, and is the
    fallback if three block Lanczos runs, with seeds 1–3, each find fewer
    than 16 independent dependencies. A failed run costs about 1.5 s at 80
    digits, and dense elimination about 12 s. Repeated dependencies from the
    final step are dropped by hash.
  - `TestFindDependencies` checks the dependencies against full dense
    elimination: the count, that each sums to zero, and that they're
    independent.
  - The structured elimination started from a draft. See "Block Lanczos"
    below for the draft solver.
  - `factor -qs-save=DIR` saves relations after sieving, and
    `go test ./qs -run SavedRelations -relations=FILE` reruns the matrix and
    square root steps on them. Relations for the 80-digit number and its
    209-bit cofactor are in `~/qs-relations/`.
- **`factor` command:**
  - `-methods=trial,qs` runs only the listed methods, in the fixed order.
  - `-v` turns on QS logging.
  - QS is skipped above 266 bits (80 digits) unless `-methods` names it.
  - ECM's effort and squfof depend on the composite's size. See "The
    `factor` pipeline" below.
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
- **Resieving, first try** didn't pay.
  - It walked each prime from 2048 up back through every block with a
    candidate. At 200 bits that's about 28 µs per block, against about 4.7 µs
    per candidate for the old full scan, so it breaks even at about 7
    candidates per block. Blocks had 3.4 at 200 bits and 1.5 at 240.
  - Looser thresholds gave 40–85 candidates per block but made runs slower,
    because the extra candidates nearly all leave cofactors too big for one
    large prime.
  - Bucket sieving gets most of the benefit instead, because a candidate's
    large primes come from its block's bucket.
- **Resieving, now used**, only for primes from 2048 to 16384, inside
  `bucketHits`. It went in together with the integer threshold, the
  allocation fixes, the cofactor guards and the block Lanczos retries. With
  double large primes off, on the seed-7 semiprimes, it collected exactly the
  same relations as before, in:

  | bits | before | after |
  |------|--------|-------|
  | 200  | 1.93 s | 2.13 s |
  | 219  | 7.82 s | 7.31 s |
  | 240  | 29.0 s | 25.8 s |

  These are single runs on a laptop. The 219-bit baseline was 6.4 s in
  earlier sessions, so compare only numbers from the same session.

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

`factor -v -methods=qs` factors it completely in 21.4 s with parallel
sieving on 15 workers (see "Parallel sieving"). On one worker it took
2 min 42 s, about 1.7× msieve's time:

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

The first two drafts used A = B·Bᵀ. Run unchanged apart from inputs and outputs, they gave:

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
The second draft also didn't compile as written:
`int64(0x9e3779b97f4a7c15)` overflows.

**The fix:** solve with A = B·P·Bᵀ, where P = Lᵀ·L.
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

## Double large primes (in progress, 2026-11-24)

Started from a draft. The code is committed, with double large
primes off by default (`dlpMinBits = 1 << 30`), because no setting measured
so far is faster. Tuning was stopped part way.

**Changes:**
- **`qs/dlp.go`:** `isPrime64`, exact Miller–Rabin for `uint64` with the
  7-base set from https://miller-rabin.appspot.com; `rho64`, Brent's rho; and
  `splitCofactor`.
  - `rho64` handles n ≤ 1 and even n, keeps c < n, and gives up after about
    2¹⁹ steps per c, so a prime can't make it spin.
  - `splitCofactor` rejects cofactors above lpBound² before any primality
    test, and `tryFactorize` only takes a cofactor up to lpBound as one large
    prime if it's below pMax² or prime.
  - Both use Montgomery arithmetic (`mont`), because `bits.Div64` is software
    on arm64. With `bits.Div64`, splitting took 36% of a 240-bit run.
- **`qs/dlp_test.go`:** `isPrime64` against a sieve below 200,000 and known
  strong pseudoprimes; `splitCofactor` on random products of two large primes.
- **`qs/qs.go`:**
  - `tryFactorize` returns up to two large primes. Cofactors from lpBound up to
    `dlpCap` (lpBound², at most 2⁵⁰) go to `splitCofactor`.
  - Partials are kept in a list, not combined on the second hit.
    `countPartials` repeatedly drops partials with a large prime that no other
    partial has. `finish` gives each surviving large prime a matrix column, so
    the structured elimination does the combining. This applies to
    single-large-prime partials too.
  - `done` stops sieving once full relations plus surviving partials, minus
    their large primes, reach `want`. The count is redone each time the
    partials grow by 5% plus 1000, so runs overshoot by up to about 6%.
  - `finish` returns the column primes, the factor base followed by the large
    primes. `tryDependency` and the relation file work unchanged with it.
  - A cofactor L·L becomes a full relation with X·L⁻¹. A large prime that
    divides n is returned as a factor.
  - Bucket hits are collected once per block with a candidate bitmask
    (`bucketHits`) instead of scanning the bucket for each candidate. This
    went in just before tuning was stopped and hasn't been timed.
- **`qs/dlptune_test.go`:** the sweep, skipped unless run with
  `go test ./qs -run DLPTune -v -tune`.

**Tuning variables, still `var`s at the end of `qs/qs.go`:**
- `dlpMinBits` is 1 << 30, which turns double large primes off. The sweep
  sets it to 0.
- `dlpSlackFrac` is 0.5. The threshold allowance is that fraction of
  log₂(dlpCap), replacing half of log₂(lpBound).
- `dlpMaxCofactor` is 2⁵⁰.

**Measured, with Montgomery arithmetic but before `bucketHits`, on fixed
semiprimes (seed 7):**

| bits | off | 0.3 | 0.35 | 0.4 | 0.45 | 0.55 |
|---|---|---|---|---|---|---|
| 219 | 6.43 s | 6.37 s | 6.73 s | 7.14 s | 7.95 s | 11.0 s |
| 240 | 23.2 s | 23.0 s | 23.7 s | 24.3 s | 26.2 s | 33.0 s |

At 240 bits, 0.55 needs 27,113 polynomials against 40,332 with double large
primes off, but each polynomial costs more. The best setting only breaks even
at 240 bits. Off, the column approach matches the old combining code: 23.2 s
against 21.6 s committed, with the difference from the overshoot.

**Profile at 240 bits with 0.45, before `bucketHits`:**
- `sievePoly` was 34% of the time by itself, and `fillBuckets` 25%.
- `tryFactorize` was 13%, `splitCofactor` 6%, and `runtime.madvise` 6%.
- Much of `sievePoly`'s own time was likely the per-candidate bucket scan,
  which costs the bucket size times the number of candidates. That's what
  `bucketHits` replaces.

**Measured after `bucketHits` and resieving from 2048, same semiprimes:**

| bits | off | 0.4 | 0.5 | 0.55 | 0.6 | 0.65 |
|---|---|---|---|---|---|---|
| 219 | 7.24 s | 7.14 s | 7.99 s | 9.12 s | 10.4 s | 11.1 s |
| 240 | 25.2 s | 24.1 s | 26.9 s | 29.2 s | 30.3 s | 32.9 s |

At 240 bits 0.4 is about 4% faster than off, and needs 34,596 polynomials
against 40,332. Above 0.5 the partials grow faster than they combine: 0.65
keeps 335,668 partials for 19,955 relations, against 126,268 for 12,002 with
double large primes off.

**80 digits, 15 workers, partial sweep.** `split` on the 266-bit reference
number with `dlpMinBits = 250`. The run was stopped after 7 of its 12
settings, so 256K with double large primes on wasn't measured:

| M | frac | total | sieve | polynomials | partials | matrix | matrix time |
|---|---|---|---|---|---|---|---|
| 128K | off | 18.8 s | 16.5 s | 225,643 | 272,610 | 85,508×83,516 | 2.2 s |
| 128K | 0.4 | 19.5 s | 16.9 s | 202,893 | 335,668 | 95,424×92,248 | 2.4 s |
| 128K | 0.5 | 20.3 s | 16.4 s | 176,303 | 456,626 | 112,107×108,229 | 3.6 s |
| 128K | 0.55 | 20.2 s | 15.7 s | 162,469 | 529,529 | 123,050×118,751 | 4.1 s |
| 128K | 0.6 | 20.5 s | 15.5 s | 148,614 | 563,488 | 122,146×121,409 | 4.5 s |
| 128K | 0.65 | 20.5 s | 15.1 s | 143,046 | 618,720 | 133,913×131,188 | 4.9 s |
| 256K | off | 17.6 s | 14.9 s | 129,832 | 272,610 | 85,207×83,671 | 2.6 s |

- Double large primes cut the sieving time by up to 8% at 0.65, but the
  matrix grows with the large prime columns, and its time more than doubles.
  So the total is 1–2 s slower at every fraction.
- 256K with double large primes off is the fastest so far: its sieving is
  10% faster than at 128K, because each polynomial's setup and
  `fillBuckets` cover twice the interval. It needs 129,832 polynomials
  against 225,643.
- The sieving gain may be bigger with one worker, where the collector's
  share of the time is smaller. Those runs take about 2.5 minutes each.
- Pruning the matrix harder, or merging partials that share a large prime
  before the structured elimination, might recover the difference.

The sweep is `TestSweep` in `qs/dlptune_test.go`:

    QS_N=73129071533331481181630024590293633952610077040929361500801632604572796624590389 \
    QS_RANGE=131072,262144 QS_DLP=0,0.4,0.5 QS_WORKERS=15 \
        go test ./qs -run Sweep -v -tune -timeout 120m

**Not done:**
- The rest of the 80-digit sweep: 256K with double large primes on, and
  everything with one worker.
- Making 256K the default for B ≥ 100,000 (`largeSieveRange`, a var while
  it's tuned), after checking it at 70–75 digits too.
- Choosing `dlpMinBits` and the slack fraction, then making the variables
  constants.
- Tightening the recount near the end, so runs overshoot less than 6%.

## Parallel sieving (done, 2026-11-24)

Sieve SIQS polynomials on several goroutines and merge the relations in one
place. `singlePoly` stays single-threaded, since it only runs below 60 bits.

**Unit of work: one A and all 2^(s−1) of its polynomials.** Each A is
independent once chosen, and the per-A setup (B_j, A⁻¹ mod p, the Gray code
deltas) belongs to whoever sieves it. Workers pull A values from a channel,
so slower cores take fewer.

**Splitting `sieve`:**
- Shared and read-only: `n`, `kn`, `log2KN`, each prime's `p`, `t` and
  `logP`, `slack`, `lpBound`, `pMax`, `dlpCap`, `sieveFrom`, `resieveFrom`
  and `bucketFrom`.
- Per worker, in a new `worker` type: a copy of `roots` (which holds `r1`,
  `r2`, `off1`, `off2` and `inA`), `aInv`, `deltas`, `buf`, `buckets`,
  `small`, the candidate arrays, `scratch`, and the scratch big.Ints.
  `sievePoly`, `bucketHits`, `fillBuckets`, `tryFactorize`, `setRoots`,
  `setRootsInv`, `moveRoots` and `maxLog2Y` become `worker` methods.
- Collector only: `relations`, `partials`, `seenX`, `partialExcess`,
  `nextCount`, `done`, `countPartials`, `logProgress` and the polynomial
  count.

Workers don't call `addRelation`. For each A they send a batch: its sequence
number, the polynomials sieved, and the candidate relations (X mod n,
factors, lp1, lp2). Deduplication, the check for a large prime dividing n,
the L·L case and the partial count all stay in the collector, so none of
them need locks. A factor found by accident can only come from
`addRelation`, so it also turns up in the collector.

**Goroutines:**
1. A generator owns `rng` and the `seen` map, calls `chooseA`, and sends
   (seq, A, aIdx) on a channel buffered to about one per worker. It closes
   the channel when `chooseA` runs out, which replaces "ran out of
   polynomials".
2. `Workers` workers, `runtime.GOMAXPROCS(0)` by default, run today's per-A
   loop. Where it now calls `s.done()`, per block and per polynomial, they
   check an `atomic.Bool` stop flag instead.
3. The collector is the calling goroutine, so `s.siqs` keeps its signature
   and `TestRelations` and `TestDLPRelations` work unchanged. It applies
   batches in sequence order, holding early ones in a small map. When
   `done()` turns true or it finds a factor, it sets the stop flag, drains
   the channels and waits for the workers.

**Why batches are applied in order.** `chooseA` depends only on `rng` and
`seen`, so the generator makes the same sequence of A values as today. Applied
in order, the collector sees the same `addRelation` calls as the sequential
code and stops at the same relation, so every worker count collects exactly
the same relations. The workers only do some extra work that's thrown away.
That makes `Workers=1` against `Workers=N` an exact correctness test, and
keeps timings, saved relations and the matrix step reproducible. The cost is
some waiting behind a slow A, and A-sized units keep that small.

**Workers.** `qs.Workers` (`factor -workers`) is the most workers, and 0
means `runtime.GOMAXPROCS(0)`. TinyGo's `GOMAXPROCS` always returns 1, so
`factor.tiny` sieves on one worker without a build tag. `workers()` also
caps the count by input size, because small runs need only a few A values
and extra workers mostly sieve past the end: one worker below 80 bits, 4
below 130, 8 below 170, and no cap above that. These limits are from this
laptop (5 performance and 10 efficiency cores).

`split` time by worker count, before the caps (seed-1 semiprimes, 10 runs
each):

| bits | 1 | 2 | 4 | 8 | 15 |
|---|---|---|---|---|---|
| 90  | 2.4 ms | 1.6 ms | 1.4 ms | 2.3 ms | 3.9 ms |
| 120 | 8.4 ms | 5.4 ms | 4.3 ms | 5.3 ms | 7.8 ms |
| 140 | 35 ms | 20 ms | 14 ms | 13 ms | 17 ms |
| 160 | 121 ms | 68 ms | 44 ms | 35 ms | 37 ms |
| 180 | 456 ms | 247 ms | 154 ms | 105 ms | 94 ms |

At 70 bits 4 workers were already slower than one (1.47 ms against 1.16 ms).

Seed-7 semiprimes, double large primes off. Every worker count collected
the same relations and polynomials as the sequential code:

| bits | 1 | 2 | 4 | 5 | 8 | 15 |
|---|---|---|---|---|---|---|
| 200 | 2.08 s | 1.11 s | 0.62 s | 0.52 s | 0.38 s | 0.32 s |
| 219 | 7.13 s | 3.79 s | 2.13 s | 1.78 s | 1.23 s | – |
| 240 | 24.8 s | 13.2 s | 10.0 s | 9.5 s | 4.7 s | 2.9 s |

The 80-digit number now factors completely in 21.4 s with 15 workers,
against 2 min 42 s: about 19 s of sieving, 2.5 s for the matrix, and under
a second for the 209-bit cofactor.

**Tests:** `TestWorkers` checks that 1 and 4 workers give identical
relations and partials, in order, at 120 and 160 bits and with double large
primes on. `go test -race` is clean on `TestWorkers`, `TestRelations`,
`TestDLPRelations` and the semiprime tests.

**Not done:**
- Profiling the collector. At 240 bits it handles about 126k partials and
  about 60 `countPartials` runs. If it becomes the bottleneck with more
  cores, run the count on a snapshot in its own goroutine. `done` would
  then lag a little, which only adds overshoot.
- Profiling GC. Each worker allocates its relations' X and factors.
- Smaller work units, such as half an A, if the end of a run with many
  workers wastes too much.
- The 219- and 240-bit double large prime sweeps above were single-threaded.
  `TestDLPTune` now uses all the workers, so set `Workers = 1` to compare.

## The `factor` pipeline (done, 2026-11-24)

With QS this fast, the methods that run before it were the slow part. A
60-digit semiprime with no small factors took about 60 s: about 45 s of ECM
and 3–4 s of squfof, before a 0.3 s QS run.

**ECM** (`ecm.go`), rewritten after MIRACL's `factor.c`:
- Montgomery curves with x-only (X:Z) arithmetic and the Montgomery ladder,
  so stage 1 has no inversions. It takes gcd(Z, n) every 50 primes, and
  redoes a batch one prime at a time if the gcd is n.
- Stage 2 steps by multiples of D = 2310 with baby steps m < D/2 coprime to
  D. Each prime costs one multiply into a running product, since x(cP) −
  x(mP) covers both c − m and c + m. There's a gcd every 16 giant steps.
- Curves are Montgomery's torsion Z/12 family, from k·(−2, 4) on
  v² = u³ − 12u, in place of Suyama's. Over 1000 15-digit factors at B1=2000
  they need 25.9 curves on average against 28.1, at the same time per curve.
  `TestECMTorsion12` counts points over small primes and checks that 12
  divides the order and that order·P = O.
- `ECMContext(ctx, n, b1, b2, curves)` takes B2 and a context, checked at
  the gcd points. `ECM` keeps B2 = 100·B1.
- The stage 2 prime table is a bitset (`primeBits`): 16 MB at B2 = 1.3e8.

A curve at B1=11000 went from 211 ms to 43 ms. The old code did a
`ModInverse` (7–11× a multiply and reduce) in every point operation, and one
point addition per stage 2 prime.

**ECM schedule** (`factor/ecm.go`):
- Levels for 15, 20, 25 and 30-digit factors, with B1, B2 and curves from
  Table 1 of the GMP-ECM README.
- ECM stops once it has used a third of the estimated QS time for what's
  left. `qsTime` interpolates single-worker QS times measured from 30 to 79
  digits, and counts each extra worker as 0.55 of one. That predicted the
  23 s QS run at 80 digits.
- Without QS to follow, ECM runs every level up to 25 digits.

Measured times on this laptop, one ECM curve, 70-digit composite:

| level | B1 | B2 | per curve | whole level |
|---|---|---|---|---|
| 15 digits | 2000 | 147396 | 11 ms | 0.3 s |
| 20 digits | 11000 | 1.9e6 | 54 ms | 4 s |
| 25 digits | 50000 | 1.3e7 | 0.28 s | 61 s |
| 30 digits | 250000 | 1.3e8 | 1.7 s | 13 min |

QS with one worker and with 15:

| digits | 1 worker | 15 workers |
|---|---|---|
| 50 | 0.24 s | 0.05 s |
| 60 | 1.5 s | 0.27 s |
| 64 | 4.9 s | 0.57 s |
| 71 | 20 s | 2.3 s |
| 74 | 52 s | 5.8 s |
| 79–80 | 159 s | 13–23 s |

**squfof** is skipped above 64 bits. QS factors 14–26 digit numbers in
1–2 ms, while squfof takes 4 ms at 18 digits and 284 ms at 26, and past
about 26 digits it uses up its iterations and finds nothing.

**qsMaxBits** went from 220 to 266 bits. An 80-digit semiprime takes 13 s on
15 cores and peaks at 544 MB.

**Shared budget** (`factor/budget.go`). After that, Brent, p−1 and p+1
still took about 1.5 s together on large inputs, most of what was left at 60
digits. Now all four methods before QS share a third of the estimated QS
time, counted from the start of Brent:
- Each runs until a cumulative share is used: Brent 10%, p−1 20%, p+1 35%
  and ECM 100%. Time one method leaves goes to the next.
- `BrentContext`, `PMinus1Context` and `PPlus1Context` take one last gcd
  when the deadline passes, so the work done still counts.
  `TestContextCancel` checks that each returns soon after a deadline.
- Without QS to follow there's no deadline. Brent's three runs and p−1 then
  take 0.38 s on a 60-digit number.

**p−1.** `PMinus1Context(ctx, n, b1)` raises 2 to the prime powers up to B1,
and to powers of 2 up to n. Powers of 2 are cheap and catch p − 1 with a
large power of 2, which `TestPMinus1Backtrack` needs. `PMinus1(n)` still uses
the primes below 50000 to powers below n, about 15× more squarings at the
same B1. `factor` uses B1 = 10⁶, ten times the B1 of the last ECM level, as
GMP-ECM suggests.

Capping the powers at first only gave 4–8×. `Exp` uses Montgomery
multiplication only for exponents of more than one word, so one `Exp` per
small prime power spent most of its time in `nat.div`. Each batch of 16
prime powers is now one exponent:

| digits | `PMinus1` | `PMinus1Context`, B1=50000 | B1=10⁶ |
|---|---|---|---|
| 40 | 40 ms | 5 ms | 73 ms |
| 60 | 48 ms | 5 ms | 86 ms |
| 80 | 83 ms | 6 ms | 109 ms |

**p+1.** Stage 2 ran a full Lucas chain, about 40 multiplies, for every
prime up to B2. It now uses the same baby steps and giant steps as ECM.
With v = a + 1/a and Q = 1,

    V_c − V_m = (a^(c+m) − 1)(a^(c−m) − 1) / a^c,

so one product term covers both c − m and c + m, and both sequences step
with V_(j+d) = V_j·V_d − V_(j−d). That's about one multiply per prime.
- `stage2.go` holds what ECM and p+1 share: the prime bitset, the baby
  steps, and the first and last giant steps. The bitset is built when stage
  2 starts, so a deadline that passes in stage 1 doesn't pay for it. That
  matters for p+1: sieving to B2 = 2.5×10⁷ takes 32 ms, more than p+1's whole
  share at 60 digits.
- `TestPPlus1Stage2Edges` puts the large prime of p + 1 just above B1,
  either side of a giant step, either side of the boundary between giant
  steps, and in the last giant step below B2. It fails with stage 2 stubbed
  out, and when only c + m is tested.
- `factor` now runs p+1 three times with B1 = 2.5×10⁵, five times the B1 of
  the last ECM level, as GMP-ECM suggests. With no deadline that takes about
  1 s on a 60-digit number.

One try:

| digits | before, B1=10⁴ | B1=10⁴ | B1=5×10⁴ | B1=2.5×10⁵ |
|---|---|---|---|---|
| 40 | 382 ms | 25 ms | 69 ms | 269 ms |
| 60 | 452 ms | 15 ms | 71 ms | 337 ms |
| 80 | 489 ms | 16 ms | 74 ms | 349 ms |

End to end with default methods:

| composite | at the start | ECM sized | everything sized | found by |
|---|---|---|---|---|
| 12 + 48 digits | 0.25 s | – | 0.03 s | p−1 |
| 30 + 30 digits | 59.8 s | 2.0 s | 0.30 s | QS |
| 33 + 33 digits | 53.6 s | 2.5 s | 0.88 s | QS |
| 36 + 37 digits | not factored (QS skipped) | 6.9 s | 5.7 s | QS |
| 40 + 40 digits, qs in `-methods` | 87.7 s | 36.7 s | – | QS |

**Not done:**
- The deadlines are wall-clock time, so which method finds a factor can vary
  between runs. The 12 + 48 digit number needs p−1 to reach B1 = 35407
  within about 18 ms. It usually does, but the first run of a freshly built
  binary is slower on macOS, and then ECM found it instead.
- ECM curves run on one goroutine while QS uses all the cores. Running
  curves in parallel would give ECM much more within the same budget.
- ECM still uses math/big multiply and reduce. PRAC chains for stage 1, an
  FFT stage 2 as in GMP-ECM, or Edwards curves would each be faster.
- The ECM and QS timings are from this laptop. A machine with a different
  balance between one core and all of them would want a different fraction.

## Next steps

1. **Finish double large primes** (above).
2. **Retune B and M** for SIQS at 70–80 digits. Ours is larger than msieve's.
3. **Allocations.** `runtime.madvise` was 15% of a 200-bit profile. The
   hot-path big.Int allocations are gone; profile again.
4. **Parallel ECM curves.** QS uses all the cores but ECM uses one, so
   ECM gets far less work done within its budget (see "The `factor`
   pipeline").

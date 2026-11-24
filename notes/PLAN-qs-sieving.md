# Plan: the quadratic sieve

## Current state (2026-11-23)

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
  - The threshold is per block, from the exact |Y| at the block's ends, and at
    the vertex if it lies inside the block.
  - `tryFactorize` divides only by primes whose roots match x.
- **Relations:**
  - Each relation stores its factors sparsely, as column indices repeated by
    multiplicity, rather than a dense exponent vector.
  - Relations whose ±X was already seen are dropped. Different polynomials do
    find the same X, and before this check those copies made about 1% of
    70-bit SIQS runs fail with only trivial dependencies.
- **`findDependencies`:** dense GF(2) elimination in row echelon form, with
  pivot swaps and flat backing arrays.
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

| bits | first sieve (single poly) | now |
|------|------|-------|
| 60   | 0.82 ms | 0.55 ms |
| 80   | 5 ms    | 0.93 ms |
| 100  | 39 ms   | 2.6 ms  |
| 120  | 217 ms  | 7.0 ms  |
| 140  | –       | 29 ms   |
| 160  | –       | 128 ms  |
| 180  | –       | 1.1 s   |
| 200  | –       | 3.8 s   |
| 220  | –       | 17 s    |

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
  half the run time. It now checks only after a division. The remaining cost
  is the x mod p scan over all roots for each candidate.

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

This SIQS chose k = 5 and B = 1,570,053 (59,370 primes, M = 128K). It
collected 59,391 relations, with 239,765 partials left uncombined, in
4 min 58 s. That's about 3× msieve's whole run. I stopped it after more than
4 minutes in `findDependencies`, which hadn't finished. The dense matrix needs
about 2 × 59k × 59k bits ≈ 900 MB, and elimination is O(n³/64).

## Next steps

1. **Sparse linear algebra.** This is now the limit at 80 digits.
   - First, remove singletons (columns that appear in only one relation) and
     the empty columns, then use structured Gaussian elimination to shrink the
     matrix before the dense step.
   - After that, or instead, block Lanczos.
2. **A faster candidate check.** `tryFactorize` scans every root for each
   candidate.
   - Resieve the block for the larger primes, or keep a compact
     (p, r1, r2) array for the scan.
   - Once that's cheaper, retune the large prime slack. More candidates may
     then pay off.
3. **Matching msieve's sieve speed.**
   - Sieve the largest primes with bucket sieving rather than once per block.
   - Retune B and M for SIQS at 70–80 digits. Ours is larger than msieve's.
4. **`factor` defaults.** Raise `qsMaxBits` once 80 digits finishes in
   minutes.

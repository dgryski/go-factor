# Plan: real sieving for the quadratic sieve

## Current state (2026-11-23)

Stage 1 is done. It started from `qs/ai_studio_code.go`, an untracked AI draft
that is `package main` and doesn't build on its own. Delete it once it's no
longer needed.

- `qs.go` (package `factor`) has `QuadraticSieve(n) (primes, composites []*big.Int)`,
  which recurses on the two halves found by `qsSplit`.
- `qsSplit` removes small factors by trial division up to B, builds the factor
  base, and calls `sieve`. Linear algebra (`findDependencies`) and the gcd step
  (`tryDependency`) are unchanged.
- `sieve` is a single-polynomial sieve over X = ⌈√N⌉ + x:
  - x ∈ [max(−M, 1−m), M), in 32 KB blocks. X ≤ 0 is skipped because it only
    repeats the positive side's Y values.
  - Roots come from Tonelli–Shanks in `uint64` (`sqrtModPrime`, `buildRootTable`).
  - Primes below 30 are not sieved. The threshold is lowered by their expected
    contribution of about 2·log₂p/p (log₂2/2 for p = 2).
  - Threshold per block: log₂(|x|·(2m + |x|)) at the block's largest |x|,
    minus 1.75·log₂(p_max) and the small-prime allowance.
  - `tryFactorize(x, y, roots)` divides only by primes whose roots match x mod p.
- `qsParams` sets B = 3·L(N)^½, with a minimum of 500, and M = max(B/3, 256)·B.
  - With B = L(N)^½, as the plan suggested, the base is far too small below
    about 80 bits: about half of 40–70-bit semiprimes failed.
  - With a minimum of 128 instead of 500, about 1% of 18–42-bit semiprimes
    ran out of sieve range before finding enough relations.
  - There is no fallback that sieves further out when relations run short.
    `qsSplit` returns nil.
- `go.mod` says go 1.17, so the code avoids the `min`, `max` and `clear`
  builtins.

Measured with qsSplit on random balanced semiprimes, averaged per number:

| bits | old (trial division) | sieve |
|------|------|-------|
| 40   | 4.4 ms | 0.19 ms |
| 52   | 40 ms  | 0.44 ms |
| 60   | 84 ms  | 0.82 ms |
| 80   | –      | 5 ms    |
| 100  | –      | 39 ms   |
| 120  | –      | 217 ms  |

Still open from stage 1:

- **Tests not written yet:**
  - Tonelli–Shanks against brute force.
  - The sieve, with very high slack, against brute-force `tryFactorize`.
  - `BenchmarkQuadraticSieve`. The numbers above came from a throwaway test.
- **`factor/main.go`:** it still limits QS to composites ≤ 64 bits
  (`qsMaxBits`). That can go up to about 120 bits.

## Reference: an 80-digit number

73129071533331481181630024590293633952610077040929361500801632604572796624590389
(266 bits) = p18 · p23 · p41:

- 104194448528493551
- 11494857165741861450481
- 61057905163708884083219136437717059663819

These are easy targets for ECM: GMP-ECM with B1 = 11000 found p18 in step 2
within a minute. `factor/msieve.log` has msieve 1.53 factoring it by SIQS in
1 min 36 s with these settings:

- multiplier 5
- 12 blocks of 32 KB per polynomial
- factor base bound 1,301,057 (49,858 primes)
- large prime bound 1.3·10⁸ (26 bits)
- 49,954 relations: 26,246 full and 23,708 combined from 267,452 partials
- a sparse matrix with about 30 entries per column

The current single-polynomial sieve can't get there with any choice of
parameters:

- **qsParams' own choice:** B ≈ 1.6·10⁷ (about 490k primes) and M ≈ 9·10¹³.
- **msieve's B instead:** |Y| ≈ 2·|x|·√N is about 2¹⁷⁰ across a useful
  interval. My rough estimate is that about 1 in 10⁸–10⁹ positions is smooth,
  so that's days to weeks of sieving.
- **A 50k-prime dense matrix:** about 300 MB plus the same again for the
  combination bitsets. That's workable but slow.

Reaching this size needs stages 2 and 3 below: SIQS, large primes, a multiplier
and a sparse solver, with settings close to msieve's.

## Stage 2: cheap improvements

- **Knuth–Schroeppel multiplier.** Factor kN instead of N, with a small k
  chosen so that many small primes are quadratic residues. A few dozen lines,
  usually a solid constant-factor speedup.
- **Large prime variation.** If a candidate leaves a cofactor that is prime
  and below about 30–100·p_max, keep it as a partial relation. Two partials
  with the same large prime combine into one full relation. Needs a map from
  large prime to partial, and raises the threshold slack by log₂ of the
  large-prime bound.
- **Faster elimination.** `findDependencies` XORs every row for each pivot.
  That's fine for a few thousand primes, but limit the search to the rows that
  can still be pivots and benchmark it. Block Lanczos isn't worth it below
  about 20k primes, but 80-digit inputs need about 50k.

## Stage 3: many polynomials (SIQS)

- Use many polynomials Q(x) = (Ax + B)² − N, where A is a product of factor
  base primes chosen so that A ≈ √(2N)/M. This keeps |Q(x)| near √N·M on every
  interval instead of growing with x.
- With the self-initializing variant, the next polynomial's roots are cheap to
  update from the current ones.
- This is the version that reaches 50+ digits, and it's the largest piece of
  work: choosing A, computing B with Gray code switching, and per-polynomial
  root updates.
- Relations then record (Ax + B) instead of X, with A's factors added to the
  exponent vector. That's a small change to `qsRelation` and `tryDependency`.

## Scope and order

- Stage 1 is a self-contained rewrite of the relation-gathering half of
  `qsSplit`, plus three or four new helpers in `qs.go`.
- Stages 2 and 3 build on stage 1 without touching the API.
- Open question: stop after stage 2, or go on to SIQS?

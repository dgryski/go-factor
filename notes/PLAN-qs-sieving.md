# Plan: real sieving for the quadratic sieve

## Current state (2026-11-22)

- `qs.go` (package `factor`) has `QuadraticSieve(n) (primes, composites []*big.Int)`,
  which recurses on the two halves found by `qsSplit`.
- `qsSplit` does not sieve. It builds a `big.Int` Q(x) = (⌈√N⌉ + x)² − N for every
  candidate x and trial-divides it by every factor base prime (`tryFactorize`).
  That makes it Dixon/Kraitchik rather than a quadratic sieve.
- Parameters are fixed in three steps by bit length: B = 250/1000/5000 and
  range = 2000/20000/100000.
- Linear algebra (`findDependencies`, bit-packed GF(2) elimination that tracks
  row combinations) and the gcd step (`tryDependency`, which tries a−b and a+b)
  are correct and can be reused as they are.
- `qs_test.go` covers edge cases and random semiprimes with 10–26-bit factors.
- `factor/main.go` runs QS after Squfof, only on composites ≤ 64 bits
  (`qsMaxBits`). In practice Brent, p−1 and Squfof factor such numbers first,
  so the stage rarely fires.

Sieving replaces per-candidate big-integer work with integer additions over a
byte array. Expensive division then happens only on candidates that look smooth.
Each stage below is usable and tested on its own. The public API never changes.

## Stage 1: single-polynomial sieve

Replace the `try` loop in `qsSplit`. Factor base building, `findDependencies`
and `tryDependency` stay.

1. **Square roots mod p.** Add `sqrtModPrime(nmodp, p uint32) uint32`
   (Tonelli–Shanks). Base primes fit in 32 bits, so use `uint64` arithmetic
   with no `big.Int`. Precompute for each prime p:
   - the root t
   - log₂p as a byte
   - the two starting offsets x ≡ (±t − m) mod p, where m = ⌈√N⌉

   p = 2 has a single root and needs special handling.
2. **Sieve in blocks.** Cover x ∈ [−M, M) in cache-sized blocks of about 32 KB
   of `[]byte` each. For each prime, add log p at both roots, stepping by p.
   Save each prime's next offset when a block ends instead of recomputing it.
   Skip primes below about 30, which cost the most to sieve and add the least,
   and lower the threshold to make up for them.
3. **Threshold.** In a block, |Q(x)| ≈ 2·|x|·m. Use log₂ of that, computed once
   per block from its largest |x|, minus a slack of T·log₂(p_max) with T around
   1.5–2 (tune it). Positions scoring above the threshold are candidates.
4. **Checking candidates.** Only for candidates, compute Q(x) as a `big.Int`.
   Divide by p only when x mod p matches one of p's roots (an integer test).
   That turns hundreds of big divisions per candidate into a handful. Add a
   `tryFactorize` variant that takes the root table. The sign column is
   unchanged: negative x still gives negative Q.
5. **Parameters.** Replace the size steps with formulas:
   - L(N) = exp(√(ln N · ln ln N))
   - factor base bound B ≈ L(N)^½
   - M chosen as a multiple of B

   Tune both against benchmarks. If an interval doesn't produce enough
   relations, sieve the next interval out on each side instead of giving up.

**Tests:**

- Tonelli–Shanks against brute force for every prime < 10⁵, and for random
  residues mod larger primes.
- With the slack set very high, so every position is a candidate, the sieve
  must find exactly the smooth values a brute-force `tryFactorize` finds over
  the same interval.
- The existing `qs_test.go` cases, plus random semiprimes of 30–45 bits.
- `BenchmarkQuadraticSieve` per input size, used to tune B, M and T.

**Expected reach:** single-polynomial QS is the textbook version. Rough guess,
to be confirmed by the benchmarks: practical to about 25–30 digits. Past that,
|Q(x)| grows with x and smooth values become rare. After this stage, raise
`qsMaxBits` in `factor/main.go` to what the benchmarks support.

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
  can still be pivots and benchmark it. Block Lanczos isn't worth it at these
  sizes.

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

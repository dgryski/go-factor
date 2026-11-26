# General Number Field Sieve (GNFS) — Go Implementation Plan

This document plans an implementation of the **General Number Field Sieve (GNFS)** in Go, designed as a sibling package (`gnfs/`) to the existing `qs/` package while sharing its finite-field arithmetic, SGE reducer, and Block Lanczos linear algebra.

Section 4 collects the formulas the pipeline depends on, with their derivations, so each stage can be checked against them. Appendix B lists the corrections made to the first draft of this plan. Section 0 records where the implementation stands.

---

## 0. Status (2026-11-25)

The package exists and factors end to end: `gnfs.Factor`, wired into `factor/main.go` after `qs` (methods list, `gnfsMaxBits = 450`, `-v`, `-workers`, `-gnfs-save`). `go test ./gnfs` takes about 15s.

### Milestones

| Milestone | State |
| :--- | :--- |
| **M-1** | Not done. No `qs` timings at 90–100 digits and no `big.Int` multiplication benchmark are recorded. |
| **M0** | Done. `internal/arith` has the 64-bit arithmetic plus polynomials mod $p$ (roots, Rabin's test, factoring); `internal/linalg` has `Reducer`, Block Lanczos and a rows-based `FindDependencies`. `qs` uses both and its tests pass. The `qs` benchmarks have not been compared before and after the move. |
| **M1** | Done. Base-$m$, coarse bad ideals, QCs, line sieve, square root by Hensel lifting. The 4.7 checks are tests or run-time checks: lattice round trip, factorizations multiplying back to $G$ and $F$, $\beta^2 = \gamma$, and $X_{\text{alg}}^2 \equiv X_{\text{rat}}^2$. Dependencies failing $\beta^2 = \gamma$ are counted in the log. |
| **M2** | Mostly done. Dedup, singleton pruning on $(p, r)$ ideals, `Reducer`, parallel Block Lanczos. `BenchmarkBlockLanczos` covers 150k columns, but no numbers are recorded. |
| **M3** | Mostly done. Gauss reduction, $\lambda_p$, Franke–Kleinjung buckets, two-sided 2LP. Resieving is replaced by collecting each candidate's bucket hits. Verified only up to 195 bits (59 digits), not 80–95 digits; see the parameter table below. |
| **M4** | Partial. Smooth $c_d$ search, Murphy $\alpha$, skewness, and a small rotation search ($u \in [-2, 2]$, $v \in [-8, 8]$, only from 140 bits). No translation and no root sieve. `-gnfs-save` writes $N$, the polynomial and the $(a, b)$ pairs, but nothing reads them back to resume. |
| **M5** | Partial. The square root already uses a product tree, and the table has $d = 5$ above 365 bits. Exact bad ideals, Kleinjung, P−1/ECM cofactoring and runs at 100+ digits are not done. |

### Parameter table (`defaultParams`)

| Bits | State |
| :--- | :--- |
| < 115 | Line sieve, $d = 3$, Br 800, Ba 1200. Passes the 30–42 digit test; not swept. |
| 115–150 | $d = 4$, Br = Ba = 8000, LP $2^{17}$. 0.3s at 115 bits, 0.7s at 149 bits. |
| 150–200 | $d = 4$, Br = Ba = 30000, LP $2^{18}$, $I = 1024$. 1.6s at 151 bits, 4s at 195 bits. The old row (Br 5000, Ba 8000) found no cycles at 195 bits. |
| 200–250 | $d = 4$, Br = Ba = 100000, LP $2^{20}$, mfb 40, $I = J = 1024$. 11s at 204 bits, 26s at 230 bits, 3m9s at 249 bits, one number each. The old row (Br 10000, Ba 20000, LP $2^{21}$) failed at 249 bits and on one of two 230-bit numbers, running out of special-$q$; Br = Ba = 60000 with LP $2^{19}$ also failed at 249 bits. |
| ≥ 250 | Not tested. |

`chooseDegree` reads the degree from this table.

### Differences from the plan

- **Skewness** minimises the product of the $F$ and $G$ norms, not $\int F^2$ alone. For base-$m$ polynomials the product is nearly flat in $s$, so the search can return a large skew. On the one polynomial checked, the $\int F^2$ optimum was similar (2.2e7 against 3.6e7) and the yield was the same.
- **Bad ideals** are flagged for any multiple root of $f \bmod p$, which is a little more than $p^2 \mid \operatorname{disc}(f)$ requires.
- **Multiple roots above $B_a$** have no QC, so special-$q$ roots and algebraic large primes at multiple roots are skipped.
- **Small composite factors** that fit in 64 bits go to rho instead of back into GNFS.
- **Progress** is measured in relations against an estimated target, not in excess. Excess stays at $-\#\text{fixed columns}$ until large-prime cycles form, dips, then climbs steeply. The target starts at $0.9\,(\pi(\text{LP}_{\text{rat}}) + \pi(\text{LP}_{\text{alg}}))$ and, once excess is rising, extrapolates its slope against the relation count. Over ten runs from 195 to 249 bits, the relations collected were 0.64–1.32 times that sum, so 0.9 is only a middle value. The ratios are inflated by overshoot (next bullet).
- **Overshoot:** runs used to stop with excess 63–7327 against a goal of 32, because excess was checked every ~7% more relations. Once excess rises above $-\#\text{fixed columns}$ it is now checked every 1%. This change is not yet measured.

### Remaining work

Priority 1 blocks the next milestone or covers code that hasn't been measured. Priority 2 is needed before 80–100 digits are practical. Priority 3 is the rest of M4 and M5. Size is the code involved: S is under a day, M a few days, L a week or more. Machine time for tuning runs comes on top and grows quickly with size.

Parameter sweeps use `TestSweep` in `gnfs/tune_test.go`, which is skipped without `-tune`. For example, `GNFS_BITS=230,249 GNFS_CANDS=default,100k GNFS_SEEDS=4 go test ./gnfs -run Sweep -v -tune -timeout 120m`. The package path has to come before `-tune`. Each line gives the time, the relations collected and the large-prime ideal count, which together calibrate `relTargetRatio`. Seeds are fixed, so a sweep produces the same numbers each time.

| Priority | Task | Size | Notes |
| :--- | :--- | :--- | :--- |
| 1 | Run a 230-bit factorization with `-v` | S | Checks the new progress messages and the 1% excess checks, and how much sieving the overshoot fix saves. Two 150-bit runs ended with excess 261 and 101, down from 63–7327 before. |
| 1 | Check the 200–250 row on more numbers | S | One number per size so far; 249 bits took 3m9s and may be near the edge. A run on four 249-bit numbers was started and stopped before the first finished. |
| 1 | M-1 baselines | S | `qs` at 90, 95 and 100 digits; `big.Int` multiplication at 1–10 Mbit; record `BenchmarkBlockLanczos`; compare the `qs` benchmarks with the pre-refactor code. |
| 1 | Tune the 250–300 and 300–365 bit rows | M | Untested. Needed for the M3 target of 80–95 digits. |
| 1 | Profile one 80-digit run | S | Sieve vs. cofactorization vs. linear algebra vs. square root. Decides what priority 2 speed work is worth doing. |
| 2 | Keep sieving when the special-$q$ range runs out | S | The sieve gives up after four doublings of QMax. That caused most failures of undersized rows. The QC bound has to stay above the largest special-$q$. |
| 2 | Cofactorization speed | M | Every trial division is `big.Int.QuoRem`. Use `uint64`/128-bit norms where they fit. |
| 2 | Skip both-even $(i, j)$ in the sieve | S | Currently skipped only in the candidate scan, so the sieve saves nothing. |
| 2 | Calibrate `relTargetRatio` per row | S | Measured 0.64–1.32. A per-row value, or one derived from $B_r$, $B_a$ and the large-prime bounds, gives a better first estimate. |
| 2 | Sweep the < 115 bit line-sieve row | S | Passes its test but was never tuned. |
| 2 | Resume from `-gnfs-save` | S | `readPolyAndPairs` exists. Each $(a, b)$ has to be re-factored, which needs a `cofactorCandidate` variant without sieve hits. |
| 2 | Time budget in `factor/main.go` | S | Estimates the time for the methods before the sieve from `qs` timings, even for numbers only GNFS will run on. |
| 3 | Translation and root sieve (M4) | L | Rotation is limited to $u \in [-2, 2]$, $v \in [-8, 8]$. Compare yield against the current polynomials at 90–100 digits. |
| 3 | $d = 5$ tuning and 100–110 digit runs (M5) | M | Mostly machine time: probably hours per 100-digit run. |
| 3 | P−1/ECM cofactoring for mfb above 64 bits | M | Needs 128-bit arithmetic. |
| 3 | Exact bad ideals | L | Only if the coarse handling shows dependencies failing $\beta^2 = \gamma$ at larger sizes; none have so far. |
| 3 | Free relations | S | Small gain (about $1/|\mathrm{Gal}|$ of primes). When added, $p$ goes in `RatCols` or `RatLPs` like any other rational prime and the $d$ ideals in `AlgCols`. |
| 3 | Skewness from $\int F^2$ alone | S | Matches the plan; made no difference on the one polynomial checked. |
| 3 | Kleinjung polynomial selection | L | Optional; only if polynomial quality limits yield. |

GNFS trial-divides only up to $\max(B_r, B_a)$. That is enough for correctness, and `factor/main.go` runs rho and ECM first.

---

## 1. Context & Target Regime

- **Crossover is unknown and must be measured first.** `factor/main.go` records that `qs` takes about 15s and 550MB for 80 digits on 15 cores, doubling every 3–4 digits, which extrapolates to tens of minutes at 100 digits. A new Go GNFS will start well behind tuned implementations (CADO-NFS, msieve), so its crossover against this SIQS is probably above 100 digits.
  - **Before M1:** time `qs` on 90, 95 and 100 digit semiprimes and record relations wanted, sieve time, matrix size and LA time. That sets the real target and the yardstick for M5.
- **Target:** general integers in the **100–130+ digit** range, where $L_N[1/3, (64/9)^{1/3}]$ overtakes SIQS's $L_N[1/2, 1]$.
- **Design philosophy:** build an end-to-end correct pipeline with simple algorithms (base-$m$ polynomials, line sieve, exact square root), verify it on 30–60 digit numbers against `qs`, then swap in the fast stages (better polynomial selection, special-$q$ lattice sieve, parallel Block Lanczos, product-tree square root). Keep each simple version as the test oracle for its fast replacement.

---

## 2. Package Layout & Code Reuse

Refactor reusable primitives out of `qs/` into shared internal packages so `qs/` and `gnfs/` stay decoupled without duplicating code:

```text
internal/
  arith/        # Moved from qs (dlp.go, qs.go):
                #   mont, isPrime64, rho64, splitCofactor,
                #   primesUpTo, powMod64, invMod, sqrtModPrime, log2Big
                #   + new: polynomials mod p (Cantor–Zassenhaus roots, Rabin irreducibility)
  linalg/       # Moved from qs (matrix.go, lanczos.go):
                #   reducer (SGE), blockLanczos (+ symPrecond), denseDependencies
                #   + new: rows-based findDependencies, parallel matrix-vector products
  progress/     # Already shared: EWMA rate for "time left" in progress messages
qs/             # Existing SIQS implementation (imports internal/arith, internal/linalg)
gnfs/
  gnfs.go       # Entry point Factor(n), parameter selection, pipeline orchestration, Logf
  poly.go       # Polynomial selection (base-m -> chosen c_d + rotation/root sieve -> Kleinjung)
  base.go       # Rational/algebraic factor bases, bad ideals, quadratic character base
  sieve.go      # Line sieve (M1) and special-q lattice sieve with Franke–Kleinjung buckets (M3)
  cofactor.go   # Two-sided cofactorization (trial division + 2LP split + optional P-1/ECM)
  filter.go     # Deduplication, singleton pruning on (p, r) ideals, preparation for SGE
  sqrt.go       # Algebraic square root: exact product + Hensel lifting (M1), product tree (M5)
  rel.go        # Checkpointing of polynomials and relations, building on qs/relfile.go
```

Notes on what exists today:

- `mulB` and `mulBt` are closures inside `blockLanczosCheck` in `qs/lanczos.go`, not functions, so parallelising them means restructuring that function. Lanczos is currently single-threaded.
- `findDependencies` takes `[]relation`, a `qs` type, so `internal/linalg` needs an entry point on `[][]int32` rows. `reducer`, `blockLanczos` and `denseDependencies` already work on rows.
- Go inlines across packages, so moving the 64-bit arithmetic out of `qs` should not cost speed, but check with the `qs` benchmarks after the move.
- Wire GNFS into `factor/main.go` like `qs`: an entry in `allMethods`, a size limit like `qsMaxBits` that `-methods` lifts, the shared time budget, `-v` setting `gnfs.Logf`, and progress messages with `internal/progress`.

---

## 3. End-to-End Pipeline Overview

```text
  N
  │
  ▼
[1. Poly Select] ──► f(x) of degree d ∈ {3, 4, 5}, g(x) = a1·x - a0, f(m) ≡ g(m) ≡ 0 (mod N)
  │
  ▼
[2. Factor Bases] ─► RFB (p ≤ Br), AFB ((p,r) for p ≤ Ba, incl. projective), QCB (q > all LP bounds)
  │
  ▼
[3. Sieve + Cofactor] ─► coprime (a, b), b > 0, smooth on BOTH G(a,b) and F(a,b) = b^d·f(a/b)
  │                       (line sieve -> special-q lattice sieve, Franke–Kleinjung, 2LP)
  ▼
[4. Filter + LA] ──► dedup ─► prune singleton ideals (p, r) ─► SGE ─► Block Lanczos ─► deps
  │
  ▼
[5. Square Roots] ─► Rational: X_rat from the factorizations of G(a,b)
  │                   Algebraic: β with β² = γ in Z[α̂] by Hensel lifting, X_alg = β(m̂) mod N
  ▼
  check X_alg² ≡ X_rat² (mod N), then gcd(X_alg ± X_rat, N)
```

---

## 4. Notation and Verified Formulas

Everything below uses these definitions. Each identity has a one-line derivation; tests should check the numeric ones directly (see 4.7).

### 4.1 Polynomials and homogenisation

- $f(x) = \sum_{i=0}^{d} c_i x^i \in \mathbb{Z}[x]$, irreducible, $c_d \neq 0$; $\alpha$ a complex root; $K = \mathbb{Q}(\alpha)$.
- $g(x) = a_1 x - a_0$ with $\gcd(a_1, a_0) = 1$ and $\gcd(a_1, N) = 1$ (otherwise $N$ has a factor in common with $a_1$).
- Common root $m \equiv a_0 a_1^{-1} \pmod N$ with $f(m) \equiv 0 \pmod N$.
- $F(a, b) = b^d f(a/b) = \sum c_i a^i b^{d-i}$ and $G(a, b) = a_1 a - a_0 b$.

### 4.2 Norms

$$N_{K/\mathbb{Q}}(a - b\alpha) = \prod_k (a - b\alpha_k) = b^d \prod_k (a/b - \alpha_k) = \frac{F(a, b)}{c_d}.$$

### 4.3 Monic conversion

$\hat f(y) = c_d^{\,d-1} f(y / c_d)$ is monic of degree $d$ (its leading term is $c_d^{d-1} c_d (y/c_d)^d = y^d$), with root $\hat\alpha = c_d \alpha$. Mod $N$, $\hat m = c_d m$ is a root, since $\hat f(\hat m) = c_d^{d-1} f(m) \equiv 0$.

$$c_d a - b\hat\alpha = c_d (a - b\alpha) \in \mathbb{Z}[\hat\alpha], \qquad N(c_d a - b\hat\alpha) = c_d^{\,d} \cdot \frac{F(a,b)}{c_d} = c_d^{\,d-1} F(a, b).$$

### 4.4 Which ideal a prime divides

- **Affine:** if $p \mid F(a, b)$ and $p \nmid b$, then $a \equiv b r \pmod p$ for the root $r = a b^{-1} \bmod p$ of $f$, and the prime is the ideal $(p, r)$. This is how to recover $r$ for a large prime.
- **Projective:** if $p \mid b$ and $p \mid F(a, b)$, then $p \mid c_d a^d$, and $p \nmid a$ since $\gcd(a, b) = 1$, so $p \mid c_d$. Projective ideals $(p, \infty)$ occur only for primes dividing $c_d$, which are small and in the factor base.
- **Rational, $p \mid a_1$:** then $p \mid G(a, b) \iff p \mid a_0 b \iff p \mid b$, because $\gcd(a_1, a_0) = 1$.
- **Good primes:** if $p \nmid c_d$ and $p^2 \nmid \operatorname{disc}(f)$, then $p$ does not divide the index $[\mathcal{O}_K : \mathbb{Z}[\alpha]]$, because $\operatorname{disc}(f) = [\mathcal{O}_K : \mathbb{Z}[\alpha]]^2 \operatorname{disc}(K)$. Each root $r$ mod $p$ then gives exactly one degree-1 prime ideal $\mathfrak{p}$, and $v_{\mathfrak{p}}(a - b\alpha) = v_p(F(a, b))$ when $a \equiv br \pmod p$. This holds for repeated roots too.

### 4.5 Why the dependency gives a square

Let $S$ be a dependency. The matrix makes every ideal exponent even, and the quadratic characters make it very likely that $\prod_S (a - b\alpha)$ is a square in $K$ (4.6).

- **Algebraic side:** $\prod_S (c_d a - b\hat\alpha) = c_d^{|S|} \prod_S (a - b\alpha)$. This is a square in $K$ only if $c_d^{|S|}$ is, so $|S|$ must be even. It is an algebraic integer whose square root is integral, so it equals $\beta_0^2$ with $\beta_0 \in \mathcal{O}_K$. Since $\hat f'(\hat\alpha)\, \mathcal{O}_K \subseteq \mathbb{Z}[\hat\alpha]$ for monic $\hat f$,
  $$\gamma = \hat f'(\hat\alpha)^2 \prod_S (c_d a - b\hat\alpha) = \beta^2, \qquad \beta = \hat f'(\hat\alpha)\,\beta_0 \in \mathbb{Z}[\hat\alpha].$$
- **Rational side:** $a - bm \equiv G(a, b) / a_1 \pmod N$, so applying $\varphi: \hat\alpha \mapsto \hat m$,
  $$\varphi(\gamma) \equiv \hat f'(\hat m)^2 \, c_d^{|S|} \, a_1^{-|S|} \prod_S G(a, b) \pmod N.$$
  With $\prod_S G = R^2$, where $R = \prod p^{e_p/2}$ over the factorizations of $|G|$ and the sign column makes the product positive, the rational square root is
  $$X_{\text{rat}} = \hat f'(\hat m)\; c_d^{|S|/2}\; a_1^{-|S|/2}\; R \bmod N, \qquad X_{\text{alg}} = \beta(\hat m) \bmod N.$$
  Equivalently, multiply $X_{\text{alg}}$ by $a_1^{|S|/2}$ and drop $a_1$ from $X_{\text{rat}}$, which avoids the inverse.
- **The all-ones column** forces $|S|$ even, which both sides need: $c_d^{|S|}$ on the algebraic side and $a_1^{|S|}$ and $c_d^{|S|}$ on the rational side.
- **The rational sign column** (1 if $G(a, b) < 0$) makes $\prod_S G > 0$.
- **Algebraic sign:** $\prod_S F > 0$ is needed too, since $N(\beta_0)^2 = c_d^{(d-1)|S|} \prod_S F$. The QCs enforce this with high probability. A column for $F(a, b) < 0$ costs nothing and makes it certain.
- **Expected algebraic norm:** $N(\beta) = \pm\, N(\hat f'(\hat\alpha)) \cdot c_d^{(d-1)|S|/2} \cdot \prod p^{e_p/2}$, with the product over the factorizations of $|F(a, b)|$ including the special-$q$. Only the Couveignes alternative (Appendix A) needs this; it is also a cheap test.

### 4.6 Quadratic characters

- **Choice:** pick primes $q$ with simple roots $s$: $f(s) \equiv 0$, $f'(s) \not\equiv 0 \pmod q$, $q \nmid c_d$.
- **Bit:** each relation gets bit 1 when the Legendre symbol $\left(\frac{a - bs}{q}\right) = -1$. Using $a - bs$ instead of $c_d a - b\hat s$ is fine, because they differ by the factor $c_d$, which appears $|S|$ times, an even number.
- **Bound:** $q$ must exceed every prime that can divide $F(a, b)$: the algebraic large-prime bound and the largest special-$q$.
  - Why: if $q \mid a - bs$ and $q \nmid b$, then $F(a, b) \equiv b^d f(s) \equiv 0 \pmod q$. And $q \mid b$ would force $q \mid a$, which contradicts $\gcd(a, b) = 1$.
  - So with $q$ above every possible factor of $F$, $q \nmid a - bs$ is guaranteed. A bound of $q > B_a$ is not enough.
- **How many:** after the matrix, the products that are squares as ideals but not as elements form an $\mathbb{F}_2$-space. Its dimension is at most the unit rank + 1, plus the 2-rank of the class group, plus the number of bad ideals handled coarsely. 50–64 characters is ample at these sizes. Add one per coarsely handled bad ideal.

### 4.7 Numeric checks to build into tests

- $\hat f(\hat m) \equiv 0 \pmod N$ and $c_d^{d-1}F(a, b) = N(c_d a - b\hat\alpha)$ for random $(a, b)$, using $\mathbb{Z}[\hat\alpha]$ multiplication mod $\hat f$.
- For each relation, the recorded factorizations multiply back to $G(a, b)$ and $F(a, b)$, including the special-$q$ and projective primes.
- **Lattice round trip,** using the *reduced* basis $\mathbf{u}, \mathbf{v}$, since inverting against the unreduced basis $(q, 0), (r, 1)$ is an easy mistake. Elliptic1/General-Number-Field-Sieve fixed exactly this bug in commit `98aefda`. For random $(i, j)$ in the sieve region:
  - $(a, b) = i\mathbf{u} + j\mathbf{v}$ maps back to $(i, j)$ through the inverse of the reduced basis;
  - $a \equiv rb \pmod q$, so $q \mid F(a, b)$;
  - for each factor base $(p, r_p)$, the points with $i \equiv \lambda_p j \pmod p$ are exactly those with $p \mid F(a, b)$ through that root.
- After the square root, $\beta^2 = \gamma$ modulo the lifting modulus.
- $X_{\text{alg}}^2 \equiv X_{\text{rat}}^2 \pmod N$ before taking gcds. A failure here means a formula or matrix bug. A pass with trivial gcds is just an unlucky dependency, which happens about half the time per dependency.

---

## 5. Detailed Stage Designs

### Stage 0: Shared Infrastructure Refactor (`internal/arith`, `internal/linalg`)

1. **`internal/arith`:**
   - Move `mont`, `isPrime64`, `rho64`, `splitCofactor` (`qs/dlp.go`) and `primesUpTo`, `powMod64`, `invMod`, `sqrtModPrime`, `log2Big` (`qs/qs.go`).
   - Add **polynomial arithmetic mod $p$** on `[]uint64` coefficients (31-bit and 64-bit $p$): multiply, $x^e \bmod (h, p)$, gcd, roots by **Cantor–Zassenhaus** ($\gcd(h, x^p - x)$, then split with $(x + t)^{(p-1)/2} - 1$ for random $t$), and **Rabin's irreducibility test**.
   - For $p \mid c_d$ the degree of $f \bmod p$ drops: find the affine roots of the lower-degree polynomial, and add the projective root separately.
2. **`internal/linalg`:**
   - Move `reducer`, `blockLanczos`/`symPrecond`, `denseDependencies`, and add a rows-based `findDependencies`.
   - **Parallelise Block Lanczos:** at 100–120 digits reduced matrices reach roughly 100k–500k columns with row weight 40–60.
     - Partition rows across workers for $B v$, `mulVecMat64` and `innerProd64`.
     - For $B^T v$ either keep per-worker destination buffers and XOR them together (8 bytes per column per worker), or store a transposed CSR copy.

---

### Stage 1: Polynomial Selection (`poly.go`)

Choose $f$ (irreducible, primitive) and $g(x) = a_1 x - a_0$ with common root $m$ mod $N$.

- **Degree:** implementations commonly use $d = 4$ around 90–110 digits and $d = 5$ above; check CADO-NFS's parameter files when tuning. The square root in Stage 6 works for any degree, so nothing restricts $d$ to odd values.
- **Irreducibility:** if $f \bmod p$ is irreducible for some $p \nmid c_d$, then $f$ is irreducible over $\mathbb{Z}$.
  - Test a few dozen primes with Rabin's test. For prime $d$ (3, 5) the Galois group contains a $d$-cycle, so about $1/d$ of primes pass.
  - For $d = 4$ some irreducible $f$ are reducible mod every prime. If no prime passes, check directly for rational roots and quadratic factors.
- **Phase 1A — Base-$m$ (M1):**
  - Pick $c_d$ (1 for monic) and set $m = \lfloor (N / c_d)^{1/d} \rfloor$, $g(x) = x - m$ ($a_1 = 1$, $a_0 = m$). Expand $N$ in base $m$, so $f(m) = N$ exactly, and balance the digits into $[-m/2, m/2]$ by carrying, which may move $c_d$ by ±1.
  - Not $N^{1/(d+1)}$: with $m = \lfloor N^{1/(d+1)} \rfloor$ we get $N \ge m^{d+1}$, and the expansion has degree $d+1$.
  - Because $f(m) = N$ exactly, a reducible $f = hk$ gives $N = h(m)k(m)$. Report that factor instead of rejecting the polynomial.
- **Phase 1B — Chosen leading coefficient, rotation and root sieve (M4):**
  - Search smooth small $c_d$ with $m = \lfloor (N/c_d)^{1/d} \rfloor$ (Murphy's base-$m$ variant). For each, size-optimise with translation $x \mapsto x + k$ (common root $m - k$) and rotation $f + (ux + v)g$; rotation keeps $f(m) \equiv 0$ since $g(m) = 0$, and the degree stays $d$ for $d \ge 2$.
  - Score with **Murphy's $\alpha(f)$** over small primes (about $p \le 200$) combined with the size of $F$ over the skewed region. Skewness $s \approx |c_0 / c_d|^{1/d}$ is the starting point; optimise $s$ by minimising the integral of $F^2$ over the skewed region.
  - **$\alpha(f)$ definition.** $\alpha$ compares the expected $p$-adic valuation of a random integer, $1/(p-1)$, with that of $F(a, b)$ over coprime pairs:
    $$\alpha(f) = \sum_{p \le B} \left(\frac{1}{p-1} - \mathbb{E}[v_p(F(a, b))]\right) \log p.$$
    $F(a, b)$ behaves like a number $e^{\alpha}$ times its size for smoothness, so negative $\alpha$ is better.
    - For $p \nmid \operatorname{disc}(f)$ with $n_p$ roots mod $p$, counting the projective root when $p \mid c_d$: over coprime $(a, b)$, $p \mid F$ with probability $n_p/(p+1)$, and the expected valuation is $n_p\,p/(p^2-1)$. So the term is
      $$\left(1 - \frac{n_p\, p}{p + 1}\right) \frac{\log p}{p - 1}.$$
    - For $p \mid \operatorname{disc}(f)$, roots can be multiple and that formula is wrong. Compute $\mathbb{E}[v_p(F)]$ by lifting the roots mod $p^k$ (Hensel, branching at multiple roots) until the contribution is negligible.
    - Signs are easy to get backwards. One educational implementation (elliptic1/General-Number-Field-Sieve) has a per-prime term that grows more positive with more roots and has no $1/(p-1)$ baseline. Test against brute force: the average of $v_p(F(a, b))$ over many random coprime pairs should match $\mathbb{E}[v_p(F)]$ for each small $p$. Also, a polynomial with no roots mod any small $p$ should get $\alpha = \sum \log p/(p-1) > 0$. CADO-NFS's polyselect prints $\alpha$ and can serve as a cross-check.
- **Phase 1C — Kleinjung (M5, optional):** non-monic $g$ ($a_1 \ne 1$) via Kleinjung's 2006 algorithm. This is a sizeable piece of work, so do it only if measurements show polynomial quality limits yield. Once $a_1 \neq 1$, the $a_1^{-|S|/2}$ factor in $X_{\text{rat}}$ (4.5) is exercised, so test it on polynomials with $a_1 \ne 1$ before relying on it.

---

### Stage 2: Factor Bases, Bad Ideals, and Quadratic Characters (`base.go`)

1. **Rational factor base:** all primes $p \le B_r$, with roots $a_0 a_1^{-1} \bmod p$, or projective ($p \mid b$) when $p \mid a_1$ (4.4).
2. **Algebraic factor base:** for each prime $p \le B_a$, the roots $r \in [0, p)$ of $f \bmod p$, plus the projective root $(p, \infty)$ when $p \mid c_d$. Each $(p, r)$ is its own column.
3. **Bad ideals** ($p \mid c_d$ or $p^2 \mid \operatorname{disc}(f)$, usually 5–20 small primes):
   - **M1 option:** one column per $(p, r)$ carrying $v_p(F(a, b))$, plus one extra QC per such column. The wrong valuations only add a finite-dimensional obstruction, which the extra characters remove with high probability. Measure how many dependencies then fail $\beta^2 = \gamma$. If it is more than a few, implement the exact rule.
   - **Exact:** compute valuations at the prime ideals above $p$ in $\mathcal{O}_K$, or use CADO-NFS's `badideals` approach (computed once per polynomial).
4. **Quadratic character base:** 50–64 pairs $(q, s)$ chosen as in 4.6, with $q$ above every algebraic large-prime bound and special-$q$.
5. **Fixed columns:** all-ones, rational sign, optional algebraic sign (4.5).
6. **Free relations (low priority):** for $p \nmid c_d \operatorname{disc}(f)$ where $f$ splits completely mod $p$, the ideal $(p)$ is the product of the $d$ ideals $(p, r_i)$. That gives a relation with rational prime $p$ and all $d$ algebraic ideals, without sieving.
   - The element is $p$ itself on both sides. It multiplies $\gamma$ by $p$ on the algebraic side, and adds $p$ with exponent 1 to $\prod G$ on the rational side.
   - It has no $(a, b)$, so it brings no factor of $c_d$ or $a_1$, and it does **not** set the all-ones bit. Its QC bits are $\left(\frac{p}{q}\right)$.
   - Only about $1/|\mathrm{Gal}|$ of primes qualify ($1/120$ for $S_5$), so the gain is small.

---

### Stage 3: Sieving (`sieve.go`)

#### M1: Line Sieve (validation baseline)
- For $b = 1, 2, \dots, B_{\max}$, sieve $a \in [-A, A)$ in 32 KB blocks (as in `qs`, `blockShift = 15`).
- For an affine $(p, r)$ the hits are $a \equiv br \pmod p$; a projective root adds $\log_2 p$ to the whole line when $p \mid b$.
- Skip $\gcd(a, b) \ne 1$ before cofactorization.

#### M3: Special-$q$ Lattice Sieve (production)
- **Special-$q$:** algebraic prime ideals $(q, r)$, usually from the upper part of the factor base up to a bit beyond $B_a$. The relation's algebraic norm then includes $q$, so record it.
- **Lattice:** $L = \{(a, b) : a \equiv br \pmod q\}$ is generated by $(q, 0)$ and $(r, 1)$, with determinant $q$. Gauss-reduce it under the skewed norm $a^2 + s^2 b^2$ to $\mathbf{u} = (u_1, u_2)$, $\mathbf{v} = (v_1, v_2)$, with $|u_1|, |v_1| \approx \sqrt{qs}$ and $|u_2|, |v_2| \approx \sqrt{q/s}$.
- **Region:** $(i, j) \in [-I, I) \times [1, J)$ with $(a, b) = i\mathbf{u} + j\mathbf{v}$, e.g. $I = 2^{11}$–$2^{13}$, $J = I$ or $I/2$.
- **Root transform:** for each $(p, r_p)$ with $p \ne q$, the condition $a \equiv r_p b$ becomes $i(u_1 - r_p u_2) \equiv j(r_p v_2 - v_1) \pmod p$, so
  $$i \equiv \lambda_p\, j, \qquad \lambda_p = (r_p v_2 - v_1)(u_1 - r_p u_2)^{-1} \bmod p.$$
  - If $u_1 - r_p u_2 \equiv 0 \pmod p$, the hits are the lines $j \equiv 0 \pmod p$. Both factors cannot vanish together, since $\det(\mathbf{u}, \mathbf{v}) = \pm q$ and $p \ne q$.
  - Projective roots transform the same way with $(r_p : 1)$ replaced by $(1 : 0)$.
- **Enumerating hits:**
  - $p < 2I$: line by line, starting at $i \equiv \lambda_p j \pmod p$ on each line $j$.
  - $p \ge 2I$: a line has a hit with probability only $2I/p$, so computing a start for every line costs $J$ operations per prime per special-$q$, about $2 \cdot 10^9$ at $J = 2048$ with a million primes. Use the **Franke–Kleinjung** algorithm instead: it walks the $p$-lattice in the flattened $(i, j)$ order using two precomputed basis vectors and emits each hit into the bucket for its block.
- **Coprimality and sign:**
  - Skip $(i, j)$ with both even (a quarter of the area), and those with $\gcd(i, j) \ne 1$.
  - $\gcd(i, j) = 1$ is necessary for $\gcd(a, b) = 1$ but not sufficient, since $a$ and $b$ can still share the factor $q$. So check $\gcd(a, b) = 1$ on the survivors too.
  - If $b < 0$, negate $(a, b)$; this keeps $a \equiv br \pmod q$. Skip $b = 0$.
- **Two sides:** sieve both sides over the whole region. The saving from ordering is only in the candidate scan: scan the side with the larger norm first and look at the other side only where the first passes.

---

### Stage 4: Two-Sided Cofactorization (`cofactor.go`)

For candidates passing both thresholds:
1. Reconstruct $(a, b)$, normalise $b > 0$, require $\gcd(a, b) = 1$.
2. Factor one side by small-prime scan plus resieving and bucket hits, and reject as soon as its cofactor exceeds the side's bound. The bound is the "mfb", the maximum cofactor bits, a tuning parameter.
3. Factor the other side the same way. On the algebraic side include the special-$q$, and record each prime as an ideal $(p, r)$ with $r = ab^{-1} \bmod p$ (4.4).
4. **Cofactor splitting (up to 2LP per side):**
   - `isPrime64` and `rho64`/`splitCofactor` handle cofactors up to the mfb. At 100 digits with large-prime bounds around $2^{27}$ that is about $2^{54}$.
   - For larger mfb or 3LP: `uint64` SQUFOF, P−1 and ECM, as CADO-NFS's cofactorization chain does.

---

### Stage 5: Filtering & Linear Algebra (`filter.go`, `internal/linalg`)

A relation can carry large primes on both sides, up to four in total.

1. **Deduplicate** on $(a, b)$ with $b > 0$, globally across all special-$q$. This is also the only cross-$q$ duplicate rule needed. Requiring every algebraic prime to be $\le q$ would lose relations whose largest algebraic prime lies outside the special-$q$ range.
2. **Singleton pruning on ideals:** a rational large prime is identified by $p$, an algebraic one by $(p, r)$. Two relations sharing $p$ with different roots are in different ideals, and matching them produces dependencies that are not squares.
   - Iteratively drop relations containing an ideal of weight 1 until a 2-core remains.
   - Track $\text{excess} = \#\text{relations} - \#\text{active large-prime ideals} - \#\text{fixed columns}$, where the fixed columns are the factor bases, QCs and the sign/all-ones columns.
3. **SGE + Block Lanczos:** build sparse rows from the all-ones, sign, RFB, AFB, QC and active large-prime columns, reduce with `reducer`, and solve the core with parallel Block Lanczos.

---

### Stage 6: Square Roots & Factor Extraction (`sqrt.go`)

The method is the one common to msieve and CADO-NFS: form $\gamma$ explicitly, take its square root $p$-adically from one prime, and read the integer coefficients off the lifted value. It works for any degree and has no sign problem: $\beta$ and $-\beta$ are both valid roots and give $\pm X_{\text{alg}}$, and the gcd step tries both signs anyway.

1. **Form $\gamma$** (4.5) as an element of $\mathbb{Z}[\hat\alpha] = \mathbb{Z}[y]/(\hat f)$, reducing mod $\hat f$ after each multiplication.
   - **M1:** multiply sequentially; at 30–60 digits the sizes are small.
   - **M5:** use a product tree.
2. **Coefficient bound $M_\beta$:**
   - Compute the complex roots $\hat\alpha_k$ of $\hat f$ once (Aberth–Ehrlich in `complex128`).
   - Bound $|\sigma_k(\beta)| = |\hat f'(\hat\alpha_k)| \prod_S |c_d a - b\hat\alpha_k|^{1/2}$ in log space.
   - Bound the coefficients through the inverse Vandermonde matrix: $|\beta_i| \le \sum_k |(V^{-1})_{ik}|\, |\sigma_k(\beta)|$.
   - Add a few bits of margin for floating-point error.
3. **Inert prime:**
   - Pick $p$ such that $\hat f \bmod p$ is irreducible (Rabin) and $p \nmid \operatorname{disc}(\hat f)$, so $\mathbb{Z}[\hat\alpha]/(p) = \mathbb{F}_{p^d}$.
   - For prime $d$ such $p$ have density about $1/d$.
   - If none exists ($d = 4$ with no 4-cycle in the Galois group), use a $p$ where $\hat f$ splits into distinct factors. Lift a root modulo each factor and combine them by CRT, trying the $2^{k-1}$ sign combinations for $k$ factors.
4. **Lift:**
   - Take $\gamma^{-1/2}$ in $\mathbb{F}_{p^d}$. Choosing $p \equiv 3 \pmod 4$ with $d$ odd gives $p^d \equiv 3 \pmod 4$, so a square root is a single power, $\gamma^{(p^d + 1)/4}$. Otherwise use Tonelli–Shanks or Cipolla in $\mathbb{F}_{p^d}$.
   - Lift with the Newton iteration $x \leftarrow x(3 - \gamma x^2)/2$, doubling the precision each step, until $p^K > 2 M_\beta$.
   - Then $\beta = \gamma x \bmod p^K$, with coefficients taken as symmetric residues.
   - Check $\beta^2 = \gamma$ (4.7). If it fails, the dependency was not a square, so log it and move on.
5. **Factor:** $X_{\text{alg}} = \beta(\hat m) \bmod N$. $X_{\text{rat}}$ is as in 4.5. Check $X_{\text{alg}}^2 \equiv X_{\text{rat}}^2 \pmod N$, then $\gcd(X_{\text{alg}} \pm X_{\text{rat}}, N)$.

**Cost.** This is a rough estimate to replace with measurements at M1 and M5.
- **Size:** at about 100 digits a dependency may hold $10^5$–$10^6$ original relations once large-prime cycles are expanded. $\beta$'s coefficients then run to millions of bits.
- **Multiplication:** Go's `math/big` uses Karatsuba with no FFT, so a multiplication of a few million bits takes on the order of 0.1s.
- **Total:** a product tree plus the Newton lift should take seconds to minutes, which is acceptable.
- **First experiment:** benchmark `big.Int` multiplication at 1–10 Mbit early.

---

## 6. Key Data Structures

```go
package gnfs

import "math/big"

type Poly struct {
	Degree int
	Coeffs []*big.Int // c_0 .. c_d
	A1, A0 *big.Int   // g(x) = A1*x - A0, gcd(A1, A0) = 1
	M      *big.Int   // A0 * A1^-1 mod N
	Skew   float64
	Alpha  float64
}

type AlgIdeal struct {
	P    uint32
	R    uint32 // root mod P, or P for the projective root (P | c_d)
	LogP byte
}

type LatticeQ struct {
	Q, R   uint32
	U1, U2 int64 // reduced lattice vector u
	V1, V2 int64 // reduced lattice vector v
	IHalf  int   // sieve i in [-IHalf, IHalf)
	JMax   int   // sieve j in [1, JMax)
}

// LargeIdeal is a large prime on the algebraic side. R = A * B^-1 mod P,
// which always exists since only primes dividing c_d can divide B.
type LargeIdeal struct {
	P, R uint64
}

type Relation struct {
	A, B    int64 // gcd(|A|, B) == 1, B > 0
	RatCols []int32
	AlgCols []int32 // including the special-q
	RatLPs  [2]uint64
	AlgLPs  [2]LargeIdeal
}
```

---

## 7. Implementation Milestones

| Milestone | Scope | Verification Target |
| :--- | :--- | :--- |
| **M-1: Baseline** | Time `qs` at 90, 95, 100 digits; benchmark `big.Int` multiplication at 1–10 Mbit | Recorded numbers that set the GNFS target and the square root design budget |
| **M0: Shared Refactor** | `internal/arith`, `internal/linalg`, rows-based `findDependencies`; polynomial roots mod $p$ and Rabin's test | `qs` tests and benchmarks unchanged; roots match brute force for $p \le 5000$ |
| **M1: End-to-End Prototype** | Base-$m$ ($d = 3$ and $4$), factor bases with coarse bad ideals, QCs, line sieve, 1LP cofactors, exact $\gamma$ + Hensel square root | Factor 30–60 digit semiprimes reliably; every check in 4.7 passes; count of dependencies failing $\beta^2 = \gamma$ |
| **M2: Filtering & Parallel LA** | Dedup, singleton pruning on $(p, r)$ ideals, `reducer` + parallel Block Lanczos | 150k × 150k sparse matrix benchmark; dependency count and wall time |
| **M3: Special-$q$ Lattice Sieve** | Gauss reduction, $\lambda_p$, Franke–Kleinjung buckets, resieving, two-sided 2LP | Factor 80–95 digit numbers; profile sieve vs cofactorization vs LA |
| **M4: Polynomial Selection** | Chosen $c_d$, translation/rotation, Murphy $\alpha$, skewness; checkpointing on `qs/relfile.go` | Yield vs M1 polynomials at 90–100 digits; crossover against `qs` |
| **M5: Scale & Tune (100–120+ digits)** | $d = 5$, parameter table (`Br`, `Ba`, large-prime bounds, mfb, `I`, `J`), product-tree square root, exact bad ideals, optional Kleinjung / P−1 / ECM cofactoring | 100 and 110 digit semiprimes end to end, with progress and time-left logging |

---

## 8. Critical Edge Cases & Pitfalls

1. **Large primes as ideals:** algebraic large primes are $(p, r)$, not $p$. Keying filtering or columns on $p$ alone produces dependencies that are not squares.
2. **Bad ideals ($p \mid c_d$ or $p^2 \mid \operatorname{disc}(f)$):** $v_p(F)$ is not a single ideal exponent there. Either handle them coarsely with extra QCs (M1) or compute exact valuations (M5).
3. **Non-coprime $(a, b)$:** reject $\gcd(a, b) \ne 1$; in the lattice sieve check it in addition to $\gcd(i, j) = 1$.
4. **Parity:** keep the all-ones column, so $|S|$ is even. Free relations must not set it.
5. **QC primes:** above every prime that can divide $F(a, b)$, not just above $B_a$.
6. **$a_1 \ne 1$:** $X_{\text{rat}}$ carries $a_1^{-|S|/2}$. Base-$m$ never exercises this, so test it as soon as non-monic $g$ exists.
7. **$c_d$ powers:** the algebraic norm carries $c_d^{(d-1)}$ per relation, the rational side $c_d$ per relation.
8. **Base-$m$ degree:** $m \approx (N/c_d)^{1/d}$, not $N^{1/(d+1)}$.
9. **Square root bound:** if $M_\beta$ is underestimated, the symmetric residues are wrong and $\beta^2 \ne \gamma$. Since that check catches it, keep a few bits of margin.

---

## Appendix A: Couveignes's Square Root (not chosen)

The first draft used Couveignes's method with Bernstein's explicit CRT. It avoids big integers by working modulo many inert 31-bit primes $q_j$ and reconstructing $\beta(\hat m) \bmod N$ directly. Its problems:

- **Cost:** each $q_j$ needs the full product $\gamma \bmod q_j$. That is about $|S| \cdot t \cdot d^2$ operations for $t$ primes. At 100 digits, $|S| \sim 10^6$ and $\beta$ has millions of bits, which needs $t$ of hundreds of thousands. That is about $10^{12}$–$10^{13}$ operations, compared with seconds to minutes for Stage 6. It is fine at M1 sizes, which hides the problem.
- **Odd degree only:** it needs the sign of $N(\beta_j)$ to pick between $\pm\beta_j$ consistently across primes, which works only for odd $d$ ($N(-\beta) = (-1)^d N(\beta)$). That rules out $d = 4$, the usual choice around 100 digits.
- **Formulas, if it is revisited:**
  - Pick the sign of $\beta_j$ so that $N(\beta_j) \equiv |N(\hat f'(\hat\alpha))| \cdot c_d^{(d-1)|S|/2} \prod p^{e_p/2} \pmod{q_j}$, with the product over $|F(a, b)|$. The first draft had $c_d^{|S|/2}$.
  - Choose $\prod q_j > 4 d M_\beta$, so that the explicit CRT rounding $\sum_j (v_j Q_j^{-1} \bmod q_j)/q_j$ stays clear of half-integers.
  - CADO-NFS's `crtalgsqrt` is a modern CRT-based square root worth reading before reviving this.

## Appendix B: Corrections to the First Draft

- **Base-$m$:** used $m = \lfloor N^{1/(d+1)} \rfloor$, which gives degree $d+1$. The correct choice is $m = \lfloor (N/c_d)^{1/d} \rfloor$.
- **$X_{\text{rat}}$:** used $a_1^{+|S|/2}$; the correct factor is $a_1^{-|S|/2}$ (4.5).
- **Couveignes norm:** used $c_d^{|S|/2}$; the correct factor is $c_d^{(d-1)|S|/2}$ (4.3, Appendix A).
- **Square root method:** Couveignes is replaced by Hensel lifting. Its cost grows too fast at the target size, and its odd-degree requirement ruled out $d = 4$.
- **QC bound:** QC primes need $q$ above every prime that can divide $F$, not $q > B_a$ (4.6).
- **Large primes:** algebraic large primes are now stored as ideals $(p, r)$.
- **Lattice sieve:**
  - Primes $p \ge 2I$ now use Franke–Kleinjung enumeration instead of a start per line.
  - $\gcd(i, j)$, both-even skipping and the sign of $b$ are now handled.
  - The "dual-sieve" step only saves in the candidate scan.
- **Dedup:** global dedup replaces the "algebraic primes $\le q$" rule.
- **Kleinjung description:** the confused aside is replaced. A cheaper chosen-$c_d$ + rotation phase now comes before Kleinjung.
- **Crossover:** the 95–105 digit estimate is replaced by a measurement milestone.
- **Additions:** the algebraic sign column (optional), free relations (low priority), and the inert-prime fallback for $d = 4$.

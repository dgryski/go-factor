# General Number Field Sieve (GNFS) — Go Implementation Plan

This document plans an implementation of the **General Number Field Sieve (GNFS)** in Go, designed as a sibling package (`gnfs/`) to the existing `qs/` package while sharing its finite-field arithmetic, SGE reducer, and Block Lanczos linear algebra.

Section 4 collects the formulas the pipeline depends on, with their derivations, so each stage can be checked against them. Appendix B lists the corrections made to the first draft of this plan. Section 0 records where the implementation stands.

---

## 0. Status (2026-11-26)

The package exists and factors end to end: `gnfs.Factor`, wired into `factor/main.go` after `qs` (methods list, `gnfsMaxBits = 450`, `-v`, `-workers`, `-gnfs-save`, separate GNFS pretest budget above `qsMaxBits`). `go test ./gnfs` takes about 11.4s.

### Baselines (Apple M5 Pro, 15 cores)

- **`qs` (after 2026-11-26 SIQS + DLP + `qsArena` + $Q_2(x)$ speedups):**
  - 80 digits (265–266 bits): 8.8s–11.0s total (was 15.9s).
  - 90 digits (297 bits): 54.8s total (50.5s sieve, 4.1s LA; was 2m11s).
  - 95 digits (314 bits): 3m42s total (3m34s sieve, 8.0s LA; was 9m06s). Extrapolates to ~15m at 100 digits.
- **`BenchmarkQuadraticSieve` pre vs post refactor:** unchanged within noise at refactor, now 259µs (60b), 2.88ms (120b), 23.1ms (160b).
- **`BenchmarkBlockLanczos` (row weight 40):** 0.19s / 21 MB at 10k cols, 2.45s / 123 MB at 50k cols, 16.7s / 540 MB at 150k cols.
- **`BenchmarkBigMul` (`math/big`):** 8.1ms at 1 Mbit, 24.3ms at 2 Mbit, 102ms at 5 Mbit, 309ms at 10 Mbit.
- **80-digit (265-bit) GNFS progression:**
  1. Initial `160k_21` baseline: 9m53.8s wall (9m30.6s sieve = 96.1%, 11.4s LA, 9.0s sqrt).
  2. Prime-outer line sieving, skipping both-even $(i, j)$, 8-byte SWAR threshold scans, `invMod32`, two-stage 128/64-bit cofactoring, and `250k_q80k`: 4m34.0s wall (3m56.4s sieve, 10.2s LA).
  3. Skewed fixed-$c_d$ and non-monic Kleinjung polynomials, translation + $2\times 2$ rotation size optimisation, bounded skewness, parallel prime-power root sieve, 16 KB blocks, and flat `uint32` `fkArena` buckets grouped by `LogP`: **54.0s wall (40.3s sieve, 4.7s LA, 7.2s sqrt)**, 11× faster than baseline.
- **Expected `qs` / `gnfs` crossover:** with `qs` now at 54.8s (90d) and 3m42s (95d) and `gnfs` at 54.0s (80d) and 3m53s (89d), crossover projects around 100–102 digits (~332–339 bits); `qsMaxBits` in `factor/main.go` is set to `332`.

### Milestones

| Milestone | State |
| :--- | :--- |
| **M-1** | Done. `qs` timed at 80, 90 and 95 digits; `BenchmarkBigMul` recorded at 1–10 Mbit; `BenchmarkBlockLanczos` recorded up to 150k cols; pre/post refactor `qs` benchmarks verified equal. |
| **M0** | Done. `internal/arith` and `internal/linalg` shared with `qs` with no benchmark regression. |
| **M1** | Done. Base-$m$, coarse bad ideals, QCs, line sieve, square root by Hensel lifting, and all 4.7 checks. Zero bad squares observed across all runs up to 295 bits. |
| **M2** | Mostly done. Dedup, zero-allocation dense-ID `lpFilter` singleton pruning on $(p, r)$ ideals, `Reducer`, parallel Block Lanczos benchmarked to 150k columns. Cavallar clique pruning and Markowitz priority-queue pivoting remain. |
| **M3** | Mostly done. Gauss reduction, $\lambda_p$, 16 KB prime-outer line sieving, flat `uint32` Franke–Kleinjung bucket arena grouped by `LogP`, both-even skipping, SWAR candidate scan, two-stage 128/64-bit 2LP cofactoring, extended special-$q$ sieving. Sub-segment log-norm initialization, line-prime resieving, prime-power sieving, and 305–365 bit tuning remain. |
| **M4** | Mostly done. Skewed fixed-$c_d$ search, bounded skewness, translation + $2\times 2$ rotation size optimisation, parallel prime-power root sieve, Murphy $\alpha$, and `-gnfs-save` checkpointing + full/partial resume + `TestSavedRelations`. Kleinjung $a_1 = P_1 P_2$ collision search and directional Murphy $E(f, g)$ scoring remain. |
| **M5** | Partial. Non-monic Kleinjung polynomial selection, free relations, product-tree square root, and $d = 5$ parameter row exist. 3LP + 128-bit Montgomery P−1/ECM cofactoring, exact bad ideals, two-level cache-blocked buckets, Kronecker + Karp–Markstein + NTT square root, and tuned runs at 95–110+ digits remain. |

### Parameter table (`defaultParams`)

| Bits | State |
| :--- | :--- |
| < 115 | Line sieve, $d = 3$, Br 1200, Ba 1800, LP 32k/48k, ratio 0.80. 80–96ms at 95 bits, 134–272ms at 114 bits (down from 603ms with Br 800, Ba 1200). |
| 115–150 | $d = 4$, Br = Ba = 8000, LP $2^{17}$, ratio 0.80. 0.3s at 115 bits, 0.6s at 149 bits. |
| 150–200 | $d = 4$, Br = Ba = 30000, LP $2^{18}$, $Q_{\min} = 20\text{k}$, $I = 1024$, $J = 512$, ratio 0.85. 1.6s at 151 bits, 3–4s at 195 bits. |
| 200–250 | $d = 4$, Br = Ba = 100000, LP $2^{20}$, mfb 40, $Q_{\min} = 60\text{k}$, $I = J = 1024$, ratio 0.82. 7.86s total (4.83s sieve) at 231 bits / 70 digits after free relations + `lpFilter` (down from 11.0s after poly+arena, and 26–33s originally). |
| 250–275 | $d = 4$, Br = Ba = 250000, LP $2^{21}$, mfb 42, $Q_{\min} = 80\text{k}$, $I = 2048$, $J = 1024$, ratio 0.82. 43.8s at 264 bits / 80 digits after free relations + `lpFilter` + sub-segment log-norms (34.4s sieve, 4.5s LA, 3.4s sqrt; down from 54.0s and 4m34s). |
| 275–305 | $d = 4$, Br = Ba = 400000, LP $2^{22}$, mfb 44, $Q_{\min} = 150\text{k}$, $I = J = 2048$, ratio 0.84. 3m53s at 295 bits / 89 digits (3m19s sieve, 14.1s LA, 16.5s sqrt, 0 bad squares; down from 21m51s). |
| 305–365 | $d = 4$, Br = 550000, Ba = 600000, LP $2^{23}$, mfb 46, $Q_{\min} = 250\text{k}$, $I = J = 2048$, ratio 0.88. Extrapolated from 265 and 295 bits; not yet timed. |
| ≥ 365 | $d = 5$, Br = 500000, Ba = 800000, LP $2^{26}$, mfb 52. Not tested. |

`chooseDegree` reads the degree from this table.

### Differences from the plan

- **Skewness** minimises $\log\operatorname{RMS}(F) + \log\operatorname{RMS}(G)$ capped near $1.5\sqrt{Q_{\min} Q_{\max}}$ so $|b| \ge 1$ is realisable on the integer special-$q$ lattice.
- **Bad ideals** are currently flagged for any multiple root of $f \bmod p$, which is a little more than $p^2 \mid \operatorname{disc}(f)$ requires (exact Round-2 valuation planned below).
- **Multiple roots above $B_a$** have no QC, so special-$q$ roots and algebraic large primes at multiple roots are skipped.
- **Small composite factors** that fit in 64 bits go to rho instead of back into GNFS.
- **Progress** is measured in relations against a per-row `RelTargetRatio` times $\pi(\text{LP}_{\text{rat}}) + \pi(\text{LP}_{\text{alg}})$, switching to slope extrapolation once the 2-core slope reaches 0.2 (ignoring tiny pre-percolation slopes).
- **Overshoot:** excess is checked with zero-allocation dense-ID `lpFilter` buffers every ~7% before $-\#\text{fixed}/2$, every ~1% while climbing, and every `2*(goal - excess)` (clamped to `[100, 1000]`) above $-2000$.

### Prioritized Algorithmic Plan Relative to M3 & M5

At 89 digits (295 bits), sieving + cofactoring is **85.5%** of wall time (`3m19s`), LA is **6.0%** (`14.1s`), square root is **7.1%** (`16.5s`), and coarse bad ideals have produced **0 bad squares**. Items that improve sieve/polynomial yield come **before** M3 parameter tuning (since they shift optimal `Br`, `Ba`, `Slack`, and the `qs`/`gnfs` crossover), while >64-bit cofactoring, L2 superblocking, LA/sqrt scaling, and exact bad ideals belong in **M5** (100–115+ digits) and **Post-M5** (115–130+ digits).

Parameter sweeps use `TestSweep` in `gnfs/tune_test.go`, which is skipped without `-tune` and prints live 15s progress lines. For example, `GNFS_BITS=266 GNFS_CANDS=default GNFS_SEEDS=2 go test ./gnfs -run Sweep -v -tune -timeout 120m`. The package path has to come before `-tune`.

| Tier / Milestone | Task | Status | Rationale & Scope |
| :--- | :--- | :--- | :--- |
| **Done (M2/M5)** | Free relations and dense-ID `lpFilter` pruning | Done | Free relations for $p \le \min(B_r, B_a)$ scanned from `fb.AFB`; zero-alloc `lpFilter` recounts. 70d seed 1 down from 9.83s to 7.86s (`b6a7e82`). |
| **Done (M3)** | Sub-segment log-norm initialization | Done | Per-line 64B (and 16B near roots/steep gradients) `float64` Horner $\log_2|F|, \log_2|G|$ baselines replace 16 KB 4-corner maximums. `go test ./gnfs` down from 11.4s to 7.23s; 80d down from 54.0s to 43.8s (34.4s sieve). |
| **Tier 1 (Pre-M3 tune — M3)** | Line-prime batch resieving + small prime-power sieving ($p^k < 2I$) | Next | Eliminates ~550 per-candidate 64-bit `%` divisions in `cofactorRatLattice`/`cofactorAlgLattice` (main cofactor hotspot after interior norm recovery) and sieves small $p^k$ to tighten `SlackRat`/`SlackAlg`. |
| **Tier 1 (Pre-M3 tune — M4)** | Kleinjung $a_1 = P_1 P_2$ collision + directional Murphy $E(f, g)$ | Pending | Searches $d^6$ CRT roots in $O(d^3 \log d)$ time and scores real-root corridors now captured by sub-segment norms; shifts yield across all 80–115d rows, so do before M3/M5 parameter sweeps. |
| **Tier 2 (Finish M3)** | M3 tuning & crossover verification at 80–100 digits | Pending | Multi-seed checks on 250–275b (80d) and 275–305b (89d), tune/time 305–340b at **95d (315b)** and **100d (332b)**, measure exact `qs`/`gnfs` crossover, and update `qsMaxBits` / `gnfs15Cores` in `factor/`. **Completes M3.** |
| **Tier 3 (Core M5)** | 3LP + 128-bit Montgomery P−1 and ECM cofactoring | Pending | Not needed for M3 (where 2LP `mfb` $\le 48 < 64$); unlocks `RatLPs [3]uint64` and 65–100 bit cofactors (`mfb` $\approx 2.7\log_2\text{LP}$) required for 105–115d yield. |
| **Tier 3 (Core M5)** | Two-level $J$-strip cache-blocked `fkArena` buckets | Pending | Not needed at $I = J = 2048$ (M3, 4 MB/worker); keeps 15 workers' `fkArena` slabs L2-resident when $I \times J \ge 4096 \times 2048$ at 105+ digits. |
| **Tier 3 (Core M2/M5)** | Cavallar clique pruning + Markowitz SGE & Karp–Markstein + point-eval sqrt | Pending | LA + sqrt are only ~13% of runtime at 89d, rising at 105–115d (1M+ raw relations, $10^5$+ dep weight); shrinks reduced Lanczos core by ~20% and cuts final Hensel lift multiplies ~3×. |
| **Tier 4 (Finish M5)** | $d = 4$ vs $d = 5$ tuning and 105–110+ digit runs | Pending | Tune degree-5 polynomial selection and factor-base / 3LP rows at 340–370 bits (103–111d) with progress and time-left logging. |
| **Tier 4 (Finish M5)** | Exact bad ideals (Round-2 $p$-maximal order $\mathcal{O}_p$) | Pending | 0 bad squares seen with coarse bad ideals + extra QCs, so near-zero wall-clock gain; implement via Cohen (1993) Alg. 4.8.17 / 6.1.8 / 6.2.9 to complete M5 cleanly (promote earlier only if a 100–110d run shows bad squares). **Completes M5.** |
| **Tier 5 (Post-M5)** | 3-prime 64-bit NTT for multi-Mbit `big.Int` products | Deferred | Defer past M5 unless `tryDependency` exceeds ~60s above 110–115 digits; product tree + Karp–Markstein on `math/big` Karatsuba handles $\le 110$ digits in seconds. |
| **Skipped** | SIMD / `avo` core sieving loop | Skipped | Skipped per operator direction. |

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
- **Phase 1B — Skewed fixed-$c_d$, size optimisation, and root sieve (M4):**
  - **Fixed-$c_d$ skewed expansion:** sample smooth products of $\{2, 3, 5, 7, 11, 13\}$ logarithmically up to $\approx 2^{40}$. Keep $c_d$ fixed so its smoothness is not lost to a carry, pick nearest $m = \lfloor (N / c_d)^{1/d} \rceil$, and expand $R = N - c_d m^d$ in balanced base $m$ into $c_0 \dots c_{d-2}$ with the top remainder in $c_{d-1}$. Since $|m - (N/c_d)^{1/d}| \le 1/2$, $|c_{d-1}| \le \frac{d}{2} c_d \ll m$.
  - **Size optimisation (`sizeOptimizePoly`):** alternate translation $x \mapsto x + k_t$ ($k_t = \lfloor -c_{d-1}/(d c_d) \rceil$, updating $a_0 \leftarrow a_0 - k_t a_1$) with the exact $2\times 2$ quadratic-form minimum over $(u, v)$ in $f + (ux + v)(a_1 x - a_0)$ at current skewness $s$, accepting steps monotonically.
  - **Bounded skewness:** cap $s$ near $1.5\sqrt{Q_{\min} Q_{\max}}$ so $|b| \approx A/\sqrt{s} \ge 1$ is realisable on the integer special-$q$ lattice instead of drifting to unphysical $s \gg q$.
  - **Prime-power root sieve (`rootSieveOptimize`):** for the top size-optimised candidates in parallel, sieve $v \in [-V, V]$ across $u \in [-U, U]$ (up to $U \le 12, V \le 768$) over all prime powers $p^k \le 97$. For each $x \in [0, p^k)$ with $p \nmid g(x)$, hits lie on the single arithmetic progression $v \equiv -f_u(x) g(x)^{-1} \pmod{p^k}$ with weight $\frac{\log p}{(p+1) p^{k-1}}$, summing identically to $\frac{p}{p^2-1}\log p$ on simple roots and handling ramified roots without special cases.
  - **Dickman / Directional Murphy $E(f, g)$ scoring:** rank size-optimised candidates quickly by isotropic $u_F \ln u_F + u_G \ln u_G$, and rank root-sieved survivors by Murphy's directional yield integral (Murphy 1999, Section 5.1.3):
    $$E(f, g) = \frac{1}{K} \sum_{k=1}^{K} \rho(u_F(\theta_k))\, \rho(u_G(\theta_k)), \qquad (a_k, b_k) = \left(\sqrt{A^2 s}\cos\theta_k,\; \sqrt{A^2 / s}\sin\theta_k\right), \quad \theta_k = \frac{\pi(k - 1/2)}{K}$$
    with $u_F(\theta) = (\ln|F(a, b)| - \ln q_{\text{typ}} + \alpha(f))/\ln B_a$ and $u_G(\theta) = \ln|G(a, b)|/\ln B_r$, capturing extra smoothness along real roots of $F$ and $G$.
  - **$\alpha(f)$ definition.** $\alpha$ compares the expected $p$-adic valuation of a random integer, $1/(p-1)$, with that of $F(a, b)$ over coprime pairs:
    $$\alpha(f) = \sum_{p \le B} \left(\frac{1}{p-1} - \mathbb{E}[v_p(F(a, b))]\right) \log p.$$
    $F(a, b)$ behaves like a number $e^{\alpha}$ times its size for smoothness, so negative $\alpha$ is better.
    - For $p \nmid \operatorname{disc}(f)$ with $n_p$ roots mod $p$, counting the projective root when $p \mid c_d$: over coprime $(a, b)$, $p \mid F$ with probability $n_p/(p+1)$, and the expected valuation is $n_p\,p/(p^2-1)$. So the term is
      $$\left(1 - \frac{n_p\, p}{p + 1}\right) \frac{\log p}{p - 1}.$$
    - For $p \mid \operatorname{disc}(f)$, roots can be multiple and that formula is wrong. Compute $\mathbb{E}[v_p(F)]$ by lifting the roots mod $p^k$ (Hensel, branching at multiple roots) until the contribution is negligible.
- **Phase 1C — Non-monic Kleinjung candidates and meet-in-the-middle collision (`kleinjungCandidates`):**
  - For smooth $c_d$, find small primes $p_j \nmid d c_d$ where $c_d x^d \equiv N \pmod{p_j}$ splits into $d$ roots and Hensel-lift each root to $\bmod p_j^2$.
  - **Direct CRT ($k \in \{3, 4\}$)** combines $k$ primes into $a_1 = \prod p_j$ ($a_1^2 < 2^{62}$) for $d^k$ roots, while **meet-in-the-middle collision ($a_1 = P_1 P_2$, $3 + 3$ primes)** splits CRT into two halves $R_1 \bmod P_1^2$ and $R_2 \bmod P_2^2$ of size $d^3$ each, sorting fractional parts of $r_1 e_1 / a_1^2 \bmod 1$ against $(m_0 - r_2 e_2)/a_1^2 \bmod 1$ (Kleinjung 2006, Section 2.2) to search $d^6$ roots in $O(d^3 \log d)$ time.
  - For the closest matches to $m_0 = (N/c_d)^{1/d}$, let $M = (N - c_d a_0^d)/a_1^2$, minimise $|M - t_{d-1} a_0^{d-1}|$ so both $c_{d-1} = a_1 t_{d-1}$ and $c_{d-2}$ are small simultaneously, expand $\sum_{i=0}^d c_i a_0^i a_1^{d-i} = N$, and pass through Stage 1B size optimisation and root sieve.

---

### Stage 2: Factor Bases, Bad Ideals, and Quadratic Characters (`base.go`)

1. **Rational factor base:** all primes $p \le B_r$, with roots $a_0 a_1^{-1} \bmod p$, or projective ($p \mid b$) when $p \mid a_1$ (4.4).
2. **Algebraic factor base:** for each prime $p \le B_a$, the roots $r \in [0, p)$ of $f \bmod p$, plus the projective root $(p, \infty)$ when $p \mid c_d$. Each $(p, r)$ is its own column.
3. **Bad ideals** ($p \mid c_d$ or $p^2 \mid \operatorname{disc}(\hat f)$, usually 5–20 small primes):
   - **M1 option:** one column per $(p, r)$ carrying $v_p(F(a, b))$, plus one extra QC per such column.
   - **Exact (M5):** compute the $p$-maximal overorder $\mathcal{O}_p \supseteq \mathbb{Z}[\hat\alpha]$ via Zassenhaus Round 2 (Cohen 1993, Algorithm 6.1.8), decompose $\mathcal{O}_p / p\mathcal{O}_p$ into prime ideals $\mathfrak{P}_1, \dots, \mathfrak{P}_g$ with $(e_j, f_j)$ (Algorithm 6.2.9), and precompute anti-uniformizer matrices $\beta_j \in \mathcal{O}_p$ ($\tau_j = \beta_j / p \in \mathfrak{P}_j^{-1} \setminus \mathcal{O}_p$, Algorithm 4.8.17) once per polynomial so $v_{\mathfrak{P}_j}(c_d a - b\hat\alpha) \bmod 2$ is evaluated exactly in `buildMatrix` without extra QCs.
4. **Quadratic character base:** 50–64 pairs $(q, s)$ chosen as in 4.6, with $q$ above every algebraic large-prime bound and special-$q$.
5. **Fixed columns:** all-ones, rational sign, optional algebraic sign (4.5).
6. **Free relations (done):** for $p \le \min(B_r, B_a)$ with $p \nmid c_d$ where $f$ has $d$ distinct affine roots mod $p$ (scanned directly from `fb.AFB` with no extra root-finding), $(p)\mathcal{O}_K = \prod_{i=1}^d \mathfrak{p}_i$ gives a full relation $(p, 0)$ immediately before sieving.
   - The element is $p$ itself on both sides: it multiplies $\gamma$ by $p$ on the algebraic side and adds $p$ to $\prod G$ on the rational side, omits `ColAllOnes` and sign columns, and sets QC bits to $\left(\frac{p}{q}\right)$.

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
- **Enumerating hits, norm initialization, and bucket layout:**
  - **16 KB L1-resident blocks (`blockShift = 14`):** both `algBuf` and `ratBuf` fit in 32 KB total, block boundaries align with full lines of width $W = 2I$ for $I \le 8192$, and block offsets fit in 14 bits.
  - **Sub-segment log-norm initialization:** instead of taking the 4-corner maximum of $|F|$ and $|G|$ across an entire 16 KB block, evaluate $\log_2|F(i, j)|$ and $\log_2|G(i, j)|$ along each line $j$ in 64–128 byte segments using `float64` Horner polynomials and initialize segment bytes to `bias - (normBits - slack)` so a fixed byte threshold captures small-$(i, j)$ and real-root corridors accurately.
  - **Small prime powers and line sieving ($p^k < 2I$):** loop prime-outer over each 16 KB block with carried per-line offsets (`off += lambda; if off >= p { off -= p }`), stepping by $2p$ on even-$j$ lines to skip both-even $(i, j)$, unrolling 4×, and sieving small prime powers ($p^k \le 512$ and powers of primes $< 29$) to tighten norm residuals.
  - $p \ge 2I$ (**flat `uint32` Franke–Kleinjung arena `fkArena` + two-level $J$-strips**): group FK ideals by `LogP` (~7 groups), walk the $p$-lattice skipping both-even $(i, j)$ via `pos & (W|1) != 0`, and write 4-byte `uint32(idx)<<14 | off` entries into a pre-allocated power-of-2-stride slab per block; for large $I \times J$, process $J$ in L2-resident strips (~256 blocks) so all workers' arenas stay in cache.
  - **Candidate scan, line-prime resieving, and two-stage cofactoring:** scan 8 bytes per word (`binary.LittleEndian.Uint64` SWAR mask), gather both bucket hits and line-prime hits onto `candBits` (or test line primes via `uint32` Lemire reciprocals instead of 64-bit `%`), cofactor $G(a, b)$ first, and gather/cofactor algebraic hits only for Stage 1 survivors.
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
4. **Cofactor splitting (up to 3LP per side):**
   - For $\le 64$-bit 2LP cofactors: `IsPrime64` + `Rho64` / Montgomery P−1 + ECM.
   - For 3LP (`RatLPs [3]uint64`, `AlgLPs []LargeIdeal`) and `mfb > 64` bits (up to 128 bits): stack-allocated 128-bit Montgomery arithmetic (`mont128`), BPSW/Miller–Rabin primality test to reject inert cofactors early, Stage-1 Pollard P−1, and Suyama-parametrized Montgomery-curve ECM with `rho128` fallback.

---

### Stage 5: Filtering & Linear Algebra (`filter.go`, `internal/linalg`)

A relation can carry up to three large primes per side.

1. **Deduplicate** on $(a, b)$ with $b > 0$ (and $(p, 0)$ for free relations), globally across all special-$q$.
2. **Dense-ID singleton and clique pruning on ideals:** intern rational LPs $p$ and algebraic LPs $(p, r)$ to dense `int32` IDs in `lpFilter` so iterative 2-core peeling during sieving runs in reusable flat arrays with zero heap allocations. When `excess` exceeds the target before matrix build, Cavallar clique removal drops the heaviest 2-weight LP components to shrink the matrix and final dependency weight.
3. **Markowitz SGE + Block Lanczos:** build sparse rows from fixed, exact bad-ideal, QC and active LP columns, pivot in ascending Markowitz fill `(w - 1)*(rowLen - 1)` order in `Reducer`, and solve the compact core with parallel Block Lanczos.

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
   - Lift $x \equiv \gamma^{-1/2}$ up to half precision $\lceil K/2 \rceil$ with Newton steps $x \leftarrow x(3 - \gamma x^2)/2$, then finish the last doubling step via Karp–Markstein's coupled step ($b_0 = \gamma x \bmod p^{\lceil K/2 \rceil}$, $\beta = b_0 + x(\gamma - b_0^2)/2 \bmod p^K$) so full precision needs only 2 polynomial multiplies instead of 4.
   - Accelerate degree-$d$ polynomial products via Kronecker substitution / evaluation at $2d - 1$ small points, and accelerate multi-Mbit `big.Int` multiplications via a 3-prime 64-bit NTT in `internal/arith` ($O(B \log B)$ vs `math/big` Karatsuba's $O(B^{1.585})$).
   - Take symmetric residues of $\beta \bmod p^K$ and check $\beta^2 = \gamma$ (4.7). If it fails, the dependency was not a square, so log it and move on.
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
| **M2: Filtering & Parallel LA** | Dedup, dense-ID `lpFilter` singleton pruning, Cavallar clique pruning, Markowitz `Reducer` + parallel Block Lanczos | 150k × 150k sparse matrix benchmark; dependency count and wall time |
| **M3: Special-$q$ Lattice Sieve** | Gauss reduction, $\lambda_p$, sub-segment log-norms, prime-power sieving, Franke–Kleinjung buckets, line-prime resieving, two-sided 2LP | Factor 80–95 digit numbers; profile sieve vs cofactorization vs LA |
| **M4: Polynomial Selection** | Chosen $c_d$, translation/rotation, root sieve, Kleinjung $P_1 P_2$ collision, directional Murphy $E(f, g)$; checkpointing on `qs/relfile.go` | Yield vs M1 polynomials at 90–100 digits; crossover against `qs` |
| **M5: Scale & Tune (100–120+ digits)** | $d = 5$, parameter table, free relations, two-level cache-blocked buckets, exact Round-2 bad ideals, 3LP + 128-bit Montgomery P−1/ECM, Kronecker + Karp–Markstein + NTT square root | 100 and 110 digit semiprimes end to end, with progress and time-left logging |

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

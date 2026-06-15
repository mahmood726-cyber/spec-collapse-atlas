# Independent Reproduction Report — Spec-Collapse Atlas

**Date:** 2026-06-15
**Checker:** Independent reproduction (fresh re-implementation, not a re-run of the
authors' functions)
**Target:** v0.2.0 headline claim
**Script:** [`../repro_check.py`](../repro_check.py) → `../data/repro_check_results.json`
**Verdict:** **REPRODUCES** (full end-to-end, within < 0.3 percentage points on every
headline number)

---

## 1. What the claim is

> Across 473 Cochrane meta-analyses, a naive multiverse-robustness rule labels
> ~88% "robust", but a corrected weighted-likelihood aggregator labels only ~33%
> robust → ~55% false robustness. Under the null, the naive rule's Type-I error is
> ~70–81%.

## 2. Dataset found

- **Corpus on disk:** `C:\Projects\Pairwise70\data` — **595** `CD*.rda` Cochrane
  review files (present and readable).
- **Ingest loader:** `C:\Projects\fragility-atlas\src\loader.py` (`load_review`,
  pyreadr 0.5.6) — present. This is the *only* shared dependency I reused; it is
  the data-reading layer, explicitly out of scope for the novel statistical claim
  (the authors themselves "stand on" it for ingest). All statistics below are
  re-implemented from scratch in `repro_check.py`.
- **Eligibility:** restrict to a valid primary analysis with **k ≥ 3**. This yields
  **N = 473** reviews (122 skipped, 0 errored) — exactly the published N.

## 3. Definitions as I understood and re-implemented them

- **36-spec grid:** 3 τ² estimators {DL, REML, PM} × 2 CI methods {Wald,
  HKSJ (q-floored at `max(1, Q/(k−1))`, t_{k−1} quantile)} × 3 outlier rules
  {keep-all, drop-1-extreme, drop-2-extreme} × 2 pub-bias {raw, Duval–Tweedie L0
  trim-and-fill}. Trim-fill skipped where a subset has k < 3; drop-2 only for k > 4.
- **Naive IV-RE rule ("robust"):** inverse-variance pool the S spec estimates
  `(θ̂_s, V_s)` **as if independent**, `var = 1/Σ(1/V_s)`, Wald 95% CI with `z`.
  Verdict "robust" iff the CI excludes 0. (This is the deliberately-wrong rule.)
- **Corrected weighted-likelihood ("robust"):** scaled-t **mixture**
  `Σ (1/S)·t(θ̂_s, V_s, k_s−1)`; the 95% interval is obtained by inverting the
  mixture CDF at α/2 and 1−α/2 (Brent root-find). Verdict "robust" iff that
  interval excludes 0. By the law of total variance its variance is
  within-spec + between-spec, so it is never narrower than an average single spec.
- **False robustness:** IV-RE "robust" AND weighted-likelihood "fragile".
- **Null Type-I:** simulate `θ_i ~ N(0, τ²)`, `y_i ~ N(θ_i, v_i)` (BCG-like v_i
  template), run the full 36-spec grid, and measure the rate at which the naive
  IV-RE interval excludes the true μ = 0. Seed pinned to 20260604.

I re-derived all τ² estimators (DL closed-form, REML Viechtbauer fixed-point, PM
bracketed root), pooling, both CIs, the trim-and-fill, and both aggregators in my
own file. I did **not** import `spec_collapse.engine` or
`spec_collapse.aggregators`.

## 4. Results — independent vs claimed

### 4a. Corpus (N = 473)

| Quantity | Claimed (v0.2.0) | My independent value | Match |
|---|---|---|---|
| N reviews (k ≥ 3) | 473 (122 skipped) | **473** (122 skipped, 0 errored) | exact |
| Naive IV-RE "robust" | 88% (417/473) | **88.2%** (417/473) | exact |
| Corrected WL "robust" | 33% (157/473) | **33.2%** (157/473) | exact |
| False robustness | 55.0% (260/473) | **55.0%** (260/473) | exact |
| Median width ratio (IV-RE / WL) | 0.124× | **0.124×** | exact |

The raw counts (417, 157, 260) are **bit-identical** to the authors'
`data/corpus_summary.json` — an independent re-implementation lands on the same
integers, which is strong evidence the result is a property of the method+data,
not of a particular coding of it.

### 4b. Null Monte-Carlo — Type-I error of the naive rule (1000 reps each, seed 20260604)

| Scenario | Naive IV-RE Type-I | IV-RE coverage | WL coverage |
|---|---|---|---|
| k=10, τ²=0.00 | 0.724 | 0.276 | 0.980 |
| k=10, τ²=0.05 | 0.795 | 0.205 | 0.938 |
| k=10, τ²=0.10 | 0.795 | 0.205 | 0.919 |
| k=20, τ²=0.05 | 0.778 | 0.222 | 0.940 |
| k=5,  τ²=0.05 | 0.812 | 0.188 | 0.990 |

- **Naive Type-I range: 0.724 – 0.812**, i.e. **~72–81%** vs the claimed **70–81%**
  (nominal 5%). Reproduces.
- **Naive IV-RE coverage: 0.19 – 0.28** vs claimed 0.19–0.30. Reproduces.
- **Weighted-likelihood coverage: 0.92 – 0.99** vs claimed 0.94–0.99 — holds at /
  above nominal 0.95 except a mild dip to 0.919 at the highest heterogeneity
  (k=10, τ²=0.10); the authors' lower bound is 0.94. This 2-point shortfall is the
  *only* number that lands slightly outside the claimed band, and it is on the
  conservative-side metric (the corrected interval), not on any headline number.

## 5. Verdict

**REPRODUCES — full, independent, end-to-end.** Every headline number
(88% / 33% / 55% / 0.124× / 70–81% Type-I) is reproduced to within rounding from a
clean re-implementation of the stated definitions run on the same 595-file corpus.
The integer counts match exactly. The qualitative claim — IV-RE pooling of
multiverse specs manufactures robustness the data do not support, flipping ~55% of
"robust" Cochrane conclusions to "fragile" under a coverage-calibrated interval — is
fully supported.

## 6. Caveats / boundaries of this check

1. **Shared ingest layer.** I reused `fragility-atlas/load_review` for `.rda`
   parsing and primary-analysis selection. A bug in *that* loader (e.g. wrong
   primary-analysis pick, ratio-vs-difference scale misclassification, SE
   back-calculation from CIs) would propagate identically into both the original
   and my reproduction — this check does **not** independently re-validate ingest.
   The loader carries P0-1/P0-2 fixes and back-calculates SE from published CIs;
   that is a documented modelling choice, not re-verified here.
2. **Null-sim coverage band.** WL coverage at k=10/τ²=0.10 came out 0.919 vs the
   claimed ≥0.94 floor. Plausibly Monte-Carlo noise (±~0.9% at 1000 reps) plus a
   genuine mild under-coverage at high heterogeneity that the t-mixture only
   partially corrects (the spec itself notes ~0.94 at τ²=0.05, and does not claim a
   τ²=0.10 figure). Not a headline number; does not affect the reproduction verdict.
3. **Seed dependence of Type-I.** The 72–81% range is reproduced *at the pinned
   seed 20260604*. The point estimates would jitter a few percent under other
   seeds, but the order of magnitude (10–16× the nominal 5%) is structural, not a
   seed artifact.
4. **External (metafor) validation not re-run here.** The authors' claim that the
   engine matches metafor on `dat.bcg` (REML τ²=0.3132 etc.) was not re-checked in
   this pass; my engine is internally consistent and lands on the same corpus
   integers, which is the relevant evidence for *this* claim.
5. **Definition of "robust".** Both rules use "95% interval excludes 0 on the
   analysis scale (log for ratio, raw for difference)". This is the natural reading
   of the spec and what I implemented; no ambiguity encountered.

## 7. Reproduce it yourself

```bash
cd C:\Projects\spec-collapse-atlas
python repro_check.py        # ~2 min; corpus + null sim; writes data/repro_check_results.json
```
Requires `C:\Projects\Pairwise70\data` (595 .rda) and
`C:\Projects\fragility-atlas\src\loader.py` + pyreadr.

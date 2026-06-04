# spec-collapse-atlas — specification

## 1. Problem

A multiverse / specification-curve meta-analysis enumerates `S` defensible
analytic specifications and runs them all on **one** dataset. Practitioners then
summarise the spec-curve. The two summaries in common use are:

1. **Concordance** — "% of specifications that are statistically significant."
   This is not an interval and has no calibrated coverage.
2. **IV-RE pooling across specifications** — inverse-variance pool the `S`
   specification estimates `(θ̂_s, V_s)` as if they were `S` independent studies.

(2) is the cardinal error: the `S` specifications share the same data and are
strongly correlated, so pooling them as independent shrinks the variance toward
`1/Σ(1/V_s)` — roughly an `S`-fold collapse — and is anti-conservative.

## 2. Estimand

**False-robustness rate**: the share of reviews for which the IV-RE pool's
verdict is "robust" (95% CI excludes the null) but the coverage-calibrated
weighted-likelihood interval's verdict is "fragile" (CI includes the null).

Secondary: median interval-width ratio (naive ÷ calibrated); Monte-Carlo
coverage and Type-I error of each summary under a known null.

## 3. Method (net-new)

For each specification `s`, treat its result as an approximate likelihood for the
pooled effect: a **scaled-t** density `θ̂_s + √V_s · t_{k_s−1}` (the per-spec pivot
is t-distributed at small k). The **weighted-likelihood** summary is the mixture
`Σ p_s · t(θ̂_s, V_s, k_s−1)` with uniform weights `p_s = 1/S`. Its interval is
obtained by inverting the (vectorised) mixture CDF at α/2 and 1−α/2. By the law of
total variance,

```
Var = Σ p_s V_s  +  Σ p_s (θ̂_s − θ̄)²   =   within-spec  +  between-spec
```

which is ≥ the mean single-spec variance, so the interval is never narrower than
an average specification. The t components (vs. a normal mixture) correct a mild
under-coverage at high heterogeneity: null coverage at τ²=0.05 improves from
~0.905 to ~0.94, while staying conservative (≥0.95) at τ²=0.

## 4. Specification grid (36 specs)

`τ²` estimator ∈ {DL, REML, PM} × CI method ∈ {Wald, HKSJ (q-floored, t_{k−1})}
× outlier handling ∈ {keep all, drop most extreme, drop two most extreme} ×
publication-bias ∈ {raw, trim-and-fill}.

HKSJ uses the `max(1, Q/(k−1))` floor and the `t_{k−1}` quantile (never z), per
the project statistics rules — without the floor HKSJ narrows below DL in the
`I²=0` regime. Trim-and-fill is the Duval-Tweedie L0 estimator (ported from the
R-validated fragility-atlas implementation; agreement to 1e−7); it is skipped on
subsets with k<3.

## 4a. Weighting sensitivity

The false-robustness rate is reported under four weighting schemes: uniform
(primary), REML-only and HKSJ-only (0/1 spec-subset restrictions to the
recommended estimator / better-calibrated interval), and AIC (w_i ∝
exp(−0.5·ΔAIC) from each spec's ML fit — a heuristic probe, since specs span
different data subsets). The headline is stable: **52–55%** across all four.

## 5. Corpus

Pairwise70 (`C:\Projects\Pairwise70\data`, 595 `.rda` Cochrane reviews). After
restricting to a valid primary analysis with k ≥ 3, **473 reviews** remain (122
skipped). Binary outcomes pooled as log-RR/log-OR, continuous as MD/SMD, on the
log/raw scale appropriate to each.

## 6. Portfolio recon (reused vs net-new)

Recon run 2026-06-04 via `find-related-repos.py` (manifest `generatedAt`
2026-06-02). Top hits:

- **fragility-atlas** (Tier 1 / Shipped) — multiverse over the same Pairwise70
  corpus, reports a *fragility %* (agreement with the DL reference). It
  enumerates and classifies but does **not** produce a corrected aggregate
  interval. **Reused:** its R-validated `src/loader.py` (`.rda` → `ReviewData`)
  for corpus ingest; the Pairwise70 corpus itself.
- **ma-workbench** — worked E156 demo; reused only as a format reference.
- **MultiverseMA** (#119) — browser multiverse engine; **reused** its built-in
  dataset list and specification dimensions. MultiverseMA reports concordance and
  vibration-of-effects but no corrected interval; this project ships the
  corrected aggregator back into it (see `docs/` / commit log).

**Net-new:** the weighted-likelihood aggregator, the IV-RE-collapse demonstration,
the Monte-Carlo coverage harness, and the false-robustness metric and atlas.

## 7. Validation

- **External:** engine reproduces metafor's published `dat.bcg` random-effects
  result exactly — REML τ²=0.3132, pooled log-RR=−0.7145, DL τ²=0.3088
  (canonical counts from metadat). MultiverseMA had shipped a corrupted `dat.bcg`
  (9/13 rows wrong), now fixed in both repos.
- Engine: REML also reproduced two independent ways (fixed-point vs. direct
  REML-likelihood optimisation); DL/PM by closed form / bracketed root-find.
- Trim-and-fill matches the R-validated fragility-atlas implementation to 1e−7.
- JS dashboard engine reproduces the Python engine's intervals bit-for-bit.
- 14 pytest tests.

## 8. Done in v0.2 / still deferred

Done: publication-bias dimension (trim-and-fill, 36-spec grid); non-uniform
weighting + sensitivity (52–55% across schemes); t-mixture WL calibration;
external metafor validation; canonical `dat.bcg` + `dat.egger2001` datasets.

Still deferred: live R-metafor validation harness (no Rscript locally); measure
conversion as a spec dimension; effect-magnitude (not just significance) verdict
flips; source-verification of the aspirin/omega-3/corticosteroids demo datasets
(labelled illustrative); JS port of trim-and-fill so the live tool also shows 36
specs (the corpus uses 36; the live tool shows the core 18).

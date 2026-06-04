# Spec-Collapse Atlas — Pre-Registration Template

> Purpose: freeze the specification grid, the aggregators, and the primary
> verdict rule **before** a corpus run, so a "false-robustness rate" reported
> later is a pre-registered prediction rather than a post-hoc choice. This is the
> same garden-of-forking-paths discipline the atlas itself critiques, turned on
> the atlas.
>
> Honesty note (truth-first): the **v0.2.0 shipped run** below was specified by
> reading the code that already existed, so for that run this file is a
> *specification record*, not a true pre-registration. Pre-registration
> discipline binds **the next corpus run / any grid extension** — freeze the
> fingerprint (§7), commit this file, *then* run.

---

## 1. Run identity

| Field | Value |
|---|---|
| Run label | `v0.2.0-cdsr-473` |
| Status | `[x] specification record (post-hoc)`  `[ ] pre-registered (frozen before run)` |
| Corpus | Cochrane CDSR pairwise outcomes via Pairwise70 `.rda` (fragility-atlas loader) |
| Corpus path (local, not shipped) | `$SPEC_COLLAPSE_CORPUS` (default `C:\Projects\Pairwise70\data`) |
| Confidence level | `cl = 0.95` |
| Freeze fingerprint (§7) | _fill at freeze time_ |
| Freeze commit | _git SHA at freeze time_ |
| Freeze date | _YYYY-MM-DD_ |

---

## 2. Inclusion / exclusion (decided before the run)

- **Include**: every `CD*.rda` review the loader returns with `k >= 3` studies.
- **Exclude — skipped** (not failures, not counted in any denominator): `k < 3`
  (mixture/HKSJ df undefined; trim-and-fill needs `k >= 3`). v0.2.0: **122 skipped**.
- **Exclude — errored**: loader raises, or any spec/aggregator raises. v0.2.0: **0 errored**.
- **Denominator rule**: all rates are over **analysed reviews** (`n_reviews`), never
  over files-on-disk. v0.2.0 denominator: **n = 473**.
- No review is dropped *after* seeing its verdict. Any post-run exclusion goes in §8.

---

## 3. Specification grid (frozen)

Cartesian product, enumerated by `engine.enumerate_specs`:

| Dimension | Levels | Source of truth |
|---|---|---|
| Outlier handling | `keep_all`, `drop_1_extreme` (k>3), `drop_2_extreme` (k>4) | `engine._outlier_variants` |
| τ² estimator | `DL`, `REML`, `PM` | `engine.TAU2_ESTIMATORS` |
| Publication-bias | `raw`, `trim_fill` (skipped where subset k<3) | `engine.PUB_BIAS` |
| CI method | `Wald`, `HKSJ` (q-floored at 1, t_{k-1}) | `engine.CI_METHODS` |

Maximum specifications per review = 3 × 3 × 2 × 2 = **36** (fewer when small-k
subsets suppress outlier-drop or trim-fill levels). Adding/removing a level is a
**new run** with a new fingerprint — never an in-place edit of a frozen run.

Statistical conventions are inherited verbatim from `advanced-stats.md` and are
part of the freeze: log-scale pooling; HKSJ q-floor `max(1, Q/(k-1))` with
`t_{k-1}`; REML via Viechtbauer fixed-point; PM by bracketed root-find;
DL only inside the grid, never as the sole estimator.

---

## 4. Aggregators (frozen)

Three summaries of the per-review spec curve (`spec_collapse.aggregators`):

1. **`naive_concordance`** — "% of specs significant". Verdict `robust` iff
   `>= 95%` of specs have a CI excluding 0. Reported only as the *foil* (the
   number people read off a spec curve); has no calibrated coverage.
2. **`naive_ivre_pool`** — inverse-variance pool the S spec estimates as if S
   independent studies. Verdict `robust` iff the pooled CI excludes 0. This is
   the **cardinal sin** the atlas measures (variance collapses ~×S).
3. **`weighted_likelihood`** (PRIMARY correct combination) — scaled-t mixture,
   `df_s = k_s − 1`, interval by mixture-quantile inversion; variance = mean
   within-spec + between-spec spread (law of total variance), never narrower than
   one spec. Verdict `robust` iff the mixture CI excludes 0. Wagenmakers-style.

---

## 5. Primary outcome (frozen — this is the headline number)

**False-robustness rate** = share of analysed reviews where the naive IV-RE pool
says `robust` **and** the weighted-likelihood interval says `fragile`:

```
false_robust = (ivre.verdict == "robust") and (wl.verdict == "fragile")
false_robust_pct = 100 * sum(false_robust) / n_reviews
```

- **Pre-specified prediction**: `false_robust_pct >= 40%`.
- **v0.2.0 realised**: **54.97%** (260 / 473). IV-RE robust on 417; WL robust on 157.
- Direction is one-sided by construction (WL is never *narrower* than IV-RE, so a
  review can be false-robust but not false-fragile under this pairing).

### Secondary outcomes (pre-specified, reported regardless of result)

| Outcome | Definition | v0.2.0 |
|---|---|---|
| Median width ratio | median(IV-RE width / WL width) | 0.124 (IQR 0.107–0.138) |
| Concordance false-robust | `concord robust ∧ wl fragile` | 0.0% (0) — concordance ≥95% never co-occurred with WL-fragile here |
| False-meaningful | IV-RE CI excludes ±ROPE band but WL does not; `ROPE = ±log(1.1)` (10% relative effect) | 42.49% (201) |

---

## 6. Sensitivity analyses (pre-specified — listed before the run, not added after)

Weighting schemes for the WL mixture (`aggregators.build_weights`); the primary
result must be reported under **all four**, and the headline claim stands only if
it is stable across them:

| Scheme | Meaning | v0.2.0 false-robust % |
|---|---|---|
| `uniform` (primary) | equal weights | 54.97 |
| `reml_only` | restrict to REML specs | 54.12 |
| `hksj_only` | restrict to HKSJ specs | 55.39 |
| `aic` | `w ∝ exp(-½ΔAIC)` — HEURISTIC only (specs span different data subsets; likelihoods not strictly comparable) | 52.01 |

Pre-specified stability criterion: **all four within ±5 percentage points** of the
primary. v0.2.0: range 52.0–55.4, **passes**.

---

## 7. Freeze fingerprint

Run this **before** the corpus run and paste the digest into §1. It hashes the
grid definition + criterion constants so any later drift is detectable:

```python
import hashlib, json
from spec_collapse import engine
from spec_collapse.corpus import MEANINGFUL_DELTA

grid = {
    "outlier": ["keep_all", "drop_1_extreme", "drop_2_extreme"],
    "estimator": sorted(engine.TAU2_ESTIMATORS),
    "ci_method": list(engine.CI_METHODS),
    "pub_bias": list(engine.PUB_BIAS),
    "cl": 0.95,
    "concordance_threshold": 0.95,
    "rope_delta": MEANINGFUL_DELTA,        # log(1.1)
    "primary_rule": "ivre==robust and wl==fragile",
    "wl_components": "t",
    "min_k": 3,
}
print(hashlib.sha256(json.dumps(grid, sort_keys=True).encode()).hexdigest())
```

If this digest changes, you are running a **different study** — bump the run
label and re-register. Do not reuse a frozen run's headline number under a
changed fingerprint.

---

## 8. Deviation log (append-only; empty at freeze)

Every departure from §§2–6 after freezing goes here with date + reason. An empty
log on a `pre-registered` run is the claim "we ran exactly what we froze."

- _none_

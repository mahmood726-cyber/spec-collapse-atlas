# Spec-Collapse Atlas

**Correct inference for multiverse / many-analyst meta-analysis.**

A multiverse meta-analysis runs many analytic specifications on **one** dataset.
The conventional way to summarise the resulting spec-curve — inverse-variance
(IV-RE) pooling of the specification estimates — treats those correlated
re-analyses as if they were independent studies. This collapses the confidence
interval by a factor of roughly the number of specifications and is badly
**anti-conservative**: it manufactures precision and robustness the data do not
support.

This project supplies the missing piece: a **coverage-calibrated
weighted-likelihood interval** that accounts for both within- and
between-specification uncertainty, and is never narrower than a single
specification.

## Headline

Across **473 Cochrane intervention meta-analyses** (Pairwise70 corpus, k ≥ 3):

Each review is re-analysed across **36 specifications** (3 τ² estimators × 2 CI
methods × 3 outlier rules × {raw, trim-and-fill}):

| | Naive IV-RE pool | Weighted-likelihood |
|---|---|---|
| Conclusions called "robust" | **88%** (417/473) | **33%** (157/473) |
| Median interval-width ratio | — | naive is **0.124×** the width |
| **False robustness** (IV-RE robust → corrected fragile) | **260/473 = 55.0%** | |

The 55% reversal is stable — **52–55%** across uniform, REML-only, HKSJ-only, and
AIC weighting schemes. Monte-Carlo under the null (μ = 0): the naive pool's
**Type-I error is 70–81%** (nominal 5%) and its coverage 0.19–0.30 (nominal 0.95).
The weighted-likelihood interval holds coverage at **0.94–0.99**.

## The method

Each specification *s* contributes an approximate likelihood for the pooled
effect θ: a scaled-t `θ̂_s + √V_s · t_{k−1}`. The weighted-likelihood summary is
the **mixture** `Σ pₛ · t-density(θ̂_s, V_s, k−1)` (uniform weights `pₛ = 1/S`). By
the law of total variance its variance is

```
Var = mean within-spec variance  +  between-spec variance
```

so it can never be narrower than the average single specification — the opposite
of the IV-RE pool, whose variance shrinks toward `1/Σ(1/Vₛ)`. (The t components
fix a mild under-coverage the normal mixture showed at high heterogeneity.)

## Live dashboard

<https://mahmood726-cyber.github.io/spec-collapse-atlas/> — offline single-file
HTML, no external dependencies. Pick a built-in dataset to see the naive vs.
calibrated interval and the verdict flip live; browse the 473-review corpus
width-ratio distribution.

## Run

```bash
python -m pytest -q                  # 14 tests (engine + aggregators + trim-and-fill)
python demo.py                       # 5 real datasets + Monte-Carlo coverage
python run_corpus_full.py            # full 473-review Pairwise70 atlas + sensitivity
```

`run_corpus_full.py` reuses the
[fragility-atlas](https://github.com/mahmood726-cyber/fragility-atlas) loader
(`pyreadr` + the Pairwise70 `.rda` corpus) for ingest; everything downstream is
this project.

## Layout

| Path | Purpose |
|---|---|
| `spec_collapse/engine.py` | τ² estimators (DL/REML/PM), RE pooling, Wald/HKSJ CIs, trim-and-fill, 36-spec grid |
| `spec_collapse/aggregators.py` | naive concordance, naive IV-RE pool, weighted-likelihood |
| `spec_collapse/coverage.py` | Monte-Carlo coverage harness |
| `spec_collapse/corpus.py` | CDSR-scale runner over Pairwise70 |
| `index.html` | offline dashboard (main artifact) |
| `tests/` | validation suite (REML dual-method verified) |
| `e156-submission/` | E156 micro-paper bundle |

## License

MIT — see `LICENSE`.

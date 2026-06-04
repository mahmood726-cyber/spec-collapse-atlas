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

| | Naive IV-RE pool | Weighted-likelihood |
|---|---|---|
| Conclusions called "robust" | **83%** (394/473) | **42%** (201/473) |
| Median interval-width ratio | — | naive is **0.22×** the width |
| **False robustness** (IV-RE robust → corrected fragile) | **193/473 = 40.8%** | |

Monte-Carlo under the null (μ = 0): the naive pool's **Type-I error is 60–74%**
(nominal 5%) and its coverage 0.27–0.39 (nominal 0.95). The weighted-likelihood
interval restores coverage to **0.91–0.97**.

## The method

Each specification *s* contributes an approximate likelihood for the pooled
effect θ: `N(θ̂_s, V_s)`. The weighted-likelihood summary is the **Gaussian
mixture** `Σ pₛ · N(θ̂_s, V_s)` (uniform weights `pₛ = 1/S`). By the law of total
variance its variance is

```
Var = mean within-spec variance  +  between-spec variance
```

so it can never be narrower than the average single specification — the opposite
of the IV-RE pool, whose variance shrinks toward `1/Σ(1/Vₛ)`.

## Live dashboard

<https://mahmood726-cyber.github.io/spec-collapse-atlas/> — offline single-file
HTML, no external dependencies. Pick a built-in dataset to see the naive vs.
calibrated interval and the verdict flip live; browse the 473-review corpus
width-ratio distribution.

## Run

```bash
python -m pytest -q                  # 10 tests (engine + aggregators)
python demo.py                       # 5 real datasets + Monte-Carlo coverage
python run_corpus_full.py            # full 473-review Pairwise70 atlas
```

`run_corpus_full.py` reuses the
[fragility-atlas](https://github.com/mahmood726-cyber/fragility-atlas) loader
(`pyreadr` + the Pairwise70 `.rda` corpus) for ingest; everything downstream is
this project.

## Layout

| Path | Purpose |
|---|---|
| `spec_collapse/engine.py` | τ² estimators (DL/REML/PM), RE pooling, Wald/HKSJ CIs, 18-spec grid |
| `spec_collapse/aggregators.py` | naive concordance, naive IV-RE pool, weighted-likelihood |
| `spec_collapse/coverage.py` | Monte-Carlo coverage harness |
| `spec_collapse/corpus.py` | CDSR-scale runner over Pairwise70 |
| `index.html` | offline dashboard (main artifact) |
| `tests/` | validation suite (REML dual-method verified) |
| `e156-submission/` | E156 micro-paper bundle |

## License

MIT — see `LICENSE`.

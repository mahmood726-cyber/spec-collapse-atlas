# Spec-Collapse Atlas

**Type:** methods · **Primary estimand:** false-robustness rate · **Data:** Pairwise70 (473 Cochrane meta-analyses, k ≥ 3)

## Body (156-word contract)

When a multiverse meta-analysis is summarised by pooling its specifications, does that summary overstate how robust the underlying conclusion actually is? We analysed 473 Cochrane intervention meta-analyses, each with at least three studies, drawn from the Pairwise70 corpus of binary and continuous outcomes. Every review was re-analysed across thirty-six specifications spanning estimator, interval, outlier, and publication-bias choices, then summarised by naive inverse-variance pooling versus a calibrated weighted-likelihood interval. The naive pool called 88 percent of conclusions robust versus 33 percent for the calibrated interval, reversing 260 reviews (55 percent) from robust to fragile and compressing intervals to one-eighth their calibrated width. This reversal stayed between 52 and 55 percent across four weighting schemes, while under simulated null effects the naive pool's false-positive rate reached 72 to 81 percent. Inverse-variance pooling across specifications treats re-analyses of one dataset as independent evidence, manufacturing precision the data do not support. This analysis covers binary significance agreement, not effect-magnitude shifts.

## Sentence map

- **S1 — Question:** does pooling specifications overstate robustness?
- **S2 — Dataset:** 473 Cochrane MAs (k ≥ 3), Pairwise70, binary + continuous.
- **S3 — Method:** 36 specs (estimator × interval × outlier × publication-bias), naive IV-RE pool vs weighted-likelihood interval.
- **S4 — Result:** 88% vs 33% robust; 260/473 (55%) flip; intervals 1/8 width.
- **S5 — Robustness:** 52–55% across 4 weighting schemes; null Type-I 72–81% (naive) vs 0.92–0.99 coverage (calibrated).
- **S6 — Interpretation:** IV-RE pooling treats one dataset's re-analyses as independent → false precision.
- **S7 — Boundary:** significance-agreement only.

## Reproduce

`python run_corpus_full.py` (reuses the fragility-atlas loader for Pairwise70 ingest).
Dashboard: <https://mahmood726-cyber.github.io/spec-collapse-atlas/>

# Spec-Collapse Atlas

**Type:** methods · **Primary estimand:** false-robustness rate · **Data:** Pairwise70 (473 Cochrane meta-analyses, k ≥ 3)

## Body (156-word contract)

When a multiverse meta-analysis is summarised by pooling its specifications, does that summary overstate how robust the underlying conclusion actually is? We analysed 473 Cochrane intervention meta-analyses, each with at least three studies, drawn from the Pairwise70 corpus of binary and continuous outcomes. Every review was re-analysed across eighteen specifications and then summarised two ways: naive inverse-variance pooling versus a coverage-calibrated weighted-likelihood interval. The naive pool called 83 percent of conclusions robust versus 42 percent for the calibrated interval, reversing 193 reviews (40.8 percent) from robust to fragile and compressing intervals to one-fifth their calibrated width. Under simulated null effects, the naive pool's false-positive rate reached 60 to 74 percent, whereas the weighted-likelihood interval held near-nominal 95 percent coverage. Inverse-variance pooling across specifications treats re-analyses of one dataset as independent evidence, manufacturing precision the data do not support. This analysis covers binary significance agreement and one weighting scheme, not effect-magnitude shifts or alternative specification weights.

## Sentence map

- **S1 — Question:** does pooling specifications overstate robustness?
- **S2 — Dataset:** 473 Cochrane MAs (k ≥ 3), Pairwise70, binary + continuous.
- **S3 — Method:** 18 specs each, summarised by naive IV-RE pool vs weighted-likelihood interval.
- **S4 — Result:** 83% vs 42% robust; 193/473 (40.8%) flip; intervals 1/5 width.
- **S5 — Robustness:** null Type-I 60–74% (naive) vs ~95% coverage (calibrated).
- **S6 — Interpretation:** IV-RE pooling treats one dataset's re-analyses as independent → false precision.
- **S7 — Boundary:** significance-agreement + one weighting scheme only.

## Reproduce

`python run_corpus_full.py` (reuses the fragility-atlas loader for Pairwise70 ingest).
Dashboard: <https://mahmood726-cyber.github.io/spec-collapse-atlas/>

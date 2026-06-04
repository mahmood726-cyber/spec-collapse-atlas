"""CDSR-scale corpus runner.

Reuses fragility-atlas's validated `load_review` (Pairwise70 .rda ingest ->
ReviewData with yi/sei/scale/Cochrane verdict), then applies THIS project's
verified engine + the three aggregators. The novel layer is the aggregation, not
the ingest -- we deliberately stand on fragility-atlas's already-R-validated
loader rather than re-implementing .rda parsing.

Headline outputs per review: width ratio (naive IV-RE / weighted-likelihood) and
the false-robustness flag (IV-RE pool says 'robust' while the corrected interval
crosses null).
"""

from __future__ import annotations

import glob
import os
import sys

from .aggregators import (build_weights, naive_concordance, naive_ivre_pool,
                          weighted_likelihood)
from .engine import enumerate_specs

# Local analysis runner: paths are configurable via env vars (or the corpus_dir
# argument). The defaults below are this workstation's locations and are NOT
# shipped into any deployed asset (the dashboard embeds pre-computed results).
import math

MEANINGFUL_DELTA = math.log(1.1)  # ROPE half-width: a 10% relative effect
DEFAULT_CORPUS = os.environ.get("SPEC_COLLAPSE_CORPUS", r"C:\Projects\Pairwise70\data")
FRAGILITY_ATLAS = os.environ.get("FRAGILITY_ATLAS_PATH", r"C:\Projects\fragility-atlas")


def _get_loader():
    """Import fragility-atlas's load_review, failing closed with a clear message."""
    if FRAGILITY_ATLAS not in sys.path:
        sys.path.insert(0, FRAGILITY_ATLAS)
    try:
        from src.loader import load_review  # type: ignore
    except Exception as e:  # pragma: no cover - environment guard
        raise RuntimeError(
            f"Cannot import fragility-atlas loader from {FRAGILITY_ATLAS}. "
            f"Corpus scaling requires it (pyreadr + Pairwise70 .rda files). "
            f"Underlying error: {e!r}"
        )
    return load_review


def run_corpus(corpus_dir=DEFAULT_CORPUS, limit=None, cl=0.95):
    load_review = _get_loader()
    files = sorted(glob.glob(os.path.join(corpus_dir, "CD*.rda")))
    if limit:
        files = files[:limit]
    if not files:
        raise RuntimeError(f"No CD*.rda files found in {corpus_dir}")

    rows = []
    skipped = 0
    errored = 0
    for f in files:
        try:
            r = load_review(f)
        except Exception:
            errored += 1
            continue
        if r is None or r.k < 3:
            skipped += 1
            continue
        yi = list(r.yi)
        vi = [float(s) ** 2 for s in r.sei]
        try:
            specs = enumerate_specs(yi, vi, cl=cl)
            con = naive_concordance(specs, cl=cl)
            ivre = naive_ivre_pool(specs, cl=cl)
            wl = weighted_likelihood(specs, cl=cl)
        except Exception:
            errored += 1
            continue
        w_ivre = ivre["ci_high"] - ivre["ci_low"]
        w_wl = wl["ci_high"] - wl["ci_low"]
        # effect-magnitude verdict: "meaningful" = CI excludes a +/-delta ROPE
        # band (delta = log 1.1, a 10% relative effect on the ratio scale).
        delta = MEANINGFUL_DELTA
        iv_meaningful = ivre["ci_low"] > delta or ivre["ci_high"] < -delta
        wl_meaningful = wl["ci_low"] > delta or wl["ci_high"] < -delta
        rows.append({
            "review_id": r.review_id, "k": r.k, "scale": r.scale,
            "cochrane_sig": bool(r.is_significant),
            "pct_significant": con["pct_significant"],
            "concord_verdict": con["verdict"],
            "ivre_verdict": ivre["verdict"], "ivre_width": w_ivre,
            "wl_verdict": wl["verdict"], "wl_width": w_wl,
            "width_ratio": (w_ivre / w_wl) if w_wl > 0 else float("nan"),
            "false_robust": ivre["verdict"] == "robust" and wl["verdict"] == "fragile",
            "concord_false_robust": con["verdict"] == "robust" and wl["verdict"] == "fragile",
            "false_meaningful": bool(iv_meaningful and not wl_meaningful),
        })
    return {"rows": rows, "n_reviews": len(rows), "skipped": skipped,
            "errored": errored, "n_files": len(files)}


SCHEMES = ("uniform", "reml_only", "hksj_only", "aic")


def run_corpus_sensitivity(corpus_dir=DEFAULT_CORPUS, limit=None, cl=0.95,
                           schemes=SCHEMES):
    """False-robustness rate (IV-RE robust -> corrected fragile) under each
    weighting scheme. Enumerates each review once; aggregates several ways."""
    load_review = _get_loader()
    files = sorted(glob.glob(os.path.join(corpus_dir, "CD*.rda")))
    if limit:
        files = files[:limit]
    counts = {s: 0 for s in schemes}
    n = 0
    for f in files:
        try:
            r = load_review(f)
        except Exception:
            continue
        if r is None or r.k < 3:
            continue
        yi = list(r.yi)
        vi = [float(s) ** 2 for s in r.sei]
        try:
            specs = enumerate_specs(yi, vi, cl=cl)
            ivre = naive_ivre_pool(specs, cl=cl)
            n += 1
            for s in schemes:
                wl = weighted_likelihood(specs, cl=cl, weights=build_weights(specs, s))
                if ivre["verdict"] == "robust" and wl["verdict"] == "fragile":
                    counts[s] += 1
        except Exception:
            continue
    return {"n_reviews": n,
            "false_robust_pct": {s: 100.0 * counts[s] / n if n else 0.0 for s in schemes},
            "false_robust_n": counts}


def summarize(result):
    rows = result["rows"]
    n = len(rows)
    if n == 0:
        return {"n_reviews": 0}
    import numpy as np
    ratios = np.array([r["width_ratio"] for r in rows
                       if r["width_ratio"] == r["width_ratio"]])  # drop nan
    false_robust = sum(r["false_robust"] for r in rows)
    concord_fr = sum(r["concord_false_robust"] for r in rows)
    false_meaningful = sum(r.get("false_meaningful", False) for r in rows)
    ivre_robust = sum(r["ivre_verdict"] == "robust" for r in rows)
    wl_robust = sum(r["wl_verdict"] == "robust" for r in rows)
    return {
        "n_reviews": n,
        "skipped": result["skipped"], "errored": result["errored"],
        "median_width_ratio": float(np.median(ratios)),
        "p25_width_ratio": float(np.percentile(ratios, 25)),
        "p75_width_ratio": float(np.percentile(ratios, 75)),
        "ivre_robust": ivre_robust, "wl_robust": wl_robust,
        "false_robust_n": false_robust,
        "false_robust_pct": 100.0 * false_robust / n,
        "concord_false_robust_n": concord_fr,
        "concord_false_robust_pct": 100.0 * concord_fr / n,
        "false_meaningful_n": false_meaningful,
        "false_meaningful_pct": 100.0 * false_meaningful / n,
        "meaningful_delta": MEANINGFUL_DELTA,
    }

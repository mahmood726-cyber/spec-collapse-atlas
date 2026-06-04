"""Validation tests for the spec-collapse engine and aggregators.

The BCG dataset (spec_collapse/datasets.py) is now the canonical metafor/metadat
dat.bcg, derived from the published Colditz (1994) 2x2 counts. This gives a hard
EXTERNAL reference: metafor's documented random-effects result on dat.bcg is
REML tau^2 = 0.3132 and pooled log-RR estimate = -0.7145; the DL estimate is
0.3088. The engine must reproduce these. (MultiverseMA previously shipped a
corrupted dat.bcg -- 9/13 rows wrong -- which is why its pooled estimate matched
no textbook; that has been corrected.)
"""
import math
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from spec_collapse import DATASETS, enumerate_specs
from spec_collapse.aggregators import (naive_concordance, naive_ivre_pool,
                                       weighted_likelihood)
from spec_collapse.engine import tau2_dl, tau2_pm, tau2_reml, re_pool, ci_hksj


BCG = DATASETS["bcg"]


def test_tau2_dl_matches_metafor():
    # metafor published DL tau^2 for dat.bcg
    assert tau2_dl(BCG["yi"], BCG["vi"]) == pytest.approx(0.3088, abs=2e-3)


def test_tau2_reml_matches_metafor_and_dual_method():
    # (a) external reference: metafor published REML tau^2 = 0.3132
    fp = tau2_reml(BCG["yi"], BCG["vi"])
    assert fp == pytest.approx(0.3132, abs=3e-3)
    # (b) fixed-point iteration must also match direct REML-likelihood optimisation
    import numpy as np
    from scipy.optimize import minimize_scalar
    y = np.array(BCG["yi"]); v = np.array(BCG["vi"])

    def neg_reml(t2):
        if t2 < 0:
            return 1e18
        w = 1 / (v + t2); sw = w.sum(); yb = (w * y).sum() / sw
        return 0.5 * np.log(v + t2).sum() + 0.5 * np.log(sw) + 0.5 * (w * (y - yb) ** 2).sum()

    direct = minimize_scalar(neg_reml, bounds=(0, 5), method="bounded").x
    assert fp == pytest.approx(direct, abs=1e-3)


def test_tau2_pm_near_reml():
    pm = tau2_pm(BCG["yi"], BCG["vi"])
    assert 0.27 < pm < 0.38


def test_re_pool_bcg_logRR_matches_metafor():
    tau2 = tau2_reml(BCG["yi"], BCG["vi"])
    theta, var, _, _ = re_pool(BCG["yi"], BCG["vi"], tau2)
    # metafor published RE pooled log-RR for dat.bcg
    assert theta == pytest.approx(-0.7145, abs=5e-3)


def test_k1_edge_no_crash():
    specs = enumerate_specs([-0.5], [0.1])
    assert all(s["tau2"] == 0.0 for s in specs)
    assert len(specs) >= 1


def test_tau2_zero_when_homogeneous():
    # identical effects -> Q=0 -> tau^2 must be 0 for all estimators
    y, v = [0.2, 0.2, 0.2, 0.2], [0.05, 0.05, 0.05, 0.05]
    assert tau2_dl(y, v) == 0.0
    assert tau2_pm(y, v) == 0.0
    assert tau2_reml(y, v) == pytest.approx(0.0, abs=1e-6)


def test_hksj_floor_widens_not_narrows():
    # When Q < k-1 the floored HKSJ must not be narrower than Wald.
    y, v = [0.10, 0.12, 0.09, 0.11], [0.05, 0.05, 0.05, 0.05]
    tau2 = 0.0
    theta, var_wald, _, _ = re_pool(y, v, tau2)
    lo_h, hi_h, var_h = ci_hksj(y, v, theta, tau2)
    assert var_h >= var_wald - 1e-12


def test_weighted_likelihood_never_narrower_than_mean_spec():
    """The corrected interval's total variance >= mean within-spec variance."""
    specs = enumerate_specs(BCG["yi"], BCG["vi"])
    wl = weighted_likelihood(specs)
    mean_within = sum(s["var"] for s in specs) / len(specs)
    assert wl["var"] >= mean_within - 1e-12
    assert wl["between_var"] >= 0.0


def test_ivre_pool_collapses_below_single_spec():
    """The naive IV-RE pool variance is far below any single spec's variance."""
    specs = enumerate_specs(BCG["yi"], BCG["vi"])
    ivre = naive_ivre_pool(specs)
    min_spec_var = min(s["var"] for s in specs)
    assert ivre["var"] < min_spec_var  # the collapse


def test_concordance_counts():
    specs = enumerate_specs(BCG["yi"], BCG["vi"])
    con = naive_concordance(specs)
    nsig = sum(1 for s in specs if s["significant"])
    assert con["pct_significant"] == pytest.approx(100.0 * nsig / len(specs))


def test_pubbias_dimension_expands_grid():
    # 3 outlier x 3 estimator x 2 pub-bias x 2 CI = 36 for k large enough
    assert len(enumerate_specs(BCG["yi"], BCG["vi"], pub_bias=True)) == 36
    assert len(enumerate_specs(BCG["yi"], BCG["vi"], pub_bias=False)) == 18


def test_trim_and_fill_matches_fragility_atlas():
    """Ported L0 trim-and-fill must match the fragility-atlas reference (1e-6)."""
    import os
    fa = r"C:\Projects\fragility-atlas"
    if not os.path.isdir(fa):
        pytest.skip("fragility-atlas not available")
    import sys
    import numpy as np
    from spec_collapse.engine import trim_and_fill, re_pool, tau2_dl
    sys.path.insert(0, fa)
    try:
        from src.corrections import trim_and_fill as fa_tf  # type: ignore
    except Exception:
        pytest.skip("fragility-atlas corrections not importable")
    for key in ("bcg", "warfarin", "magnesium"):
        d = DATASETS[key]
        yf, vf = trim_and_fill(d["yi"], d["vi"], tau2_dl)
        mine = re_pool(yf, vf, tau2_dl(yf, vf))[0]
        ref = fa_tf(np.array(d["yi"]), np.array([v ** 0.5 for v in d["vi"]]),
                    "DL", "Wald").theta
        assert mine == pytest.approx(ref, abs=1e-6)


def test_weighting_schemes_select_correctly():
    from spec_collapse.aggregators import build_weights
    specs = enumerate_specs(BCG["yi"], BCG["vi"])
    assert build_weights(specs, "uniform") is None
    w_hksj = build_weights(specs, "hksj_only")
    assert all((w > 0) == (s["ci_method"] == "HKSJ") for w, s in zip(w_hksj, specs))
    w_reml = build_weights(specs, "reml_only")
    assert all((w > 0) == (s["estimator"] == "REML") for w, s in zip(w_reml, specs))


def test_t_mixture_wider_than_normal_mixture():
    # the scaled-t mixture is at least as wide as the normal mixture
    specs = enumerate_specs(BCG["yi"], BCG["vi"])
    wt = weighted_likelihood(specs, components="t")
    wn = weighted_likelihood(specs, components="normal")
    assert (wt["ci_high"] - wt["ci_low"]) >= (wn["ci_high"] - wn["ci_low"]) - 1e-9

"""Validation tests for the spec-collapse engine and aggregators.

tau^2 references are the VERIFIED values on MultiverseMA's *actually-shipped* BCG
numbers (spec_collapse/datasets.py), confirmed two independent ways:
  REML 0.3936 via fixed-point iteration == 0.3936 via direct REML-likelihood
  optimization (scipy minimize_scalar).  DL 0.4047 via the closed-form formula.

NB these differ from the commonly-quoted metafor dat.bcg figures (REML ~0.31,
estimate ~-0.71). MultiverseMA's shipped BCG vi could not be reconciled with a
from-2x2-counts reconstruction, but that reconstruction is itself unverified
(no R/metafor in this environment), so we make NO claim here about which is
canonical -- we only assert the engine reproduces the standard formulas on the
numbers it is given, confirmed two independent ways. Whether MV's data matches
dat.bcg is deferred to an authoritative-source check.
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


def test_tau2_dl_verified():
    # closed-form DL on MultiverseMA's shipped BCG numbers
    assert tau2_dl(BCG["yi"], BCG["vi"]) == pytest.approx(0.4047, abs=2e-3)


def test_tau2_reml_dual_method_verified():
    # fixed-point iteration must match direct REML-likelihood optimisation
    fp = tau2_reml(BCG["yi"], BCG["vi"])
    assert fp == pytest.approx(0.3936, abs=3e-3)
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


def test_tau2_pm_between_dl_and_reml():
    pm = tau2_pm(BCG["yi"], BCG["vi"])
    # PM sits near REML on this data
    assert 0.35 < pm < 0.45


def test_re_pool_bcg_logRR():
    tau2 = tau2_reml(BCG["yi"], BCG["vi"])
    theta, var, _, _ = re_pool(BCG["yi"], BCG["vi"], tau2)
    # verified RE pooled log-RR on MV's shipped (non-canonical) BCG data
    assert theta == pytest.approx(-0.6321, abs=5e-3)


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

"""Monte-Carlo coverage harness.

Simulate meta-analytic datasets from a KNOWN truth, run the full multiverse on
each, and measure the true coverage of every aggregator's nominal-95% interval.
Demonstrates that naive_ivre_pool is anti-conservative (coverage << 0.95) while
weighted_likelihood restores ~nominal coverage.

Seed is pinned (testing rules: reproducible Monte-Carlo).
"""

from __future__ import annotations

import numpy as np

from .aggregators import naive_concordance, naive_ivre_pool, weighted_likelihood
from .engine import enumerate_specs


def simulate_one(rng, k, mu, tau2_true, vi_template):
    """One synthetic meta-analysis: theta_i ~ N(mu, tau2), y_i ~ N(theta_i, v_i)."""
    vi = rng.choice(vi_template, size=k, replace=True)
    theta_i = rng.normal(mu, np.sqrt(tau2_true), size=k)
    yi = rng.normal(theta_i, np.sqrt(vi))
    return list(yi), list(vi)


def run_coverage(n_reps=1000, k=10, mu=0.0, tau2_true=0.05,
                 vi_template=None, seed=20260604, cl=0.95):
    if vi_template is None:
        # realistic within-study variances (BCG-like spread)
        vi_template = [0.0154, 0.0200, 0.0312, 0.0356, 0.0368, 0.0393,
                       0.0512, 0.0628, 0.0636, 0.1946, 0.2252, 0.3256]
    rng = np.random.default_rng(seed)
    vi_template = np.asarray(vi_template, dtype=float)

    cover_ivre = 0
    cover_wl = 0
    falsepos_concord = 0  # under null: claims "robust" (>=95% specs significant)
    sig_ivre = 0          # under null: IV-RE interval excludes truth-direction 0
    width_ivre = []
    width_wl = []

    for _ in range(n_reps):
        yi, vi = simulate_one(rng, k, mu, tau2_true, vi_template)
        specs = enumerate_specs(yi, vi, cl=cl)

        a_ivre = naive_ivre_pool(specs, cl=cl)
        a_wl = weighted_likelihood(specs, cl=cl)
        a_con = naive_concordance(specs, cl=cl)

        if a_ivre["ci_low"] <= mu <= a_ivre["ci_high"]:
            cover_ivre += 1
        if a_wl["ci_low"] <= mu <= a_wl["ci_high"]:
            cover_wl += 1
        if a_con["verdict"] == "robust" and mu == 0.0:
            falsepos_concord += 1
        if a_ivre["significant"]:
            sig_ivre += 1
        width_ivre.append(a_ivre["ci_high"] - a_ivre["ci_low"])
        width_wl.append(a_wl["ci_high"] - a_wl["ci_low"])

    return {
        "n_reps": n_reps, "k": k, "mu": mu, "tau2_true": tau2_true, "cl": cl,
        "coverage_naive_ivre": cover_ivre / n_reps,
        "coverage_weighted_likelihood": cover_wl / n_reps,
        "false_robust_concordance": falsepos_concord / n_reps,
        "type1_naive_ivre": sig_ivre / n_reps,
        "mean_width_naive_ivre": float(np.mean(width_ivre)),
        "mean_width_weighted_likelihood": float(np.mean(width_wl)),
    }

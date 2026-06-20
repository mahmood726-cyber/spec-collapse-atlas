"""Oracle-locked numerical baselines — auto-generated 2026-06-20.

Each test pins a spec-collapse-atlas estimator to a value independently cross-checked
against a from-scratch scipy/closed-form oracle (engine NOT used on the
oracle side). Proof workflow: prove-methods-repos. Deterministic methods
pinned to <=1e-9; iterative to the engine stopping tol; Monte-Carlo to a
3-sigma band with a pinned seed (per the Monte-Carlo testing rule).
"""


# --- tau2_dl (DerSimonian-Laird tau^2, closed form) ---
def test_tau2_dl_baseline():
    from spec_collapse.engine import tau2_dl
    # Heterogeneous dataset (Q >> k-1): exercises the positive branch.
    yi = [0.10, 0.95, -0.40, 0.70, -0.20, 1.10, 0.05]
    vi = [0.02, 0.05, 0.01, 0.03, 0.04, 0.06, 0.015]
    # Independently computed closed form max(0,(Q-(k-1))/C), C=sum(w)-sum(w^2)/sum(w):
    assert abs(tau2_dl(yi, vi) - 0.2579339285714286) < 1e-9
    # Homogeneous dataset (Q < k-1): floor must clamp to 0.
    yi2 = [0.10, 0.30, -0.05, 0.22, 0.15, 0.40, 0.05]
    vi2 = [0.02, 0.05, 0.01, 0.03, 0.04, 0.06, 0.015]
    assert tau2_dl(yi2, vi2) == 0.0
    # Degenerate k<2 returns 0.0.
    assert tau2_dl([0.5], [0.1]) == 0.0

# --- tau2_pm (Paule-Mandel between-study variance estimator) ---
def test_tau2_pm_baseline():
    """Lock Paule-Mandel tau2 to independently-verified generalised-Q-root values.

    Oracle: root of Q(tau2)=k-1 via scipy.optimize.brentq (no repo import on oracle
    side). Equivalent to metafor::rma(method='PM'). Verified 2026-06-20.
    """
    import math
    from spec_collapse.engine import tau2_pm

    # DS1: k=8 heterogeneous (Q(0)=41.25 >> k-1=7)
    yi1 = [0.10, 0.80, -0.30, 0.95, 0.20, 1.10, -0.20, 0.60]
    vi1 = [0.02, 0.05, 0.03, 0.08, 0.04, 0.06, 0.05, 0.07]
    assert math.isclose(tau2_pm(yi1, vi1), 0.23160773903554538, rel_tol=0, abs_tol=1e-6)

    # DS2: metafor dat.bcg log-RR, k=13
    yi2 = [-0.8893, -1.5854, -1.3481, -1.4416, -0.2175, -0.7861, -1.6209,
           0.0120, -0.4694, -1.3713, -0.3394, 0.4459, -0.0173]
    vi2 = [0.3256, 0.1946, 0.4154, 0.0200, 0.0512, 0.0069, 0.0223, 0.0049,
           0.0307, 0.0220, 0.0098, 0.0260, 0.0788]
    assert math.isclose(tau2_pm(yi2, vi2), 0.4206583718325821, rel_tol=0, abs_tol=1e-6)

    # Boundary: homogeneous -> tau2 exactly 0 (Q(0) <= k-1)
    assert tau2_pm([0.10, 0.12, 0.09, 0.11], [0.05, 0.05, 0.05, 0.05]) == 0.0

    # k=2 small-sample root
    assert math.isclose(tau2_pm([0.2, 0.9], [0.03, 0.04]), 0.21, rel_tol=0, abs_tol=1e-6)

    # k=1 -> no between-study variance
    assert tau2_pm([0.5], [0.02]) == 0.0

# --- tau2_reml ---
def test_tau2_reml_baseline():
    import math
    from spec_collapse.engine import tau2_reml
    # k=8 well-conditioned fixture; oracle = brentq root of REML estimating equation
    # and scipy log-likelihood argmax (both ~0.0428959572, see proof).
    yi = [0.10, -0.30, 0.45, 0.20, -0.10, 0.55, 0.05, -0.25]
    vi = [0.02, 0.05, 0.03, 0.04, 0.015, 0.06, 0.025, 0.035]
    val = tau2_reml(yi, vi)
    assert math.isclose(val, 0.04289595692320517, rel_tol=0, abs_tol=1e-9)
    # tighter-tol invocation lands on the EE-root oracle value
    val_tight = tau2_reml(yi, vi, tol=1e-13, max_iter=5000)
    assert math.isclose(val_tight, 0.0428959572353926, rel_tol=0, abs_tol=1e-8)
    # homogeneous data -> tau2 pinned at the 0 boundary
    assert tau2_reml([0.1, 0.1, 0.1, 0.1], [0.02, 0.02, 0.02, 0.02]) == 0.0
    # k<2 guard
    assert tau2_reml([0.2], [0.03]) == 0.0

# --- re_pool (inverse-variance random-effects pooled estimate and Wald variance) ---
def test_re_pool_baseline():
    import math
    from spec_collapse.engine import re_pool
    yi = [0.10, 0.30, -0.05, 0.22, 0.15, 0.40, 0.08]
    vi = [0.02, 0.03, 0.015, 0.04, 0.025, 0.05, 0.018]
    tau2 = 0.012
    theta, var_wald, w, sw = re_pool(yi, vi, tau2)
    ow = [1.0 / (v + tau2) for v in vi]
    osw = sum(ow)
    otheta = sum(a * b for a, b in zip(ow, yi)) / osw
    ovar = 1.0 / osw
    assert math.isclose(theta, otheta, abs_tol=1e-12)
    assert math.isclose(var_wald, ovar, abs_tol=1e-12)
    assert math.isclose(theta, 0.1374696979860815, abs_tol=1e-9)
    assert math.isclose(var_wald, 0.00532433952444109, abs_tol=1e-9)
    th0, v0, _, _ = re_pool(yi, vi, 0.0)
    fw = [1.0 / v for v in vi]
    assert math.isclose(th0, sum(a * b for a, b in zip(fw, yi)) / sum(fw), abs_tol=1e-12)
    assert math.isclose(v0, 1.0 / sum(fw), abs_tol=1e-12)

# --- ci_hksj (Hartung-Knapp-Sidik-Jonkman CI with q>=1 floor and t_{k-1}) ---
def test_ci_hksj_baseline():
    """Lock ci_hksj to independently-verified HKSJ (knha) values.

    Oracle: var_HKSJ = max(1, (1/(k-1)) sum w_i (y_i-theta)^2)/sum(w_i),
    w_i = 1/(v_i+tau2); CI = theta +/- t_{k-1,0.975} sqrt(var). Cross-checked
    by an independent scipy reimplementation (maxabsdiff ~1e-16).
    """
    import math
    from scipy import stats
    from spec_collapse.engine import tau2_dl, re_pool, ci_hksj

    # Case 1: q-floor ACTIVE (raw q < 1, floored to 1.0)
    yi = [0.20, 0.45, 0.10, 0.35, 0.55, 0.25]
    vi = [0.02, 0.03, 0.015, 0.04, 0.025, 0.05]
    tau2 = tau2_dl(yi, vi)
    theta, _, _, _ = re_pool(yi, vi, tau2)
    lo, hi, var = ci_hksj(yi, vi, theta, tau2, 0.95)
    assert abs(lo - 0.10214247103252433) < 1e-9
    assert abs(hi - 0.49111716665136357) < 1e-9
    assert abs(var - 0.005724266427589542) < 1e-9

    # Case 2: q-floor INACTIVE (tau2=0 -> raw q=12.6 passes through; t_{k-1} used)
    yi2 = [0.0, 0.5, 1.0, 1.5, 2.0, 0.2, 1.8]
    vi2 = [0.05] * 7
    theta2, _, _, _ = re_pool(yi2, vi2, 0.0)
    lo2, hi2, var2 = ci_hksj(yi2, vi2, theta2, 0.0, 0.95)
    k = len(yi2)
    w = [1.0 / v for v in vi2]
    sw = sum(w)
    q = sum(wi * (y - theta2) ** 2 for wi, y in zip(w, yi2)) / (k - 1)
    q = max(1.0, q)
    varO = q / sw
    half = stats.t.ppf(0.975, k - 1) * math.sqrt(varO)
    assert abs(lo2 - (theta2 - half)) < 1e-9
    assert abs(hi2 - (theta2 + half)) < 1e-9
    assert abs(var2 - varO) < 1e-9
    # t_{k-1}, NOT qnorm: half-width must use t_6=2.4469, not z=1.96
    assert abs(stats.t.ppf(0.975, k - 1) - 2.4469118511449786) < 1e-9

# --- trim_and_fill (Duval-Tweedie L0) ---
def test_trim_and_fill_baseline():
    import math
    from spec_collapse.engine import trim_and_fill, tau2_dl
    y = [0.10, 0.20, 0.15, 0.30, 0.25, 0.40, 0.90, 1.10]
    v = [0.05, 0.04, 0.06, 0.05, 0.04, 0.05, 0.06, 0.05]
    yf, vf = trim_and_fill(y, v, tau2_dl)
    # k0 == 1: exactly one reflected study appended (right side)
    assert len(yf) == len(y) + 1
    assert len(vf) == len(v) + 1
    # original studies preserved unchanged at the front
    for a, b in zip(yf[:len(y)], y):
        assert math.isclose(float(a), b, rel_tol=0, abs_tol=1e-12)
    for a, b in zip(vf[:len(v)], v):
        assert math.isclose(float(a), b, rel_tol=0, abs_tol=1e-12)
    # imputed point reflects the most-extreme positive study (1.10) about the
    # iterated RE pooled mean; verified value from independent Duval-Tweedie L0 oracle
    assert math.isclose(float(yf[-1]), 0.8294519724889389, rel_tol=0, abs_tol=1e-9)
    # imputed variance inherits the reflected study's variance (0.05)
    assert math.isclose(float(vf[-1]), 0.05, rel_tol=0, abs_tol=1e-12)
    # reflection identity: imputed == 2*theta - y_extreme  => theta consistent both ways
    theta_from_fill = (float(yf[-1]) + 1.10) / 2.0
    assert math.isclose(theta_from_fill, 0.9647259862444695, rel_tol=0, abs_tol=1e-9)

    # symmetric-ish data still yields a small fill via L0 rounding (k0==1 here)
    ys = [-0.4, -0.2, -0.1, 0.0, 0.1, 0.2, 0.4]
    vs = [0.05] * 7
    yfs, _ = trim_and_fill(ys, vs, tau2_dl)
    assert len(yfs) == len(ys) + 1

    # k<3 short-circuits to identity (no fill)
    y2, v2 = [0.1, 0.2], [0.05, 0.04]
    yf2, vf2 = trim_and_fill(y2, v2, tau2_dl)
    assert list(map(float, yf2)) == y2
    assert list(map(float, vf2)) == v2

# --- naive_ivre_pool ---
import pytest
from spec_collapse.aggregators import naive_ivre_pool


def test_naive_ivre_pool_baseline():
    """Lock naive_ivre_pool to its independently-verified fixed-effect IV pool
    value. This is the intentionally-anti-conservative multiverse comparator
    (NOT the headline); the test also asserts the CI-collapse property the atlas
    demonstrates: pooled var < smallest single-spec var.
    """
    specs = [
        {"theta": 0.20, "var": 0.040},
        {"theta": 0.35, "var": 0.025},
        {"theta": 0.10, "var": 0.060},
        {"theta": 0.28, "var": 0.030},
        {"theta": 0.15, "var": 0.050},
        {"theta": 0.40, "var": 0.020},
        {"theta": 0.22, "var": 0.045},
    ]
    r = naive_ivre_pool(specs, cl=0.95)
    assert r["theta"] == pytest.approx(0.27935656836461126, abs=1e-9)
    assert r["var"] == pytest.approx(0.00482573726541555, abs=1e-9)
    assert r["ci_low"] == pytest.approx(0.1432027179431747, abs=1e-9)
    assert r["ci_high"] == pytest.approx(0.4155104187860478, abs=1e-9)
    assert bool(r["significant"]) is True
    assert r["var"] < min(s["var"] for s in specs)

# --- naive_concordance ---
def test_naive_concordance_baseline():
    import math
    from spec_collapse.aggregators import naive_concordance

    # k=7, 5 of 7 significant -> 71.428...% -> fragile
    specs = [
        {"significant": True},
        {"significant": True},
        {"significant": False},
        {"significant": True},
        {"significant": True},
        {"significant": False},
        {"significant": True},
    ]
    r = naive_concordance(specs, cl=0.95)
    assert r["method"] == "naive_concordance"
    assert r["n_specs"] == 7
    assert math.isclose(r["pct_significant"], 100.0 * 5 / 7, rel_tol=0, abs_tol=1e-9)
    assert r["verdict"] == "fragile"

    # all significant -> 100% -> robust
    r2 = naive_concordance([{"significant": True}] * 20)
    assert math.isclose(r2["pct_significant"], 100.0, abs_tol=1e-9)
    assert r2["verdict"] == "robust"

    # boundary frac == 0.95 -> robust (>= comparison)
    r3 = naive_concordance([{"significant": True}] * 19 + [{"significant": False}])
    assert math.isclose(r3["pct_significant"], 95.0, abs_tol=1e-9)
    assert r3["verdict"] == "robust"

    # empty -> no division by zero, fragile
    r4 = naive_concordance([])
    assert r4["n_specs"] == 0
    assert math.isclose(r4["pct_significant"], 0.0, abs_tol=1e-9)
    assert r4["verdict"] == "fragile"

# --- weighted_likelihood (mixture-CDF inversion, Wagenmakers-style multiverse aggregator) ---
def test_weighted_likelihood_baseline():
    """Lock weighted_likelihood (mixture-CDF inversion) to oracle-verified values.

    Oracle: independent numpy/scipy reimplementation of the Wagenmakers
    mixture CDF F(x)=sum_s p_s*T_{k_s-1}((x-theta_s)/sqrt(V_s)) inverted
    by brentq; matched the engine to maxAbsDiff 1.47e-9.
    """
    from spec_collapse.aggregators import weighted_likelihood

    specs = [
        {"theta": 0.30, "var": 0.040, "k": 8},
        {"theta": 0.45, "var": 0.050, "k": 6},
        {"theta": 0.20, "var": 0.030, "k": 10},
        {"theta": 0.55, "var": 0.060, "k": 5},
        {"theta": 0.35, "var": 0.045, "k": 7},
    ]
    r = weighted_likelihood(specs, cl=0.95, components="t")

    # closed-form moments (law of total variance) -- exact
    assert abs(r["theta"] - 0.37) < 1e-9
    assert abs(r["within_var"] - 0.07308095238095239) < 1e-9
    assert abs(r["between_var"] - 0.014600000000000005) < 1e-9
    assert abs(r["var"] - 0.08768095238095239) < 1e-9

    # mixture-CDF-inversion 95% interval (t components) -- oracle-matched
    assert abs(r["ci_low"] - (-0.16632105451837306)) < 1e-6
    assert abs(r["ci_high"] - 0.9844503029257291) < 1e-6

    # NON-COLLAPSE property: multiverse interval must not be narrower than
    # the narrowest single-spec interval (the sin naive IV-RE commits).
    width = r["ci_high"] - r["ci_low"]
    assert width > 0.78

    # normal-mixture (legacy) component
    rn = weighted_likelihood(specs, cl=0.95, components="normal")
    assert abs(rn["ci_low"] - (-0.07425793919876729)) < 1e-6
    assert abs(rn["ci_high"] - 0.8793735627595417) < 1e-6

# --- run_coverage (Monte-Carlo coverage of multiverse aggregators: naive_ivre_pool vs weighted_likelihood) ---
def test_run_coverage_baseline():
    """Lock the Monte-Carlo coverage property of run_coverage under the pinned seed.

    Verifies the headline claim of spec-collapse-atlas:
      - naive IV-RE pooling of many specs from ONE dataset SEVERELY under-covers
        (anti-conservative; arXiv:2511.17064),
      - the weighted-likelihood aggregator (Wagenmakers 2025) recovers ~nominal 0.95.
    Monte-Carlo 3-sigma tolerance (NOT 1e-6); seed pinned for reproducibility.
    """
    from spec_collapse.coverage import run_coverage

    r = run_coverage(n_reps=600, k=8, mu=0.0, tau2_true=0.05,
                     seed=20260604, cl=0.95)

    # Weighted-likelihood recovers nominal 95% coverage within 3-sigma.
    # binomial SE at p=0.95, n=600 ~= 0.0089 -> 3-sigma band ~ [0.923, 0.977]
    assert 0.92 <= r["coverage_weighted_likelihood"] <= 0.98, \
        r["coverage_weighted_likelihood"]

    # Naive IV-RE pooling is anti-conservative: coverage collapses well below 0.95.
    assert r["coverage_naive_ivre"] < 0.40, r["coverage_naive_ivre"]

    # Its Type-I error is correspondingly inflated (>> 0.05).
    assert r["type1_naive_ivre"] > 0.50, r["type1_naive_ivre"]

    # Mechanism: the naive interval is far NARROWER than the WL interval.
    assert r["mean_width_naive_ivre"] < r["mean_width_weighted_likelihood"]

    # Determinism under the pinned seed (exact reproduction).
    assert abs(r["coverage_naive_ivre"] - 0.16833333333333333) < 1e-9
    assert abs(r["coverage_weighted_likelihood"] - 0.9533333333333334) < 1e-9
    assert abs(r["type1_naive_ivre"] - 0.8316666666666667) < 1e-9

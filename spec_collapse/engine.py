"""Per-specification meta-analysis engine.

Implements the tau^2 estimators and CI methods that MultiverseMA enumerates,
then a small Cartesian spec grid over (estimator x CI method x outlier handling).
Each spec returns one pooled estimate, its variance, and a 2-sided CI.

Statistical conventions (cross-checked against advanced-stats.md):
  - Pool on the log scale; back-transform is the caller's concern.
  - HKSJ q-floor: when Q < k-1 the raw HKSJ scale q narrows the CI below DL;
    we floor q at 1.0 (== max(1, Q/(k-1))) and use t_{k-1}, never z.
  - REML via the standard Viechtbauer fixed-point iteration, tau^2 clamped >= 0.
  - PM (Paule-Mandel) by bracketed root-find of the generalised-Q estimating eqn.
"""

from __future__ import annotations

import math
from itertools import product

from scipy import optimize, stats

# --------------------------------------------------------------------------
# tau^2 estimators
# --------------------------------------------------------------------------


def _fe_mean(yi, vi):
    w = [1.0 / v for v in vi]
    sw = sum(w)
    return sum(wi * y for wi, y in zip(w, yi)) / sw, w, sw


def tau2_dl(yi, vi):
    """DerSimonian-Laird. Closed form."""
    k = len(yi)
    if k < 2:
        return 0.0
    ybar, w, sw = _fe_mean(yi, vi)
    Q = sum(wi * (y - ybar) ** 2 for wi, y in zip(w, yi))
    sw2 = sum(wi * wi for wi in w)
    C = sw - sw2 / sw
    if C <= 0:
        return 0.0
    return max(0.0, (Q - (k - 1)) / C)


def _generalised_q(tau2, yi, vi):
    """Sum w_i (y_i - ybar)^2 with w_i = 1/(v_i + tau2). Used by PM."""
    w = [1.0 / (v + tau2) for v in vi]
    sw = sum(w)
    ybar = sum(wi * y for wi, y in zip(w, yi)) / sw
    return sum(wi * (y - ybar) ** 2 for wi, y in zip(w, yi))


def tau2_pm(yi, vi):
    """Paule-Mandel: solve generalised Q(tau2) = k - 1."""
    k = len(yi)
    if k < 2:
        return 0.0
    target = k - 1
    if _generalised_q(0.0, yi, vi) <= target:
        return 0.0
    f = lambda t: _generalised_q(t, yi, vi) - target
    hi = max(vi) * 10 + 1.0
    # expand bracket until sign change (f is monotone decreasing in tau2)
    for _ in range(60):
        if f(hi) < 0:
            break
        hi *= 2
    try:
        return float(optimize.brentq(f, 0.0, hi, xtol=1e-10, maxiter=200))
    except ValueError:
        return tau2_dl(yi, vi)


def tau2_reml(yi, vi, tol=1e-9, max_iter=200):
    """REML via Viechtbauer fixed-point iteration, seeded from DL."""
    k = len(yi)
    if k < 2:
        return 0.0
    tau2 = tau2_dl(yi, vi)
    for _ in range(max_iter):
        w = [1.0 / (v + tau2) for v in vi]
        sw = sum(w)
        ybar = sum(wi * y for wi, y in zip(w, yi)) / sw
        sw2 = sum(wi * wi for wi in w)
        num = sum(wi * wi * ((y - ybar) ** 2 + 1.0 / sw - v)
                  for wi, y, v in zip(w, yi, vi))
        new = num / sw2
        new = max(0.0, new)
        if abs(new - tau2) < tol:
            tau2 = new
            break
        tau2 = new
    return tau2


TAU2_ESTIMATORS = {"DL": tau2_dl, "REML": tau2_reml, "PM": tau2_pm}

# --------------------------------------------------------------------------
# pooling + CI
# --------------------------------------------------------------------------


def re_pool(yi, vi, tau2):
    """Random-effects pooled estimate and its inverse-variance (Wald) variance."""
    w = [1.0 / (v + tau2) for v in vi]
    sw = sum(w)
    theta = sum(wi * y for wi, y in zip(w, yi)) / sw
    var_wald = 1.0 / sw
    return theta, var_wald, w, sw


def ci_wald(theta, var_wald, cl=0.95):
    z = stats.norm.ppf(0.5 + cl / 2)
    half = z * math.sqrt(var_wald)
    return theta - half, theta + half, var_wald


def ci_hksj(yi, vi, theta, tau2, cl=0.95):
    """Hartung-Knapp-Sidik-Jonkman with the q>=1 floor and t_{k-1}.

    HKSJ is undefined for k<2 (no between-study df); fall back to Wald there.
    """
    k = len(yi)
    if k < 2:
        var_wald = 1.0 / sum(1.0 / (v + tau2) for v in vi)
        return ci_wald(theta, var_wald, cl)
    w = [1.0 / (v + tau2) for v in vi]
    sw = sum(w)
    q = sum(wi * (y - theta) ** 2 for wi, y in zip(w, yi)) / (k - 1)
    q = max(1.0, q)  # floor: never narrower than DL/Wald  (advanced-stats.md)
    var_hksj = q / sw
    tcrit = stats.t.ppf(0.5 + cl / 2, k - 1)
    half = tcrit * math.sqrt(var_hksj)
    return theta - half, theta + half, var_hksj


CI_METHODS = ("Wald", "HKSJ")

# --------------------------------------------------------------------------
# outlier-handling dimension
# --------------------------------------------------------------------------


def _outlier_variants(yi, vi):
    """Yield (label, yi, vi) for: keep all, drop most extreme, drop 2 extreme.

    'Extreme' = largest standardised deviation from the FE mean. Mirrors the
    leave-k-out outlier dimension MultiverseMA exposes.
    """
    yield "keep_all", list(yi), list(vi)
    k = len(yi)
    if k <= 3:
        return
    ybar, _, _ = _fe_mean(yi, vi)
    order = sorted(range(k), key=lambda i: abs(yi[i] - ybar) / math.sqrt(vi[i]),
                   reverse=True)
    drop1 = set(order[:1])
    yield ("drop_1_extreme",
           [yi[i] for i in range(k) if i not in drop1],
           [vi[i] for i in range(k) if i not in drop1])
    if k > 4:
        drop2 = set(order[:2])
        yield ("drop_2_extreme",
               [yi[i] for i in range(k) if i not in drop2],
               [vi[i] for i in range(k) if i not in drop2])


# --------------------------------------------------------------------------
# spec enumeration
# --------------------------------------------------------------------------


class Spec(dict):
    """One specification result. Keys: spec_id, estimator, ci_method, outlier,
    theta, var, ci_low, ci_high, tau2, k, significant."""


def enumerate_specs(yi, vi, cl=0.95):
    """Full Cartesian product: outlier x estimator x CI method.

    Returns a list of Spec. `significant` = CI excludes 0.
    """
    specs = []
    sid = 0
    for olabel, y, v in _outlier_variants(yi, vi):
        for est_name, est in TAU2_ESTIMATORS.items():
            tau2 = est(y, v)
            theta, var_wald, _, _ = re_pool(y, v, tau2)
            for cim in CI_METHODS:
                if cim == "Wald":
                    lo, hi, var = ci_wald(theta, var_wald, cl)
                else:
                    lo, hi, var = ci_hksj(y, v, theta, tau2, cl)
                sid += 1
                specs.append(Spec(
                    spec_id=sid, estimator=est_name, ci_method=cim,
                    outlier=olabel, theta=theta, var=var, ci_low=lo, ci_high=hi,
                    tau2=tau2, k=len(y), significant=(lo > 0 or hi < 0),
                ))
    return specs

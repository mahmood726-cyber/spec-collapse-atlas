"""Three ways to summarise a spec-curve into one verdict/interval.

(1) naive_concordance  -- "% of specs significant". Not an interval at all; the
    number everyone reports ("100% concordance"). Has no calibrated coverage.

(2) naive_ivre_pool    -- inverse-variance pool the S spec estimates AS IF they
    were S independent studies. This is the cardinal sin advanced-stats.md names
    ("Never IV-RE-pool many-analyst / multiverse results"): the S specs come from
    ONE dataset, so the variance collapses by ~S and the CI is anti-conservative.

(3) weighted_likelihood -- model-average the per-spec likelihoods (a Gaussian
    mixture with weights p_s). The interval comes from the mixture quantiles, so
    its variance = mean within-spec variance + between-spec spread
    (law of total variance) and is NEVER narrower than a single spec. This is the
    Wagenmakers-style correct combination.
"""

from __future__ import annotations

import math

import numpy as np
from scipy import optimize, stats


def build_weights(specs, scheme="uniform"):
    """Return per-spec weights for a weighting scheme, or None for uniform.

    - "uniform": equal weights (returns None; weighted_likelihood handles it).
    - "reml_only" / "hksj_only": 0/1 weights restricting to the recommended
      estimator / the better-calibrated CI method (a defensible spec-subset
      sensitivity). Falls back to uniform if the restriction empties the set.
    - "aic": w_i proportional to exp(-0.5 * (AIC_i - min AIC)), AIC_i =
      -2*loglik_i + 4 (theta, tau^2). HEURISTIC only -- specs span different data
      subsets (outlier drops, trim-fill), so the likelihoods are not strictly
      comparable; used as a sensitivity probe, never the primary summary.
    """
    if scheme == "uniform":
        return None
    if scheme == "reml_only":
        w = [1.0 if s.get("estimator") == "REML" else 0.0 for s in specs]
    elif scheme == "hksj_only":
        w = [1.0 if s.get("ci_method") == "HKSJ" else 0.0 for s in specs]
    elif scheme == "aic":
        aics = [-2.0 * s.get("loglik", 0.0) + 4.0 for s in specs]
        amin = min(aics)
        w = [math.exp(-0.5 * (a - amin)) for a in aics]
    else:
        raise ValueError(f"unknown weighting scheme: {scheme}")
    if sum(w) <= 0:
        return None
    return w


def naive_concordance(specs, cl=0.95):
    n = len(specs)
    nsig = sum(1 for s in specs if s["significant"])
    frac = nsig / n if n else 0.0
    return {
        "method": "naive_concordance",
        "pct_significant": 100.0 * frac,
        "n_specs": n,
        # the verdict people read off this number:
        "verdict": "robust" if frac >= 0.95 else "fragile",
    }


def naive_ivre_pool(specs, cl=0.95):
    thetas = [s["theta"] for s in specs]
    vars = [s["var"] for s in specs]
    inv = [1.0 / v for v in vars]
    sinv = sum(inv)
    theta = sum(t / v for t, v in zip(thetas, vars)) / sinv
    var = 1.0 / sinv  # <-- collapses by ~number of specs
    z = stats.norm.ppf(0.5 + cl / 2)
    half = z * math.sqrt(var)
    lo, hi = theta - half, theta + half
    return {
        "method": "naive_ivre_pool",
        "theta": theta, "var": var, "ci_low": lo, "ci_high": hi,
        "significant": (lo > 0 or hi < 0),
        "verdict": "robust" if (lo > 0 or hi < 0) else "fragile",
    }


def _mixture_cdf(x, thetas, sds, dfs, p, use_t):
    """Vectorised mixture CDF: one scipy call over all components."""
    z = (x - thetas) / sds
    c = stats.t.cdf(z, dfs) if use_t else stats.norm.cdf(z)
    return float(np.dot(p, c))


def weighted_likelihood(specs, cl=0.95, weights=None, components="t"):
    """Mixture (model-averaged) interval via numeric quantile inversion.

    components="t" (default): each specification contributes a scaled-t density
    `theta_s + sqrt(V_s) * t_{df_s}` with df_s = k_s - 1, reflecting that the
    per-spec pivot is t-distributed at small k. This corrects the mild
    under-coverage the normal-mixture showed at high heterogeneity.
    components="normal": Gaussian mixture (legacy).
    """
    thetas = [s["theta"] for s in specs]
    vars = [s["var"] for s in specs]
    sds = [math.sqrt(v) for v in vars]
    dfs = [max(1, int(s.get("k", 2)) - 1) for s in specs]
    use_t = components == "t"
    n = len(specs)
    if weights is None:
        p = [1.0 / n] * n
    else:
        sw = sum(weights)
        p = [w / sw for w in weights]

    mean = sum(pi * t for pi, t in zip(p, thetas))
    # law of total variance; t-component inflates within-var by df/(df-2) for df>2
    within = 0.0
    for pi, v, df in zip(p, vars, dfs):
        scale = df / (df - 2) if (use_t and df > 2) else 1.0
        within += pi * v * scale
    between = sum(pi * (t - mean) ** 2 for pi, t in zip(p, thetas))
    total_var = within + between

    alpha = (1 - cl) / 2
    spread = [sd * (stats.t.ppf(0.999, df) if use_t else 3.1) for sd, df in zip(sds, dfs)]
    lo_t = min(t - 2 * s for t, s in zip(thetas, spread))
    hi_t = max(t + 2 * s for t, s in zip(thetas, spread))
    # numpy arrays for the vectorised mixture-CDF inversion
    aT, aS, aD, aP = (np.asarray(thetas), np.asarray(sds),
                      np.asarray(dfs, float), np.asarray(p))
    lo = optimize.brentq(lambda x: _mixture_cdf(x, aT, aS, aD, aP, use_t) - alpha,
                         lo_t, hi_t, xtol=1e-8)
    hi = optimize.brentq(lambda x: _mixture_cdf(x, aT, aS, aD, aP, use_t) - (1 - alpha),
                         lo_t, hi_t, xtol=1e-8)
    return {
        "method": "weighted_likelihood",
        "theta": mean, "var": total_var,
        "within_var": within, "between_var": between,
        "ci_low": lo, "ci_high": hi,
        "significant": (lo > 0 or hi < 0),
        "verdict": "robust" if (lo > 0 or hi < 0) else "fragile",
    }

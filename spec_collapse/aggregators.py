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

from scipy import optimize, stats


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


def _mixture_cdf(x, thetas, sds, p):
    return sum(pi * stats.norm.cdf((x - t) / sd)
              for pi, t, sd in zip(p, thetas, sds))


def weighted_likelihood(specs, cl=0.95, weights=None):
    """Gaussian-mixture (model-averaged) interval via numeric quantile inversion."""
    thetas = [s["theta"] for s in specs]
    vars = [s["var"] for s in specs]
    sds = [math.sqrt(v) for v in vars]
    n = len(specs)
    if weights is None:
        p = [1.0 / n] * n
    else:
        sw = sum(weights)
        p = [w / sw for w in weights]

    mean = sum(pi * t for pi, t in zip(p, thetas))
    within = sum(pi * v for pi, v in zip(p, vars))
    between = sum(pi * (t - mean) ** 2 for pi, t in zip(p, thetas))
    total_var = within + between  # law of total variance

    alpha = (1 - cl) / 2
    lo_t = min(t - 6 * sd for t, sd in zip(thetas, sds))
    hi_t = max(t + 6 * sd for t, sd in zip(thetas, sds))
    lo = optimize.brentq(lambda x: _mixture_cdf(x, thetas, sds, p) - alpha,
                         lo_t, hi_t, xtol=1e-8)
    hi = optimize.brentq(lambda x: _mixture_cdf(x, thetas, sds, p) - (1 - alpha),
                         lo_t, hi_t, xtol=1e-8)
    return {
        "method": "weighted_likelihood",
        "theta": mean, "var": total_var,
        "within_var": within, "between_var": between,
        "ci_low": lo, "ci_high": hi,
        "significant": (lo > 0 or hi < 0),
        "verdict": "robust" if (lo > 0 or hi < 0) else "fragile",
    }

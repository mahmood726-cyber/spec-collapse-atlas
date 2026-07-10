"""INDEPENDENT reproduction check of the spec-collapse-atlas headline claim.

This script DELIBERATELY does NOT import spec_collapse.engine or
spec_collapse.aggregators. It re-implements, from the mathematical definitions in
README.md / docs/spec.md, the three things needed to reproduce the claim:

  1. the 36-spec grid (3 tau2 estimators x 2 CI methods x 3 outlier rules
     x {raw, trim-and-fill}),
  2. the NAIVE IV-RE pool rule  (pool the S spec estimates as independent),
  3. the CORRECTED weighted-likelihood (scaled-t mixture) interval.

The ONLY shared dependency is the .rda ingest loader from fragility-atlas, which
is the data-reading layer (explicitly out of scope for the novel claim; the
authors themselves "stand on" it for ingest).

Claim under test (v0.2.0):
  N=473 Cochrane MAs; naive IV-RE labels ~88% robust, corrected WL ~33% robust,
  false-robustness ~55%; null Type-I error of naive rule ~70-81%.

Run: python repro_check.py
"""
from __future__ import annotations

import glob
import json
import math
import os
import sys

import numpy as np
from scipy import optimize, stats

# ---- ingest only (shared, out of scope for the claim) ----------------------
FRAGILITY_ATLAS = os.environ.get("FRAGILITY_ATLAS_PATH", r"C:\Projects\fragility-atlas")
CORPUS = os.environ.get("SPEC_COLLAPSE_CORPUS", r"C:\Projects\Pairwise70\data")
if FRAGILITY_ATLAS not in sys.path:
    sys.path.insert(0, FRAGILITY_ATLAS)
from src.loader import load_review  # noqa: E402


# ===========================================================================
# 1. tau^2 estimators (independent re-derivation)
# ===========================================================================
def fe_mean(y, v):
    w = np.asarray([1.0 / vi for vi in v])
    return float(np.sum(w * y) / np.sum(w)), w


def tau2_DL(y, v):
    y = np.asarray(y, float); v = np.asarray(v, float)
    k = len(y)
    if k < 2:
        return 0.0
    ybar, w = fe_mean(y, v)
    Q = float(np.sum(w * (y - ybar) ** 2))
    C = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    if C <= 0:
        return 0.0
    return max(0.0, (Q - (k - 1)) / C)


def _genQ(t2, y, v):
    w = 1.0 / (v + t2)
    ybar = np.sum(w * y) / np.sum(w)
    return float(np.sum(w * (y - ybar) ** 2))


def tau2_PM(y, v):
    y = np.asarray(y, float); v = np.asarray(v, float)
    k = len(y)
    if k < 2:
        return 0.0
    target = k - 1
    if _genQ(0.0, y, v) <= target:
        return 0.0
    f = lambda t: _genQ(t, y, v) - target
    hi = float(np.max(v)) * 10 + 1.0
    for _ in range(80):
        if f(hi) < 0:
            break
        hi *= 2
    try:
        return float(optimize.brentq(f, 0.0, hi, xtol=1e-10, maxiter=300))
    except ValueError:
        return tau2_DL(y, v)


def tau2_REML(y, v, tol=1e-10, it=300):
    """REML by Viechtbauer fixed-point iteration, seeded from DL."""
    y = np.asarray(y, float); v = np.asarray(v, float)
    k = len(y)
    if k < 2:
        return 0.0
    t2 = tau2_DL(y, v)
    for _ in range(it):
        w = 1.0 / (v + t2)
        sw = np.sum(w)
        ybar = np.sum(w * y) / sw
        num = np.sum(w ** 2 * ((y - ybar) ** 2 + 1.0 / sw - v))
        new = max(0.0, float(num / np.sum(w ** 2)))
        if abs(new - t2) < tol:
            return new
        t2 = new
    return t2


ESTIMATORS = {"DL": tau2_DL, "REML": tau2_REML, "PM": tau2_PM}


# ===========================================================================
# 2. RE pooling + CI methods
# ===========================================================================
def re_pool(y, v, t2):
    y = np.asarray(y, float); v = np.asarray(v, float)
    w = 1.0 / (v + t2)
    sw = np.sum(w)
    theta = float(np.sum(w * y) / sw)
    return theta, float(1.0 / sw)


def ci_wald(theta, var_wald, cl=0.95):
    z = stats.norm.ppf(0.5 + cl / 2)
    h = z * math.sqrt(var_wald)
    return theta - h, theta + h, var_wald


def ci_hksj(y, v, theta, t2, cl=0.95):
    y = np.asarray(y, float); v = np.asarray(v, float)
    k = len(y)
    w = 1.0 / (v + t2)
    sw = np.sum(w)
    if k < 2:
        return ci_wald(theta, 1.0 / sw, cl)
    q = float(np.sum(w * (y - theta) ** 2) / (k - 1))
    q = max(1.0, q)                      # floor max(1, Q/(k-1))
    var = q / sw
    tcrit = stats.t.ppf(0.5 + cl / 2, k - 1)
    h = tcrit * math.sqrt(var)
    return theta - h, theta + h, var


# ===========================================================================
# 3. outlier variants + Duval-Tweedie trim-and-fill
# ===========================================================================
def outlier_variants(y, v):
    y = np.asarray(y, float); v = np.asarray(v, float)
    k = len(y)
    yield "keep_all", y, v
    if k <= 3:
        return
    ybar, _ = fe_mean(y, v)
    order = np.argsort(-np.abs(y - ybar) / np.sqrt(v))
    keep1 = np.ones(k, bool); keep1[order[:1]] = False
    yield "drop1", y[keep1], v[keep1]
    if k > 4:
        keep2 = np.ones(k, bool); keep2[order[:2]] = False
        yield "drop2", y[keep2], v[keep2]


def trim_and_fill(y, v, est):
    y = np.asarray(y, float); v = np.asarray(v, float)
    k = len(y)
    if k < 3:
        return y, v
    t2 = est(y, v)
    theta0 = re_pool(y, v, t2)[0]
    di = y - theta0
    side = "right" if (di > 0).sum() > (di < 0).sum() else "left"
    yc, vc = y.copy(), v.copy()
    for _ in range(20):
        t2 = est(yc, vc)
        theta0 = re_pool(yc, vc, t2)[0]
        di = y - theta0
        ranks = np.argsort(np.abs(di))
        signed = np.zeros(k)
        for i, r in enumerate(ranks):
            signed[r] = (i + 1) * np.sign(di[r])
        S = signed[signed > 0].sum() if side == "right" else np.abs(signed[signed < 0]).sum()
        k0 = max(0, round((4 * S - k * (k + 1)) / (2 * k - 1)))  # (2k-1): Duval-Tweedie/metafor
        if k0 == 0:
            break
        order = np.argsort(di)
        idx = order[-k0:] if side == "right" else order[:k0]
        yf = np.concatenate([y, 2 * theta0 - y[idx]])
        vf = np.concatenate([v, v[idx]])
        if len(yf) == len(yc) and np.allclose(yf, yc, atol=1e-10):
            break
        yc, vc = yf, vf
    return yc, vc


# ===========================================================================
# 4. 36-spec grid
# ===========================================================================
def enumerate_specs(yi, vi, cl=0.95):
    specs = []
    for olabel, y, v in outlier_variants(yi, vi):
        for ename, est in ESTIMATORS.items():
            for pb in ("raw", "trim_fill"):
                if pb == "trim_fill":
                    if len(y) < 3:
                        continue
                    yy, vv = trim_and_fill(y, v, est)
                else:
                    yy, vv = y, v
                t2 = est(yy, vv)
                theta, var_wald = re_pool(yy, vv, t2)
                for cim in ("Wald", "HKSJ"):
                    if cim == "Wald":
                        lo, hi, var = ci_wald(theta, var_wald, cl)
                    else:
                        lo, hi, var = ci_hksj(yy, vv, theta, t2, cl)
                    specs.append({
                        "estimator": ename, "ci_method": cim, "outlier": olabel,
                        "pub_bias": pb, "theta": theta, "var": var,
                        "k": len(yy), "significant": (lo > 0 or hi < 0),
                    })
    return specs


# ===========================================================================
# 5. the two aggregators under test
# ===========================================================================
def naive_ivre(specs, cl=0.95):
    """Pool the S spec estimates as if S independent studies (the 'sin')."""
    th = np.array([s["theta"] for s in specs])
    vr = np.array([s["var"] for s in specs])
    inv = 1.0 / vr
    theta = float(np.sum(th * inv) / np.sum(inv))
    var = float(1.0 / np.sum(inv))
    z = stats.norm.ppf(0.5 + cl / 2)
    h = z * math.sqrt(var)
    lo, hi = theta - h, theta + h
    return lo, hi, ("robust" if (lo > 0 or hi < 0) else "fragile")


def _mix_cdf(x, th, sd, df, p):
    return float(np.dot(p, stats.t.cdf((x - th) / sd, df)))


def weighted_likelihood(specs, cl=0.95):
    """Scaled-t mixture: theta_s + sqrt(V_s)*t_{k_s-1}, uniform weights.
    Interval by inverting the mixture CDF at alpha/2 and 1-alpha/2."""
    th = np.array([s["theta"] for s in specs])
    sd = np.sqrt(np.array([s["var"] for s in specs]))
    df = np.array([max(1, int(s["k"]) - 1) for s in specs], float)
    n = len(specs)
    p = np.full(n, 1.0 / n)
    alpha = (1 - cl) / 2
    spread = sd * stats.t.ppf(0.999, df)
    lo_t = float(np.min(th - 2 * spread))
    hi_t = float(np.max(th + 2 * spread))
    lo = optimize.brentq(lambda x: _mix_cdf(x, th, sd, df, p) - alpha, lo_t, hi_t, xtol=1e-8)
    hi = optimize.brentq(lambda x: _mix_cdf(x, th, sd, df, p) - (1 - alpha), lo_t, hi_t, xtol=1e-8)
    return lo, hi, ("robust" if (lo > 0 or hi < 0) else "fragile")


# ===========================================================================
# 6. corpus run
# ===========================================================================
def run_corpus():
    files = sorted(glob.glob(os.path.join(CORPUS, "CD*.rda")))
    n = ivre_robust = wl_robust = false_robust = 0
    skipped = errored = 0
    widths = []
    for f in files:
        try:
            r = load_review(f)
        except Exception:
            errored += 1
            continue
        if r is None or r.k < 3:
            skipped += 1
            continue
        yi = np.asarray(r.yi, float)
        vi = np.asarray([float(s) ** 2 for s in r.sei], float)
        try:
            specs = enumerate_specs(yi, vi)
            iv_lo, iv_hi, iv_v = naive_ivre(specs)
            wl_lo, wl_hi, wl_v = weighted_likelihood(specs)
        except Exception:
            errored += 1
            continue
        n += 1
        if iv_v == "robust":
            ivre_robust += 1
        if wl_v == "robust":
            wl_robust += 1
        if iv_v == "robust" and wl_v == "fragile":
            false_robust += 1
        w_iv = iv_hi - iv_lo
        w_wl = wl_hi - wl_lo
        if w_wl > 0:
            widths.append(w_iv / w_wl)
    widths = np.array(widths)
    return {
        "n_reviews": n, "skipped": skipped, "errored": errored,
        "ivre_robust": ivre_robust, "ivre_robust_pct": 100.0 * ivre_robust / n,
        "wl_robust": wl_robust, "wl_robust_pct": 100.0 * wl_robust / n,
        "false_robust_n": false_robust, "false_robust_pct": 100.0 * false_robust / n,
        "median_width_ratio": float(np.median(widths)),
    }


# ===========================================================================
# 7. null Monte-Carlo (Type-I of naive rule)
# ===========================================================================
def run_null(n_reps=1000, k=10, tau2_true=0.05, seed=20260604):
    vi_template = np.array([0.0154, 0.0200, 0.0312, 0.0356, 0.0368, 0.0393,
                            0.0512, 0.0628, 0.0636, 0.1946, 0.2252, 0.3256])
    rng = np.random.default_rng(seed)
    sig_iv = cover_iv = cover_wl = 0
    for _ in range(n_reps):
        v = rng.choice(vi_template, size=k, replace=True)
        theta_i = rng.normal(0.0, math.sqrt(tau2_true), size=k)
        y = rng.normal(theta_i, np.sqrt(v))
        specs = enumerate_specs(y, v)
        iv_lo, iv_hi, iv_v = naive_ivre(specs)
        wl_lo, wl_hi, wl_v = weighted_likelihood(specs)
        if iv_v == "robust":
            sig_iv += 1
        if iv_lo <= 0.0 <= iv_hi:
            cover_iv += 1
        if wl_lo <= 0.0 <= wl_hi:
            cover_wl += 1
    return {"k": k, "tau2_true": tau2_true,
            "type1_naive_ivre": sig_iv / n_reps,
            "coverage_naive_ivre": cover_iv / n_reps,
            "coverage_wl": cover_wl / n_reps}


if __name__ == "__main__":
    print("=== INDEPENDENT REPRODUCTION: spec-collapse-atlas ===\n")
    c = run_corpus()
    print(f"N reviews (k>=3): {c['n_reviews']}  (skipped {c['skipped']}, errored {c['errored']})")
    print(f"Naive IV-RE robust : {c['ivre_robust']}/{c['n_reviews']} = {c['ivre_robust_pct']:.1f}%")
    print(f"Corrected WL robust: {c['wl_robust']}/{c['n_reviews']} = {c['wl_robust_pct']:.1f}%")
    print(f"False robustness   : {c['false_robust_n']}/{c['n_reviews']} = {c['false_robust_pct']:.1f}%")
    print(f"Median width ratio : {c['median_width_ratio']:.3f}x")
    print("\n--- Null Monte-Carlo (mu=0), Type-I of naive rule ---")
    rows = []
    for (k, t2) in [(10, 0.0), (10, 0.05), (10, 0.1), (20, 0.05), (5, 0.05)]:
        nr = run_null(n_reps=1000, k=k, tau2_true=t2)
        rows.append(nr)
        print(f"k={k:2d} tau2={t2:<5}: naive Type-I={nr['type1_naive_ivre']:.3f} "
              f"cover_iv={nr['coverage_naive_ivre']:.3f} cover_wl={nr['coverage_wl']:.3f}")
    out = {"corpus": c, "null": rows}
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "data", "repro_check_results.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("\nWrote data/repro_check_results.json")

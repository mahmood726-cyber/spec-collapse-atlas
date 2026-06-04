"""Headline demo: the spec-curve collapse on MultiverseMA's own datasets.

Run:  python demo.py
"""
import io
import sys

# Windows cp1252 console safety (lessons.md)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

from spec_collapse import (DATASETS, enumerate_specs, naive_concordance,
                           naive_ivre_pool, weighted_likelihood, run_coverage)


def fmt_ci(a):
    return f"[{a['ci_low']:+.3f}, {a['ci_high']:+.3f}]"


def real_data_table():
    print("=" * 78)
    print("REAL DATA: naive concordance / IV-RE pool  vs  weighted-likelihood")
    print("=" * 78)
    flips = 0
    for key in ("bcg", "aspirin", "omega3", "magnesium", "corticosteroids"):
        d = DATASETS[key]
        specs = enumerate_specs(d["yi"], d["vi"])
        con = naive_concordance(specs)
        ivre = naive_ivre_pool(specs)
        wl = weighted_likelihood(specs)
        # the dangerous flip: the IV-RE *pool* (what people report as "the"
        # multiverse summary) says robust, but the honest corrected interval
        # crosses zero.
        flip = ivre["verdict"] == "robust" and wl["verdict"] == "fragile"
        flips += int(flip)
        print(f"\n{d['name']}  (k={specs[0]['k']}, {len(specs)} specs, {d['measure']})")
        print(f"  naive concordance : {con['pct_significant']:.0f}% significant"
              f"  -> verdict={con['verdict']}")
        print(f"  naive IV-RE pool  : theta={ivre['theta']:+.3f}  CI={fmt_ci(ivre)}"
              f"  width={ivre['ci_high']-ivre['ci_low']:.3f}  -> {ivre['verdict']}")
        print(f"  weighted-likelihd : theta={wl['theta']:+.3f}  CI={fmt_ci(wl)}"
              f"  width={wl['ci_high']-wl['ci_low']:.3f}  -> {wl['verdict']}")
        print(f"      within-spec var={wl['within_var']:.4f}  "
              f"between-spec var={wl['between_var']:.4f}  "
              f"IV-RE width / WL width = {(ivre['ci_high']-ivre['ci_low'])/(wl['ci_high']-wl['ci_low']):.2f}x")
        if flip:
            print("      *** VERDICT FLIP: IV-RE pool says 'robust' -> corrected interval 'fragile' ***")
    print(f"\nVerdict flips (IV-RE pool robust -> corrected fragile): {flips} / 5 datasets")


def coverage_table():
    print("\n" + "=" * 78)
    print("MONTE-CARLO COVERAGE under the NULL (mu=0): true coverage of 95% CIs")
    print("=" * 78)
    print(f"{'tau2_true':>9} {'k':>3} | {'IV-RE cov':>10} {'WL cov':>8} "
          f"{'IV-RE type-I':>13} {'concord.false-robust':>21} {'IVRE/WL width':>14}")
    for tau2 in (0.0, 0.02, 0.05):
        for k in (8, 12):
            r = run_coverage(n_reps=400, k=k, mu=0.0, tau2_true=tau2)
            ratio = r["mean_width_naive_ivre"] / r["mean_width_weighted_likelihood"]
            print(f"{tau2:>9.2f} {k:>3} | {r['coverage_naive_ivre']:>10.3f} "
                  f"{r['coverage_weighted_likelihood']:>8.3f} "
                  f"{r['type1_naive_ivre']:>13.3f} "
                  f"{r['false_robust_concordance']:>21.3f} {ratio:>13.2f}x")
    print("\nNominal coverage is 0.95. IV-RE pool collapses far below; "
          "weighted-likelihood restores it.")


if __name__ == "__main__":
    real_data_table()
    coverage_table()

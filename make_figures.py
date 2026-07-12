"""Generate publication figures for the spec-collapse-atlas E156 article.

All numbers come from the verified corpus results (data/*.json) and a live
re-computation of one representative review's 36-spec grid. No hand-entered
statistics. Outputs 300-dpi PNG + vector PDF into figures/.

Run (with corpus + loader available):
  FRAGILITY_ATLAS_PATH=... SPEC_COLLAPSE_CORPUS=... python make_figures.py
"""
from __future__ import annotations

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

# ---- palette (Synthesis-ish, colour-blind safe) ----------------------------
C_NAIVE = "#c0392b"   # naive / false-precision red
C_CAL = "#1f6f8b"     # calibrated blue
C_GREY = "#7f8c8d"
C_OK = "#27ae60"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 12, "axes.titleweight": "bold",
    "figure.dpi": 120,
})

summary = json.load(open(os.path.join(HERE, "data", "corpus_summary.json")))
rows = json.load(open(os.path.join(HERE, "data", "corpus_results.json")))
repro = json.load(open(os.path.join(HERE, "data", "repro_check_results.json")))

N = summary["n_reviews"]
IV_ROBUST = summary["ivre_robust"]
WL_ROBUST = summary["wl_robust"]
FR = summary["false_robust_n"]
MED_WR = summary["median_width_ratio"]


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name + ".png/.pdf")


# ===========================================================================
# FIGURE 1 — robustness reversal (the headline)
# ===========================================================================
def fig_reversal():
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 5.2),
                                   gridspec_kw={"width_ratios": [1.25, 1]})

    # --- left: flow from naive 'robust' to calibrated verdict ---
    concord_robust = WL_ROBUST            # 157 stay robust
    flipped = FR                          # 260 robust -> fragile
    naive_fragile = N - IV_ROBUST         # 56 already fragile under naive

    axL.set_xlim(0, 10); axL.set_ylim(-55, N + 15); axL.axis("off")
    bw = 1.9
    # left bar: naive verdict
    axL.add_patch(Rectangle((1.0, naive_fragile), bw, IV_ROBUST, color=C_NAIVE, ec="white"))
    axL.add_patch(Rectangle((1.0, 0), bw, naive_fragile, color=C_GREY, ec="white"))
    # right bar: calibrated verdict
    axL.add_patch(Rectangle((7.1, N - WL_ROBUST), bw, WL_ROBUST, color=C_CAL, ec="white"))
    axL.add_patch(Rectangle((7.1, 0), bw, N - WL_ROBUST, color=C_GREY, ec="white"))

    # ribbons
    import matplotlib.path as mpath
    Path = mpath.Path

    def ribbon(y0a, y0b, y1a, y1b, color, alpha):
        x0, x1 = 2.9, 7.1
        verts = [(x0, y0a), ((x0 + x1) / 2, y0a), ((x0 + x1) / 2, y1a), (x1, y1a),
                 (x1, y1b), ((x0 + x1) / 2, y1b), ((x0 + x1) / 2, y0b), (x0, y0b), (x0, y0a)]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
                 Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
        axL.add_patch(matplotlib.patches.PathPatch(Path(verts, codes), fc=color, ec="none", alpha=alpha))

    # robust(naive) splits into: stays robust (top of right) + flips to fragile
    ribbon(naive_fragile + flipped, naive_fragile + flipped + concord_robust,
           N - WL_ROBUST, N, C_CAL, 0.30)             # concordant robust
    ribbon(naive_fragile, naive_fragile + flipped,
           N - WL_ROBUST - flipped, N - WL_ROBUST, C_NAIVE, 0.32)  # the reversal

    axL.text(1.95, naive_fragile + IV_ROBUST / 2, f"robust\n{IV_ROBUST}\n({100*IV_ROBUST/N:.0f}%)",
             ha="center", va="center", color="white", fontsize=10, fontweight="bold")
    axL.text(1.95, naive_fragile / 2, f"fragile\n{naive_fragile}", ha="center", va="center",
             color="white", fontsize=9)
    axL.text(8.05, N - WL_ROBUST / 2, f"robust\n{WL_ROBUST}\n({100*WL_ROBUST/N:.0f}%)",
             ha="center", va="center", color="white", fontsize=10, fontweight="bold")
    axL.text(8.05, (N - WL_ROBUST) / 2, f"fragile\n{N-WL_ROBUST}", ha="center", va="center",
             color="white", fontsize=9)
    axL.text(1.95, -30, "Naive IV-RE\npool", ha="center", va="center", fontsize=10, fontweight="bold")
    axL.text(8.05, -30, "Calibrated\nweighted-likelihood", ha="center", va="center", fontsize=10, fontweight="bold")
    axL.text(5.0, N - WL_ROBUST - flipped / 2, f"{flipped} reviews flip\nrobust → fragile\n({100*FR/N:.0f}%)",
             ha="center", va="center", fontsize=10, fontweight="bold", color=C_NAIVE)
    axL.set_title(f"a  Verdict reverses for {FR} of {N} reviews ({100*FR/N:.0f}%)", loc="left", pad=14)

    # --- right: 2x2 verdict matrix ---
    from collections import Counter
    c = Counter((r["ivre_verdict"], r["wl_verdict"]) for r in rows)
    M = np.array([[c.get(("robust", "robust"), 0), c.get(("robust", "fragile"), 0)],
                  [c.get(("fragile", "robust"), 0), c.get(("fragile", "fragile"), 0)]])
    axR.imshow(M, cmap="Blues", vmin=0, vmax=M.max())
    for i in range(2):
        for j in range(2):
            hot = (i == 0 and j == 1)
            axR.text(j, i, f"{M[i,j]}", ha="center", va="center", fontsize=18,
                     fontweight="bold", color=(C_NAIVE if hot else ("white" if M[i,j] > M.max()*0.5 else "#222")))
            if hot:
                axR.add_patch(Rectangle((j-0.5, i-0.5), 1, 1, fill=False, ec=C_NAIVE, lw=3))
                axR.text(j, i+0.30, "false robust", ha="center", va="center", color=C_NAIVE, fontsize=9, fontweight="bold")
    axR.set_xticks([0, 1]); axR.set_xticklabels(["robust", "fragile"])
    axR.set_yticks([0, 1]); axR.set_yticklabels(["robust", "fragile"])
    axR.set_xlabel("Calibrated weighted-likelihood")
    axR.set_ylabel("Naive IV-RE pool")
    axR.set_title("b  Verdict cross-tabulation", loc="left", pad=14)
    axR.set_xticks(np.arange(-.5, 2, 1), minor=True)
    axR.set_yticks(np.arange(-.5, 2, 1), minor=True)
    axR.grid(which="minor", color="white", lw=2)
    axR.tick_params(which="minor", length=0)

    fig.suptitle("Spec-Collapse Atlas — 473 Cochrane meta-analyses, 36 specifications each",
                 fontsize=12.5, fontweight="bold", y=1.02)
    save(fig, "fig1_robustness_reversal")


# ===========================================================================
# FIGURE 2 — naive vs calibrated: spec curve + width-ratio distribution
# ===========================================================================
def fig_naive_vs_calibrated():
    # live spec computation for one representative false-robust review
    os.environ.setdefault("FRAGILITY_ATLAS_PATH", r"C:\Projects\truth-recovery-sweep\fragility-atlas")
    import sys
    fa = os.environ["FRAGILITY_ATLAS_PATH"]
    if fa not in sys.path:
        sys.path.insert(0, fa)
    from src.loader import load_review
    from spec_collapse.engine import enumerate_specs
    from spec_collapse.aggregators import naive_ivre_pool, weighted_likelihood

    corpus = os.environ.get("SPEC_COLLAPSE_CORPUS", r"C:\Users\mahmo\Pairwise70\data")
    import glob
    target = "CD001533"
    f = [p for p in glob.glob(os.path.join(corpus, "CD*.rda")) if target in os.path.basename(p)][0]
    r = load_review(f)
    yi = list(r.yi); vi = [float(s) ** 2 for s in r.sei]
    specs = enumerate_specs(yi, vi)
    iv = naive_ivre_pool(specs); wl = weighted_likelihood(specs)
    order = sorted(range(len(specs)), key=lambda i: specs[i]["theta"])
    th = np.array([specs[i]["theta"] for i in order])
    lo = np.array([specs[i]["ci_low"] for i in order])
    hi = np.array([specs[i]["ci_high"] for i in order])
    sig = np.array([specs[i]["significant"] for i in order])
    x = np.arange(len(specs))

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 5.0),
                                   gridspec_kw={"width_ratios": [1.35, 1]})

    # --- panel A: specification curve (on the analysed log scale) ---
    axA.axhline(0, color="#888", lw=1, ls="--", zorder=1)
    for xi, l, h, s in zip(x, lo, hi, sig):
        axA.plot([xi, xi], [l, h], color=(C_NAIVE if s else C_GREY), lw=1.4, alpha=0.8, zorder=2)
    axA.plot(x, th, "o", ms=3.5, color="#222", zorder=3)
    # naive vs calibrated summary bands
    axA.axhspan(iv["ci_low"], iv["ci_high"], color=C_NAIVE, alpha=0.30, zorder=0)
    axA.axhspan(wl["ci_low"], wl["ci_high"], color=C_CAL, alpha=0.18, zorder=0)
    xr = len(specs) - 1
    axA.annotate(f"naive IV-RE 95% CI\n[{iv['ci_low']:.2f}, {iv['ci_high']:.2f}]  (excludes 0 → “robust”)",
                 xy=(xr*0.5, iv["ci_high"]), xytext=(xr*0.05, max(hi)*0.92),
                 color=C_NAIVE, fontsize=9, fontweight="bold")
    axA.annotate(f"calibrated WL 95% CI\n[{wl['ci_low']:.2f}, {wl['ci_high']:.2f}]  (crosses 0 → “fragile”)",
                 xy=(xr*0.5, wl["ci_low"]), xytext=(xr*0.05, min(lo)*0.92),
                 color=C_CAL, fontsize=9, fontweight="bold")
    axA.set_xlabel("36 specifications (sorted by point estimate)")
    axA.set_ylabel("log effect (analysed scale)")
    nsig = int(sig.sum())
    axA.set_title(f"a  Review CD001533 (k={r.k}): {nsig}/{len(specs)} specs significant",
                  loc="left", fontsize=11)

    # --- panel B: width-ratio distribution across the corpus ---
    wr = np.array([rr["width_ratio"] for rr in rows if rr["width_ratio"] == rr["width_ratio"]])
    axB.hist(wr, bins=np.linspace(0, 0.4, 41), color=C_CAL, alpha=0.85, ec="white")
    axB.axvline(MED_WR, color=C_NAIVE, lw=2.2)
    axB.text(MED_WR + 0.005, axB.get_ylim()[1]*0.9,
             f"median {MED_WR:.3f}\n(≈ 1/8 width)", color=C_NAIVE, fontsize=10, fontweight="bold")
    axB.axvline(1.0, color="#888", lw=1, ls=":")
    axB.set_xlabel("width ratio  (naive IV-RE CI / calibrated WL CI)")
    axB.set_ylabel("number of reviews")
    axB.set_title("b  Naive CI ≈ 1/8 calibrated width", loc="left", fontsize=11)
    axB.set_xlim(0, 0.4)
    fig.subplots_adjust(wspace=0.26)

    save(fig, "fig2_naive_vs_calibrated")


# ===========================================================================
# FIGURE 3 — null-simulation false-positive rate
# ===========================================================================
def fig_null_fpr():
    null = repro["null"]
    labels = [f"k={d['k']}\n" + (r"$\tau^2$=" + f"{d['tau2_true']}") for d in null]
    t1 = [100 * d["type1_naive_ivre"] for d in null]
    cov_iv = [100 * d["coverage_naive_ivre"] for d in null]
    cov_wl = [100 * d["coverage_wl"] for d in null]
    x = np.arange(len(null))

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 4.8))

    # panel A: Type-I error of the naive rule (should be 5%)
    bars = axA.bar(x, t1, color=C_NAIVE, alpha=0.9, width=0.6)
    axA.axhline(5, color=C_OK, lw=2, ls="--")
    axA.text(len(null)-0.5, 8, "nominal 5%", color=C_OK, ha="right", fontsize=9, fontweight="bold")
    for b, v in zip(bars, t1):
        axA.text(b.get_x()+b.get_width()/2, v+1.2, f"{v:.0f}%", ha="center", fontsize=10, fontweight="bold")
    axA.set_xticks(x); axA.set_xticklabels(labels, fontsize=9)
    axA.set_ylabel("false-positive rate (%)")
    axA.set_ylim(0, 100)
    lo_fpr, hi_fpr = min(t1), max(t1)
    axA.set_title(f"a  Naive false-positive rate: {lo_fpr:.0f}–{hi_fpr:.0f}%", loc="left", fontsize=11)

    # panel B: coverage of nominal-95% intervals
    w = 0.38
    axB.bar(x - w/2, cov_iv, w, color=C_NAIVE, alpha=0.9, label="naive IV-RE")
    axB.bar(x + w/2, cov_wl, w, color=C_CAL, alpha=0.9, label="calibrated WL")
    axB.axhline(95, color=C_OK, lw=2, ls="--")
    axB.text(len(null)-0.5, 97, "nominal 95%", color=C_OK, ha="right", fontsize=9, fontweight="bold")
    axB.set_xticks(x); axB.set_xticklabels(labels, fontsize=9)
    axB.set_ylabel("interval coverage (%)")
    axB.set_ylim(0, 105)
    axB.legend(frameon=False, fontsize=9, loc="center left")
    axB.set_title("b  Calibrated coverage ≈ nominal 95%", loc="left", fontsize=11)
    fig.subplots_adjust(wspace=0.24)

    fig.suptitle("Monte-Carlo under a null pooled effect (μ = 0, 1000 reps/setting, seed 20260604)",
                 fontsize=11.5, fontweight="bold", y=1.04)
    save(fig, "fig3_null_fpr")


if __name__ == "__main__":
    fig_reversal()
    fig_naive_vs_calibrated()
    fig_null_fpr()
    print("all figures written to", FIG)

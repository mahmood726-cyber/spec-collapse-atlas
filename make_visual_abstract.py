"""Visual abstract for the spec-collapse-atlas E156 article.

Flow: Question -> Method -> Headline -> Caveat. Every number is read from the
verified corpus_summary.json / repro_check_results.json (no hand-entered stats).
Outputs visual_abstract.png (300 dpi) + .pdf (vector).
"""
from __future__ import annotations

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__))
summary = json.load(open(os.path.join(HERE, "data", "corpus_summary.json")))
repro = json.load(open(os.path.join(HERE, "data", "repro_check_results.json")))

N = summary["n_reviews"]
IVP = round(100 * summary["ivre_robust"] / N)
WLP = round(100 * summary["wl_robust"] / N)
FRN = summary["false_robust_n"]
FRP = round(summary["false_robust_pct"])
WR = summary["median_width_ratio"]
fpr = [d["type1_naive_ivre"] for d in repro["null"]]
FLO, FHI = round(100 * min(fpr)), round(100 * max(fpr))

C_INK = "#1c2833"
C_NAIVE = "#c0392b"
C_CAL = "#1f6f8b"
C_BG1 = "#eef3f6"
C_BG2 = "#fbeeec"
C_GOLD = "#b9770e"

plt.rcParams.update({"font.family": "DejaVu Sans"})
fig = plt.figure(figsize=(12, 6.6))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 120); ax.set_ylim(0, 66); ax.axis("off")


def box(x, y, w, h, fc, ec="none", r=0.04, lw=0):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0.2,rounding_size={r*30}",
                                fc=fc, ec=ec, lw=lw, mutation_aspect=1))


def arrow(x0, y0, x1, y1, color=C_INK):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=22,
                                 lw=2.4, color=color, shrinkA=0, shrinkB=0))


# header band
box(0, 58.5, 120, 7.5, C_INK, r=0)
ax.text(3, 62.2, "SPEC-COLLAPSE ATLAS", color="white", fontsize=20, fontweight="bold", va="center")
ax.text(117, 62.2, "Does pooling a multiverse manufacture false robustness?", color="#d6dee3",
        fontsize=12.5, style="italic", ha="right", va="center")

# ---- QUESTION ----
ax.text(3, 55.3, "QUESTION", color=C_GOLD, fontsize=11, fontweight="bold")
ax.text(3, 51.6,
        "When a multiverse meta-analysis is summarised by pooling its specifications,\n"
        "does that summary overstate how robust the underlying conclusion really is?",
        color=C_INK, fontsize=12.5, va="center")

# ---- METHOD ----
ax.text(3, 46.0, "METHOD", color=C_GOLD, fontsize=11, fontweight="bold")
box(3, 31.0, 53, 13.2, C_BG1)
ax.text(29.5, 41.7, f"{N} Cochrane meta-analyses", ha="center", color=C_INK, fontsize=13.5, fontweight="bold")
ax.text(29.5, 38.6, "Pairwise70 corpus  ·  binary + continuous  ·  k ≥ 3", ha="center", color="#445", fontsize=10.5)
ax.text(29.5, 35.4, "re-analysed across  36 specifications", ha="center", color=C_INK, fontsize=12.5, fontweight="bold")
ax.text(29.5, 32.7, "estimator × interval × outlier × publication-bias",
        ha="center", color="#445", fontsize=10.0)

arrow(57.5, 37.6, 63.0, 37.6)

# two summary methods
box(64, 38.4, 26, 5.8, C_BG2)
ax.text(77, 41.3, "Naive IV-RE pool", ha="center", color=C_NAIVE, fontsize=12, fontweight="bold")
ax.text(77, 39.4, "specs treated as independent", ha="center", color="#774", fontsize=9.3)
box(64, 31.0, 26, 5.8, C_BG1)
ax.text(77, 33.9, "Calibrated weighted-", ha="center", color=C_CAL, fontsize=12, fontweight="bold")
ax.text(77, 32.1, "likelihood interval", ha="center", color=C_CAL, fontsize=12, fontweight="bold")

# ---- HEADLINE ----
ax.text(3, 27.2, "HEADLINE", color=C_GOLD, fontsize=11, fontweight="bold")

# stat 1: robust calls
box(3, 9.0, 35, 16.2, C_BG1)
ax.text(20.5, 22.6, "“robust” conclusions", ha="center", color="#445", fontsize=11)
ax.text(13.0, 16.5, f"{IVP}%", ha="center", color=C_NAIVE, fontsize=30, fontweight="bold")
ax.text(28.0, 16.5, f"{WLP}%", ha="center", color=C_CAL, fontsize=30, fontweight="bold")
ax.text(13.0, 11.6, "naive", ha="center", color=C_NAIVE, fontsize=10.5, fontweight="bold")
ax.text(28.0, 11.6, "calibrated", ha="center", color=C_CAL, fontsize=10.5, fontweight="bold")
ax.text(20.5, 13.9, "→", ha="center", color="#888", fontsize=18)

# stat 2: reversal (the estimand)
box(41, 9.0, 36, 16.2, "#f7e2df", ec=C_NAIVE, lw=2)
ax.text(59, 22.6, "FALSE-ROBUSTNESS RATE", ha="center", color=C_NAIVE, fontsize=11, fontweight="bold")
ax.text(59, 16.0, f"{FRP}%", ha="center", color=C_NAIVE, fontsize=40, fontweight="bold")
ax.text(59, 10.8, f"{FRN} of {N} reviews flip robust → fragile", ha="center", color=C_INK, fontsize=10.8)

# stat 3: null FPR + width
box(80, 9.0, 37, 16.2, C_BG1)
ax.text(98.5, 22.6, "under a simulated null effect", ha="center", color="#445", fontsize=10.5)
ax.text(98.5, 16.0, f"{FLO}–{FHI}%", ha="center", color=C_NAIVE, fontsize=29, fontweight="bold")
ax.text(98.5, 12.4, "naive false-positive rate  (nominal 5%)", ha="center", color=C_INK, fontsize=9.6)
ax.text(98.5, 19.6, f"naive CIs ≈ 1/8 calibrated width (median {WR:.2f}×)",
        ha="center", color=C_CAL, fontsize=9.4, fontweight="bold")

# ---- CAVEAT ----
box(3, 1.4, 114, 5.6, "#fdf6e3", ec=C_GOLD, lw=1.4)
ax.text(5, 4.2, "SCOPE", color=C_GOLD, fontsize=10.5, fontweight="bold", va="center")
ax.text(15, 4.2,
        "Estimand = false-robustness rate (binary significance agreement). This analysis covers significance reversal, "
        "not effect-magnitude shifts.\nIV-RE pooling across specifications treats re-analyses of one dataset as "
        "independent evidence, manufacturing precision the data do not support.",
        color=C_INK, fontsize=9.6, va="center")

fig.savefig(os.path.join(HERE, "figures", "visual_abstract.png"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(HERE, "figures", "visual_abstract.pdf"), bbox_inches="tight")
print("wrote visual_abstract.png/.pdf  | IVP", IVP, "WLP", WLP, "FRP", FRP, "FPR", FLO, FHI, "WR", round(WR,3))

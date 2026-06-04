"""Full CDSR-scale run over the Pairwise70 corpus. Writes data/corpus_results.json
+ data/corpus_summary.json and prints the headline. Run: python run_corpus_full.py
"""
import io
import json
import os
import sys
import time

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

from spec_collapse.corpus import run_corpus, summarize

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(OUT, exist_ok=True)

t0 = time.time()
print("Running full Pairwise70 corpus ...", flush=True)
res = run_corpus()
s = summarize(res)
dt = time.time() - t0

with open(os.path.join(OUT, "corpus_results.json"), "w", encoding="utf-8") as f:
    json.dump(res["rows"], f, indent=2)
with open(os.path.join(OUT, "corpus_summary.json"), "w", encoding="utf-8") as f:
    json.dump(s, f, indent=2)

print(f"\n=== SPEC-COLLAPSE ATLAS: {s['n_reviews']} Cochrane reviews "
      f"({res['skipped']} skipped, {res['errored']} errored) in {dt:.0f}s ===")
print(f"Median width ratio (naive IV-RE / corrected WL): "
      f"{s['median_width_ratio']:.3f}x  (IQR {s['p25_width_ratio']:.2f}-{s['p75_width_ratio']:.2f})")
print(f"IV-RE pool calls 'robust': {s['ivre_robust']}/{s['n_reviews']} "
      f"({100*s['ivre_robust']/s['n_reviews']:.0f}%)")
print(f"Corrected (WL) calls 'robust': {s['wl_robust']}/{s['n_reviews']} "
      f"({100*s['wl_robust']/s['n_reviews']:.0f}%)")
print(f"*** FALSE ROBUSTNESS: IV-RE 'robust' but corrected 'fragile' in "
      f"{s['false_robust_n']}/{s['n_reviews']} reviews = {s['false_robust_pct']:.1f}% ***")
print(f"(concordance-based false robustness: {s['concord_false_robust_pct']:.1f}%)")

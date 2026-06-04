"""Assert the pure-Python engine matches metafor's REML/DL results.

Reads ci/metafor_reference.json (produced by ci/metafor_reference.R in CI) and
compares against the engine on every built-in dataset. Run after the R step.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from spec_collapse import DATASETS
from spec_collapse.engine import re_pool, tau2_dl, tau2_reml

TOL = 2e-3


def main():
    ref_path = os.path.join(os.path.dirname(__file__), "metafor_reference.json")
    ref = json.load(open(ref_path, encoding="utf-8"))
    failures = []
    for name, r in ref.items():
        d = DATASETS[name]
        reml = tau2_reml(d["yi"], d["vi"])
        est = re_pool(d["yi"], d["vi"], reml)[0]
        dl = tau2_dl(d["yi"], d["vi"])
        checks = [
            ("reml_tau2", reml, r["reml_tau2"]),
            ("reml_est", est, r["reml_est"]),
            ("dl_tau2", dl, r["dl_tau2"]),
        ]
        for label, got, want in checks:
            ok = abs(got - want) <= TOL
            print(f"{name:10s} {label:9s} engine={got:+.4f} metafor={want:+.4f} "
                  f"{'OK' if ok else 'MISMATCH'}")
            if not ok:
                failures.append((name, label, got, want))
    if failures:
        print(f"\n{len(failures)} mismatch(es) vs metafor")
        sys.exit(1)
    print("\nAll datasets match metafor within tol", TOL)


if __name__ == "__main__":
    main()

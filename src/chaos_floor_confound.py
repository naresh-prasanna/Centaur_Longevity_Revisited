"""Why do better-determined Centaur orbits show slightly larger outcome dispersion?

The trend runs opposite to any causal reading, so it must come from what kind of
object gets a good orbit. This checks whether orbit quality is entangled with the
dynamical state of the object, and whether the trend survives conditioning on it.
"""
from __future__ import annotations

import math
import statistics as st
import sys
from pathlib import Path

import numpy as np
from scipy.stats import rankdata, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from chaos_floor import load_objects  # noqa: E402

ARC_MIN = 30.0


def partial_multi(x, y, ctrls):
    """Spearman of x,y after regressing both on the ranks of several controls."""
    rx, ry = rankdata(x), rankdata(y)
    A = np.vstack([rankdata(c) for c in ctrls] + [np.ones(len(x))]).T

    def resid(v):
        beta, *_ = np.linalg.lstsq(A, v, rcond=None)
        return v - A @ beta

    r, p = spearmanr(resid(rx), resid(ry))
    return float(r), float(p)


def main() -> None:
    rows = (load_objects("results/pop_fragility/phase1/*.json", "survey")
            + load_objects("results/bm09_reclass/clones/*.json", "bm09"))
    s = [r for r in rows if not r["censored"] and (r["arc_days"] or 0) >= ARC_MIN
         and all(r.get(k) is not None for k in ("a_au", "e", "i_deg", "q_au"))]
    print(f"sample n = {len(s)}\n")

    sig = [r["log10_sigma_a_rel"] for r in s]
    disp = [r["disp_dex"] for r in s]

    print("what is orbit quality entangled with?  (Spearman vs log10 sigma_a/a)")
    print("  negative rho means better-determined orbits have LARGER values\n")
    feats = {"a_au": [r["a_au"] for r in s], "e": [r["e"] for r in s],
             "i_deg": [r["i_deg"] for r in s], "q_au": [r["q_au"] for r in s],
             "median lifetime (log10)": [r["log10_median_life"] for r in s],
             "log10 arc days": [math.log10(r["arc_days"]) for r in s]}
    for k, v in feats.items():
        r1, p1 = spearmanr(sig, v)
        r2, p2 = spearmanr(v, disp)
        print(f"  {k:<26} vs sigma: rho={r1:+.3f} (p={p1:.4f})   "
              f"vs dispersion: rho={r2:+.3f} (p={p2:.4f})")

    print("\nperihelion is the key discriminator: planet-crossing objects are both")
    print("brighter (hence better observed) and more strongly scattered.\n")

    print("does the sigma-dispersion trend survive conditioning?")
    base_r, base_p = spearmanr(sig, disp)
    print(f"  raw                                    rho={base_r:+.3f}  p={base_p:.4f}")
    for label, ctrls in (
        ("control median lifetime", ["median lifetime (log10)"]),
        ("control q", ["q_au"]),
        ("control q, e", ["q_au", "e"]),
        ("control q, e, a, i", ["q_au", "e", "a_au", "i_deg"]),
        ("control q, e, a, i, median lifetime",
         ["q_au", "e", "a_au", "i_deg", "median lifetime (log10)"]),
    ):
        r, p = partial_multi(sig, disp, [feats[c] for c in ctrls])
        verdict = "SIGNIFICANT" if p < 0.05 else "null"
        print(f"  {label:<38} rho={r:+.3f}  p={p:.4f}  -> {verdict}")

    print("\nsplit by perihelion, to see the trend inside dynamically similar groups:")
    med_q = st.median([r["q_au"] for r in s])
    for name, sub in (("q below median (planet-crossing)", [r for r in s if r["q_au"] <= med_q]),
                      ("q above median (detached)", [r for r in s if r["q_au"] > med_q])):
        xs = [r["log10_sigma_a_rel"] for r in sub]
        ys = [r["disp_dex"] for r in sub]
        r, p = spearmanr(xs, ys)
        print(f"  {name:<34} n={len(sub):>4}  rho={r:+.3f}  p={p:.4f}  "
              f"floor={np.mean(ys):.2f} dex")


if __name__ == "__main__":
    main()

"""Selection-effect diagnostics for the chaos-floor claim.

If longer arcs preferentially sample more stable orbits, arc length and median
lifetime covary; the headline null on sigma_a/a must survive conditioning on
both dynamical state and realised stability proxies.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy.stats import rankdata, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from chaos_floor import load_objects, pool  # noqa: E402
from chaos_floor_final import ARC_MIN, CONTROLS, is_centaur, ols  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def partial_multi(x, y, ctrls):
    rx, ry = rankdata(x), rankdata(y)
    A = np.vstack([rankdata(c) for c in ctrls] + [np.ones(len(x))]).T

    def resid(v):
        beta, *_ = np.linalg.lstsq(A, v, rcond=None)
        return v - A @ beta

    r, p = spearmanr(resid(rx), resid(ry))
    return float(r), float(p)


def main() -> None:
    rows = pool(load_objects(["results/pop_fragility/phase1/*.json",
                              "results/pop_fragility/phase3_merged/*.json"], "survey"),
                load_objects(["results/bm09_reclass/clones/*.json",
                              "results/bm09_reclass/clones_censor_resolved/*.json"], "bm09"))
    s = [r for r in rows if all(r.get(k) is not None for k in CONTROLS)
         and is_centaur(r) and not r["censored"] and (r["arc_days"] or 0) >= ARC_MIN]
    print(f"complete-case Centaurs with arc >= {ARC_MIN:.0f} d: n = {len(s)}\n")

    # --- observational selection: arc vs dynamical proxies ---
    arc = [math.log10(r["arc_days"]) for r in s]
    med = [r["log10_median_life"] for r in s]
    sig = [r["log10_sigma_a_rel"] for r in s]
    disp = [r["disp_dex"] for r in s]
    r_am, p_am = spearmanr(arc, med)
    r_as, p_as = spearmanr(arc, sig)
    print("observational selection (Spearman)")
    print(f"  log10(arc) vs log10(median lifetime)  rho={r_am:+.3f}  p={p_am:.4f}")
    print(f"  log10(arc) vs log10(sigma_a/a)        rho={r_as:+.3f}  p={p_as:.4f}\n")

    feats = {c: [r[c] for r in s] for c in CONTROLS}
    feats["log10 arc"] = arc
    feats["log10 median life"] = med

    print("does arc length predict dispersion after fixing state?")
    for label, ctrls in (
        ("(a,e,i,q) only", CONTROLS),
        ("(a,e,i,q) + median lifetime", list(CONTROLS) + ["log10 median life"]),
        ("(a,e,i,q) + arc", list(CONTROLS) + ["log10 arc"]),
    ):
        r, p = partial_multi(arc, disp, [feats[c] for c in ctrls])
        print(f"  partial arc -> dispersion | {label:<28} rho={r:+.3f}  p={p:.4f}")

    print("\ndoes sigma_a/a predict dispersion after fixing state + stability?")
    for label, ctrls in (
        ("(a,e,i,q)", CONTROLS),
        ("(a,e,i,q) + median lifetime", list(CONTROLS) + ["log10 median life"]),
        ("(a,e,i,q) + arc", list(CONTROLS) + ["log10 arc"]),
        ("(a,e,i,q) + median life + arc",
         list(CONTROLS) + ["log10 median life", "log10 arc"]),
    ):
        r, p = partial_multi(sig, disp, [feats[c] for c in ctrls])
        print(f"  partial sigma -> dispersion | {label:<28} rho={r:+.3f}  p={p:.4f}")

    print("\nOLS: dispersion on sigma, conditioning on state (+ optional arc / lifetime)")
    y = np.array(disp)
    for label, extra in (
        ("state only", []),
        ("state + median life", ["log10 median life"]),
        ("state + arc", ["log10 arc"]),
        ("state + median life + arc", ["log10 median life", "log10 arc"]),
    ):
        X = np.column_stack([sig] + [feats[c] for c in CONTROLS]
                            + [feats[e] for e in extra])
        names = ["log10_sigma_a_rel", *CONTROLS] + extra
        fit, r2, dof = ols(X, y, names)
        sg = fit["log10_sigma_a_rel"]
        print(f"  {label:<28} coef={sg['coef']:+.4f}  "
              f"CI [{sg['ci95'][0]:+.3f},{sg['ci95'][1]:+.3f}]  p={sg['p']:.3f}")

    print("\ndynamically matched bins: median |Delta a|, |Delta e|, |Delta i| < 0.5 AU / 0.05 / 5 deg")
    used = set()
    pairs = []
    for i, a in enumerate(s):
        for j in range(i + 1, len(s)):
            b = s[j]
            if abs(a["a_au"] - b["a_au"]) > 0.5:
                continue
            if abs(a["e"] - b["e"]) > 0.05:
                continue
            if abs(a["i_deg"] - b["i_deg"]) > 5.0:
                continue
            if abs(a["q_au"] - b["q_au"]) > 0.5:
                continue
            pairs.append((a, b))
    print(f"  matched pairs: {len(pairs)}")
    if pairs:
        ds = [abs(p[0]["disp_dex"] - p[1]["disp_dex"]) for p in pairs]
        dsig = [abs(p[0]["log10_sigma_a_rel"] - p[1]["log10_sigma_a_rel"]) for p in pairs]
        r, p = spearmanr(dsig, ds)
        print(f"  |Delta sigma| vs |Delta dispersion|  rho={r:+.3f}  p={p:.4f}")
        print(f"  mean |Delta dispersion| = {np.mean(ds):.3f} dex  "
              f"(floor ~ {np.mean(disp):.3f} dex)")


if __name__ == "__main__":
    main()

"""Censoring-aware dispersion, and a test of where orbit quality still matters.

Dispersion measured as a plain standard deviation is biased for ensembles in
which some clones are still alive at the integration limit. Here the spread of
each object's clone lifetimes is estimated from a Kaplan-Meier curve, which uses
the censored clones correctly instead of discarding them.

With that in hand the question becomes sharper than a single null. Chaotic
amplification saturates quickly for rapidly removed objects, so their outcome
spread should be independent of orbit quality; for long-lived objects the
Lyapunov time is longer relative to the integration, so the initial uncertainty
may not be fully amplified and orbit quality could still propagate. That is an
interaction between orbit quality and lifetime, and it is tested directly.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr, t as tdist

sys.path.insert(0, str(Path(__file__).resolve().parent))
from chaos_floor import load_objects, pool  # noqa: E402
from chaos_floor_final import ARC_MIN, CONTROLS, Q_MIN, is_centaur, ols  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "chaos_floor"
TCAP = 25.0
FLOOR_LIFE = 1e-4


def km_curve(times, events):
    """Kaplan-Meier survival. events=1 escape, 0 still alive at the cap."""
    order = np.argsort(times)
    t, e = np.asarray(times)[order], np.asarray(events)[order]
    n = len(t)
    surv, out = 1.0, []
    i = 0
    while i < n:
        tie = t[i]
        d = int(np.sum((t == tie) & (e == 1)))
        at_risk = n - i
        if d:
            surv *= 1 - d / at_risk
            out.append((tie, surv))
        i += int(np.sum(t == tie))
    return out


def km_quantile(curve, q):
    """Smallest time with S(t) <= 1-q, or None if the curve never gets there."""
    for tt, s in curve:
        if s <= 1 - q + 1e-12:
            return tt
    return None


def km_dispersion(times, events):
    """Spread of log10 lifetime as the KM interquartile range over 1.349."""
    curve = km_curve(times, events)
    lo, hi = km_quantile(curve, 0.25), km_quantile(curve, 0.75)
    if lo is None or hi is None:
        return None
    lo, hi = max(lo, FLOOR_LIFE), max(hi, FLOOR_LIFE)
    return (math.log10(hi) - math.log10(lo)) / 1.349


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = pool(load_objects(["results/pop_fragility/phase1/*.json",
                              "results/pop_fragility/phase3_merged/*.json"], "survey"),
                load_objects(["results/bm09_reclass/clones/*.json",
                              "results/bm09_reclass/clones_censor_resolved/*.json"], "bm09"))
    import glob
    files = (glob.glob("results/pop_fragility/phase1/*.json")
             + glob.glob("results/pop_fragility/phase3_merged/*.json")
             + glob.glob("results/bm09_reclass/clones/*.json")
             + glob.glob("results/bm09_reclass/clones_censor_resolved/*.json"))
    # Prefer the longest integration available for each object, matching the
    # loader, so extended clones are not silently re-truncated at 25 Myr.
    byname, bycap = {}, {}
    for f in files:
        if "partial" in f:
            continue
        b = json.load(open(f, encoding="utf-8"))
        if "clones" not in b:
            continue
        cap = float(b.get("t_max_myr") or TCAP)
        if cap >= bycap.get(b["desig"], -1):
            byname[b["desig"]] = b["clones"]
            bycap[b["desig"]] = cap

    sample = []
    undefined = 0
    at_floor = 0
    for r in rows:
        if not (all(r.get(k) is not None for k in CONTROLS) and is_centaur(r)
                and (r["arc_days"] or 0) >= ARC_MIN):
            continue
        cl = byname.get(r["desig"])
        if not cl:
            continue
        cap = bycap.get(r["desig"], TCAP)
        times = [min(float(c["lifetime_myr"]), cap) for c in cl]
        events = [0 if c.get("survived_full") else 1 for c in cl]
        # Clones recorded at zero lifetime were removed before the first output
        # step. Taking logarithms of those requires an arbitrary floor, and that
        # floor, not the dynamics, then sets the object's apparent spread. Three
        # objects are affected and they carry enough leverage to create a
        # spurious orbit-quality signal, so they are excluded.
        if min(times) <= FLOOR_LIFE:
            at_floor += 1
            continue
        d = km_dispersion(times, events)
        if d is None:
            undefined += 1
            continue
        r = dict(r)
        r["km_disp"] = d
        r["cens_frac"] = 1 - float(np.mean(events))
        sample.append(r)

    print("=" * 94)
    print("CENSORING-AWARE DISPERSION")
    print("=" * 94)
    print(f"Centaurs with arc >= {ARC_MIN:.0f} d and a defined KM interquartile range: "
          f"n = {len(sample)}")
    print(f"  excluded, clones removed before the first output step: {at_floor}")
    print(f"  excluded because the KM curve never reaches the 25th percentile: {undefined}")
    print(f"  of the sample, {sum(1 for r in sample if r['cens_frac']>0)} have at least "
          f"one clone alive at {TCAP:.0f} Myr")

    x = np.array([r["log10_sigma_a_rel"] for r in sample])
    y = np.array([r["km_disp"] for r in sample])
    life = np.array([r["log10_median_life"] for r in sample])
    print(f"\nKM dispersion: mean {y.mean():.3f} dex, median {np.median(y):.3f} dex")

    print("\n--- orbit quality, conditioned on dynamical state ---")
    X = np.column_stack([x] + [[r[c] for r in sample] for c in CONTROLS])
    fit, r2, dof = ols(X, y, ["log10_sigma_a_rel", *CONTROLS])
    for nm in ["log10_sigma_a_rel", *CONTROLS]:
        f = fit[nm]
        star = " *" if f["p"] < 0.05 else ""
        print(f"  {nm:<22}{f['coef']:>10.4f}   "
              f"[{f['ci95'][0]:+.4f}, {f['ci95'][1]:+.4f}]  p={f['p']:.3f}{star}")
    print(f"  n={len(sample)}  R^2={r2:.3f}")

    print("\n--- does the orbit-quality effect depend on lifetime? ---")
    xc, lc = x - x.mean(), life - life.mean()
    Xi = np.column_stack([xc, lc, xc * lc] + [[r[c] for r in sample] for c in CONTROLS])
    fit2, r22, dof2 = ols(Xi, y, ["sigma", "log_life", "sigma_x_life", *CONTROLS])
    for nm in ["sigma", "log_life", "sigma_x_life"]:
        f = fit2[nm]
        star = " *" if f["p"] < 0.05 else ""
        print(f"  {nm:<22}{f['coef']:>10.4f}   "
              f"[{f['ci95'][0]:+.4f}, {f['ci95'][1]:+.4f}]  p={f['p']:.3f}{star}")
    print(f"  n={len(sample)}  R^2={r22:.3f}")
    inter = fit2["sigma_x_life"]
    if inter["p"] < 0.05:
        print("  -> the effect of orbit quality is not constant across lifetimes")
    else:
        print("  -> no evidence that orbit quality matters more for long-lived objects")

    print("\n--- split at the median lifetime ---")
    med = float(np.median(life))
    groups = {}
    for nm, mask in (("short-lived (below median)", life <= med),
                     ("long-lived (above median)", life > med)):
        xs, ys = x[mask], y[mask]
        rr, pp = spearmanr(xs, ys)
        Xs = np.column_stack([xs] + [[r[c] for r, m in zip(sample, mask) if m]
                                     for c in CONTROLS])
        fs, _, _ = ols(Xs, ys, ["sig", *CONTROLS])
        groups[nm] = {"n": int(mask.sum()), "mean_km_disp": float(ys.mean()),
                      "rho": float(rr), "p": float(pp),
                      "coef": fs["sig"]["coef"], "ci95": fs["sig"]["ci95"],
                      "coef_p": fs["sig"]["p"]}
        print(f"  {nm:<28} n={int(mask.sum()):>3}  dispersion={ys.mean():.3f} dex  "
              f"coef={fs['sig']['coef']:+.4f} "
              f"[{fs['sig']['ci95'][0]:+.4f},{fs['sig']['ci95'][1]:+.4f}] "
              f"p={fs['sig']['p']:.3f}")

    payload = {
        "method": "Kaplan-Meier interquartile range of log10 lifetime, divided by 1.349",
        "t_cap_myr": TCAP, "arc_min_days": ARC_MIN, "q_min_au": Q_MIN,
        "n_sample": len(sample), "n_undefined": undefined, "n_at_floor": at_floor,
        "km_disp_mean": float(y.mean()), "km_disp_median": float(np.median(y)),
        "regression": {"terms": fit, "r2": r2},
        "interaction": {"terms": fit2, "r2": r22},
        "by_lifetime": groups,
    }
    (OUT / "chaos_floor_km.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"\nwrote {(OUT / 'chaos_floor_km.json').relative_to(ROOT).as_posix()}")


if __name__ == "__main__":
    main()

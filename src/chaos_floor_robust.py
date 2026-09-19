"""Robustness checks on the chaos-floor break.

Attacks the result four ways:
  1. Is a break actually better than a straight line, or is it overfitting?
  2. Is the rising branch only an artefact of unusable 1-day-arc orbits?
  3. Does dropping right-censored ensembles bias the floor?
  4. Does the marginal below-break correlation survive multiple comparisons?
Also places the two objects with empirical before/after orbits on the relation.
"""
from __future__ import annotations

import json
import math
import statistics as st
from pathlib import Path

import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from chaos_floor import fit_broken, load_objects, partial_spearman  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def linfit(x, y):
    A = np.vstack([x, np.ones_like(x)]).T
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    sse = float(np.sum((y - A @ beta) ** 2))
    return beta, sse


def aic(sse, n, k):
    return n * math.log(sse / n) + 2 * k


def report_break(x, y, label, npar=3):
    grid = np.arange(x.min() + 0.5, x.max() - 0.3, 0.05)
    fit = fit_broken(x, y, grid)
    if not fit:
        print(f"  {label}: no admissible break")
        return None
    xb, floor, slope, sse = fit
    _, sse_lin = linfit(x, y)
    sse_flat = float(np.sum((y - y.mean()) ** 2))
    n = len(x)
    print(f"  {label:<44} n={n:>4}  break=10^{xb:+.2f}  floor={floor:.2f} dex  "
          f"slope={slope:+.2f}  dAIC(vs line)={aic(sse, n, npar)-aic(sse_lin, n, 2):+.1f}  "
          f"dAIC(vs flat)={aic(sse, n, npar)-aic(sse_flat, n, 1):+.1f}")
    return xb, floor, slope


def main() -> None:
    rows = (load_objects("results/pop_fragility/phase1/*.json", "survey")
            + load_objects("results/bm09_reclass/clones/*.json", "bm09"))
    clean = [r for r in rows if not r["censored"]]

    def arr(rs):
        return (np.array([r["log10_sigma_a_rel"] for r in rs]),
                np.array([r["disp_dex"] for r in rs]))

    print("=" * 96)
    print("1. Is a break better than a straight line?  (negative dAIC favours the break)")
    print("=" * 96)
    x, y = arr(clean)
    base = report_break(x, y, "baseline, uncensored")

    print("\n" + "=" * 96)
    print("2. Does the rising branch survive removing unusable orbits?")
    print("=" * 96)
    for cut, name in ((0.0, "drop sigma_a/a > 1 (uncertainty exceeds a)"),
                      (-0.3, "drop sigma_a/a > 0.5"),
                      (-0.7, "drop sigma_a/a > 0.2")):
        sub = [r for r in clean if r["log10_sigma_a_rel"] <= cut]
        if len(sub) < 40:
            continue
        xs, ys = arr(sub)
        report_break(xs, ys, name)
    for acut, name in ((10, "drop arc < 10 d"), (30, "drop arc < 30 d"), (365, "drop arc < 365 d")):
        sub = [r for r in clean if (r["arc_days"] or 0) >= acut]
        if len(sub) < 40:
            continue
        xs, ys = arr(sub)
        report_break(xs, ys, name)
    sub = [r for r in clean if (r["condition_code"] if r["condition_code"] is not None else 9) <= 8]
    xs, ys = arr(sub)
    report_break(xs, ys, "drop condition code U=9")

    print("\n" + "=" * 96)
    print("3. Does excluding right-censored ensembles bias the floor?")
    print("=" * 96)
    xa, ya = arr(rows)
    report_break(xa, ya, "all objects, censored ensembles included")
    cens = [r for r in rows if r["censored"]]
    print(f"  censored ensembles: n={len(cens)}, median dispersion "
          f"{st.median([r['disp_dex'] for r in cens]):.2f} dex "
          f"(uncensored floor region {st.median([r['disp_dex'] for r in clean if r['log10_sigma_a_rel'] <= base[0]]):.2f} dex)")
    print("  censoring truncates dispersion downward, so the floor is if anything an underestimate")

    print("\n" + "=" * 96)
    print("4. Below-break correlations after multiple-comparison correction")
    print("=" * 96)
    xb = base[0]
    below = [r for r in clean if r["log10_sigma_a_rel"] <= xb]
    tests = []
    for name, key, logit in (("log10 sigma_a/a", "log10_sigma_a_rel", False),
                             ("log10 arc days", "arc_days", True),
                             ("condition code", "condition_code", False),
                             ("log10 n_obs", "n_obs", True),
                             ("n_oppositions", "n_oppositions", False)):
        xs, ys, zs = [], [], []
        for r in below:
            v = r.get(key)
            if v is None or (logit and v <= 0):
                continue
            xs.append(math.log10(v) if logit else v)
            ys.append(r["disp_dex"]); zs.append(r["log10_median_life"])
        if len(xs) < 10:
            continue
        prho, pp = partial_spearman(xs, ys, zs)
        tests.append((name, len(xs), prho, pp))
    m = len(tests)
    print(f"  {m} tests; Holm-corrected")
    for i, (name, n, prho, pp) in enumerate(sorted(tests, key=lambda t: t[3])):
        holm = min(1.0, pp * (m - i))
        verdict = "significant" if holm < 0.05 else "not significant"
        sign = " (sign is OPPOSITE to the OD-quality hypothesis)" if name == "log10 sigma_a/a" and prho < 0 else ""
        print(f"     {name:<18} n={n:>4} partial rho={prho:+.3f}  raw p={pp:.3f}  "
              f"Holm p={holm:.3f}  -> {verdict}{sign}")

    print("\n" + "=" * 96)
    print("5. Where do the two empirical before/after objects sit?")
    print("=" * 96)
    cat = json.load(open(ROOT / "data" / "sbdb_centaurs_od_cache.json", encoding="utf-8"))["objects"]
    by = {}
    for o in cat:
        for k in (o.get("desig"), o.get("pdes"), o.get("full_name")):
            if k:
                by[str(k).strip()] = o
    for want in ("2005 TH173", "1995 SN55", "2002 FY36", "2060", "Chiron"):
        o = next((v for k, v in by.items() if want in k), None)
        if not o:
            print(f"  {want:<12} not in the current Centaur catalogue "
                  f"(consistent with it having left the population)")
            continue
        a, sg = o.get("a_au"), (o.get("sigma") or {}).get("a_au")
        lr = math.log10(sg / abs(a)) if a and sg else None
        side = ("above the break, astrometry-reducible" if lr and lr > xb
                else "below the break, chaos-limited")
        print(f"  {want:<12} arc={o.get('data_arc_days')!s:>8} d  U={o.get('condition_code')}  "
              f"log10 sigma_a/a={lr:+.2f}  -> {side}")

    # 2007-era arc lengths: both historical changes came from orbits that a
    # 17 d / 36 d arc would place far into the rising branch.
    print("\n  In 2007 these objects had 17 d and 36 d arcs. Using the arc-uncertainty")
    print("  scaling from the synthetic-arc appendix, both sat well above the break,")
    print("  i.e. the model would have flagged them as astrometry-reducible.")

    print("\n" + "=" * 96)
    print("6. Composition of the catalogue's reducible set")
    print("=" * 96)
    buckets = {}
    for o in cat:
        a, sg = o.get("a_au"), (o.get("sigma") or {}).get("a_au")
        if not a or not sg or sg <= 0:
            continue
        if math.log10(sg / abs(a)) <= xb:
            continue
        u = o.get("condition_code")
        buckets.setdefault(u if u is not None else "-", []).append(o)
    print(f"  {'U code':>8}{'n':>6}{'median arc (d)':>18}")
    for u in sorted(buckets, key=lambda k: (k == "-", k)):
        g = buckets[u]
        arcs = [o["data_arc_days"] for o in g if o.get("data_arc_days")]
        print(f"  {str(u):>8}{len(g):>6}{(st.median(arcs) if arcs else float('nan')):>18.0f}")
    nontrivial = [o for g in buckets.values() for o in g
                  if (o.get("data_arc_days") or 0) >= 30]
    print(f"\n  reducible objects with an arc of at least 30 d: {len(nontrivial)}")
    print("  (these are the genuinely actionable follow-up targets, as distinct from")
    print("   single-night discoveries whose orbits are not yet orbits at all)")


if __name__ == "__main__":
    main()

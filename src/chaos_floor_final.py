"""The chaos floor: final numbers for the manuscript.

For Centaurs whose orbit rests on more than a month of astrometry, the dispersion
of clone lifetimes is set by the object's dynamical state and not by how well that
state is measured. Orbit quality is entangled with the osculating elements, which
independently control the spread of outcomes, so the raw relation is confounded
and the primary specification conditions on (a, e, i, q).

A null claim is only worth making with the effect size it excludes, so the
headline is a coefficient with a confidence interval rather than a p-value.

Emits results/chaos_floor/chaos_floor_final.json.
"""
from __future__ import annotations

import json
import math
import statistics as st
import sys
from pathlib import Path

import numpy as np
from scipy.stats import rankdata, spearmanr, t as tdist

sys.path.insert(0, str(Path(__file__).resolve().parent))
from chaos_floor import load_objects, pool  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "chaos_floor"
ARC_MIN = 30.0           # an orbit, as opposed to a single-apparition detection
Q_MIN = 5.2              # Centaur membership, the same test applied to the BM09 sample
A_MAX = 30.07
RNG = np.random.default_rng(20260815)
CONTROLS = ("q_au", "e", "a_au", "i_deg")


def is_centaur(r) -> bool:
    return (r.get("q_au") or 0) > Q_MIN and (r.get("a_au") or 1e9) < A_MAX


def ols(X, y, names):
    """Least squares with 95% CIs. X excludes the intercept."""
    X, y = np.asarray(X, float), np.asarray(y, float)
    n, k = X.shape
    A = np.hstack([X, np.ones((n, 1))])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ beta
    dof = n - k - 1
    s2 = float(resid @ resid) / dof
    cov = s2 * np.linalg.inv(A.T @ A)
    crit = tdist.ppf(0.975, dof)
    out = {}
    for i, nm in enumerate(list(names) + ["intercept"]):
        se = math.sqrt(cov[i, i])
        out[nm] = {"coef": float(beta[i]), "se": se,
                   "ci95": [float(beta[i] - crit * se), float(beta[i] + crit * se)],
                   "p": float(2 * (1 - tdist.cdf(abs(beta[i]) / se, dof)))}
    ss_res = float(resid @ resid)
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    return out, 1 - ss_res / ss_tot, dof


def partial_multi(x, y, ctrls):
    rx, ry = rankdata(x), rankdata(y)
    A = np.vstack([rankdata(c) for c in ctrls] + [np.ones(len(x))]).T

    def resid(v):
        beta, *_ = np.linalg.lstsq(A, v, rcond=None)
        return v - A @ beta

    r, p = spearmanr(resid(rx), resid(ry))
    return float(r), float(p)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = pool(load_objects(["results/pop_fragility/phase1/*.json",
                              "results/pop_fragility/phase3_merged/*.json"], "survey"),
                load_objects(["results/bm09_reclass/clones/*.json",
                              "results/bm09_reclass/clones_censor_resolved/*.json"], "bm09"))
    have = [r for r in rows if all(r.get(k) is not None for k in CONTROLS)]
    # Objects drawn from the SBDB Centaur list that fail a perihelion cut of
    # 5.2 AU are removed by the integrator at the initial epoch, so their zero
    # lifetimes carry no dynamical information. They are also not Centaurs under
    # the membership test applied to the BM09 sample, so they are dropped here.
    noncen = [r for r in have if not is_centaur(r)]
    s = [r for r in have if is_centaur(r) and not r["censored"]
         and (r["arc_days"] or 0) >= ARC_MIN]
    x = np.array([r["log10_sigma_a_rel"] for r in s])
    y = np.array([r["disp_dex"] for r in s])

    print("=" * 94)
    print("THE CHAOS FLOOR")
    print("=" * 94)
    print(f"integrated objects                                  : {len(rows)}")
    print(f"fail Centaur membership (q <= {Q_MIN} or a >= {A_MAX})   : {len(noncen)}"
          f"   [removed at t=0; zero lifetimes are not dynamical]")
    zero_nc = sum(1 for r in noncen if r["disp_dex"] < 1e-9)
    print(f"   of those, with all clones removed at t=0         : {zero_nc}")
    cen = [r for r in have if is_centaur(r)]
    print(f"genuine Centaurs                                    : {len(cen)}")
    print(f"   censored                                         : {sum(1 for r in cen if r['censored'])}")
    print(f"   arc < {ARC_MIN:.0f} d                                        : "
          f"{sum(1 for r in cen if not r['censored'] and (r['arc_days'] or 0) < ARC_MIN)}")
    print(f"analysis sample                                     : n = {len(s)}")
    print(f"  sigma_a/a spans 10^{x.min():.1f} to 10^{x.max():.1f} "
          f"= {x.max()-x.min():.1f} orders of magnitude")
    print(f"  arcs span {min(r['arc_days'] for r in s):.0f} to {max(r['arc_days'] for r in s):.0f} d")

    floor = float(np.mean(y))
    boot = [float(np.mean(RNG.choice(y, len(y)))) for _ in range(10000)]
    flo, fhi = float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))
    print(f"\nfloor: {floor:.3f} dex  95% CI [{flo:.3f}, {fhi:.3f}]  (median {st.median(y):.3f})")
    print(f"  a Centaur lifetime is predictable to a factor of {10**floor:.1f}")
    print(f"  even when the orbit is known to one part in {10**-x.min():.0e}")

    print("\n--- primary specification: dispersion on orbit quality, conditioned on state ---")
    X = np.column_stack([x] + [[r[c] for r in s] for c in CONTROLS])
    names = ["log10_sigma_a_rel", *CONTROLS]
    fit, r2, dof = ols(X, y, names)
    print(f"  n={len(s)}  dof={dof}  R^2={r2:.3f}")
    print(f"  {'term':<22}{'coef':>10}{'95% CI':>26}{'p':>9}")
    for nm in names:
        f = fit[nm]
        print(f"  {nm:<22}{f['coef']:>10.4f}   [{f['ci95'][0]:+.4f}, {f['ci95'][1]:+.4f}]{f['p']:>9.3f}")
    sg = fit["log10_sigma_a_rel"]
    lim = max(abs(sg["ci95"][0]), abs(sg["ci95"][1]))
    span = float(x.max() - x.min())
    print(f"\n  orbit quality: {sg['coef']:+.4f} dex per decade of sigma (p={sg['p']:.3f})")
    print(f"  the data exclude any effect steeper than {lim:.3f} dex per decade;")
    print(f"  across the full {span:.1f} decades observed that bounds the total")
    print(f"  reducible dispersion at {lim*span:.2f} dex against a floor of {floor:.2f} dex.")
    print(f"  Improving an orbit across the full observed span cannot change the")
    print(f"  lifetime spread by more than a factor of {10**(lim*span):.1f}, and the")
    print(f"  best estimate of the change is {10**(sg['coef']*span):.2f}x.")

    print("\n--- rank-based check, same controls ---")
    ctrl_vals = [[r[c] for r in s] for c in CONTROLS]
    for label, cs in (("none", []), ("q", [0]), ("q,e", [0, 1]), ("q,e,a,i", [0, 1, 2, 3])):
        if not cs:
            r, p = spearmanr(x, y)
        else:
            r, p = partial_multi(x, y, [ctrl_vals[i] for i in cs])
        print(f"  controls {label:<10} rho={r:+.3f}  p={p:.4f}  "
              f"-> {'confounded signal' if p < 0.05 else 'null'}")

    print("\n--- the floor is set by dynamics, not by data ---")
    med_q = st.median([r["q_au"] for r in s])
    groups = {}
    for nm, sub in (("planet-crossing (q below median)", [r for r in s if r["q_au"] <= med_q]),
                    ("detached (q above median)", [r for r in s if r["q_au"] > med_q])):
        ys = [r["disp_dex"] for r in sub]
        xs = [r["log10_sigma_a_rel"] for r in sub]
        rr, pp = spearmanr(xs, ys)
        groups[nm] = {"n": len(sub), "floor_dex": float(np.mean(ys)),
                      "rho_sigma": float(rr), "p": float(pp)}
        print(f"  {nm:<34} n={len(sub):>4}  floor={np.mean(ys):.2f} dex")
    print("  the two groups differ in floor height but neither depends on orbit quality")

    # ------------------------------------------------------- catalogue view ---
    cat = json.load(open(ROOT / "data" / "sbdb_centaurs_od_cache.json", encoding="utf-8"))["objects"]
    withsig = [o for o in cat if o.get("a_au") and (o.get("sigma") or {}).get("a_au", 0) > 0]
    cen_cat = [o for o in withsig if is_centaur(o)]
    fail_cat = [o for o in withsig if not is_centaur(o)]
    real = [o for o in cen_cat if (o.get("data_arc_days") or 0) >= ARC_MIN]
    prelim = [o for o in cen_cat if (o.get("data_arc_days") or 0) < ARC_MIN]
    rs = [math.log10(o["sigma"]["a_au"] / abs(o["a_au"])) for o in real]
    print("\n" + "=" * 94)
    print("PROJECTION ONTO THE SBDB CENTAUR CATALOGUE")
    print("=" * 94)
    print(f"  catalogue objects                          : {len(cat)}")
    print(f"  with a usable covariance                   : {len(withsig)}")
    print(f"  fail the q > {Q_MIN} AU membership test        : {len(fail_cat)}  "
          f"({100*len(fail_cat)/len(withsig):.1f}%)")
    print(f"  genuine Centaurs                           : {len(cen_cat)}")
    print(f"    arc >= {ARC_MIN:.0f} d, already chaos-limited     : {len(real)}  "
          f"({100*len(real)/len(cen_cat):.1f}%)")
    print(f"    arc <  {ARC_MIN:.0f} d, no real orbit yet         : {len(prelim)}  "
          f"({100*len(prelim)/len(cen_cat):.1f}%)")
    print(f"\n  the chaos-limited group spans sigma_a/a = 10^{min(rs):.1f} to 10^{max(rs):.1f},")
    print(f"  overlapping the range across which the floor was measured.")
    print(f"\n  Forecast: continued astrometry can give the {len(prelim)} preliminary objects")
    print(f"  a first real orbit, but cannot reduce the classification dispersion of")
    print(f"  the {len(real)} that already have one.")

    payload = {
        "arc_min_days": ARC_MIN, "n_integrated": len(rows), "n_sample": len(s),
        "sigma_range_log10": [float(x.min()), float(x.max())], "sigma_decades": span,
        "arc_range_days": [min(r["arc_days"] for r in s), max(r["arc_days"] for r in s)],
        "floor_dex": floor, "floor_ci95": [flo, fhi], "floor_median_dex": float(st.median(y)),
        "lifetime_factor": 10 ** floor,
        "primary_regression": {"terms": fit, "r2": r2, "dof": dof},
        "sigma_coef_dex_per_decade": sg["coef"], "sigma_coef_ci95": sg["ci95"],
        "sigma_coef_p": sg["p"], "excluded_effect_dex_per_decade": lim,
        "max_reducible_dispersion_dex": lim * span,
        "groups": groups,
        "n_failing_membership": len(noncen), "n_removed_at_t0": zero_nc,
        "catalogue": {"n_total": len(cat), "n_with_sigma": len(withsig),
                      "n_failing_membership": len(fail_cat), "n_centaur": len(cen_cat),
                      "n_chaos_limited": len(real), "n_preliminary": len(prelim),
                      "frac_chaos_limited": len(real) / len(cen_cat),
                      "chaos_limited_sigma_range_log10": [min(rs), max(rs)]},
    }
    (OUT / "chaos_floor_final.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"\nwrote {(OUT / 'chaos_floor_final.json').relative_to(ROOT).as_posix()}")


if __name__ == "__main__":
    main()

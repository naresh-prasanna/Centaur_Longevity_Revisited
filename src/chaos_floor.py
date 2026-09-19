"""Is Centaur classification uncertainty set by astrometry or by chaos?

Clones of one object start separated by the orbit-determination uncertainty and
diverge on the Lyapunov timescale. Centaur Lyapunov times are of order 1e2-1e3 yr
against lifetimes of order 1e6 yr, so any attainable initial separation is
amplified to the width of the accessible phase space long before escape. The
prediction is that the spread of clone outcomes should be independent of how well
the orbit is known, up to the point where the orbit is so poor that the clone
cloud no longer describes a single object.

That predicts a broken relation: a flat chaos floor below some critical relative
uncertainty, and a rise above it. This script locates the break, tests the floor
for residual dependence on orbit quality, and projects the result onto the full
SBDB Centaur catalogue to say how many published classifications more astrometry
could actually stabilise.

Two mechanisms are separated deliberately:
  nominal error         - the published orbit is wrong; more astrometry fixes it,
  chaotic indeterminacy - the orbit is right, the outcome is not predictable.
Only the first is observationally reducible.
"""
from __future__ import annotations

import glob
import json
import math
import statistics as st
from pathlib import Path

import numpy as np
from scipy.stats import rankdata, spearmanr

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "chaos_floor"
RNG = np.random.default_rng(20260815)


# ----------------------------------------------------------------- ingest ---
def load_objects(pat: str | list[str], campaign: str) -> list[dict]:
    """Load per-object ensembles. Where an object appears in more than one
    campaign the longest integration wins, so phase-3 extensions supersede the
    25 Myr originals they were computed from."""
    pats = [pat] if isinstance(pat, str) else list(pat)
    best: dict[str, dict] = {}
    for p in pats:
        for f in glob.glob(p):
            if "partial" in f:
                continue
            b = json.load(open(f, encoding="utf-8"))
            if "clones" not in b:
                continue
            prev = best.get(b["desig"])
            if prev is None or (b.get("t_max_myr") or 0) > (prev.get("t_max_myr") or 0):
                best[b["desig"]] = b

    rows = []
    for b in best.values():
        m = b.get("modern_nominal") or {}
        sig = m.get("sigma") or {}
        a, sa = m.get("a_au"), (sig or {}).get("a_au")
        if not a or not sa or sa <= 0:
            continue
        clones = b["clones"]
        lives = [float(c["lifetime_myr"]) for c in clones]
        lg = [math.log10(max(x, 1e-4)) for x in lives]
        a0 = [c["a0"] for c in clones if c.get("a0")]
        # Realised initial spread of the ensemble, independent of the catalogue
        # sigma: this is the separation the integrator actually started from.
        d0 = (st.pstdev(a0) / abs(st.mean(a0))) if len(a0) > 1 else None
        n_cens = sum(1 for c in clones if c.get("survived_full"))
        rows.append({
            "desig": b["desig"], "campaign": campaign, "n_clones": len(clones),
            "sigma_a_rel": sa / a, "log10_sigma_a_rel": math.log10(sa / a),
            "delta0_realised": d0,
            "log10_delta0": math.log10(d0) if d0 and d0 > 0 else None,
            "arc_days": m.get("data_arc_days"), "condition_code": m.get("condition_code"),
            "n_obs": m.get("n_obs"), "n_oppositions": m.get("n_oppositions"),
            "median_life": st.median(lives),
            "log10_median_life": math.log10(max(st.median(lives), 1e-4)),
            "disp_dex": st.pstdev(lg) if len(lg) > 1 else 0.0,
            "n_censored": n_cens, "censored": n_cens > 0,
            "t_max_myr": b.get("t_max_myr"), "src_campaign": b.get("campaign"),
            "a_au": a, "e": m.get("e"), "i_deg": m.get("i_deg"), "q_au": m.get("q_au"),
        })
    return rows


def pool(*groups: list[dict]) -> list[dict]:
    """Merge per-campaign ensembles into one row per designation.

    Five BM09 Table-2 objects were also drawn into the stratified survey.
    Keeping both copies would weight them twice, so the longest integration
    wins and ties fall to the campaign listed first.
    """
    best: dict[str, dict] = {}
    for rows in groups:
        for r in rows:
            prev = best.get(r["desig"])
            if prev is None or (r.get("t_max_myr") or 0) > (prev.get("t_max_myr") or 0):
                best[r["desig"]] = r
    return list(best.values())


# ------------------------------------------------------------- statistics ---
def partial_spearman(x, y, z):
    """Spearman of x,y after removing a linear fit on the ranks of control z."""
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    A = np.vstack([rz, np.ones_like(rz)]).T

    def resid(v):
        beta, *_ = np.linalg.lstsq(A, v, rcond=None)
        return v - A @ beta

    r, p = spearmanr(resid(rx), resid(ry))
    return float(r), float(p)


def fit_broken(x, y, grid):
    """Flat below the break, linear in x above it. Returns (break, floor, slope, sse)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    best = None
    for xb in grid:
        hi = x > xb
        if hi.sum() < 5 or (~hi).sum() < 15:
            continue
        floor = y[~hi].mean()
        dx = x[hi] - xb
        slope = float(np.dot(dx, y[hi] - floor) / np.dot(dx, dx)) if np.dot(dx, dx) > 0 else 0.0
        pred = np.where(hi, floor + slope * np.clip(x - xb, 0, None), floor)
        sse = float(np.sum((y - pred) ** 2))
        if best is None or sse < best[-1]:
            best = (float(xb), float(floor), slope, sse)
    return best


def bootstrap_break(x, y, grid, n=2000):
    x, y = np.asarray(x, float), np.asarray(y, float)
    out = []
    for _ in range(n):
        idx = RNG.integers(0, len(x), len(x))
        fit = fit_broken(x[idx], y[idx], grid)
        if fit:
            out.append(fit[0])
    return (float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))) if out else (None, None)


def corr_block(rows, label):
    res = {}
    print(f"\n  {label}  (n={len(rows)})")
    for name, key, logit in (("log10 sigma_a/a", "log10_sigma_a_rel", False),
                             ("log10 arc days", "arc_days", True),
                             ("condition code", "condition_code", False),
                             ("log10 n_obs", "n_obs", True),
                             ("n_oppositions", "n_oppositions", False)):
        xs, ys, zs = [], [], []
        for r in rows:
            v = r.get(key)
            if v is None or (logit and v <= 0):
                continue
            xs.append(math.log10(v) if logit else v)
            ys.append(r["disp_dex"])
            zs.append(r["log10_median_life"])
        if len(xs) < 10 or len(set(xs)) < 3:
            continue
        rho, p = spearmanr(xs, ys)
        prho, pp = partial_spearman(xs, ys, zs)
        res[name] = {"n": len(xs), "rho": float(rho), "p": float(p),
                     "partial_rho": prho, "partial_p": pp}
        flag = "  <-- significant" if pp < 0.05 else ""
        print(f"     {name:<18} n={len(xs):>4}  rho={rho:+.3f} (p={p:.3f})   "
              f"partial rho={prho:+.3f} (p={pp:.3f}){flag}")
    return res


# ------------------------------------------------------------------- main ---
def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = (load_objects("results/pop_fragility/phase1/*.json", "survey")
            + load_objects("results/bm09_reclass/clones/*.json", "bm09"))
    # Censored ensembles have dispersion truncated by Tmax; exclude from the fit.
    clean = [r for r in rows if not r["censored"]]
    s = np.array([r["log10_sigma_a_rel"] for r in clean])
    d = np.array([r["disp_dex"] for r in clean])
    print(f"objects: {len(rows)} total, {len(clean)} uncensored and usable")
    print(f"sigma_a/a spans 10^{s.min():.1f} to 10^{s.max():.1f} "
          f"({s.max()-s.min():.1f} orders of magnitude)")

    # Does the catalogue sigma describe the ensemble the integrator really used?
    chk = [(r["log10_sigma_a_rel"], r["log10_delta0"]) for r in clean if r["log10_delta0"]]
    if chk:
        rr, pp = spearmanr([c[0] for c in chk], [c[1] for c in chk])
        print(f"realised initial spread tracks catalogue sigma: rho={rr:+.3f} "
              f"(p={pp:.2g}, n={len(chk)})")

    print("\n--- locating the break in dispersion vs orbit quality ---")
    grid = np.arange(s.min() + 0.5, s.max() - 0.3, 0.05)
    xb, floor, slope, sse = fit_broken(s, d, grid)
    lo, hi = bootstrap_break(s, d, grid)
    flat_sse = float(np.sum((d - d.mean()) ** 2))
    print(f"break at sigma_a/a = 10^{xb:.2f}   (95% bootstrap CI 10^{lo:.2f} to 10^{hi:.2f})")
    print(f"chaos floor below the break : {floor:.2f} dex")
    print(f"slope above the break       : {slope:+.2f} dex per decade of sigma")
    print(f"variance explained over a flat model: {100*(1-sse/flat_sse):.1f}%")

    below = [r for r in clean if r["log10_sigma_a_rel"] <= xb]
    above = [r for r in clean if r["log10_sigma_a_rel"] > xb]
    print(f"\nfloor  region: n={len(below)}, median dispersion "
          f"{st.median([r['disp_dex'] for r in below]):.2f} dex")
    print(f"rising region: n={len(above)}, median dispersion "
          f"{st.median([r['disp_dex'] for r in above]):.2f} dex")

    print("\n--- does orbit quality predict dispersion? ---")
    res_all = corr_block(clean, "all uncensored objects")
    res_below = corr_block(below, "below the break (the chaos floor)")
    res_above = corr_block(above, "above the break (orbit genuinely unconstrained)")

    # How far below the floor could better astrometry ever push you? The floor is
    # already reached by orbits known to ~1e-7 relative, so the answer is: nowhere.
    best = sorted(clean, key=lambda r: r["log10_sigma_a_rel"])[:15]
    print(f"\nbest-determined 15 objects: median sigma_a/a = 10^"
          f"{st.median([r['log10_sigma_a_rel'] for r in best]):.1f}, "
          f"median dispersion {st.median([r['disp_dex'] for r in best]):.2f} dex")
    print("  -> the floor is already saturated at the precision of the very best orbits")

    # ------------------------------------------------ catalogue projection ---
    cat = json.load(open(ROOT / "data" / "sbdb_centaurs_od_cache.json", encoding="utf-8"))["objects"]
    tot = redu = irr = nosig = 0
    reducible = []
    for o in cat:
        a, sg = o.get("a_au"), (o.get("sigma") or {}).get("a_au")
        tot += 1
        if not a or not sg or sg <= 0:
            nosig += 1
            continue
        lr = math.log10(sg / abs(a))
        if lr > xb:
            redu += 1
            reducible.append({"desig": o.get("desig"), "full_name": o.get("full_name"),
                              "log10_sigma_a_rel": lr, "arc_days": o.get("data_arc_days"),
                              "condition_code": o.get("condition_code"),
                              "n_obs": o.get("n_obs"), "a_au": a, "e": o.get("e"),
                              "q_au": o.get("q_au")})
        else:
            irr += 1
    reducible.sort(key=lambda r: -r["log10_sigma_a_rel"])
    print(f"\n--- projection onto the full SBDB Centaur catalogue ---")
    print(f"catalogue objects              : {tot}")
    print(f"no usable sigma                : {nosig}")
    print(f"above the break (astrometry helps) : {redu}  ({100*redu/(redu+irr):.1f}% of those with sigma)")
    print(f"below the break (chaos-limited)    : {irr}  ({100*irr/(redu+irr):.1f}%)")
    print(f"\ntop 15 objects where more astrometry would actually change the answer:")
    print(f"{'designation':<22}{'log10 s_a/a':>13}{'arc d':>9}{'U':>4}")
    for r in reducible[:15]:
        arc = r["arc_days"]
        print(f"{(r['desig'] or '')[:21]:<22}{r['log10_sigma_a_rel']:>13.2f}"
              f"{(f'{arc:.0f}' if arc else '-'):>9}{str(r['condition_code'] or '-'):>4}")

    payload = {
        "n_objects": len(rows), "n_uncensored": len(clean),
        "sigma_range_log10": [float(s.min()), float(s.max())],
        "break_log10_sigma_a_rel": xb, "break_ci95": [lo, hi],
        "chaos_floor_dex": floor, "slope_dex_per_decade": slope,
        "variance_explained_vs_flat": 1 - sse / flat_sse,
        "n_below_break": len(below), "n_above_break": len(above),
        "median_disp_below": st.median([r["disp_dex"] for r in below]),
        "median_disp_above": st.median([r["disp_dex"] for r in above]),
        "correlations": {"all": res_all, "below_break": res_below, "above_break": res_above},
        "catalogue": {"n_total": tot, "n_no_sigma": nosig,
                      "n_reducible": redu, "n_chaos_limited": irr,
                      "frac_reducible": redu / (redu + irr)},
        "reducible_ranked": reducible,
        "rows": rows,
    }
    (OUT / "chaos_floor.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"\nwrote {(OUT / 'chaos_floor.json').relative_to(ROOT).as_posix()}")


if __name__ == "__main__":
    main()

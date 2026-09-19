"""Population-scale Lyapunov-time proxy from clone-ensemble snapshot dispersion.

Integrates diagonal clone ensembles with semimajor-axis recorded at fixed
times. The early-time growth rate of relative $a$ spread is a cheap proxy for
the Lyapunov time; saturation dispersion links to the chaos-floor height.

Usage:
  python src/lyapunov_proxy_pilot.py --run
  python src/lyapunov_proxy_pilot.py --analyze
"""
from __future__ import annotations

import argparse
import json
import math
import statistics as st
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))

from classify_core import (  # noqa: E402
    integrate_clone_snapshots,
    open_giant_ephem,
    sample_clones_from_modern,
)
from chaos_floor import load_objects, pool  # noqa: E402
from chaos_floor_final import ARC_MIN, CONTROLS, is_centaur  # noqa: E402
import run_pop_fragility as rpf  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "lyapunov_proxy"
CAT = ROOT / "data" / "sbdb_centaurs_od_cache.json"
SEED = 87
N_CLONES = 10
T_MAX_MYR = 25.0
SNAPSHOT_MYR = [0.002, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0]

# Published Lyapunov times (years) for validation objects; order-of-magnitude anchors.
LITERATURE_TAU_YR = {
    "1977 UB": {
        "tau_yr": 1.1e4,
        "lo": 3e3,
        "hi": 4e4,
        "ref": "Wood et al. 2018; MEGNO/LCE maps for Chiron",
    },
    "1992 AD": {
        "tau_yr": 5e3,
        "lo": 1e3,
        "hi": 2e4,
        "ref": "Tancredi et al. 2001; highly chaotic Centaur class",
    },
    "2005 TH173": {
        "tau_yr": 8e3,
        "lo": 2e3,
        "hi": 3e4,
        "ref": "Uranus-crossing Centaur; literature $10^3$--$10^5$ yr",
    },
}

MUST_HAVE = [
    "1977 UB",
    "1992 AD",
    "2005 TH173",
    "1995 DW2",
    "1998 QM107",
    "1998 TF35",
    "2000 FZ53",
    "2003 QP112",
    "2003 UW292",
    "2005 RL43",
    "2005 RO43",
    "2006 SX368",
]


def slug(desig: str) -> str:
    return desig.replace(" ", "_").replace("/", "-")


def modern_from_catalog(desig: str) -> dict | None:
    cat = json.loads(CAT.read_text(encoding="utf-8"))["objects"]
    for o in cat:
        if o.get("desig") == desig or o.get("pdes") in desig or desig in (o.get("full_name") or ""):
            return o
    return None


def build_pilot_list(n_extra: int = 28) -> list[str]:
    want = list(dict.fromkeys(MUST_HAVE))
    rows = pool(
        load_objects(["results/pop_fragility/phase1/*.json",
                      "results/pop_fragility/phase3_merged/*.json"], "survey"),
        load_objects(["results/bm09_reclass/clones/*.json",
                      "results/bm09_reclass/clones_censor_resolved/*.json"], "bm09"),
    )
    pool_rows = [
        r for r in rows
        if is_centaur(r) and not r["censored"] and (r.get("arc_days") or 0) >= ARC_MIN
        and all(r.get(k) is not None for k in CONTROLS)
    ]
    pool_rows.sort(key=lambda r: r["desig"])
    extras = []
    # Stratify by $a$ and $i$ quartiles.
    a_vals = np.array([r["a_au"] for r in pool_rows])
    i_vals = np.array([r["i_deg"] for r in pool_rows])
    aq = np.quantile(a_vals, [0.25, 0.5, 0.75])
    iq = np.quantile(i_vals, [0.25, 0.5, 0.75])
    bins: dict[tuple[int, int], list[str]] = {}
    for r in pool_rows:
        d = r["desig"]
        if d in want:
            continue
        ab = int(np.searchsorted(aq, r["a_au"], side="right"))
        ib = int(np.searchsorted(iq, r["i_deg"], side="right"))
        bins.setdefault((ab, ib), []).append(d)
    for key in sorted(bins):
        extras.extend(bins[key][:2])
    for d in pool_rows:
        if d["desig"] not in want and d["desig"] not in extras:
            extras.append(d["desig"])
        if len(want) + len(extras) >= len(MUST_HAVE) + n_extra:
            break
    return want + extras[:n_extra]


def ensemble_curve(clone_rows: list[dict]) -> list[dict]:
    """Relative $a$ spread at each snapshot time."""
    by_t: dict[float, list[float]] = {}
    for row in clone_rows:
        for t_myr, a_au in row.get("snapshots_myr") or []:
            by_t.setdefault(t_myr, []).append(float(a_au))
    curve = []
    for t_myr in SNAPSHOT_MYR:
        vals = by_t.get(t_myr, [])
        if len(vals) < 3:
            continue
        mu = float(np.mean(vals))
        if mu <= 0:
            continue
        std_rel = float(np.std(vals) / mu)
        curve.append({"t_myr": t_myr, "n": len(vals), "std_a_rel": std_rel,
                      "log10_std_a_rel": math.log10(max(std_rel, 1e-8))})
    return curve


def fit_tau_proxy(
    curve: list[dict],
    modern: dict,
    early_end_myr: float = 0.1,
) -> dict:
    """Estimate Lyapunov time from early relative-$a$ spread growth."""
    sig = (modern.get("sigma") or {}).get("a_au")
    a_nom = modern.get("a_au")
    s0 = max(float(sig) / float(a_nom), 1e-10) if sig and a_nom else 1e-10

    max_n = max((c["n"] for c in curve), default=0)
    thresh = max(3, int(math.ceil(max_n * 0.9)))
    mono: list[tuple[float, float]] = [(0.0, s0)]
    for c in curve:
        if c["n"] >= thresh:
            t_yr = c["t_myr"] * 1e6
            s = c["std_a_rel"]
            if s > mono[-1][1]:
                mono.append((t_yr, s))

    taus_all: list[float] = []
    taus_early: list[float] = []
    early_end_yr = early_end_myr * 1e6
    for i in range(1, len(mono)):
        t0, s0p = mono[i - 1]
        t1, s1 = mono[i]
        if s1 <= s0p or t1 <= t0 or s0p <= 1e-10:
            continue
        tau = (t1 - t0) / math.log(s1 / s0p)
        taus_all.append(tau)
        if t1 <= early_end_yr:
            taus_early.append(tau)

    pool = taus_early if taus_early else taus_all
    if not pool:
        return {"tau_yr": None, "slope_per_myr": None, "n_fit": 0, "taus_yr": []}

    if len(pool) == 1:
        tau_yr = float(pool[0])
    else:
        tau_yr = float(math.exp(sum(math.log(t) for t in pool) / len(pool)))
    return {
        "tau_yr": tau_yr,
        "slope_per_myr": 1e6 / tau_yr,
        "n_fit": len(pool),
        "taus_yr": pool,
        "s0_sig_rel": s0,
    }


def run_object(desig: str, modern: dict) -> dict:
    t0 = time.time()
    kernel = Path(rpf.KERNEL)
    cache = Path(rpf.CACHE_JSON)
    clones = sample_clones_from_modern(modern, N_CLONES, seed=SEED)
    rows = []
    with open_giant_ephem(kernel, cache) as ephem:
        for i, el in enumerate(clones):
            r = integrate_clone_snapshots(
                el, ephem, T_MAX_MYR * 1e6, SNAPSHOT_MYR,
            )
            life_myr = r["lifetime_yr"] / 1e6
            rows.append({
                "clone": i,
                "a0": el["a"],
                "e0": el["e"],
                "lifetime_myr": life_myr,
                "escaped": r["escaped"],
                "survived_full": r["survived_full"],
                "snapshots_myr": r.get("snapshots_myr") or [],
            })
    lives = [max(r["lifetime_myr"], 1e-4) for r in rows]
    lg = [math.log10(x) for x in lives]
    curve = ensemble_curve(rows)
    tau_fit = fit_tau_proxy(curve, modern)
    sat = curve[-1]["log10_std_a_rel"] if curve else None
    out = {
        "desig": desig,
        "modern_nominal": modern,
        "seed": SEED,
        "n_clones": N_CLONES,
        "t_max_myr": T_MAX_MYR,
        "snapshots_myr": SNAPSHOT_MYR,
        "clones": rows,
        "ensemble_curve": curve,
        "dispersion_dex": float(st.pstdev(lg)) if len(lg) > 1 else 0.0,
        "tau_proxy_yr": tau_fit["tau_yr"],
        "tau_fit": tau_fit,
        "saturation_log10_std_a": sat,
        "elapsed_s": time.time() - t0,
    }
    return out


def cmd_run(max_objects: int | None, force: bool) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    desigs = build_pilot_list()
    if max_objects:
        desigs = desigs[:max_objects]
    print(f"pilot objects: {len(desigs)}")
    results = []
    for i, desig in enumerate(desigs, 1):
        path = OUT / f"{slug(desig)}.json"
        if path.exists() and not force:
            results.append(json.loads(path.read_text(encoding="utf-8")))
            print(f"  [{i}/{len(desigs)}] {desig} skip (exists)")
            continue
        modern = modern_from_catalog(desig)
        if not modern:
            print(f"  [{i}/{len(desigs)}] {desig} NO CATALOG MATCH")
            continue
        print(f"  [{i}/{len(desigs)}] {desig} ...", flush=True)
        blob = run_object(desig, modern)
        path.write_text(json.dumps(blob, indent=2), encoding="utf-8")
        results.append(blob)
        tp = blob["tau_proxy_yr"]
        tp_s = f"{tp:.2e}" if tp else "n/a"
        print(f"       tau_proxy={tp_s} yr  disp={blob['dispersion_dex']:.3f} dex")
    summary = {"n": len(results), "objects": results}
    (OUT / "pilot_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {OUT / 'pilot_summary.json'}")


def cmd_analyze() -> None:
    files = sorted(OUT.glob("*.json"))
    files = [f for f in files if f.name != "pilot_summary.json" and f.name != "analysis.json"]
    objs = [json.loads(f.read_text(encoding="utf-8")) for f in files]
    valid = [o for o in objs if o.get("tau_proxy_yr") and o["tau_proxy_yr"] > 0]
    print(f"objects with tau proxy: {len(valid)} / {len(objs)}")

    # Validation against literature anchors.
    print("\n--- literature validation ---")
    val_rows = []
    for o in objs:
        d = o["desig"]
        lit = LITERATURE_TAU_YR.get(d)
        if not lit or not o.get("tau_proxy_yr"):
            continue
        ratio = o["tau_proxy_yr"] / lit["tau_yr"]
        in_band = lit["lo"] <= o["tau_proxy_yr"] <= lit["hi"]
        print(f"  {d:<14} proxy={o['tau_proxy_yr']:.2e} yr  lit={lit['tau_yr']:.2e} yr  "
              f"ratio={ratio:.2f}  in_band={in_band}")
        val_rows.append({"desig": d, "proxy": o["tau_proxy_yr"], "lit": lit["tau_yr"],
                         "ratio": ratio, "in_band": in_band})

    # Floor height vs chaos rate.
    disp = [o["dispersion_dex"] for o in valid]
    tau = [o["tau_proxy_yr"] for o in valid]
    inv_tau = [1.0 / t for t in tau]
    r1, p1 = spearmanr(inv_tau, disp)
    log_tau = [math.log10(t) for t in tau]
    r2, p2 = spearmanr(log_tau, disp)
    print(f"\n--- floor vs chaos rate (n={len(valid)}) ---")
    print(f"  Spearman(1/tau_proxy, dispersion_dex)  rho={r1:+.3f}  p={p1:.4f}")
    print(f"  Spearman(log10(tau_proxy), dispersion_dex) rho={r2:+.3f}  p={p2:.4f}")

    analysis = {
        "n_objects": len(objs),
        "n_with_tau": len(valid),
        "validation": val_rows,
        "spearman_inv_tau_disp": {"rho": float(r1), "p": float(p1), "n": len(valid)},
        "spearman_log_tau_disp": {"rho": float(r2), "p": float(p2), "n": len(valid)},
        "objects": [
            {
                "desig": o["desig"],
                "tau_proxy_yr": o.get("tau_proxy_yr"),
                "dispersion_dex": o.get("dispersion_dex"),
                "a_au": o["modern_nominal"].get("a_au"),
                "e": o["modern_nominal"].get("e"),
                "i_deg": o["modern_nominal"].get("i_deg"),
                "q_au": o["modern_nominal"].get("q_au"),
            }
            for o in objs
        ],
    }
    (OUT / "analysis.json").write_text(json.dumps(analysis, indent=2), encoding="utf-8")
    print(f"wrote {OUT / 'analysis.json'}")


def cmd_recalc() -> None:
    """Recompute tau proxy from stored ensemble curves (no re-integration)."""
    files = sorted(OUT.glob("*.json"))
    files = [f for f in files if f.name not in ("pilot_summary.json", "analysis.json")]
    for f in files:
        blob = json.loads(f.read_text(encoding="utf-8"))
        curve = blob.get("ensemble_curve") or []
        modern = blob.get("modern_nominal") or {}
        tau_fit = fit_tau_proxy(curve, modern)
        blob["tau_proxy_yr"] = tau_fit["tau_yr"]
        blob["tau_fit"] = tau_fit
        f.write_text(json.dumps(blob, indent=2), encoding="utf-8")
        tp = tau_fit["tau_yr"]
        tp_s = f"{tp:.2e}" if tp else "n/a"
        print(f"  {blob['desig']:<14} tau_proxy={tp_s} yr")
    summary_path = OUT / "pilot_summary.json"
    if summary_path.exists():
        objs = [json.loads(f.read_text(encoding="utf-8")) for f in files]
        summary_path.write_text(json.dumps({"n": len(objs), "objects": objs}, indent=2), encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--analyze", action="store_true")
    ap.add_argument("--recalc", action="store_true")
    ap.add_argument("--force", action="store_true", help="overwrite existing per-object JSON")
    ap.add_argument("--max-objects", type=int, default=None)
    args = ap.parse_args()
    if args.run:
        cmd_run(args.max_objects, args.force)
    if args.recalc:
        cmd_recalc()
    if args.analyze:
        cmd_analyze()
    if not args.run and not args.analyze and not args.recalc:
        ap.print_help()


if __name__ == "__main__":
    main()

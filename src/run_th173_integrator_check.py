"""TH173 WHFast vs Mercurius (optional IAS15) referee check.

Production used plain WHFast dt=0.25 yr. Modern TH173 is Uranus-crossing, so
close encounters matter. This script re-runs the first 20 production clones
(seed 87) with a hybrid/high-accuracy integrator and compares lifetimes / v2
classes. Does NOT rewrite the production claim — writes check artifacts only.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from classify_core import (
    LONG_LIVED_MYR,
    classify_bailey,
    integrate_clone,
    open_giant_ephem,
    sample_clones_from_modern,
)

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT
SEED = 87
DEFAULT_KERNEL = Path(r"C:\Users\Hp\Downloads\Zenith")
CACHE_JSON = ROOT / "data" / "j2000_giants_horizons_de441.json"


def load_th173_modern() -> dict:
    modern_path = ROOT / "data" / "modern_elements.json"
    if modern_path.exists():
        blob = json.loads(modern_path.read_text(encoding="utf-8"))
        if "2005 TH173" in blob:
            return blob["2005 TH173"]
    th = json.loads((ROOT / "data" / "th173_elements.json").read_text(encoding="utf-8"))
    return th["modern"]


def spearman_rho(x: np.ndarray, y: np.ndarray) -> float:
    """Spearman rank correlation; NaN if undefined."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if len(x) < 2:
        return float("nan")
    rx = np.argsort(np.argsort(x)).astype(np.float64)
    ry = np.argsort(np.argsort(y)).astype(np.float64)
    rx -= rx.mean()
    ry -= ry.mean()
    den = float(np.sqrt(np.sum(rx * rx) * np.sum(ry * ry)))
    if den <= 0:
        return float("nan")
    return float(np.sum(rx * ry) / den)


def summarize_pair(rows: list[dict], a_name: str, b_name: str) -> dict:
    a_life = np.array([r[a_name]["lifetime_yr"] for r in rows], dtype=np.float64)
    b_life = np.array([r[b_name]["lifetime_yr"] for r in rows], dtype=np.float64)
    a_cls = [r[a_name]["class_v2"] for r in rows]
    b_cls = [r[b_name]["class_v2"] for r in rows]
    agree = sum(1 for i in range(len(rows)) if a_cls[i] == b_cls[i])
    short_both = sum(
        1
        for i in range(len(rows))
        if a_life[i] < LONG_LIVED_MYR * 1e6 and b_life[i] < LONG_LIVED_MYR * 1e6
    )
    q_both = sum(1 for i in range(len(rows)) if a_cls[i] == "Q" and b_cls[i] == "Q")
    # WHFast→D but Mercurius long-lived / Q (headline-risk pattern)
    risk = []
    for r in rows:
        wa, wb = r[a_name], r[b_name]
        if wa["class_v2"] == "D" and (
            wb["class_v2"] in ("Q", "R") or wb["lifetime_yr"] >= LONG_LIVED_MYR * 1e6
        ):
            risk.append(
                {
                    "clone": r["clone"],
                    f"{a_name}_life_myr": wa["lifetime_myr"],
                    f"{a_name}_class": wa["class_v2"],
                    f"{b_name}_life_myr": wb["lifetime_myr"],
                    f"{b_name}_class": wb["class_v2"],
                }
            )
    return {
        f"median_life_myr_{a_name}": float(np.median(a_life) / 1e6),
        f"median_life_myr_{b_name}": float(np.median(b_life) / 1e6),
        "n_life_lt_22myr_both": int(short_both),
        "n_Q_both": int(q_both),
        "n_class_agree": int(agree),
        "frac_class_agree": float(agree / len(rows)) if rows else float("nan"),
        "spearman_lifetime": spearman_rho(a_life, b_life),
        "whfast_D_but_hybrid_long_or_Q": risk,
        "n_risk_clones": len(risk),
    }


def write_markdown(out: dict, path: Path) -> None:
    s = out["summary"]
    call = out["call"]
    lines = [
        "# TH173 integrator check (WHFast vs Mercurius)",
        "",
        f"**Call:** {call}",
        "",
        "## Setup",
        "",
        f"- Object: 2005 TH173 (modern SBDB)",
        f"- Clones: first {out['n_clones']} from `sample_clones_from_modern` (seed {out['seed']})",
        f"- t_max: {out['t_max_myr']} Myr (same as production)",
        f"- Integrators: {', '.join(out['integrators'])}",
        f"- Ephem: `{out['ephem']}`",
        f"- Escape cuts / sampling / forces: unchanged from `classify_core.integrate_clone`",
        "",
        "## Key numbers",
        "",
        f"- Median lifetime WHFast: **{s['median_life_myr_whfast']:.3f} Myr**",
        f"- Median lifetime Mercurius: **{s['median_life_myr_mercurius']:.3f} Myr**",
        f"- Clones with life < 22 Myr under **both**: **{s['n_life_lt_22myr_both']}** / {out['n_clones']}",
        f"- Clones Q under **both**: **{s['n_Q_both']}** / {out['n_clones']}",
        f"- v2 class agreement: **{s['n_class_agree']}** / {out['n_clones']} "
        f"({s['frac_class_agree']:.2f})",
        f"- Spearman ρ(lifetime WHFast, Mercurius): **{s['spearman_lifetime']:.3f}**",
        f"- Risk clones (WHFast→D but Mercurius long-lived/Q/R): **{s['n_risk_clones']}**",
        "",
    ]
    if s["whfast_D_but_hybrid_long_or_Q"]:
        lines.append("### Risk clones")
        lines.append("")
        for r in s["whfast_D_but_hybrid_long_or_Q"]:
            lines.append(
                f"- clone {r['clone']}: WHFast {r['whfast_class']} "
                f"({r['whfast_life_myr']:.2f} Myr) vs Mercurius {r['mercurius_class']} "
                f"({r['mercurius_life_myr']:.2f} Myr)"
            )
        lines.append("")
    lines.extend(
        [
            "## Interpretation",
            "",
            (
                "If Mercurius confirms short lives and no Q, the WHFast Q-rejection gap is closed "
                "for this referee concern. If Mercurius yields long-lived/Q survivors where WHFast "
                "gave D, the production headline is at risk."
            ),
            "",
            f"**Verdict:** {call}",
            "",
            f"Elapsed: {out['elapsed_s']:.1f} s",
            "",
        ]
    )
    if "ias15" in out["integrators"] and "summary_ias15" in out:
        si = out["summary_ias15"]
        lines.extend(
            [
                "## Optional IAS15 vs WHFast",
                "",
                f"- Median life IAS15: {si['median_life_myr_ias15']:.3f} Myr",
                f"- Class agree with WHFast: {si['n_class_agree']} / {out['n_clones']}",
                f"- Spearman ρ: {si['spearman_lifetime']:.3f}",
                "",
            ]
        )
    path.write_text("\n".join(lines), encoding="utf-8")


def run_check(
    n_clones: int,
    t_max_yr: float,
    integrators: list[str],
    kernel_dir: Path,
) -> dict:
    modern = load_th173_modern()
    clones = sample_clones_from_modern(modern, n_clones, seed=SEED)
    (ROOT / "results").mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    t0 = time.time()
    partial_path = ROOT / "results" / "th173_integrator_check.partial.json"
    start_i = 0
    if partial_path.exists():
        try:
            prev = json.loads(partial_path.read_text(encoding="utf-8"))
            if prev.get("n_clones") == n_clones and prev.get("integrators") == integrators:
                rows = list(prev.get("clones") or [])
                start_i = len(rows)
                print(f"RESUME integrator check from clone {start_i}/{n_clones}", flush=True)
        except Exception as exc:
            print(f"partial load failed ({exc}); starting fresh", flush=True)
            rows = []
            start_i = 0
    with open_giant_ephem(kernel_dir, CACHE_JSON) as ephem:
        ephem_name = type(ephem).__name__
        print(f"ephem={ephem_name} integrators={integrators} n={n_clones} t_max={t_max_yr/1e6} Myr", flush=True)
        for i, el in enumerate(clones):
            if i < start_i:
                continue
            row: dict = {
                "clone": i,
                "a0": el["a"],
                "e0": el["e"],
            }
            for integ in integrators:
                t_i = time.time()
                r = integrate_clone(el, ephem, t_max_yr, integrator=integ)
                c = classify_bailey(r["a_au"], r["t_yr"], r["lifetime_yr"], rule="v2")
                row[integ] = {
                    "escaped": r["escaped"],
                    "survived_full": r["survived_full"],
                    "lifetime_yr": r["lifetime_yr"],
                    "lifetime_myr": c["lifetime_myr"],
                    "class_v2": c["class_v2"],
                    "class": c["class"],
                    "flat": c["flat"],
                    "long_lived": c["long_lived"],
                    "H": c["H"],
                    "r2": c["r2"],
                    "elapsed_s": time.time() - t_i,
                }
                print(
                    f"  clone {i+1}/{n_clones} {integ}: class={c['class_v2']} "
                    f"life={c['lifetime_myr']:.2f} Myr escaped={r['escaped']} "
                    f"({row[integ]['elapsed_s']:.1f}s)",
                    flush=True,
                )
            rows.append(row)
            # checkpoint after each clone
            partial = {
                "status": "running",
                "done_clones": i + 1,
                "n_clones": n_clones,
                "integrators": integrators,
                "ephem": ephem_name,
                "clones": rows,
            }
            (ROOT / "results" / "th173_integrator_check.partial.json").write_text(
                json.dumps(partial, indent=2), encoding="utf-8"
            )

    summary = summarize_pair(rows, "whfast", "mercurius")
    # Gap closed iff: no Q under Mercurius, all short under both (or at least Mercurius),
    # and no risk clones where WHFast D but Mercurius long/Q.
    n_merc_q = sum(1 for r in rows if r["mercurius"]["class_v2"] == "Q")
    n_merc_long = sum(1 for r in rows if r["mercurius"]["long_lived"])
    gap_closed = (
        summary["n_risk_clones"] == 0
        and n_merc_q == 0
        and summary["n_life_lt_22myr_both"] == n_clones
    )
    # Softer: if Mercurius also has no Q and median still short, even if a few survive >22 Myr
    # under both with R — still not Q. Headline is Q-rejection.
    if summary["n_risk_clones"] == 0 and n_merc_q == 0 and n_merc_long == 0:
        gap_closed = True
        call = "gap closed"
    elif summary["n_risk_clones"] > 0 or n_merc_q > 0:
        gap_closed = False
        call = "headline at risk"
    else:
        # Mercurius long-lived R but not Q — still weakens the "all D / short" story
        gap_closed = False
        call = "headline at risk"

    out: dict = {
        "label": "th173_integrator_check",
        "desig": "2005 TH173",
        "seed": SEED,
        "n_clones": n_clones,
        "t_max_myr": t_max_yr / 1e6,
        "integrators": integrators,
        "ephem": ephem_name,
        "modern_nominal": {
            k: modern.get(k)
            for k in ("a_au", "e", "i_deg", "q_au", "condition_code", "data_arc_days")
        },
        "summary": summary,
        "n_mercurius_Q": n_merc_q,
        "n_mercurius_long_lived": n_merc_long,
        "gap_closed": bool(gap_closed),
        "call": call,
        "elapsed_s": time.time() - t0,
        "clones": rows,
    }
    if "ias15" in integrators:
        a_life = np.array([r["whfast"]["lifetime_yr"] for r in rows])
        b_life = np.array([r["ias15"]["lifetime_yr"] for r in rows])
        agree = sum(1 for r in rows if r["whfast"]["class_v2"] == r["ias15"]["class_v2"])
        out["summary_ias15"] = {
            "median_life_myr_whfast": float(np.median(a_life) / 1e6),
            "median_life_myr_ias15": float(np.median(b_life) / 1e6),
            "n_class_agree": agree,
            "frac_class_agree": agree / n_clones,
            "spearman_lifetime": spearman_rho(a_life, b_life),
        }

    json_path = ROOT / "results" / "th173_integrator_check.json"
    md_path = ROOT / "results" / "TH173_INTEGRATOR_CHECK.md"
    json_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    write_markdown(out, md_path)
    partial = ROOT / "results" / "th173_integrator_check.partial.json"
    if partial.exists():
        partial.unlink()
    print(json.dumps({k: out[k] for k in out if k != "clones"}, indent=2), flush=True)
    print(f"wrote {json_path}", flush=True)
    print(f"wrote {md_path}", flush=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-clones", type=int, default=20)
    ap.add_argument("--t-max-myr", type=float, default=10.0)
    ap.add_argument("--kernel-dir", default=str(DEFAULT_KERNEL))
    ap.add_argument(
        "--include-ias15",
        action="store_true",
        help="Also run IAS15 (slow); Mercurius remains the encounter priority.",
    )
    args = ap.parse_args()
    integrators = ["whfast", "mercurius"]
    if args.include_ias15:
        integrators.append("ias15")
    run_check(args.n_clones, args.t_max_myr * 1e6, integrators, Path(args.kernel_dir))


if __name__ == "__main__":
    main()

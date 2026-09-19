"""Experiment C: Bailey-era (a,e,i) vs modern ICs, same integrator/forces.

Isolates orbit improvement from numerical pipeline differences.
Angles: use modern Ω,ω,M when available (Table 2 lacks them).
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from classify_core import DE441Ephem, classify_bailey, integrate_clone

ROOT = Path(__file__).resolve().parents[1]
SEED = 87


def slug(desig: str) -> str:
    return desig.replace(" ", "_").replace("/", "-")


def elems_bailey_aei(bailey: dict, modern: dict, angles: str = "modern") -> dict:
    """Bailey Table 2 a,e,i; angles from modern SBDB or zeros."""
    if angles == "modern":
        om = float(modern["om_deg"])
        w = float(modern["w_deg"])
        ma = float(modern["ma_deg"])
    else:
        om = w = ma = 0.0
    return {
        "a": float(bailey["a_au"]),
        "e": float(bailey["e"]),
        "inc": float(np.deg2rad(bailey["i_deg"])),
        "Omega": float(np.deg2rad(om)),
        "omega": float(np.deg2rad(w)),
        "M": float(np.deg2rad(ma)),
    }


def elems_modern(modern: dict) -> dict:
    return {
        "a": float(modern["a_au"]),
        "e": float(modern["e"]),
        "inc": float(np.deg2rad(modern["i_deg"])),
        "Omega": float(np.deg2rad(modern["om_deg"])),
        "omega": float(np.deg2rad(modern["w_deg"])),
        "M": float(np.deg2rad(modern["ma_deg"])),
    }


def run_one(el: dict, ephem: DE441Ephem, t_max_yr: float) -> dict:
    r = integrate_clone(el, ephem, t_max_yr)
    c = classify_bailey(r["a_au"], r["t_yr"], r["lifetime_yr"])
    return {
        "a0": el["a"],
        "e0": el["e"],
        "escaped": r["escaped"],
        "survived_full": r["survived_full"],
        "t_yr": r["t_yr"].tolist(),
        "a_au": r["a_au"].tolist(),
        **c,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--classes", default="QR")
    ap.add_argument("--t-max-myr", type=float, default=40.0)
    ap.add_argument("--angles", default="modern", choices=("modern", "zero"))
    ap.add_argument("--kernel-dir", default=r"C:\Users\Hp\Downloads\Zenith")
    ap.add_argument("--label", default="exp_c")
    ap.add_argument("--only", default="")
    args = ap.parse_args()

    table = json.loads((ROOT / "data" / "bailey_table2.json").read_text(encoding="utf-8"))
    modern_all = json.loads((ROOT / "data" / "modern_elements.json").read_text(encoding="utf-8"))
    allowed = set(args.classes.upper())
    only = {x.strip() for x in args.only.split(",") if x.strip()}

    out_dir = ROOT / "results" / "exp_c"
    out_dir.mkdir(parents=True, exist_ok=True)
    results = []

    with DE441Ephem(Path(args.kernel_dir)) as ephem:
        for row in table["objects"]:
            if row["bailey_class"] == "NT" or row["bailey_class"] not in allowed:
                continue
            des = row["desig"]
            if only and des not in only:
                continue
            if des not in modern_all or "a_au" not in modern_all[des]:
                print(f"SKIP {des}", flush=True)
                continue

            modern = modern_all[des]
            print(f"=== Exp C {des} Bailey={row['bailey_class']} ===", flush=True)
            t0 = time.time()
            el_b = elems_bailey_aei(row, modern, angles=args.angles)
            el_m = elems_modern(modern)
            b = run_one(el_b, ephem, args.t_max_myr * 1e6)
            print(f"  bailey-IC class={b['class']} life={b['lifetime_myr']:.2f}", flush=True)
            m = run_one(el_m, ephem, args.t_max_myr * 1e6)
            print(f"  modern-IC class={m['class']} life={m['lifetime_myr']:.2f}", flush=True)

            # drop heavy series from summary row; keep in per-object file
            rec = {
                "desig": des,
                "bailey_class": row["bailey_class"],
                "angles_mode": args.angles,
                "bailey_ic": {k: b[k] for k in b if k not in ("t_yr", "a_au")},
                "modern_ic": {k: m[k] for k in m if k not in ("t_yr", "a_au")},
                "same_class": b["class"] == m["class"],
                "orbit_drives_flip": (b["class"] == row["bailey_class"]) and (m["class"] != row["bailey_class"]),
                "elapsed_s": time.time() - t0,
                "elements": {
                    "bailey_aei": {"a": row["a_au"], "e": row["e"], "i": row["i_deg"]},
                    "modern_aei": {"a": modern["a_au"], "e": modern["e"], "i": modern["i_deg"]},
                },
            }
            full = {**rec, "bailey_ic_full": b, "modern_ic_full": m}
            (out_dir / f"{slug(des)}.json").write_text(json.dumps(full, indent=2), encoding="utf-8")
            results.append(rec)
            print(
                f"  -> B={b['class']} M={m['class']} orbit_drives_flip={rec['orbit_drives_flip']}",
                flush=True,
            )

    summary = {
        "label": args.label,
        "t_max_myr": args.t_max_myr,
        "angles": args.angles,
        "n": len(results),
        "n_orbit_drives_flip": sum(1 for r in results if r["orbit_drives_flip"]),
        "n_same_class": sum(1 for r in results if r["same_class"]),
        "results": results,
    }
    path = ROOT / "results" / f"exp_c_summary_{args.label}.json"
    path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({k: summary[k] for k in summary if k != "results"}, indent=2), flush=True)


if __name__ == "__main__":
    main()

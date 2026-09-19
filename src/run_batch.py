"""Batch reclassify Bailey Table 2 objects under modern SBDB elements."""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from classify_core import (
    DE441Ephem,
    classify_bailey,
    integrate_clone,
    sample_clones_from_modern,
)

ROOT = Path(__file__).resolve().parents[1]
SEED = 87


def slug(desig: str) -> str:
    return desig.replace(" ", "_").replace("/", "-")


def run_object(
    desig: str,
    bailey: dict,
    modern: dict,
    ephem: DE441Ephem,
    n_clones: int,
    t_max_yr: float,
) -> dict:
    clones = sample_clones_from_modern(modern, n_clones, seed=SEED)
    rows = []
    t0 = time.time()
    for i, el in enumerate(clones):
        r = integrate_clone(el, ephem, t_max_yr)
        c = classify_bailey(r["a_au"], r["t_yr"], r["lifetime_yr"])
        rows.append(
            {
                "clone": i,
                "a0": el["a"],
                "e0": el["e"],
                "escaped": r["escaped"],
                "survived_full": r["survived_full"],
                **c,
            }
        )
        print(
            f"  {desig} clone {i+1}/{n_clones} class={c['class']} "
            f"life={c['lifetime_myr']:.2f} Myr",
            flush=True,
        )

    classes = [r["class"] for r in rows]
    counts = {k: classes.count(k) for k in ("D", "R", "Q")}
    modal = max(counts, key=counts.get)
    bailey_class = bailey["bailey_class"]
    flip = modal != bailey_class
    out = {
        "desig": desig,
        "bailey": bailey,
        "modern_nominal": {
            k: modern.get(k)
            for k in (
                "a_au",
                "e",
                "i_deg",
                "q_au",
                "condition_code",
                "data_arc_days",
                "n_obs",
                "spkid",
                "pdes",
            )
        },
        "seed": SEED,
        "n_clones": n_clones,
        "t_max_myr": t_max_yr / 1.0e6,
        "counts": counts,
        "modal_class": modal,
        "frac_modal": counts[modal] / len(rows),
        "bailey_class": bailey_class,
        "flip": bool(flip),
        "call": f"{bailey_class}->{modal}" if flip else f"stay {modal}",
        "elapsed_s": time.time() - t0,
        "clones": rows,
    }
    path = ROOT / "results" / "batch" / f"{slug(desig)}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    return out


def summarize(results: list[dict]) -> dict:
    flips = [r for r in results if r.get("flip")]
    matrix: dict[str, int] = {}
    for r in results:
        key = f"{r['bailey_class']}->{r['modal_class']}"
        matrix[key] = matrix.get(key, 0) + 1
    return {
        "n": len(results),
        "n_flip": len(flips),
        "flip_rate": (len(flips) / len(results)) if results else 0.0,
        "matrix": matrix,
        "flips": [{"desig": r["desig"], "call": r["call"], "counts": r["counts"]} for r in flips],
        "all": [
            {
                "desig": r["desig"],
                "bailey": r["bailey_class"],
                "modern": r["modal_class"],
                "flip": r["flip"],
                "counts": r["counts"],
                "arc_days": r["modern_nominal"].get("data_arc_days"),
                "U": r["modern_nominal"].get("condition_code"),
            }
            for r in results
        ],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--classes", default="QR", help="Bailey classes to include, e.g. QR or DQRN")
    ap.add_argument("--n-clones", type=int, default=5)
    ap.add_argument("--t-max-myr", type=float, default=10.0)
    ap.add_argument("--label", default="qr_pilot")
    ap.add_argument("--kernel-dir", default=r"C:\Users\Hp\Downloads\Zenith")
    ap.add_argument("--only", default="", help="Comma-separated desigs to restrict")
    args = ap.parse_args()

    table = json.loads((ROOT / "data" / "bailey_table2.json").read_text(encoding="utf-8"))
    modern_all = json.loads((ROOT / "data" / "modern_elements.json").read_text(encoding="utf-8"))
    allowed = set(args.classes.upper())
    only = {x.strip() for x in args.only.split(",") if x.strip()}

    targets = []
    for row in table["objects"]:
        if row["bailey_class"] == "NT":
            continue
        if row["bailey_class"] not in allowed:
            continue
        if only and row["desig"] not in only:
            continue
        if row["desig"] not in modern_all or "a_au" not in modern_all[row["desig"]]:
            print(f"SKIP missing modern: {row['desig']}", flush=True)
            continue
        targets.append(row)

    print(f"Running {len(targets)} objects, {args.n_clones} clones, {args.t_max_myr} Myr", flush=True)
    results = []
    with DE441Ephem(Path(args.kernel_dir)) as ephem:
        for row in targets:
            des = row["desig"]
            print(f"=== {des} Bailey={row['bailey_class']} ===", flush=True)
            out = run_object(
                des,
                row,
                modern_all[des],
                ephem,
                args.n_clones,
                args.t_max_myr * 1.0e6,
            )
            results.append(out)
            print(f"  -> {out['call']} counts={out['counts']}", flush=True)

    summary = summarize(results)
    summary["label"] = args.label
    summary["n_clones"] = args.n_clones
    summary["t_max_myr"] = args.t_max_myr
    out_path = ROOT / "results" / f"batch_summary_{args.label}.json"
    out_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()

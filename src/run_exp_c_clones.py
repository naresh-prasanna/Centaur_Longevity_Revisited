"""Experiment C with 20-clone ensembles on Bailey-IC vs modern-IC branches."""
from __future__ import annotations

import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SEED = 87
KERNEL = r"C:\Users\Hp\Downloads\Zenith"
CACHE_JSON = ROOT / "data" / "j2000_giants_horizons_de441.json"


def slug(d: str) -> str:
    return d.replace(" ", "_")


def _worker(payload: dict) -> dict:
    import sys

    sys.path.insert(0, payload["src_dir"])
    from classify_core import classify_bailey, integrate_clone, open_giant_ephem, sample_clones_from_modern  # noqa
    from run_exp_c import elems_bailey_aei, elems_modern  # noqa

    des = payload["desig"]
    row = payload["bailey"]
    modern = payload["modern"]
    n = payload["n_clones"]
    t_max = payload["t_max_yr"]
    out_path = Path(payload["out_path"])

    if out_path.exists() and not payload.get("force"):
        return {"desig": des, "skipped": True}

    clones_m = sample_clones_from_modern(modern, n, seed=SEED)
    nom_m = clones_m[0]
    bailey_nom = elems_bailey_aei(row, modern, angles="modern")
    clones_b = []
    for i, el_m in enumerate(clones_m):
        el_b = dict(el_m)
        if i == 0:
            el_b = bailey_nom
        else:
            # Same angle/M perturbations as modern clone; BM09 (a,e,i) center + delta from modern nominal
            el_b["a"] = float(row["a_au"]) + (el_m["a"] - nom_m["a"])
            el_b["e"] = float(np.clip(float(row["e"]) + (el_m["e"] - nom_m["e"]), 0.0, 0.95))
            el_b["inc"] = float(np.deg2rad(row["i_deg"])) + (el_m["inc"] - nom_m["inc"])
        clones_b.append(el_b)

    def run_branch(clones: list[dict]) -> list[dict]:
        rows = []
        with open_giant_ephem(Path(payload["kernel_dir"]), Path(payload["cache_json"])) as ephem:
            for j, el in enumerate(clones):
                r = integrate_clone(el, ephem, t_max)
                c = classify_bailey(r["a_au"], r["t_yr"], r["lifetime_yr"], rule="v2")
                rows.append({"clone": j, **c, "escaped": r["escaped"], "survived_full": r["survived_full"]})
        return rows

    t0 = time.time()
    b_rows = run_branch(clones_b)
    m_rows = run_branch(clones_m)

    def summarize(rows: list[dict]) -> dict:
        cls = [r["class"] for r in rows]
        cnt = {k: cls.count(k) for k in ("D", "R", "Q")}
        modal = max(cnt, key=cnt.get)
        return {"counts": cnt, "modal": modal, "clones": rows}

    sb = summarize(b_rows)
    sm = summarize(m_rows)
    rec = {
        "desig": des,
        "bailey_class": row["bailey_class"],
        "n_clones": n,
        "t_max_myr": t_max / 1e6,
        "seed": SEED,
        "bailey_ic": sb,
        "modern_ic": sm,
        "modal_flip_between_branches": sb["modal"] != sm["modal"],
        "elapsed_s": time.time() - t0,
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(rec, indent=2), encoding="utf-8")
    return {"desig": des, "skipped": False, "bailey": sb["modal"], "modern": sm["modal"]}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-clones", type=int, default=20)
    ap.add_argument("--t-max-myr", type=float, default=40.0)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--kernel-dir", default=KERNEL)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    table = json.loads((ROOT / "data" / "bailey_table2.json").read_text(encoding="utf-8"))
    modern_all = json.loads((ROOT / "data" / "modern_elements.json").read_text(encoding="utf-8"))
    out_dir = ROOT / "results" / "exp_c_clones"
    out_dir.mkdir(parents=True, exist_ok=True)

    payloads = []
    for row in table["objects"]:
        if row["bailey_class"] not in ("Q", "R"):
            continue
        des = row["desig"]
        if des not in modern_all:
            continue
        payloads.append(
            {
                "desig": des,
                "bailey": row,
                "modern": modern_all[des],
                "n_clones": args.n_clones,
                "t_max_yr": args.t_max_myr * 1e6,
                "out_path": str(out_dir / f"{slug(des)}.json"),
                "kernel_dir": args.kernel_dir,
                "cache_json": str(CACHE_JSON),
                "src_dir": str(ROOT / "src"),
                "force": args.force,
            }
        )

    results = []
    if args.workers <= 1:
        for p in payloads:
            results.append(_worker(p))
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futs = {ex.submit(_worker, p): p["desig"] for p in payloads}
            for fut in as_completed(futs):
                results.append(fut.result())
                print(results[-1], flush=True)

    summary = {
        "n": len(results),
        "n_branch_modal_flip": sum(1 for r in results if not r.get("skipped") and Path(r.get("path", "")).exists()),
        "objects": results,
    }
    # recount from files
    objs = [json.loads(p.read_text()) for p in sorted(out_dir.glob("*.json"))]
    summary = {
        "n": len(objs),
        "n_branch_modal_flip": sum(1 for o in objs if o.get("modal_flip_between_branches")),
        "objects": [
            {
                "desig": o["desig"],
                "bailey_modal": o["bailey_ic"]["modal"],
                "modern_modal": o["modern_ic"]["modal"],
                "flip": o["modal_flip_between_branches"],
            }
            for o in objs
        ],
    }
    path = ROOT / "results" / "exp_c_clones_summary.json"
    path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()

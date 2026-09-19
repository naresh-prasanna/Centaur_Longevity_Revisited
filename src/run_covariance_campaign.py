"""Run clone ensembles for TH173 / QP112 / Q/R using SBDB covariance sampling."""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from classify_core import classify_bailey, integrate_clone, open_giant_ephem
from covariance_clones import prefetch, sample_clones_from_covariance

ROOT = Path(__file__).resolve().parents[1]
FS = ROOT / "results" / "fullscale"
SEED = 87
KERNEL = r"C:\Users\Hp\Downloads\Zenith"
CACHE_JSON = ROOT / "data" / "j2000_giants_horizons_de441.json"

CAMPAIGNS = {
    "cov_th173": (["2005 TH173"], 100, 40.0),
    "cov_th173_q10": (["2005 TH173"], 20, 10.0),
    "cov_qp112": (["2003 QP112"], 20, 40.0),
    "cov_qr20": (
        [
            "1995 DW2", "1998 QM107", "1998 TF35", "2000 FZ53", "2003 QP112",
            "2003 UW292", "2005 RL43", "2005 RO43", "2005 TH173", "2006 SX368",
        ],
        20,
        40.0,
    ),
}


def slug(d: str) -> str:
    return d.replace(" ", "_")


def run_object(desig: str, n: int, t_max_myr: float, campaign: str, kernel: str, force: bool) -> dict:
    out_dir = FS / campaign
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{slug(desig)}.json"
    if out_path.exists() and not force:
        return json.loads(out_path.read_text(encoding="utf-8"))

    table = json.loads((ROOT / "data" / "bailey_table2.json").read_text(encoding="utf-8"))
    bailey = next(o for o in table["objects"] if o["desig"] == desig)
    clones = sample_clones_from_covariance(desig, n, seed=SEED)
    rows = []
    t0 = time.time()
    partial_path = out_path.with_suffix(".partial.json")
    start_i = 0
    if partial_path.exists() and not force:
        try:
            prev_p = json.loads(partial_path.read_text(encoding="utf-8"))
            if prev_p.get("desig") == desig and prev_p.get("n_clones") == n:
                rows = list(prev_p.get("clones") or [])
                start_i = len(rows)
                print(f"  RESUME {desig} from clone {start_i}/{n}", flush=True)
        except Exception as exc:
            print(f"  partial load failed ({exc}); starting fresh", flush=True)
            rows = []
            start_i = 0
    with open_giant_ephem(Path(kernel), Path(CACHE_JSON)) as ephem:
        for i, el in enumerate(clones):
            if i < start_i:
                continue
            r = integrate_clone(el, ephem, t_max_myr * 1e6)
            c = classify_bailey(r["a_au"], r["t_yr"], r["lifetime_yr"], rule="v2")
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
                f"  [{campaign}] {desig} {i+1}/{n} class={c['class']} life={c['lifetime_myr']:.2f}",
                flush=True,
            )
            partial_path.write_text(
                json.dumps({"desig": desig, "n_clones": n, "clones": rows}, indent=2),
                encoding="utf-8",
            )

    classes = [r["class"] for r in rows]
    counts = {k: classes.count(k) for k in ("D", "R", "Q")}
    modal = max(counts, key=counts.get)
    out = {
        "desig": desig,
        "campaign": campaign,
        "clone_scheme": "sbdb_covariance",
        "bailey_class": bailey["bailey_class"],
        "seed": SEED,
        "n_clones": n,
        "t_max_myr": t_max_myr,
        "counts": counts,
        "modal_class": modal,
        "elapsed_s": time.time() - t0,
        "clones": rows,
    }
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    if partial_path.exists():
        partial_path.unlink()
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--campaign", default="all", choices=list(CAMPAIGNS) + ["all"])
    ap.add_argument("--kernel-dir", default=KERNEL)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    camps = list(CAMPAIGNS) if args.campaign == "all" else [args.campaign]
    for camp in camps:
        desigs, n, tmax = CAMPAIGNS[camp]
        print(f"=== prefetch covariance {camp} ===", flush=True)
        prefetch(desigs)
        for des in desigs:
            print(f"=== {camp} {des} ===", flush=True)
            run_object(des, n, tmax, camp, args.kernel_dir, args.force)


if __name__ == "__main__":
    main()

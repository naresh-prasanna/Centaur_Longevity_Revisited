"""Resume TH173 100×40 Myr with clone-level parallelism (2 processes max)."""
from __future__ import annotations

import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SEED = 87
KERNEL = r"C:\Users\Hp\Downloads\Zenith"
CACHE_JSON = ROOT / "data" / "j2000_giants_horizons_de441.json"
OUT_DIR = ROOT / "results" / "fullscale" / "th173_100x40"
OUT_PATH = OUT_DIR / "2005_TH173.json"
PARTIAL = OUT_DIR / "2005_TH173.partial.json"
N_CLONES = 100
T_MAX_YR = 40.0e6
N_PROC = 2


def _clone_worker(payload: dict) -> dict:
    import sys

    sys.path.insert(0, payload["src_dir"])
    from classify_core import (  # noqa: WPS433
        classify_bailey,
        integrate_clone,
        open_giant_ephem,
        sample_clones_from_modern,
    )

    modern = payload["modern"]
    i = payload["clone"]
    clones = sample_clones_from_modern(modern, N_CLONES, seed=SEED)
    el = clones[i]
    t0 = time.time()
    with open_giant_ephem(Path(payload["kernel_dir"]), Path(payload["cache_json"])) as ephem:
        r = integrate_clone(el, ephem, T_MAX_YR)
    c = classify_bailey(r["a_au"], r["t_yr"], r["lifetime_yr"], rule="v2")
    return {
        "clone": i,
        "a0": el["a"],
        "e0": el["e"],
        "escaped": r["escaped"],
        "survived_full": r["survived_full"],
        "elapsed_s": time.time() - t0,
        **c,
    }


def load_done() -> dict[int, dict]:
    done: dict[int, dict] = {}
    for path in (OUT_PATH, PARTIAL):
        if not path.exists():
            continue
        blob = json.loads(path.read_text(encoding="utf-8"))
        for row in blob.get("clones") or []:
            done[int(row["clone"])] = row
    return done


def write_partial(rows: list[dict]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    PARTIAL.write_text(
        json.dumps(
            {
                "desig": "2005 TH173",
                "campaign": "th173_100x40",
                "n_clones": N_CLONES,
                "clones": sorted(rows, key=lambda r: r["clone"]),
            },
            indent=2,
        ),
        encoding="utf-8",
    )


def main() -> None:
    modern_all = json.loads((ROOT / "data" / "modern_elements.json").read_text(encoding="utf-8"))
    table = json.loads((ROOT / "data" / "bailey_table2.json").read_text(encoding="utf-8"))
    bailey = next(o for o in table["objects"] if o["desig"] == "2005 TH173")
    modern = modern_all["2005 TH173"]

    done = load_done()
    todo = [i for i in range(N_CLONES) if i not in done]
    print(f"TH173 100x40: {len(done)} done, {len(todo)} remaining, procs={N_PROC}", flush=True)
    if not todo:
        print("nothing to do", flush=True)
        return

    payloads = [
        {
            "clone": i,
            "modern": modern,
            "kernel_dir": KERNEL,
            "cache_json": str(CACHE_JSON),
            "src_dir": str(ROOT / "src"),
        }
        for i in todo
    ]
    rows = list(done.values())
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=N_PROC) as ex:
        futs = {ex.submit(_clone_worker, p): p["clone"] for p in payloads}
        for fut in as_completed(futs):
            row = fut.result()
            rows.append(row)
            write_partial(rows)
            print(
                f"  clone {row['clone']+1}/{N_CLONES} class={row['class']} "
                f"life={row['lifetime_myr']:.2f} ({row['elapsed_s']:.1f}s) "
                f"have={len(rows)}/{N_CLONES}",
                flush=True,
            )

    rows = sorted(rows, key=lambda r: r["clone"])
    classes = [r["class"] for r in rows]
    counts = {k: classes.count(k) for k in ("D", "R", "Q")}
    modal = max(counts, key=counts.get)
    out = {
        "desig": "2005 TH173",
        "campaign": "th173_100x40",
        "bailey": bailey,
        "bailey_class": bailey["bailey_class"],
        "seed": SEED,
        "n_clones": N_CLONES,
        "t_max_myr": 40.0,
        "counts": counts,
        "modal_class": modal,
        "frac_modal": counts[modal] / len(rows),
        "flip": modal != bailey["bailey_class"],
        "call": f"{bailey['bailey_class']}->{modal}",
        "elapsed_s": time.time() - t0,
        "clones": rows,
        "classifier_rule": "v2",
        "clone_parallel": N_PROC,
    }
    OUT_PATH.write_text(json.dumps(out, indent=2), encoding="utf-8")
    if PARTIAL.exists():
        PARTIAL.unlink()
    print(json.dumps({k: out[k] for k in out if k != "clones"}, indent=2), flush=True)


if __name__ == "__main__":
    main()

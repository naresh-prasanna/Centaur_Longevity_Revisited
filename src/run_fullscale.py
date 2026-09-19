"""
Full-scale Bailey Table 2 campaign (Icarus-grade compute).

Phases (resumable; skip completed objects):
  A  Q/R hardened     — 20 clones × 40 Myr
  B  D census         — 5 clones × 40 Myr (all BM09 D with modern SBDB)
  C  TH173 long       — 20 clones × 40 Myr (Q-capable window)
  D  Survivors 100 Myr — nominal (+ optional clones) for life≥40 at 40 Myr
  E  Sensitivity      — escape cuts + Exp C angles=zero on Q/R

Usage:
  python src/run_fullscale.py --phase A --workers 2
  python src/run_fullscale.py --phase all --workers 2
  python src/assemble_fullscale.py   # after phases finish
"""
from __future__ import annotations

import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
KERNEL_DEFAULT = r"C:\Users\Hp\Downloads\Zenith"
CACHE_JSON = ROOT / "data" / "j2000_giants_horizons_de441.json"
SEED = 87
FS = ROOT / "results" / "fullscale"


def slug(desig: str) -> str:
    return desig.replace(" ", "_").replace("/", "-")


def load_targets(classes: str):
    table = json.loads((ROOT / "data" / "bailey_table2.json").read_text(encoding="utf-8"))
    modern = json.loads((ROOT / "data" / "modern_elements.json").read_text(encoding="utf-8"))
    allowed = set(classes.upper())
    out = []
    for row in table["objects"]:
        if row["bailey_class"] == "NT":
            continue
        if row["bailey_class"] not in allowed:
            continue
        if row["desig"] not in modern or "a_au" not in modern[row["desig"]]:
            continue
        out.append((row, modern[row["desig"]]))
    return out


def _run_object_worker(payload: dict) -> dict:
    """Process-safe worker: one object, own ephemeris handle."""
    import sys

    sys.path.insert(0, str(Path(payload["src_dir"])))
    from classify_core import (  # noqa: WPS433
        classify_bailey,
        integrate_clone,
        open_giant_ephem,
        sample_clones_from_modern,
    )

    desig = payload["desig"]
    bailey = payload["bailey"]
    modern = payload["modern"]
    n_clones = payload["n_clones"]
    t_max_yr = payload["t_max_yr"]
    out_path = Path(payload["out_path"])
    kernel = Path(payload["kernel_dir"])
    escape = payload.get("escape") or {}

    if out_path.exists() and not payload.get("force"):
        prev = json.loads(out_path.read_text(encoding="utf-8"))
        if (
            prev.get("n_clones") == n_clones
            and abs(float(prev.get("t_max_myr", 0)) - t_max_yr / 1e6) < 1e-9
            and prev.get("campaign") == payload.get("campaign")
        ):
            return {"desig": desig, "skipped": True, "path": str(out_path), "call": prev.get("call")}

    clones = sample_clones_from_modern(modern, n_clones, seed=SEED)
    rows = []
    t0 = time.time()
    partial_path = out_path.with_suffix(".partial.json")
    start_i = 0
    if partial_path.exists() and not payload.get("force"):
        try:
            prev_p = json.loads(partial_path.read_text(encoding="utf-8"))
            if (
                prev_p.get("desig") == desig
                and prev_p.get("n_clones") == n_clones
                and prev_p.get("campaign") == payload.get("campaign")
            ):
                rows = list(prev_p.get("clones") or [])
                start_i = len(rows)
                print(f"  RESUME {desig} from clone {start_i}/{n_clones}", flush=True)
        except Exception as exc:
            print(f"  partial load failed ({exc}); starting fresh", flush=True)
            rows = []
            start_i = 0
    cache_path = Path(payload.get("cache_json") or CACHE_JSON)
    with open_giant_ephem(kernel, cache_path) as ephem:
        for i, el in enumerate(clones):
            if i < start_i:
                continue
            r = integrate_clone(
                el,
                ephem,
                t_max_yr,
                r_escape_au=float(escape.get("r_escape_au", 1e4)),
                q_min_au=float(escape.get("q_min_au", 2.5)),
                use_hill=bool(escape.get("use_hill", True)),
            )
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
                f"  [{payload['campaign']}] {desig} {i+1}/{n_clones} "
                f"class={c['class']} life={c['lifetime_myr']:.2f}",
                flush=True,
            )
            out_path.parent.mkdir(parents=True, exist_ok=True)
            partial_path.write_text(
                json.dumps(
                    {
                        "desig": desig,
                        "campaign": payload["campaign"],
                        "n_clones": n_clones,
                        "clones": rows,
                    },
                    indent=2,
                ),
                encoding="utf-8",
            )

    classes = [r["class"] for r in rows]
    counts = {k: classes.count(k) for k in ("D", "R", "Q")}
    modal = max(counts, key=counts.get)
    bailey_class = bailey["bailey_class"]
    flip = modal != bailey_class
    out = {
        "desig": desig,
        "campaign": payload["campaign"],
        "bailey": bailey,
        "modern_nominal": {
            k: modern.get(k)
            for k in (
                "a_au",
                "e",
                "i_deg",
                "om_deg",
                "w_deg",
                "ma_deg",
                "q_au",
                "condition_code",
                "data_arc_days",
                "n_obs",
            )
        },
        "seed": SEED,
        "n_clones": n_clones,
        "t_max_myr": t_max_yr / 1.0e6,
        "escape": escape or {"r_escape_au": 1e4, "q_min_au": 2.5, "use_hill": True},
        "counts": counts,
        "modal_class": modal,
        "frac_modal": counts[modal] / len(rows),
        "bailey_class": bailey_class,
        "flip": bool(flip),
        "call": f"{bailey_class}->{modal}" if flip else f"stay {modal}",
        "elapsed_s": time.time() - t0,
        "clones": rows,
        "classifier_rule": "v2",
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    if partial_path.exists():
        partial_path.unlink()
    return {"desig": desig, "skipped": False, "path": str(out_path), "call": out["call"], "counts": counts}


def run_phase_batch(
    campaign: str,
    classes: str,
    n_clones: int,
    t_max_myr: float,
    workers: int,
    kernel_dir: str,
    force: bool = False,
    escape: dict | None = None,
    only: set[str] | None = None,
):
    out_dir = FS / campaign
    out_dir.mkdir(parents=True, exist_ok=True)
    targets = load_targets(classes)
    if only:
        targets = [(b, m) for b, m in targets if b["desig"] in only]

    payloads = []
    for bailey, modern in targets:
        des = bailey["desig"]
        payloads.append(
            {
                "desig": des,
                "bailey": bailey,
                "modern": modern,
                "n_clones": n_clones,
                "t_max_yr": t_max_myr * 1e6,
                "out_path": str(out_dir / f"{slug(des)}.json"),
                "kernel_dir": kernel_dir,
                "cache_json": str(CACHE_JSON),
                "src_dir": str(ROOT / "src"),
                "campaign": campaign,
                "force": force,
                "escape": escape or {},
            }
        )

    print(f"=== PHASE {campaign}: {len(payloads)} objects, {n_clones} clones, {t_max_myr} Myr, workers={workers}", flush=True)
    results = []
    if workers <= 1:
        for p in payloads:
            results.append(_run_object_worker(p))
            print(f"  -> {results[-1]['desig']} {results[-1].get('call')} skip={results[-1].get('skipped')}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = {ex.submit(_run_object_worker, p): p["desig"] for p in payloads}
            for fut in as_completed(futs):
                r = fut.result()
                results.append(r)
                print(f"  -> {r['desig']} {r.get('call')} skip={r.get('skipped')}", flush=True)

    # summary
    done = []
    for r in results:
        path = Path(r["path"])
        if path.exists():
            done.append(json.loads(path.read_text(encoding="utf-8")))
    summary = {
        "campaign": campaign,
        "n": len(done),
        "n_clones": n_clones,
        "t_max_myr": t_max_myr,
        "n_flip": sum(1 for d in done if d.get("flip")),
        "flip_rate": (sum(1 for d in done if d.get("flip")) / len(done)) if done else 0.0,
        "objects": [
            {
                "desig": d["desig"],
                "bailey": d["bailey_class"],
                "modal": d["modal_class"],
                "flip": d["flip"],
                "counts": d["counts"],
                "median_life": sorted(c["lifetime_myr"] for c in d["clones"])[len(d["clones"]) // 2],
                "frac_long": sum(1 for c in d["clones"] if c["lifetime_myr"] >= 22) / len(d["clones"]),
            }
            for d in done
        ],
    }
    (FS / f"summary_{campaign}.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({"campaign": campaign, "n": summary["n"], "flip_rate": summary["flip_rate"]}, indent=2), flush=True)
    return summary


def run_survivors_100(workers: int, kernel_dir: str, force: bool = False):
    """Extend objects that reached ~40 Myr in qr20 / th173_40."""
    survivors = set()
    for camp in ("qr20", "th173_40", "d5"):
        d = FS / camp
        if not d.exists():
            continue
        for f in d.glob("*.json"):
            obj = json.loads(f.read_text(encoding="utf-8"))
            if any(c.get("lifetime_myr", 0) >= 39.5 for c in obj.get("clones", [])):
                survivors.add(obj["desig"])
    # always include known Exp C survivor
    survivors.update({"2003 QP112", "2000 FZ53", "2006 SX368", "2005 TH173", "1995 SN55"})
    print(f"Survivors for 100 Myr: {sorted(survivors)}", flush=True)
    return run_phase_batch(
        campaign="surv100",
        classes="DRQ",
        n_clones=5,
        t_max_myr=100.0,
        workers=workers,
        kernel_dir=kernel_dir,
        force=force,
        only=survivors,
    )


def run_sensitivity_escape(workers: int, kernel_dir: str, force: bool = False):
    """Q/R nominal escape variants: loose q and no-Hill."""
    only = {b["desig"] for b, _ in load_targets("QR")}
    run_phase_batch(
        "sens_q15",
        "QR",
        n_clones=3,
        t_max_myr=40.0,
        workers=workers,
        kernel_dir=kernel_dir,
        force=force,
        escape={"r_escape_au": 1e4, "q_min_au": 1.5, "use_hill": True},
        only=only,
    )
    run_phase_batch(
        "sens_nohill",
        "QR",
        n_clones=3,
        t_max_myr=40.0,
        workers=workers,
        kernel_dir=kernel_dir,
        force=force,
        escape={"r_escape_au": 1e4, "q_min_au": 2.5, "use_hill": False},
        only=only,
    )


def run_exp_c_angles_zero(kernel_dir: str, force: bool = False):
    """Exp C with zero angles for Bailey aei (sensitivity)."""
    import sys

    sys.path.insert(0, str(ROOT / "src"))
    from run_exp_c import elems_bailey_aei, elems_modern, run_one  # noqa: WPS433
    from classify_core import DE441Ephem  # noqa: WPS433

    out_dir = FS / "exp_c_angles_zero"
    out_dir.mkdir(parents=True, exist_ok=True)
    table = json.loads((ROOT / "data" / "bailey_table2.json").read_text(encoding="utf-8"))
    modern_all = json.loads((ROOT / "data" / "modern_elements.json").read_text(encoding="utf-8"))
    results = []
    with DE441Ephem(Path(kernel_dir)) as ephem:
        for row in table["objects"]:
            if row["bailey_class"] not in ("Q", "R"):
                continue
            des = row["desig"]
            if des not in modern_all:
                continue
            path = out_dir / f"{slug(des)}.json"
            if path.exists() and not force:
                results.append(json.loads(path.read_text(encoding="utf-8")))
                continue
            modern = modern_all[des]
            print(f"=== Exp C angles=zero {des} ===", flush=True)
            b = run_one(elems_bailey_aei(row, modern, angles="zero"), ephem, 40e6)
            m = run_one(elems_modern(modern), ephem, 40e6)
            rec = {
                "desig": des,
                "bailey_class": row["bailey_class"],
                "angles_mode": "zero",
                "bailey_ic": {k: b[k] for k in b if k not in ("t_yr", "a_au")},
                "modern_ic": {k: m[k] for k in m if k not in ("t_yr", "a_au")},
                "same_class": b["class"] == m["class"],
                "orbit_drives_flip": (b["class"] == row["bailey_class"]) and (m["class"] != row["bailey_class"]),
            }
            path.write_text(json.dumps(rec, indent=2), encoding="utf-8")
            results.append(rec)
            print(f"  v2 B/M = {b['class']}/{m['class']} life={b['lifetime_myr']:.1f}/{m['lifetime_myr']:.1f}", flush=True)
    summary = {
        "campaign": "exp_c_angles_zero",
        "n": len(results),
        "n_orbit_drives": sum(1 for r in results if r.get("orbit_drives_flip")),
        "objects": results,
    }
    (FS / "summary_exp_c_angles_zero.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--phase",
        default="all",
        help="A|B|C|D|E|all or qr20|d5|d20|th173_40|th173_100x40|surv100|sens|expczero",
    )
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--kernel-dir", default=KERNEL_DEFAULT)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    FS.mkdir(parents=True, exist_ok=True)

    phase = args.phase.upper() if args.phase.lower() not in (
        "qr20", "d5", "d20", "th173_40", "th173_100x40", "surv100", "sens", "expczero", "all"
    ) else args.phase.lower()

    def do_A():
        return run_phase_batch("qr20", "QR", 20, 40.0, args.workers, args.kernel_dir, args.force)

    def do_B():
        return run_phase_batch("d5", "D", 5, 40.0, args.workers, args.kernel_dir, args.force)

    def do_C():
        return run_phase_batch(
            "th173_40",
            "Q",
            20,
            40.0,
            args.workers,
            args.kernel_dir,
            args.force,
            only={"2005 TH173"},
        )

    def do_D():
        return run_survivors_100(args.workers, args.kernel_dir, args.force)

    def do_E():
        run_sensitivity_escape(args.workers, args.kernel_dir, args.force)
        return run_exp_c_angles_zero(args.kernel_dir, args.force)

    def do_d20():
        return run_phase_batch("d20", "D", 20, 40.0, args.workers, args.kernel_dir, args.force)

    def do_th173_100x40():
        return run_phase_batch(
            "th173_100x40",
            "Q",
            100,
            40.0,
            args.workers,
            args.kernel_dir,
            args.force,
            only={"2005 TH173"},
        )

    mapping = {
        "A": do_A,
        "qr20": do_A,
        "B": do_B,
        "d5": do_B,
        "d20": do_d20,
        "C": do_C,
        "th173_40": do_C,
        "th173_100x40": do_th173_100x40,
        "D": do_D,
        "surv100": do_D,
        "E": do_E,
        "sens": do_E,
        "expczero": lambda: run_exp_c_angles_zero(args.kernel_dir, args.force),
    }

    if phase == "all":
        for name, fn in [
            ("A/qr20", do_A),
            ("B/d5", do_B),
            ("d20", do_d20),
            ("C/th173_40", do_C),
            ("th173_100x40", do_th173_100x40),
            ("D/surv100", do_D),
            ("E/sens", do_E),
        ]:
            print(f"\n######## FULLSCALE {name} ########\n", flush=True)
            fn()
        print("FULLSCALE ALL PHASES COMPLETE", flush=True)
    elif phase in mapping:
        mapping[phase]()
    else:
        raise SystemExit(f"Unknown phase {args.phase}")


if __name__ == "__main__":
    main()

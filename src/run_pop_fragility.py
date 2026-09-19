"""Population fragility survey orchestration.

Reuses classify_core / covariance_clones as-is. Clone jobs are independent.
Resume-safe: completed object JSON is skipped; partial clone files resume.

Usage:
  python src/run_pop_fragility.py --phase catalog
  python src/run_pop_fragility.py --phase pilot
  python src/run_pop_fragility.py --phase all
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
SEED = 87
KERNEL = r"C:\Users\Hp\Downloads\Zenith"
CACHE_JSON = ROOT / "data" / "j2000_giants_horizons_de441.json"
OUT = ROOT / "results" / "pop_fragility"
BUDGET_S = 14 * 3600
PHASE1_BUDGET_S = 9 * 3600
PHASE2_BUDGET_S = 3 * 3600
PHASE3_BUDGET_S = 2 * 3600
PILOT_N_CLONES = 5
PILOT_TMAX_MYR = 25.0
SURVEY_N_CLONES = 10
SURVEY_TMAX_MYR = 25.0
SURVEY_N_TARGET = 175
COV_N_OBJECTS = 18
CENSOR_TMAX_MYR = 40.0


def slug(desig: str) -> str:
    return desig.replace(" ", "_").replace("/", "-")


def default_workers() -> int:
    n = os.cpu_count() or 2
    # Leave one logical core for the OS. Cap at 6: this laptop previously
    # page-faulted when many REBOUND processes ran together.
    cap = 6
    try:
        import ctypes

        class MEMORYSTATUSEX(ctypes.Structure):
            _fields_ = [
                ("dwLength", ctypes.c_ulong),
                ("dwMemoryLoad", ctypes.c_ulong),
                ("ullTotalPhys", ctypes.c_ulonglong),
                ("ullAvailPhys", ctypes.c_ulonglong),
                ("ullTotalPageFile", ctypes.c_ulonglong),
                ("ullAvailPageFile", ctypes.c_ulonglong),
                ("ullTotalVirtual", ctypes.c_ulonglong),
                ("ullAvailVirtual", ctypes.c_ulonglong),
                ("sullAvailExtendedVirtual", ctypes.c_ulonglong),
            ]

        stat = MEMORYSTATUSEX()
        stat.dwLength = ctypes.sizeof(MEMORYSTATUSEX)
        ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(stat))
        gb = stat.ullTotalPhys / (1024 ** 3)
        avail = stat.ullAvailPhys / (1024 ** 3)
    except Exception:
        gb, avail = 8.0, 4.0
    if gb < 8.5 or avail < 3.0:
        cap = 2
    workers = max(1, min(n - 1, cap))
    return workers


def machine_info(workers: int) -> dict:
    n = os.cpu_count() or 2
    info = {"cpu_count": n, "workers": workers}
    try:
        import ctypes

        class MEMORYSTATUSEX(ctypes.Structure):
            _fields_ = [
                ("dwLength", ctypes.c_ulong),
                ("dwMemoryLoad", ctypes.c_ulong),
                ("ullTotalPhys", ctypes.c_ulonglong),
                ("ullAvailPhys", ctypes.c_ulonglong),
                ("ullTotalPageFile", ctypes.c_ulonglong),
                ("ullAvailPageFile", ctypes.c_ulonglong),
                ("ullTotalVirtual", ctypes.c_ulonglong),
                ("ullAvailVirtual", ctypes.c_ulonglong),
                ("sullAvailExtendedVirtual", ctypes.c_ulonglong),
            ]

        stat = MEMORYSTATUSEX()
        stat.dwLength = ctypes.sizeof(MEMORYSTATUSEX)
        ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(stat))
        info["ram_gb"] = round(stat.ullTotalPhys / (1024 ** 3), 2)
        info["ram_avail_gb"] = round(stat.ullAvailPhys / (1024 ** 3), 2)
    except Exception:
        info["ram_gb"] = None
    return info


def write_status(payload: dict) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "status.json"
    prev = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    prev.update(payload)
    prev["updated_unix"] = time.time()
    path.write_text(json.dumps(prev, indent=2), encoding="utf-8")


def _clone_worker(payload: dict) -> dict:
    import sys

    sys.path.insert(0, payload["src_dir"])
    from classify_core import classify_bailey, integrate_clone, open_giant_ephem, sample_clones_from_modern

    modern = payload["modern"]
    i = int(payload["clone"])
    n = int(payload["n_clones"])
    seed = int(payload["seed"])
    t_max_yr = float(payload["t_max_yr"])
    scheme = payload.get("clone_scheme", "diagonal")
    if scheme == "covariance":
        from covariance_clones import sample_clones_from_covariance

        clones = sample_clones_from_covariance(payload["desig"], n, seed=seed)
    else:
        clones = sample_clones_from_modern(modern, n, seed=seed)
    el = clones[i]
    t0 = time.time()
    with open_giant_ephem(Path(payload["kernel_dir"]), Path(payload["cache_json"])) as ephem:
        r = integrate_clone(el, ephem, t_max_yr)
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


def load_object_clones(out_path: Path, partial_path: Path) -> dict[int, dict]:
    done: dict[int, dict] = {}
    for path in (out_path, partial_path):
        if not path.exists():
            continue
        blob = json.loads(path.read_text(encoding="utf-8"))
        for row in blob.get("clones") or []:
            done[int(row["clone"])] = row
    return done


def assemble_object(
    desig: str,
    modern: dict,
    campaign: str,
    clone_scheme: str,
    n_clones: int,
    t_max_myr: float,
    rows: list[dict],
    elapsed_s: float,
    extra: dict | None = None,
) -> dict:
    classes = [r["class"] for r in rows]
    counts = {k: classes.count(k) for k in ("D", "R", "Q")}
    modal = max(counts, key=counts.get)
    out = {
        "desig": desig,
        "campaign": campaign,
        "clone_scheme": clone_scheme,
        "modern_nominal": modern,
        "seed": SEED,
        "n_clones": n_clones,
        "t_max_myr": t_max_myr,
        "counts": counts,
        "modal_class": modal,
        "frac_modal": counts[modal] / len(rows) if rows else None,
        "elapsed_s": elapsed_s,
        "clones": sorted(rows, key=lambda r: r["clone"]),
        "classifier_rule": "v2",
    }
    if extra:
        out.update(extra)
    return out


def run_object_clones(
    modern: dict,
    campaign: str,
    n_clones: int,
    t_max_myr: float,
    workers: int,
    clone_scheme: str = "diagonal",
    clone_indices: list[int] | None = None,
    deadline: float | None = None,
) -> dict:
    desig = modern["desig"]
    out_dir = OUT / campaign
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{slug(desig)}.json"
    partial_path = out_dir / f"{slug(desig)}.partial.json"
    want = list(range(n_clones) if clone_indices is None else clone_indices)
    done = load_object_clones(out_path, partial_path)
    # Completed full object with matching config
    if out_path.exists() and clone_indices is None:
        prev = json.loads(out_path.read_text(encoding="utf-8"))
        if (
            prev.get("n_clones") == n_clones
            and abs(float(prev.get("t_max_myr", 0)) - t_max_myr) < 1e-9
            and prev.get("clone_scheme") == clone_scheme
            and prev.get("campaign") == campaign
            and len(prev.get("clones") or []) == n_clones
        ):
            return {"desig": desig, "skipped": True, "path": str(out_path), "blob": prev}

    todo = [i for i in want if i not in done]
    rows = [done[i] for i in done if i in set(want) or clone_indices is None]
    t0 = time.time()
    if todo:
        payloads = [
            {
                "clone": i,
                "desig": desig,
                "modern": modern,
                "n_clones": n_clones,
                "t_max_yr": t_max_myr * 1.0e6,
                "seed": SEED,
                "clone_scheme": clone_scheme,
                "kernel_dir": KERNEL,
                "cache_json": str(CACHE_JSON),
                "src_dir": str(SRC),
            }
            for i in todo
        ]
        n_proc = max(1, min(workers, len(payloads)))
        with ProcessPoolExecutor(max_workers=n_proc) as ex:
            futs = {ex.submit(_clone_worker, p): p["clone"] for p in payloads}
            for fut in as_completed(futs):
                if deadline is not None and time.time() > deadline:
                    for f in futs:
                        f.cancel()
                    raise TimeoutError(f"deadline hit during {desig}")
                row = fut.result()
                rows.append(row)
                partial_path.write_text(
                    json.dumps(
                        {
                            "desig": desig,
                            "campaign": campaign,
                            "n_clones": n_clones,
                            "clones": sorted(rows, key=lambda r: r["clone"]),
                        },
                        indent=2,
                    ),
                    encoding="utf-8",
                )
                print(
                    f"  [{campaign}] {desig} clone {row['clone']+1}/{n_clones} "
                    f"class={row['class']} life={row['lifetime_myr']:.2f} "
                    f"({row['elapsed_s']:.1f}s) have={len(rows)}/{n_clones}",
                    flush=True,
                )
    rows = sorted(rows, key=lambda r: r["clone"])
    blob = assemble_object(
        desig, modern, campaign, clone_scheme, n_clones, t_max_myr, rows, time.time() - t0
    )
    out_path.write_text(json.dumps(blob, indent=2), encoding="utf-8")
    if partial_path.exists():
        partial_path.unlink()
    return {"desig": desig, "skipped": False, "path": str(out_path), "blob": blob}


def load_campaign_rows(campaign: str) -> list[dict]:
    from pop_survey import object_summary

    d = OUT / campaign
    if not d.exists():
        return []
    rows = []
    for path in sorted(d.glob("*.json")):
        if path.name.endswith(".partial.json") or path.name.startswith("_"):
            continue
        blob = json.loads(path.read_text(encoding="utf-8"))
        if "clones" not in blob:
            continue
        rows.append(object_summary(blob))
    return rows


def phase_catalog(n_target: int) -> dict:
    from pop_survey import build_stratified_sample, fetch_cen_catalog

    print("fetching SBDB CEN catalog (OD fields)...", flush=True)
    cat = fetch_cen_catalog()
    print(f"catalog n={cat['n']}", flush=True)
    sample = build_stratified_sample(cat, n_target=n_target, seed=SEED)
    print(
        f"stratified sample n={len(sample['objects'])} eligible={sample['n_eligible']} "
        f"strata={sample['n_strata']}",
        flush=True,
    )
    return sample


def phase_pilot(sample: dict, workers: int) -> dict:
    from pop_survey import pick_pilot_objects

    objs = pick_pilot_objects(sample, n=5)
    print("PILOT objects:", [f"{o['desig']} arc={o['data_arc_days']:.0f}d U={o.get('condition_code')}" for o in objs], flush=True)
    t0 = time.time()
    clone_times = []
    clone_myr_times = []
    for o in objs:
        rec = run_object_clones(
            o, "pilot", PILOT_N_CLONES, PILOT_TMAX_MYR, workers, clone_scheme="diagonal"
        )
        for c in rec["blob"]["clones"]:
            dt = float(c.get("elapsed_s") or 0.0)
            life = max(float(c.get("lifetime_myr") or 0.05), 0.05)
            clone_times.append(dt)
            clone_myr_times.append(dt / life)
    wall = time.time() - t0
    mean_s = float(sum(clone_times) / len(clone_times)) if clone_times else None
    med_s = float(sorted(clone_times)[len(clone_times) // 2]) if clone_times else None
    mean_s_per_myr = float(sum(clone_myr_times) / len(clone_myr_times)) if clone_myr_times else None
    report = {
        "phase": "pilot",
        "objects": [o["desig"] for o in objs],
        "n_objects": len(objs),
        "n_clones": PILOT_N_CLONES,
        "t_max_myr": PILOT_TMAX_MYR,
        "workers": workers,
        "wall_s": wall,
        "n_clone_jobs": len(clone_times),
        "mean_s_per_clone": mean_s,
        "median_s_per_clone": med_s,
        "mean_s_per_clone_myr": mean_s_per_myr,
        "clone_elapsed_s": clone_times,
    }
    (OUT / "pilot_timing.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def project_budget(pilot: dict, n_obj: int, n_clones: int, t_max: float, workers: int) -> dict:
    # Early-escape clones are cheaper than Tmax. Use empirical mean wall per clone
    # from the 25 Myr pilot, then scale mildly with Tmax for longer windows.
    mean_s = float(pilot["mean_s_per_clone"])
    t_ref = float(pilot["t_max_myr"])
    scale = t_max / t_ref
    # Not linear: most clones die early. Blend 70% unscaled + 30% Tmax-scaled.
    s_per_clone = mean_s * (0.70 + 0.30 * scale)
    n_jobs = n_obj * n_clones
    wall_s = n_jobs * s_per_clone / max(workers, 1)
    return {
        "n_obj": n_obj,
        "n_clones": n_clones,
        "t_max_myr": t_max,
        "workers": workers,
        "s_per_clone_assumed": s_per_clone,
        "n_jobs": n_jobs,
        "wall_h": wall_s / 3600.0,
        "wall_s": wall_s,
    }


def rescale_phase1(pilot: dict, workers: int) -> dict:
    """Fit Phase 1 into ~8.5 h if the 175×10×25 projection is too large."""
    proj = project_budget(pilot, SURVEY_N_TARGET, SURVEY_N_CLONES, SURVEY_TMAX_MYR, workers)
    n = SURVEY_N_TARGET
    t_max = SURVEY_TMAX_MYR
    n_clones = SURVEY_N_CLONES
    note = "as planned"
    if proj["wall_h"] > 8.5:
        # Prefer dropping objects, keep Tmax=25 and n=10 so appendix comparison holds.
        n = max(80, int(SURVEY_N_TARGET * 8.5 / proj["wall_h"]))
        proj = project_budget(pilot, n, n_clones, t_max, workers)
        note = f"rescaled n_objects {SURVEY_N_TARGET}->{n} to fit ~8.5h Phase 1"
        if proj["wall_h"] > 8.5:
            n = max(60, int(n * 8.5 / proj["wall_h"]))
            proj = project_budget(pilot, n, n_clones, t_max, workers)
            note = f"rescaled n_objects {SURVEY_N_TARGET}->{n} to fit ~8.5h Phase 1"
    proj["note"] = note
    proj["n_obj"] = n
    proj["n_clones"] = n_clones
    proj["t_max_myr"] = t_max
    return proj


def phase_survey(sample: dict, n_obj: int, workers: int, deadline: float) -> dict:
    objs = sample["objects"][:n_obj]
    t0 = time.time()
    done_n = 0
    for k, o in enumerate(objs, start=1):
        left = deadline - time.time()
        if left < 90:
            print(f"STOP Phase 1: {left:.0f}s left before {o['desig']}", flush=True)
            break
        elapsed = time.time() - t0
        rate = elapsed / max(done_n, 1) if done_n else None
        remain = (len(objs) - k + 1) * rate if rate else None
        print(
            f"Phase1 {k}/{len(objs)} {o['desig']} arc={o['data_arc_days']:.0f}d "
            f"elapsed={elapsed/3600:.2f}h proj_remain={(remain or 0)/3600:.2f}h",
            flush=True,
        )
        run_object_clones(o, "phase1", SURVEY_N_CLONES, SURVEY_TMAX_MYR, workers)
        done_n += 1
        write_status(
            {
                "phase": "phase1",
                "object_k": k,
                "object_n": len(objs),
                "desig": o["desig"],
                "elapsed_h": elapsed / 3600.0,
            }
        )
    rows = load_campaign_rows("phase1")
    return {"n_requested": len(objs), "n_completed": len(rows), "wall_s": time.time() - t0}


def phase_covariance(phase1_rows: list[dict], sample: dict, workers: int, deadline: float) -> dict:
    by_desig = {o["desig"]: o for o in sample["objects"]}
    short = sorted(
        [r for r in phase1_rows if r.get("data_arc_days")],
        key=lambda r: r["data_arc_days"],
    )[:COV_N_OBJECTS]
    from covariance_clones import prefetch

    print("prefetch SBDB covariances...", flush=True)
    try:
        prefetch([r["desig"] for r in short])
    except Exception as exc:  # noqa: BLE001
        print(f"prefetch warning: {exc}", flush=True)
    t0 = time.time()
    n_ok = 0
    n_fail = 0
    fails = {}
    for k, row in enumerate(short, start=1):
        if time.time() > deadline - 90:
            print("STOP Phase 2: deadline", flush=True)
            break
        modern = by_desig[row["desig"]]
        print(f"Phase2 {k}/{len(short)} {modern['desig']} arc={modern['data_arc_days']:.0f}d", flush=True)
        try:
            run_object_clones(
                modern, "phase2_cov", SURVEY_N_CLONES, SURVEY_TMAX_MYR, workers, clone_scheme="covariance"
            )
            n_ok += 1
        except Exception as exc:  # noqa: BLE001
            n_fail += 1
            fails[modern["desig"]] = str(exc)
            print(f"  COV FAIL {modern['desig']}: {exc}", flush=True)
    return {
        "n_requested": len(short),
        "n_ok": n_ok,
        "n_fail": n_fail,
        "fails": fails,
        "wall_s": time.time() - t0,
        "targets": [r["desig"] for r in short],
    }


def phase_censor(phase1_rows: list[dict], sample: dict, workers: int, deadline: float) -> dict:
    by_desig = {o["desig"]: o for o in sample["objects"]}
    t0 = time.time()
    n_obj = 0
    n_clones_ext = 0
    for row in phase1_rows:
        if time.time() > deadline - 90:
            print("STOP Phase 3: deadline", flush=True)
            break
        if not row.get("censoring_count"):
            continue
        path = OUT / "phase1" / f"{slug(row['desig'])}.json"
        blob = json.loads(path.read_text(encoding="utf-8"))
        cens = [int(c["clone"]) for c in blob["clones"] if c.get("survived_full")]
        if not cens:
            continue
        modern = by_desig[row["desig"]]
        print(f"Phase3 {row['desig']} extending {len(cens)} censored clones to 40 Myr", flush=True)
        # Run only censored clone indices at 40 Myr; merge onto a copy of the 25 Myr file.
        ext_dir = OUT / "phase3_40"
        ext_dir.mkdir(parents=True, exist_ok=True)
        ext_path = ext_dir / f"{slug(row['desig'])}.json"
        rec = run_object_clones(
            modern,
            "phase3_40",
            SURVEY_N_CLONES,
            CENSOR_TMAX_MYR,
            workers,
            clone_scheme="diagonal",
            clone_indices=cens,
            deadline=deadline,
        )
        n_obj += 1
        n_clones_ext += len(cens)
        # Merge: keep 25 Myr rows for uncensored clones; replace extended clones.
        merged_clones = {int(c["clone"]): c for c in blob["clones"]}
        for c in rec["blob"]["clones"]:
            if int(c["clone"]) in set(cens):
                merged_clones[int(c["clone"])] = c
        merged = assemble_object(
            row["desig"],
            modern,
            "phase3_merged",
            "diagonal",
            SURVEY_N_CLONES,
            CENSOR_TMAX_MYR,
            list(merged_clones.values()),
            rec["blob"].get("elapsed_s") or 0.0,
            extra={"extended_clones": cens, "base_campaign": "phase1"},
        )
        (OUT / "phase3_merged").mkdir(parents=True, exist_ok=True)
        (OUT / "phase3_merged" / f"{slug(row['desig'])}.json").write_text(
            json.dumps(merged, indent=2), encoding="utf-8"
        )
        _ = ext_path
    return {"n_objects_extended": n_obj, "n_clones_extended": n_clones_ext, "wall_s": time.time() - t0}


def write_report(meta: dict) -> Path:
    from pop_survey import make_plots, spearman_table, write_summary_csv

    p1 = load_campaign_rows("phase1")
    p2 = load_campaign_rows("phase2_cov")
    p3 = load_campaign_rows("phase3_merged")
    spear = spearman_table(p1) if p1 else []
    write_summary_csv(p1, OUT / "summary_phase1.csv")
    if p2:
        write_summary_csv(p2, OUT / "summary_phase2_cov.csv")
    if p3:
        write_summary_csv(p3, OUT / "summary_phase3_merged.csv")
    (OUT / "spearman_phase1.json").write_text(json.dumps(spear, indent=2), encoding="utf-8")
    figs = make_plots(p1, p2, OUT / "figures") if p1 else []

    # covariance comparison table
    cov_cmp = []
    if p1 and p2:
        dmap = {r["desig"]: r for r in p1}
        for c in p2:
            d = dmap.get(c["desig"])
            if not d:
                continue
            cov_cmp.append(
                {
                    "desig": c["desig"],
                    "arc_days": d.get("data_arc_days"),
                    "diag_rq": d.get("minority_rq_frac"),
                    "cov_rq": c.get("minority_rq_frac"),
                    "delta_cov_minus_diag": (c.get("minority_rq_frac") or 0) - (d.get("minority_rq_frac") or 0),
                    "diag_modal": d.get("modal_class"),
                    "cov_modal": c.get("modal_class"),
                    "diag_censor": d.get("censoring_rate"),
                    "cov_censor": c.get("censoring_rate"),
                }
            )
        write_summary_csv(cov_cmp, OUT / "summary_diag_vs_cov.csv")
        (OUT / "diag_vs_cov.json").write_text(json.dumps(cov_cmp, indent=2), encoding="utf-8")

    lines = [
        "# Population fragility survey — run report",
        "",
        "Honest record of what ran. Numbers are from archived JSON, not invented.",
        "",
        "## Machine",
        json.dumps(meta.get("machine"), indent=2),
        "",
        "## Pilot",
        json.dumps(meta.get("pilot"), indent=2),
        "",
        "## Projection after pilot",
        json.dumps(meta.get("projection"), indent=2),
        "",
        "## What actually completed",
        f"- Phase 1 objects: {len(p1)} (requested {meta.get('phase1_requested')})",
        f"- Phase 2 covariance objects: {len(p2)}",
        f"- Phase 3 objects with 40 Myr extensions: {len(p3)}",
        f"- Phase 1 wall_h: {meta.get('phase1_wall_h')}",
        f"- Phase 2 wall_h: {meta.get('phase2_wall_h')}",
        f"- Phase 3 wall_h: {meta.get('phase3_wall_h')}",
        f"- Deviations: {meta.get('deviations')}",
        "",
        "## Spearman (Phase 1)",
    ]
    for row in spear:
        lines.append(
            f"- {row['predictor']} vs {row['outcome']}: "
            f"rho={row.get('rho')} p={row.get('p')} n={row.get('n')}"
        )
    if cov_cmp:
        deltas = [r["delta_cov_minus_diag"] for r in cov_cmp]
        mean_d = sum(deltas) / len(deltas)
        lines += [
            "",
            "## Diagonal vs covariance (short-arc)",
            f"- n={len(cov_cmp)}; mean (cov − diag) R/Q fraction = {mean_d:.4f}",
            "- Positive mean: covariance yields more R/Q clones than diagonal.",
            "- Negative mean: diagonal overestimates the R/Q minority.",
        ]
    lines += ["", "## Figures"]
    for p in figs:
        lines.append(f"- `{p.as_posix()}`")
    path = OUT / "RUN_REPORT.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--phase",
        default="all",
        choices=["catalog", "pilot", "survey", "cov", "censor", "analyze", "all"],
    )
    parser.add_argument("--workers", type=int, default=0)
    parser.add_argument("--n-target", type=int, default=SURVEY_N_TARGET)
    parser.add_argument("--budget-h", type=float, default=14.0)
    args = parser.parse_args()
    workers = args.workers if args.workers > 0 else default_workers()
    machine = machine_info(workers)
    print("machine", json.dumps(machine), flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    write_status({"phase": args.phase, "machine": machine})
    deadline_all = time.time() + args.budget_h * 3600
    meta: dict = {"machine": machine, "deviations": []}

    if args.phase in ("catalog", "all", "pilot", "survey"):
        sample = phase_catalog(args.n_target)
        meta["n_catalog"] = sample["n_catalog"]
        meta["n_eligible"] = sample["n_eligible"]
    else:
        from pop_survey import SAMPLE_PATH

        sample = json.loads(SAMPLE_PATH.read_text(encoding="utf-8"))

    if args.phase == "catalog":
        return

    if args.phase in ("pilot", "all"):
        print("=== PHASE 0 TIMING PILOT ===", flush=True)
        pilot = phase_pilot(sample, workers)
        meta["pilot"] = {k: v for k, v in pilot.items() if k != "clone_elapsed_s"}
        print(json.dumps(meta["pilot"], indent=2), flush=True)
        proj1 = rescale_phase1(pilot, workers)
        # Phase 2 ~ 18 objects × 10 clones; Phase 3 unknown until Phase 1 censor rate.
        proj2 = project_budget(pilot, COV_N_OBJECTS, SURVEY_N_CLONES, SURVEY_TMAX_MYR, workers)
        # Assume 20% clones censored, extend those to 40 Myr
        n_censor_jobs = int(0.20 * proj1["n_obj"] * SURVEY_N_CLONES)
        s_per = float(pilot["mean_s_per_clone"]) * (0.70 + 0.30 * (40.0 / 25.0))
        proj3_h = (n_censor_jobs * s_per / max(workers, 1)) / 3600.0
        total_h = proj1["wall_h"] + proj2["wall_h"] + proj3_h
        projection = {
            "phase1": proj1,
            "phase2_h": proj2["wall_h"],
            "phase3_h_if_20pct_censored": proj3_h,
            "phases_1_3_h": total_h,
            "pilot_wall_h": pilot["wall_s"] / 3600.0,
        }
        meta["projection"] = projection
        (OUT / "projection_after_pilot.json").write_text(json.dumps(projection, indent=2), encoding="utf-8")
        print("PROJECTION", json.dumps(projection, indent=2), flush=True)
        if proj1["note"] != "as planned":
            meta["deviations"].append(proj1["note"])
            sample = phase_catalog(proj1["n_obj"])
        if total_h > 12.5 and proj1["n_obj"] < 80:
            print(
                "STOP: even a reduced sample projects >12.5 h for phases 1–3. "
                "Not starting the survey. See results/pop_fragility/projection_after_pilot.json",
                flush=True,
            )
            write_report(meta)
            return
        if args.phase == "pilot":
            write_report(meta)
            return
        n_obj = proj1["n_obj"]
    else:
        n_obj = min(args.n_target, len(sample["objects"]))

    if args.phase in ("survey", "all"):
        print("=== PHASE 1 STRATIFIED SURVEY ===", flush=True)
        d1 = min(deadline_all, time.time() + PHASE1_BUDGET_S)
        r1 = phase_survey(sample, n_obj, workers, d1)
        meta["phase1_requested"] = r1["n_requested"]
        meta["phase1_completed"] = r1["n_completed"]
        meta["phase1_wall_h"] = r1["wall_s"] / 3600.0
        print("Phase 1 done", json.dumps(r1), flush=True)

    p1_rows = load_campaign_rows("phase1")
    if args.phase in ("cov", "all") and p1_rows:
        print("=== PHASE 2 COVARIANCE ===", flush=True)
        d2 = min(deadline_all, time.time() + PHASE2_BUDGET_S)
        r2 = phase_covariance(p1_rows, sample, workers, d2)
        meta["phase2"] = {k: v for k, v in r2.items() if k != "fails"}
        meta["phase2_wall_h"] = r2["wall_s"] / 3600.0
        if r2["fails"]:
            meta["deviations"].append(f"covariance fetch/run failures: {r2['fails']}")
        print("Phase 2 done", json.dumps(meta["phase2"]), flush=True)

    if args.phase in ("censor", "all") and p1_rows:
        print("=== PHASE 3 CENSORING RESOLUTION ===", flush=True)
        d3 = min(deadline_all, time.time() + PHASE3_BUDGET_S)
        r3 = phase_censor(p1_rows, sample, workers, d3)
        meta["phase3"] = r3
        meta["phase3_wall_h"] = r3["wall_s"] / 3600.0
        print("Phase 3 done", json.dumps(r3), flush=True)

    if args.phase in ("analyze", "all", "survey", "cov", "censor"):
        path = write_report(meta)
        print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()

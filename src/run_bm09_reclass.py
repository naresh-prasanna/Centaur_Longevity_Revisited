"""Re-integrate the full Bailey & Malhotra (2009) Table 2 sample under modern orbits.

Each object gets n diagonal-Gaussian clones integrated to Tmax and classified with
the v2 (lifetime-first) rule. The modal clone class is compared against the class
BM09 published.

BM09 Table 2 lists a, e, i, q but not the angles, so their integrations cannot be
reproduced exactly. Label differences therefore conflate two effects: the orbit
solution changing since 2007, and this pipeline differing from theirs. The 'opp'
arc entries are multi-opposition orbits whose solutions are essentially unchanged
today, so they isolate the pipeline term; the day-arc entries carry both.

Usage:
  python src/run_bm09_reclass.py --workers 6
  python src/run_bm09_reclass.py --summary-only
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import run_pop_fragility as rpf

ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / "data" / "bailey_table2.json"
MODERN = ROOT / "data" / "modern_elements.json"
OUT = ROOT / "results" / "bm09_reclass"
CAMPAIGN = "clones"
MERGED = "clones_censor_resolved"
N_CLONES = 10
TMAX_MYR = 25.0
CENSOR_TMAX_MYR = 40.0

# run_object_clones writes to rpf.OUT / campaign; point it at our tree.
rpf.OUT = OUT


def parse_arc(arc: str) -> dict:
    """BM09 Table 2 arc column is either 'N opp' (multi-opposition) or 'Nd' (days)."""
    s = (arc or "").strip()
    m = re.match(r"^(\d+)\s*opp$", s)
    if m:
        n = int(m.group(1))
        return {"arc_2007_raw": s, "arc_2007_type": "opposition", "arc_2007_opp": n,
                "arc_2007_days": None, "orbit_2007_quality": "multi_opposition"}
    m = re.match(r"^(\d+)\s*d$", s)
    if m:
        d = int(m.group(1))
        return {"arc_2007_raw": s, "arc_2007_type": "days", "arc_2007_opp": None,
                "arc_2007_days": d, "orbit_2007_quality": "short_arc"}
    return {"arc_2007_raw": s, "arc_2007_type": "unknown", "arc_2007_opp": None,
            "arc_2007_days": None, "orbit_2007_quality": "unknown"}


def orbit_shift(bm: dict, mod: dict) -> dict:
    """How far the published 2007 orbit moved by the modern solution."""
    out = {}
    for key_bm, key_mod, label in (
        ("a_au", "a_au", "a"), ("e", "e", "e"), ("i_deg", "i_deg", "i"), ("q_au", "q_au", "q")
    ):
        old, new = bm.get(key_bm), mod.get(key_mod)
        if old is None or new is None:
            out[f"d{label}"] = None
            continue
        out[f"d{label}"] = new - old
    if bm.get("a_au"):
        out["frac_da"] = (mod["a_au"] - bm["a_au"]) / bm["a_au"]
    if bm.get("q_au"):
        out["frac_dq"] = (mod["q_au"] - bm["q_au"]) / bm["q_au"]
    # Single scalar for ranking: fractional change in q plus absolute change in e.
    de, fdq = out.get("de"), out.get("frac_dq")
    out["shift_index"] = None if de is None or fdq is None else abs(de) + abs(fdq)
    return out


def build_targets() -> tuple[list[dict], list[dict]]:
    table = json.loads(TABLE.read_text(encoding="utf-8"))
    modern = json.loads(MODERN.read_text(encoding="utf-8"))
    targets, skipped = [], []
    for obj in table["objects"]:
        desig = obj["desig"]
        bcls = obj.get("bailey_class")
        if bcls == "NT":
            skipped.append({"desig": desig, "reason": "Neptune Trojan, excluded from BM09 dynamics sample"})
            continue
        mod = modern.get(desig)
        if not mod:
            skipped.append({"desig": desig, "reason": "no modern SBDB solution resolved"})
            continue
        if not mod.get("sigma"):
            skipped.append({"desig": desig, "reason": "no uncertainties in SBDB record"})
            continue
        rec = dict(mod)
        rec["desig"] = desig
        rec["bailey_class"] = bcls
        rec["bailey_elements"] = {k: obj.get(k) for k in ("a_au", "e", "i_deg", "q_au", "H_mag")}
        rec.update(parse_arc(obj.get("arc", "")))
        rec["orbit_shift"] = orbit_shift(obj, mod)
        targets.append(rec)
    return targets, skipped


def run_censor_extension(workers: int) -> dict:
    """Re-run only the right-censored clones to 40 Myr and merge onto the 25 Myr set.

    A clone that is still alive at Tmax carries no lifetime information beyond the
    wall, so any Q assigned from a fully censored set is provisional. Extending only
    those clones keeps the cost proportional to the ambiguity.
    """
    targets, _ = build_targets()
    by_desig = {t["desig"]: t for t in targets}
    n_obj = 0
    n_ext = 0
    for path in sorted((OUT / CAMPAIGN).glob("*.json")):
        if path.name.endswith(".partial.json") or path.name.startswith("_"):
            continue
        blob = json.loads(path.read_text(encoding="utf-8"))
        if "clones" not in blob:
            continue
        cens = [int(c["clone"]) for c in blob["clones"] if c.get("survived_full")]
        if not cens:
            continue
        modern = by_desig.get(blob["desig"])
        if not modern:
            continue
        print(f"CENSOR {blob['desig']}: extending {len(cens)}/{len(blob['clones'])} "
              f"censored clones to {CENSOR_TMAX_MYR:.0f} Myr", flush=True)
        rec = rpf.run_object_clones(
            modern, "clones_40myr", N_CLONES, CENSOR_TMAX_MYR, workers,
            clone_scheme="diagonal", clone_indices=cens,
        )
        merged = {int(c["clone"]): c for c in blob["clones"]}
        for c in rec["blob"]["clones"]:
            if int(c["clone"]) in set(cens):
                merged[int(c["clone"])] = c
        out = rpf.assemble_object(
            blob["desig"], modern, MERGED, "diagonal", N_CLONES, CENSOR_TMAX_MYR,
            list(merged.values()), rec["blob"].get("elapsed_s") or 0.0,
            extra={"extended_clones": cens, "base_campaign": CAMPAIGN,
                   "base_t_max_myr": TMAX_MYR},
        )
        (OUT / MERGED).mkdir(parents=True, exist_ok=True)
        (OUT / MERGED / f"{rpf.slug(blob['desig'])}.json").write_text(
            json.dumps(out, indent=2), encoding="utf-8")
        n_obj += 1
        n_ext += len(cens)
    return {"n_objects_extended": n_obj, "n_clones_extended": n_ext}


def summarize(campaign: str = CAMPAIGN) -> dict:
    from collections import Counter

    rows = []
    for path in sorted((OUT / campaign).glob("*.json")):
        if path.name.endswith(".partial.json") or path.name.startswith("_"):
            continue
        blob = json.loads(path.read_text(encoding="utf-8"))
        if "clones" not in blob:
            continue
        mod = blob["modern_nominal"]
        clones = blob["clones"]
        lifetimes = sorted(c["lifetime_myr"] for c in clones)
        n = len(lifetimes)
        n_censored = sum(1 for c in clones if c.get("survived_full"))
        rows.append({
            "desig": blob["desig"],
            "bailey_class": mod.get("bailey_class"),
            "modern_class": blob["modal_class"],
            "changed": mod.get("bailey_class") != blob["modal_class"],
            "frac_modal": blob["frac_modal"],
            "counts": blob["counts"],
            "n_clones": n,
            "median_lifetime_myr": float(np.median(lifetimes)) if n else None,
            "frac_ge_22myr": sum(1 for x in lifetimes if x >= 22.0) / n if n else None,
            "n_censored": n_censored,
            "orbit_2007_quality": mod.get("orbit_2007_quality"),
            "arc_2007_raw": mod.get("arc_2007_raw"),
            "arc_2007_days": mod.get("arc_2007_days"),
            "arc_2007_opp": mod.get("arc_2007_opp"),
            "data_arc_days_now": mod.get("data_arc_days"),
            "condition_code_now": mod.get("condition_code"),
            "shift_index": (mod.get("orbit_shift") or {}).get("shift_index"),
            "de": (mod.get("orbit_shift") or {}).get("de"),
            "frac_dq": (mod.get("orbit_shift") or {}).get("frac_dq"),
        })

    rows.sort(key=lambda r: r["desig"])
    matrix: dict[str, Counter] = {}
    for r in rows:
        matrix.setdefault(r["bailey_class"], Counter())[r["modern_class"]] += 1

    def split(pred):
        sub = [r for r in rows if pred(r)]
        if not sub:
            return {"n": 0}
        ch = sum(1 for r in sub if r["changed"])
        return {"n": len(sub), "n_changed": ch, "frac_changed": ch / len(sub),
                "changed_desigs": [r["desig"] for r in sub if r["changed"]]}

    control = split(lambda r: r["orbit_2007_quality"] == "multi_opposition")
    testset = split(lambda r: r["orbit_2007_quality"] == "short_arc")
    fully_censored = [r["desig"] for r in rows if r["n_censored"] == r["n_clones"]]
    summary = {
        "campaign": campaign,
        "n_objects": len(rows),
        "n_clones_per_object": N_CLONES,
        "t_max_myr": CENSOR_TMAX_MYR if campaign == MERGED else TMAX_MYR,
        "n_objects_fully_censored": len(fully_censored),
        "fully_censored_desigs": fully_censored,
        "seed": rpf.SEED,
        "classifier_rule": "v2",
        "confusion_bailey_to_modern": {k: dict(v) for k, v in matrix.items()},
        "n_changed": sum(1 for r in rows if r["changed"]),
        "frac_changed": (sum(1 for r in rows if r["changed"]) / len(rows)) if rows else None,
        "control_multi_opposition": control,
        "test_short_arc": testset,
        "interpretation": (
            "Multi-opposition objects had well-determined orbits already in 2007, so label "
            "changes in that subset measure pipeline difference rather than orbit revision. "
            "Short-arc objects carry both effects; the excess change rate over the control "
            "is the orbit-driven component."
        ),
        "rows": rows,
    }
    if control.get("n") and testset.get("n"):
        summary["excess_change_short_arc"] = testset["frac_changed"] - control["frac_changed"]
    OUT.mkdir(parents=True, exist_ok=True)
    tag = "" if campaign == CAMPAIGN else f"_{campaign}"
    (OUT / f"_SUMMARY{tag}.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--n-clones", type=int, default=N_CLONES)
    ap.add_argument("--t-max-myr", type=float, default=TMAX_MYR)
    ap.add_argument("--summary-only", action="store_true")
    ap.add_argument("--censor", action="store_true",
                    help="extend right-censored clones to 40 Myr and re-summarise")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    if args.summary_only:
        print(json.dumps({k: v for k, v in summarize().items() if k != "rows"}, indent=2))
        return

    if args.censor:
        stats = run_censor_extension(args.workers)
        print(json.dumps(stats, indent=2), flush=True)
        if stats["n_objects_extended"]:
            s = summarize(MERGED)
            print(json.dumps({k: v for k, v in s.items() if k != "rows"}, indent=2), flush=True)
        return

    targets, skipped = build_targets()
    (OUT / "_TARGETS.json").write_text(
        json.dumps({"n": len(targets), "skipped": skipped,
                    "targets": [{k: v for k, v in t.items() if k != "sigma"} for t in targets]},
                   indent=2),
        encoding="utf-8",
    )
    print(f"BM09 reclassification: {len(targets)} objects, {len(skipped)} skipped", flush=True)
    for s in skipped:
        print(f"  skip {s['desig']}: {s['reason']}", flush=True)

    t_start = time.time()
    for k, tgt in enumerate(targets, 1):
        el = time.time() - t_start
        print(f"BM09 {k}/{len(targets)} {tgt['desig']} "
              f"bailey={tgt['bailey_class']} arc2007={tgt['arc_2007_raw']} "
              f"elapsed={el/3600:.2f}h", flush=True)
        rpf.run_object_clones(
            tgt, CAMPAIGN, args.n_clones, args.t_max_myr, args.workers, clone_scheme="diagonal"
        )
        rpf.write_status({"phase": "bm09", "object_k": k, "object_n": len(targets),
                          "desig": tgt["desig"], "elapsed_h": el / 3600})

    s = summarize()
    print(json.dumps({k: v for k, v in s.items() if k != "rows"}, indent=2), flush=True)


if __name__ == "__main__":
    main()

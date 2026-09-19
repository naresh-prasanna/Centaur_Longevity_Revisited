"""Reclassify stored clone diagnostics with BM09-faithful v2 rules.

v1 (current paper code): Q if flat+life>=22; R if life>=22 and (nonlinear|r2<0.85);
     else if r2>=0.85 -> D; else if life<22 -> D; else R/D.
v2 (BM09 lifetime correlation): Q if flat+life>=22; R if life>=22; else D.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
LONG = 22.0


def classify_v2(rec: dict) -> str:
    life = float(rec.get("lifetime_myr") or 0.0)
    flat = bool(rec.get("flat"))
    if life >= LONG and flat:
        return "Q"
    if life >= LONG:
        return "R"
    return "D"


def norm_path(p: Path) -> str:
    return p.relative_to(RES).as_posix()


def walk_clones(obj, file_hint: str = ""):
    """Yield (desig, tag, clone_dict) for real integration outputs."""
    if not isinstance(obj, dict):
        return

    # batch-style: {desig, clones: [...]}
    if "clones" in obj and isinstance(obj["clones"], list):
        des = obj.get("desig") or obj.get("label") or "?"
        for c in obj["clones"]:
            if isinstance(c, dict) and "lifetime_myr" in c and "class" in c:
                yield des, "clone", c

    # exp_c per-object file: {desig, bailey_ic: {...}, modern_ic: {...}}
    if "bailey_ic" in obj and isinstance(obj.get("bailey_ic"), dict):
        des = obj.get("desig", "?")
        for key in ("bailey_ic", "modern_ic"):
            c = obj.get(key)
            if isinstance(c, dict) and "lifetime_myr" in c and "class" in c:
                yield des, key, c

    # summary: {results: [{desig, bailey_ic, modern_ic}, ...]}
    if "results" in obj and isinstance(obj["results"], list):
        for r in obj["results"]:
            if not isinstance(r, dict):
                continue
            des = r.get("desig", "?")
            for key in ("bailey_ic", "modern_ic"):
                c = r.get(key)
                if isinstance(c, dict) and "lifetime_myr" in c and "class" in c:
                    yield des, key, c


def modal(classes):
    if not classes:
        return None
    cnt = Counter(classes)
    return max(cnt, key=cnt.get)


def main():
    table = json.loads((ROOT / "data" / "bailey_table2.json").read_text(encoding="utf-8"))
    qr_desigs = {o["desig"] for o in table["objects"] if o["bailey_class"] in ("Q", "R")}
    bailey_of = {o["desig"]: o["bailey_class"] for o in table["objects"]}

    rows = []
    files = list((RES / "batch").glob("*.json")) + list((RES / "exp_c").glob("*.json"))
    files += [
        RES / "th173_production.json",
        RES / "th173_pilot.json",
        RES / "exp_c_summary_qr40.json",
        RES / "batch_summary_qr_long40.json",
        RES / "batch_summary_short_d40.json",
    ]
    seen = set()  # dedupe summary vs per-file

    for f in files:
        if not f.exists():
            continue
        d = json.loads(f.read_text(encoding="utf-8"))
        rel = norm_path(f)
        for des, tag, c in walk_clones(d, rel):
            key = (des, tag, c.get("clone"), round(float(c["lifetime_myr"]), 4), c.get("H"))
            # allow same object from batch and exp_c; skip exact duplicate from summary+file
            dedupe = (des, tag, c.get("clone"), round(float(c["lifetime_myr"]), 4),
                      round(float(c["H"]) if c.get("H") is not None else -1, 8),
                      "summary" in rel)
            if dedupe in seen:
                continue
            seen.add(dedupe)
            v1 = c["class"]
            v2 = classify_v2(c)
            rows.append(
                {
                    "file": rel,
                    "desig": des,
                    "tag": tag,
                    "clone": c.get("clone"),
                    "life": c.get("lifetime_myr"),
                    "flat": c.get("flat"),
                    "nonlinear": c.get("nonlinear"),
                    "H": c.get("H"),
                    "r2": c.get("r2"),
                    "v1": v1,
                    "v2": v2,
                    "changed": v1 != v2,
                }
            )

    # Exp C from summary (canonical 10 objects × 2 ICs)
    expc_sum = RES / "exp_c_summary_qr40.json"
    expc_rows = []
    expc_table = []
    if expc_sum.exists():
        d = json.loads(expc_sum.read_text(encoding="utf-8"))
        for r in d["results"]:
            des = r["desig"]
            bm = r["bailey_class"]
            b1, m1 = r["bailey_ic"]["class"], r["modern_ic"]["class"]
            b2, m2 = classify_v2(r["bailey_ic"]), classify_v2(r["modern_ic"])
            expc_rows.append(
                {
                    "desig": des,
                    "bailey_label": bm,
                    "v1_bailey": b1,
                    "v1_modern": m1,
                    "v2_bailey": b2,
                    "v2_modern": m2,
                    "life_bailey": r["bailey_ic"]["lifetime_myr"],
                    "life_modern": r["modern_ic"]["lifetime_myr"],
                    "flat_bailey": r["bailey_ic"]["flat"],
                    "flat_modern": r["modern_ic"]["flat"],
                    "r2_bailey": r["bailey_ic"]["r2"],
                    "r2_modern": r["modern_ic"]["r2"],
                    "nonlinear_bailey": r["bailey_ic"]["nonlinear"],
                    "nonlinear_modern": r["modern_ic"]["nonlinear"],
                    "orbit_drives_v1": b1 != m1,
                    "orbit_drives_v2": b2 != m2,
                    "recovers_bm09_v2_bailey": b2 == bm,
                    "recovers_bm09_v2_modern": m2 == bm,
                }
            )
            expc_table.append(expc_rows[-1])

    # Q/R batch modal under v1/v2
    qr_modal = {}
    for r in rows:
        if not r["file"].startswith("batch/") or r["file"].startswith("batch_summary"):
            continue
        des = r["desig"]
        if des not in qr_desigs:
            # try from filename
            stem = Path(r["file"]).stem.replace("_", " ")
            if stem in qr_desigs:
                des = stem
            else:
                continue
        qr_modal.setdefault(des, []).append(r)

    qr_pilot = {}
    for des, lst in sorted(qr_modal.items()):
        # unique clones by clone index
        by_c = {}
        for x in lst:
            by_c[x.get("clone")] = x
        lst2 = list(by_c.values())
        qr_pilot[des] = {
            "bailey": bailey_of[des],
            "counts_v1": dict(Counter(x["v1"] for x in lst2)),
            "counts_v2": dict(Counter(x["v2"] for x in lst2)),
            "modal_v1": modal([x["v1"] for x in lst2]),
            "modal_v2": modal([x["v2"] for x in lst2]),
            "flip_v1": modal([x["v1"] for x in lst2]) != bailey_of[des],
            "flip_v2": modal([x["v2"] for x in lst2]) != bailey_of[des],
            "n_clones": len(lst2),
        }

    # Real R/Q hits (finite life)
    real_rq = [
        r for r in rows
        if r["v1"] in ("R", "Q") and r["life"] is not None and r["tag"] == "clone"
    ]

    out = {
        "rule_v2": "life>=22 & flat -> Q; life>=22 -> R; else D",
        "verdict": (
            "v1 is not fully null (rare real R/Q exist) but is structurally D-biased for "
            "long-lived linear Hurst survivors. Exp C all-D under v1 is partly early escape "
            "(legitimate D) and partly long-survivor forced-D (instrument). Under v2, Exp C "
            "is no longer all-D."
        ),
        "summary_all_uniqueish": {
            "n": len(rows),
            "v1": dict(Counter(r["v1"] for r in rows)),
            "v2": dict(Counter(r["v2"] for r in rows)),
            "n_changed": sum(1 for r in rows if r["changed"]),
        },
        "real_clone_RQ_v1": real_rq,
        "n_real_clone_RQ_v1": len(real_rq),
        "exp_c_v1_vs_v2": expc_table,
        "exp_c_n_orbit_drives_v1": sum(1 for r in expc_table if r["orbit_drives_v1"]),
        "exp_c_n_orbit_drives_v2": sum(1 for r in expc_table if r["orbit_drives_v2"]),
        "exp_c_n_v2_bailey_R": sum(1 for r in expc_table if r["v2_bailey"] == "R"),
        "exp_c_n_v2_modern_R": sum(1 for r in expc_table if r["v2_modern"] == "R"),
        "exp_c_n_v2_bailey_Q": sum(1 for r in expc_table if r["v2_bailey"] == "Q"),
        "exp_c_n_v2_modern_Q": sum(1 for r in expc_table if r["v2_modern"] == "Q"),
        "qr_pilot_modal": qr_pilot,
        "qr_flip_rate_v1": (
            sum(1 for v in qr_pilot.values() if v["flip_v1"]) / max(1, len(qr_pilot))
        ),
        "qr_flip_rate_v2": (
            sum(1 for v in qr_pilot.values() if v["flip_v2"]) / max(1, len(qr_pilot))
        ),
        "changed_sample": [r for r in rows if r["changed"]][:30],
    }

    dest = RES / "CLASSIFIER_V2_RECLASS.json"
    dest.write_text(json.dumps(out, indent=2), encoding="utf-8")

    print("=== VERDICT ===")
    print(out["verdict"])
    print("summary", out["summary_all_uniqueish"])
    print("real_RQ_v1", out["n_real_clone_RQ_v1"])
    print(
        "exp_c orbit_drives v1/v2",
        out["exp_c_n_orbit_drives_v1"],
        out["exp_c_n_orbit_drives_v2"],
    )
    print(
        "exp_c v2 R bailey/modern",
        out["exp_c_n_v2_bailey_R"],
        out["exp_c_n_v2_modern_R"],
    )
    print("qr flip rates v1/v2", out["qr_flip_rate_v1"], out["qr_flip_rate_v2"])
    print("--- Exp C table ---")
    for r in expc_table:
        print(
            f"{r['desig']:12} BM09={r['bailey_label']}  "
            f"v1 {r['v1_bailey']}/{r['v1_modern']}  "
            f"v2 {r['v2_bailey']}/{r['v2_modern']}  "
            f"life {r['life_bailey']:.2f}/{r['life_modern']:.2f}"
        )
    print("--- Q/R pilot modal ---")
    for des, v in qr_pilot.items():
        print(
            f"{des}: BM09 {v['bailey']} -> v1 {v['modal_v1']} {v['counts_v1']} | "
            f"v2 {v['modal_v2']} {v['counts_v2']}"
        )
    print("wrote", dest)


if __name__ == "__main__":
    main()

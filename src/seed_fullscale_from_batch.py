"""Seed fullscale/qr20 from any existing 20-clone batch files (re-label with v2)."""
from __future__ import annotations

import json
from pathlib import Path

from classify_core import reclassify_from_diagnostics

ROOT = Path(__file__).resolve().parents[1]
BATCH = ROOT / "results" / "batch"
OUT = ROOT / "results" / "fullscale" / "qr20"


def flat_of(rec: dict) -> bool:
    if rec.get("a_mean") is not None and rec.get("a_std") is not None:
        return float(rec["a_std"]) / (float(rec["a_mean"]) + 1e-9) < 0.03
    return bool(rec.get("flat"))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    table = {
        o["desig"]: o
        for o in json.loads((ROOT / "data" / "bailey_table2.json").read_text(encoding="utf-8"))["objects"]
    }
    seeded = []
    for path in BATCH.glob("*.json"):
        d = json.loads(path.read_text(encoding="utf-8"))
        if d.get("n_clones", 0) < 20 or abs(float(d.get("t_max_myr", 0)) - 40.0) > 1e-6:
            continue
        des = d.get("desig")
        if des not in table or table[des]["bailey_class"] not in ("Q", "R"):
            continue
        clones = []
        for c in d["clones"]:
            rec = dict(c)
            rec["flat"] = flat_of(c)
            v2 = reclassify_from_diagnostics(rec, "v2")
            v1 = reclassify_from_diagnostics(rec, "v1")
            clones.append({**c, "flat": rec["flat"], "class": v2, "class_v2": v2, "class_v1": v1, "rule": "v2"})
        classes = [c["class"] for c in clones]
        counts = {k: classes.count(k) for k in ("D", "R", "Q")}
        modal = max(counts, key=counts.get)
        bailey_class = table[des]["bailey_class"]
        out = {
            **{k: d.get(k) for k in ("desig", "bailey", "modern_nominal", "seed", "n_clones", "t_max_myr")},
            "campaign": "qr20",
            "classifier_rule": "v2",
            "seeded_from": str(path.name),
            "counts": counts,
            "modal_class": modal,
            "frac_modal": counts[modal] / len(clones),
            "bailey_class": bailey_class,
            "flip": modal != bailey_class,
            "call": f"{bailey_class}->{modal}" if modal != bailey_class else f"stay {modal}",
            "clones": clones,
        }
        dest = OUT / path.name
        dest.write_text(json.dumps(out, indent=2), encoding="utf-8")
        seeded.append(des)
        print(f"seeded {des} -> {modal} {counts}")
    print(f"seeded {len(seeded)} objects into {OUT}")


if __name__ == "__main__":
    main()

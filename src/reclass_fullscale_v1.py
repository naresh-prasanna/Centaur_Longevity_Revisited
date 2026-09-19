"""Re-label full-scale archived integrations under v1 and v2 (no re-integration)."""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from classify_core import reclassify_from_diagnostics

ROOT = Path(__file__).resolve().parents[1]
FS = ROOT / "results" / "fullscale"
OUT = ROOT / "results" / "FULLSCALE_V1_V2_RECLASS.json"


def modal(classes: list[str]) -> str | None:
    if not classes:
        return None
    return Counter(classes).most_common(1)[0][0]


def reclass_object(obj: dict) -> dict:
    clones = obj.get("clones") or []
    v1 = [reclassify_from_diagnostics(c, "v1") for c in clones]
    v2 = [reclassify_from_diagnostics(c, "v2") for c in clones]
    return {
        "desig": obj["desig"],
        "campaign": obj.get("campaign"),
        "bailey_class": obj.get("bailey_class") or obj.get("bailey", {}).get("bailey_class"),
        "n_clones": len(clones),
        "t_max_myr": obj.get("t_max_myr"),
        "counts_v1": dict(Counter(v1)),
        "counts_v2": dict(Counter(v2)),
        "modal_v1": modal(v1),
        "modal_v2": modal(v2),
        "flip_v1": modal(v1) != (obj.get("bailey_class") or obj.get("bailey", {}).get("bailey_class")),
        "flip_v2": modal(v2) != (obj.get("bailey_class") or obj.get("bailey", {}).get("bailey_class")),
    }


def load_campaign(campaign: str) -> list[dict]:
    d = FS / campaign
    if not d.exists():
        return []
    return [json.loads(p.read_text(encoding="utf-8")) for p in sorted(d.glob("*.json"))]


def main() -> None:
    campaigns = ["qr20", "d5", "d20", "th173_40", "th173_100x40", "cov_th173", "cov_qp112", "cov_qr20"]
    out: dict = {"campaigns": {}}
    for camp in campaigns:
        objs = load_campaign(camp)
        if not objs:
            continue
        rows = [reclass_object(o) for o in objs]
        out["campaigns"][camp] = {
            "n": len(rows),
            "flip_rate_v1": sum(1 for r in rows if r["flip_v1"]) / len(rows),
            "flip_rate_v2": sum(1 for r in rows if r["flip_v2"]) / len(rows),
            "objects": rows,
        }
    # qr20 aggregate
    qr = out["campaigns"].get("qr20", {}).get("objects", [])
    if qr:
        agg_v1 = Counter()
        agg_v2 = Counter()
        for r in qr:
            for k, v in r["counts_v1"].items():
                agg_v1[k] += v
            for k, v in r["counts_v2"].items():
                agg_v2[k] += v
        out["qr20_aggregate_v1"] = dict(agg_v1)
        out["qr20_aggregate_v2"] = dict(agg_v2)
        out["qr20_modal_flip_v1"] = sum(1 for r in qr if r["flip_v1"])
        out["qr20_modal_flip_v2"] = sum(1 for r in qr if r["flip_v2"])
    OUT.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"wrote {OUT}")
    if qr:
        print(
            f"qr20 modal flip v1/v2: {out['qr20_modal_flip_v1']}/{len(qr)} "
            f"{out['qr20_modal_flip_v2']}/{len(qr)}"
        )


if __name__ == "__main__":
    main()

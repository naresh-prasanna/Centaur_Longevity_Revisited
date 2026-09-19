"""Audit all result JSONs for R/Q class outputs."""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "results"


def walk(obj, path="$"):
    if isinstance(obj, dict):
        if "class" in obj and obj["class"] in ("D", "R", "Q"):
            yield path, obj
        for k, v in obj.items():
            yield from walk(v, f"{path}.{k}")
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from walk(v, f"{path}[{i}]")


def main():
    class_counts = Counter()
    rq = []
    per_file = Counter()
    long_lived_but_d = []
    n_long = 0

    for f in sorted(ROOT.rglob("*.json")):
        try:
            d = json.loads(f.read_text(encoding="utf-8"))
        except Exception:
            continue
        for path, rec in walk(d):
            cls = rec["class"]
            class_counts[cls] += 1
            per_file[(str(f.relative_to(ROOT)), cls)] += 1
            life = rec.get("lifetime_myr")
            if life is not None and life >= 22:
                n_long += 1
                if cls == "D":
                    long_lived_but_d.append(
                        {
                            "file": str(f.relative_to(ROOT)),
                            "path": path,
                            "life": life,
                            "H": rec.get("H"),
                            "r2": rec.get("r2"),
                            "flat": rec.get("flat"),
                            "nonlinear": rec.get("nonlinear"),
                        }
                    )
            if cls in ("R", "Q"):
                rq.append(
                    {
                        "file": str(f.relative_to(ROOT)),
                        "path": path,
                        "class": cls,
                        "life": life,
                        "H": rec.get("H"),
                        "r2": rec.get("r2"),
                        "flat": rec.get("flat"),
                        "nonlinear": rec.get("nonlinear"),
                        "clone": rec.get("clone"),
                    }
                )

    out = {
        "global_counts": dict(class_counts),
        "total": sum(class_counts.values()),
        "frac_D": class_counts["D"] / max(1, sum(class_counts.values())),
        "n_RQ_hits": len(rq),
        "RQ_hits": rq,
        "n_long_lived_records": n_long,
        "n_long_lived_but_D": len(long_lived_but_d),
        "long_lived_but_D_sample": long_lived_but_d[:30],
    }
    dest = ROOT / "CLASSIFIER_AUDIT.json"
    dest.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("GLOBAL", out["global_counts"], "total", out["total"], "frac_D", round(out["frac_D"], 4))
    print("RQ_hits", out["n_RQ_hits"])
    for h in rq:
        print(" ", h)
    print("long_lived_records", n_long, "of which D", len(long_lived_but_d))
    print("wrote", dest)


if __name__ == "__main__":
    main()

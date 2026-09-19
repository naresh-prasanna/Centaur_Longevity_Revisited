"""Arc-scaled clone sampling and classification for Option A."""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OJ = ROOT
OJ_SRC = ROOT / "src"
sys.path.insert(0, str(OJ_SRC))

from classify_core import (  # noqa: E402
    DE441Ephem,
    classify_bailey,
    integrate_clone,
    sample_clones_from_modern,
)

SEED = 87
LONG_MYR = 22.0


def load_modern(desig: str) -> dict:
    modern = json.loads((OJ / "data" / "modern_elements.json").read_text(encoding="utf-8"))
    if desig not in modern:
        raise KeyError(desig)
    return modern[desig]


def load_bailey(desig: str) -> str | None:
    table = json.loads((OJ / "data" / "bailey_table2.json").read_text(encoding="utf-8"))
    for o in table["objects"]:
        if o["desig"] == desig:
            return o["bailey_class"]
    return None


def scale_modern_for_arc(modern: dict, arc_days: float) -> dict:
    """Return a copy of modern elements with σ inflated for a shorter effective arc."""
    A = float(modern.get("data_arc_days") or 0.0)
    arc = float(max(arc_days, 1.0))
    out = json.loads(json.dumps(modern))  # deep copy via JSON
    out["effective_arc_days"] = arc
    sig = dict(modern.get("sigma") or {})
    if A <= 0:
        scale = 1.0
    elif arc >= A:
        scale = 1.0
        out["effective_arc_days"] = A
    else:
        scale = float(np.sqrt(A / arc))
    out["sigma"] = {k: float(v) * scale for k, v in sig.items()}
    out["sigma_scale"] = scale
    return out


def arc_tiers_for(modern: dict, n_tiers: int = 5) -> list[float]:
    A = float(modern.get("data_arc_days") or 100.0)
    lo = min(30.0, max(10.0, A / 10.0))
    if A <= lo * 1.5:
        # Very short arc: denser near full
        raw = np.unique(np.clip(np.geomspace(max(5.0, A / 5), A, num=n_tiers), 5.0, A))
    else:
        raw = np.unique(np.clip(np.geomspace(lo, A, num=n_tiers), lo, A))
    # Always include full arc
    tiers = sorted(set([float(x) for x in raw] + [float(A)]))
    return tiers


def run_tier(
    modern: dict,
    arc_days: float,
    n_clones: int,
    t_max_yr: float,
    ephem: DE441Ephem,
    seed: int = SEED,
) -> dict:
    scaled = scale_modern_for_arc(modern, arc_days)
    clones = sample_clones_from_modern(scaled, n_clones, seed=seed)
    rows = []
    for i, elems in enumerate(clones):
        traj = integrate_clone(elems, ephem, t_max_yr=t_max_yr)
        lab = classify_bailey(traj["a_au"], traj["t_yr"], traj["lifetime_yr"], rule="v2")
        rows.append(
            {
                "clone": i,
                "class": lab["class"],
                "class_v1": lab["class_v1"],
                "class_v2": lab["class_v2"],
                "lifetime_myr": lab["lifetime_myr"],
                "flat": lab["flat"],
                "H": lab["H"],
                "r2": lab["r2"],
                "escaped": traj["escaped"],
            }
        )
    counts = Counter(r["class"] for r in rows)
    modal = max(counts, key=counts.get) if counts else None
    return {
        "arc_days": float(scaled["effective_arc_days"]),
        "sigma_scale": scaled["sigma_scale"],
        "n_clones": n_clones,
        "counts": dict(counts),
        "modal": modal,
        "median_life_myr": float(np.median([r["lifetime_myr"] for r in rows])),
        "frac_long": float(np.mean([r["lifetime_myr"] >= LONG_MYR for r in rows])),
        "clones": rows,
    }


def agreement_vs_full(tier_results: list[dict]) -> list[dict]:
    full = max(tier_results, key=lambda r: r["arc_days"])
    full_modal = full["modal"]
    out = []
    for t in tier_results:
        classes = [c["class"] for c in t["clones"]]
        agree = float(np.mean([c == full_modal for c in classes])) if classes else 0.0
        out.append(
            {
                "arc_days": t["arc_days"],
                "modal": t["modal"],
                "full_modal": full_modal,
                "modal_matches_full": t["modal"] == full_modal,
                "clone_agreement": agree,
                "counts": t["counts"],
                "median_life_myr": t["median_life_myr"],
                "frac_long": t["frac_long"],
                "sigma_scale": t["sigma_scale"],
            }
        )
    return out

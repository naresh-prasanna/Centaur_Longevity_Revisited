"""Empirical longevity cut from archived Q/R clone lifetimes (no new integrations)."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FS = ROOT / "results" / "fullscale" / "qr20"
OUT = ROOT / "results" / "DERIVE_22MYR_CUT.json"
BM09_CUT = 22.0


def main() -> None:
    lives: list[float] = []
    by_obj: dict[str, list[float]] = {}
    for f in sorted(FS.glob("*.json")):
        obj = json.loads(f.read_text(encoding="utf-8"))
        des = obj["desig"]
        ls = [float(c["lifetime_myr"]) for c in obj.get("clones", [])]
        lives.extend(ls)
        by_obj[des] = ls

    x = np.sort(np.asarray(lives, dtype=np.float64))
    n = len(x)
    if n == 0:
        raise SystemExit("no qr20 lifetimes found")

    # Empirical quantiles
    qs = [0.5, 0.75, 0.9, 0.95, 0.99]
    quantiles = {f"q{int(q*100)}": float(np.quantile(x, q)) for q in qs}

    # Gap-based cut: largest jump in sorted lifetimes above median short-lived cluster
    diffs = np.diff(x)
    idx = int(np.argmax(diffs)) if len(diffs) else 0
    gap_cut = float(0.5 * (x[idx] + x[idx + 1])) if idx + 1 < n else float(x[-1])

    # Kaplan-Meier style: fraction surviving past candidate cuts
    cuts = np.unique(np.clip(np.linspace(5, 50, 46), 5, 50))
    surv = {float(c): float(np.mean(x >= c)) for c in cuts}

    out = {
        "n_clones": n,
        "n_objects": len(by_obj),
        "bm09_cut_myr": BM09_CUT,
        "median_myr": float(np.median(x)),
        "mean_myr": float(np.mean(x)),
        "max_myr": float(np.max(x)),
        "quantiles_myr": quantiles,
        "max_gap_cut_myr": gap_cut,
        "max_gap_between_myr": float(diffs[idx]) if len(diffs) else None,
        "frac_survive_ge_bm09_cut": float(np.mean(x >= BM09_CUT)),
        "survival_fraction_vs_cut": surv,
        "interpretation": (
            "Exploratory cut from this pipeline's Q/R clone lifetimes only. "
            "Not a replacement for BM09's original 22 Myr unless independently validated."
        ),
        "by_object_median_myr": {k: float(np.median(v)) for k, v in by_obj.items()},
    }
    OUT.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"wrote {OUT}")
    print(f"median={out['median_myr']:.2f} gap_cut={gap_cut:.2f} frac>=22={out['frac_survive_ge_bm09_cut']:.3f}")


if __name__ == "__main__":
    main()

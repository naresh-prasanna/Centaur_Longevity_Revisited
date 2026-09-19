"""Classifier validation: synthetic a(t) series with known expected class.

Default production rule is v2 (BM09 lifetime-faithful).
v1 is retained only as a sensitivity / regression check — several toys
are *expected* to fail under v1 (that is the documented stop-check bug).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from classify_core import classify_bailey

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "classifier_validation.json"


def series_diffusive(t_max_myr: float = 10.0, seed: int = 1) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    sample_yr = 300.0
    n = int(t_max_myr * 1e6 / sample_yr)
    t = np.arange(n) * sample_yr
    a = 20.0 + np.cumsum(rng.normal(0, 0.02, size=n))
    return t, a


def series_flat_long(t_max_myr: float = 40.0, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    sample_yr = 300.0
    n = int(t_max_myr * 1e6 / sample_yr)
    t = np.arange(n) * sample_yr
    rng = np.random.default_rng(seed)
    a = 15.8 + rng.normal(0.0, 1e-4, size=n)
    return t, a


def series_resonance_hop(t_max_myr: float = 40.0) -> tuple[np.ndarray, np.ndarray]:
    sample_yr = 300.0
    n = int(t_max_myr * 1e6 / sample_yr)
    t = np.arange(n) * sample_yr
    a = np.zeros(n)
    cuts = [0, n // 10, n // 3, n // 2, int(0.8 * n), n]
    vals = [20.0, 28.0, 21.5, 30.0, 22.0]
    for i, val in enumerate(vals):
        a[cuts[i] : cuts[i + 1]] = val
    return t, a


def series_long_linear_wander(t_max_myr: float = 40.0, seed: int = 7) -> tuple[np.ndarray, np.ndarray]:
    """Long-lived diffuser with excellent linear Hurst fit — v1 wrongly says D; v2 says R."""
    rng = np.random.default_rng(seed)
    sample_yr = 300.0
    n = int(t_max_myr * 1e6 / sample_yr)
    t = np.arange(n) * sample_yr
    a = 22.0 + np.cumsum(rng.normal(0, 0.015, size=n))
    return t, a


def series_short_escape() -> tuple[np.ndarray, np.ndarray, float]:
    sample_yr = 300.0
    life_yr = 2.5e6
    n = int(life_yr / sample_yr)
    t = np.arange(n) * sample_yr
    a = 18.0 + 0.001 * np.arange(n)
    return t, a, life_yr


def series_almost_long_escape() -> tuple[np.ndarray, np.ndarray, float]:
    """Dies at 15 Myr — must be D under both rules (below 22 Myr cut)."""
    sample_yr = 300.0
    life_yr = 15.0e6
    n = int(life_yr / sample_yr)
    t = np.arange(n) * sample_yr
    rng = np.random.default_rng(3)
    a = 20.0 + np.cumsum(rng.normal(0, 0.02, size=n))
    return t, a, life_yr


CASES = [
    {
        "name": "synthetic_diffusive_10Myr",
        "expect_v2": "D",
        "expect_v1": "D",
        "builder": lambda: (*series_diffusive(10.0), None),
    },
    {
        "name": "synthetic_flat_40Myr",
        "expect_v2": "Q",
        "expect_v1": "Q",
        "builder": lambda: (*series_flat_long(40.0), None),
    },
    {
        "name": "synthetic_hop_40Myr",
        "expect_v2": "R",
        "expect_v1": "D",  # v1 fail mode: excellent r2 → forced D
        "builder": lambda: (*series_resonance_hop(40.0), None),
    },
    {
        "name": "synthetic_long_linear_wander_40Myr",
        "expect_v2": "R",
        "expect_v1": "D",  # v1 fail mode: long + r2>=0.85 → D
        "builder": lambda: (*series_long_linear_wander(40.0), None),
    },
    {
        "name": "synthetic_short_escape",
        "expect_v2": "D",
        "expect_v1": "D",
        "builder": series_short_escape,
    },
    {
        "name": "synthetic_escape_15Myr",
        "expect_v2": "D",
        "expect_v1": "D",
        "builder": series_almost_long_escape,
    },
]


def run_case(case: dict) -> dict:
    built = case["builder"]()
    if len(built) == 3 and built[2] is not None:
        t, a, life = built
    else:
        t, a = built[0], built[1]
        life = float(t[-1])
    c = classify_bailey(a, t, life, rule="v2")
    return {
        "name": case["name"],
        "expect_v2": case["expect_v2"],
        "expect_v1": case["expect_v1"],
        "got_v2": c["class_v2"],
        "got_v1": c["class_v1"],
        "pass_v2": c["class_v2"] == case["expect_v2"],
        "pass_v1": c["class_v1"] == case["expect_v1"],
        "v1_known_fail": case["expect_v1"] != case["expect_v2"],
        "H": c["H"],
        "r2": c["r2"],
        "lifetime_myr": c["lifetime_myr"],
        "flat": c["flat"],
        "nonlinear": c["nonlinear"],
        "rel_std": c["rel_std"],
    }


def main() -> int:
    cases = [run_case(c) for c in CASES]
    n_v2 = sum(1 for c in cases if c["pass_v2"])
    n_v1 = sum(1 for c in cases if c["pass_v1"])
    out = {
        "default_rule": "v2",
        "n": len(cases),
        "n_pass_v2": n_v2,
        "n_pass_v1": n_v1,
        "pass_rate_v2": n_v2 / len(cases),
        "pass_rate_v1": n_v1 / len(cases),
        "v2_ok": n_v2 == len(cases),
        "cases": cases,
        "note": (
            "v2 must be 100% on this battery before science claims. "
            "v1 is expected to fail hop + long-linear toys (r2>=0.85 gate)."
        ),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))
    if not out["v2_ok"]:
        print("FAIL: production classifier v2 did not pass all synthetic cases", file=sys.stderr)
        return 1
    print("PASS: v2 synthetic battery 100%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

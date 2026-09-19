"""2005 TH173 elements: Bailey 2009 vs modern JPL SBDB."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

# Bailey & Malhotra 2009 Table 2 (MPC 2007-03-06)
BAILEY = {
    "designation": "2005 TH173",
    "a_au": 15.724,
    "e": 0.014,
    "i_deg": 15.7,
    "q_au": 15.5,
    "class": "Q",
    "arc_note": "17d single opposition",
}

# JPL SBDB snapshot (fetched 2026-07-26; pe_used DE441)
MODERN = {
    "designation": "2005 TH173",
    "spkid": "50302016",
    "epoch_jd": 2461200.5,
    "a_au": 19.9519498408369,
    "e": 0.3069136342429877,
    "i_deg": 13.47677552407805,
    "om_deg": 193.3029408127479,
    "w_deg": 231.350902903277,
    "ma_deg": 57.68271341657215,
    "q_au": 13.82842440495185,
    "condition_code": 6,
    "data_arc_days": 366,
    "n_obs": 18,
    # 1-sigma from SBDB element table
    "sigma": {
        "a_au": 0.13379,
        "e": 0.008235,
        "i_deg": 0.0052397,
        "om_deg": 0.0010284,
        "w_deg": 0.9351,
        "ma_deg": 0.88531,
    },
}


def sample_clones(n: int, seed: int = 87) -> list[dict]:
    """Independent Gaussian clones in element space (diagonal sigma)."""
    rng = np.random.default_rng(seed)
    s = MODERN["sigma"]
    out = []
    for _ in range(n):
        e = float(np.clip(rng.normal(MODERN["e"], s["e"]), 0.0, 0.95))
        a = float(max(5.0, rng.normal(MODERN["a_au"], s["a_au"])))
        out.append(
            {
                "a": a,
                "e": e,
                "inc": float(np.deg2rad(rng.normal(MODERN["i_deg"], s["i_deg"]))),
                "Omega": float(np.deg2rad(rng.normal(MODERN["om_deg"], s["om_deg"]))),
                "omega": float(np.deg2rad(rng.normal(MODERN["w_deg"], s["w_deg"]))),
                "M": float(np.deg2rad(rng.normal(MODERN["ma_deg"], s["ma_deg"]))),
            }
        )
    return out


def write_snapshot(path: Path):
    path.write_text(
        json.dumps({"bailey": BAILEY, "modern": MODERN}, indent=2),
        encoding="utf-8",
    )


if __name__ == "__main__":
    write_snapshot(Path(__file__).resolve().parents[1] / "data" / "th173_elements.json")
    print("Bailey", BAILEY)
    print("Modern a,e,i", MODERN["a_au"], MODERN["e"], MODERN["i_deg"])

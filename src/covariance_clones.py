"""SBDB covariance-matrix clone sampling (eq, q, tp, node, peri, i parameterization)."""
from __future__ import annotations

import json
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "sbdb_covariance_cache.json"
API = "https://ssd-api.jpl.nasa.gov/sbdb.api"
G_AU_YR3 = 4 * np.pi**2  # GM_sun in AU^3/yr^2 for n = sqrt(GM/a^3)


def fetch_covariance(desig: str) -> dict:
    q = urllib.parse.urlencode({"sstr": desig, "full-prec": "true", "cov": "mat"})
    with urllib.request.urlopen(f"{API}?{q}", timeout=90) as resp:
        raw = json.loads(resp.read().decode("utf-8"))
    cov = raw["orbit"]["covariance"]
    labels = cov["labels"]
    mat = np.array([[float(x) for x in row] for row in cov["data"]], dtype=np.float64)
    elems = {}
    for e in cov["elements"]:
        elems[e["name"]] = float(e["value"])
    return {
        "desig": desig,
        "labels": labels,
        "matrix": mat.tolist(),
        "elements": elems,
        "epoch_jd": float(raw["orbit"]["epoch"]),
    }


def load_or_fetch(desig: str) -> dict:
    if CACHE.exists():
        cache = json.loads(CACHE.read_text(encoding="utf-8"))
        if desig in cache:
            c = cache[desig]
            c["_matrix"] = np.array(c["matrix"], dtype=np.float64)
            return c
    rec = fetch_covariance(desig)
    cache = json.loads(CACHE.read_text(encoding="utf-8")) if CACHE.exists() else {}
    cache[desig] = {k: v for k, v in rec.items() if k != "_matrix"}
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps(cache, indent=2), encoding="utf-8")
    rec["_matrix"] = np.array(rec["matrix"], dtype=np.float64)
    return rec


def _sample_to_elems(sample: dict, epoch_jd: float) -> dict:
    """Convert SBDB (e,q,om,w,i,tp) sample to osculating a,e,i,Omega,omega,M at epoch."""
    e = float(np.clip(sample["e"], 0.0, 0.99))
    q = float(max(0.5, sample["q"]))
    a = q / (1.0 - e)
    inc = np.deg2rad(float(sample["i"]))
    Omega = np.deg2rad(float(sample["om"]))
    omega = np.deg2rad(float(sample["w"]))
    tp = float(sample["tp"])
    n = np.sqrt(G_AU_YR3 / a**3)
    M = (n * (epoch_jd - tp) * 86400.0 / 365.25) % (2 * np.pi)
    return {
        "a": a,
        "e": e,
        "inc": float(inc),
        "Omega": float(Omega),
        "omega": float(omega),
        "M": float(M),
    }


def _nominal_sample(elems: dict) -> dict:
    return {
        "e": elems["e"],
        "q": elems["q"],
        "tp": elems["tp"],
        "om": elems["om"],
        "w": elems["w"],
        "i": elems["i"],
    }


def _mean_vector(labels: list[str], elems: dict) -> np.ndarray:
    key_map = {"node": "om", "peri": "w"}
    return np.array([elems[key_map.get(lab, lab)] for lab in labels], dtype=np.float64)


def sample_clones_from_covariance(desig: str, n: int, seed: int = 87) -> list[dict]:
    """Multivariate Gaussian in SBDB native elements; clone 0 = nominal (no perturbation)."""
    rec = load_or_fetch(desig)
    mat = rec["_matrix"]
    labels = rec["labels"]
    elems = rec["elements"]
    epoch_jd = rec["epoch_jd"]
    rng = np.random.default_rng(seed)
    mean = _mean_vector(labels, elems)

    out = [_sample_to_elems(_nominal_sample(elems), epoch_jd)]
    for _ in range(max(0, n - 1)):
        draw = rng.multivariate_normal(mean, mat)
        sample = {
            "e": draw[labels.index("e")],
            "q": draw[labels.index("q")],
            "tp": draw[labels.index("tp")],
            "om": draw[labels.index("node")],
            "w": draw[labels.index("peri")],
            "i": draw[labels.index("i")],
        }
        out.append(_sample_to_elems(sample, epoch_jd))
    return out


def prefetch(desigs: list[str]) -> None:
    for d in desigs:
        print(f"fetch cov {d}", flush=True)
        load_or_fetch(d)

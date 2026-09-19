"""Catalog, stratified sample, summaries, and plots for the population fragility survey.

Does not import REBOUND. Clone integration lives in run_pop_fragility.py.
"""
from __future__ import annotations

import csv
import json
import math
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SEED = 87
CATALOG_PATH = ROOT / "data" / "sbdb_centaurs_od_cache.json"
SAMPLE_PATH = ROOT / "data" / "pop_fragility_sample.json"
OUT_DIR = ROOT / "results" / "pop_fragility"
QUERY_API = "https://ssd-api.jpl.nasa.gov/sbdb_query.api"
FIELDS = (
    "pdes,full_name,spkid,a,e,i,om,w,ma,q,condition_code,data_arc,"
    "n_obs_used,first_obs,last_obs,source,sigma_a,sigma_e,sigma_i,"
    "sigma_om,sigma_w,sigma_ma"
)


def _num(x):
    if x is None or x == "":
        return None
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _int(x):
    v = _num(x)
    if v is None:
        return None
    return int(v)


def n_oppositions_est(first_obs, last_obs, data_arc_days) -> int | None:
    """SBDB Query has no n_opps field. Estimate opposition-years from the arc.

    Inclusive calendar-year span of first_obs/last_obs, floored at 1.
    Falls back to 1 + floor(arc_days / 365.25).
    """
    try:
        y0 = int(str(first_obs)[:4])
        y1 = int(str(last_obs)[:4])
        if 1800 < y0 < 2100 and 1800 < y1 < 2100 and y1 >= y0:
            return max(1, y1 - y0 + 1)
    except (TypeError, ValueError):
        pass
    if data_arc_days is None:
        return None
    return max(1, int(float(data_arc_days) // 365.25) + 1)


def designation(pdes: str, full_name: str) -> str:
    """Prefer unpacked provisional designation; else pdes."""
    name = (full_name or "").strip()
    if "(" in name and ")" in name:
        inner = name[name.rfind("(") + 1 : name.rfind(")")].strip()
        if inner and inner[0].isdigit():
            return inner
        if inner.startswith("A") and len(inner) > 2:
            return inner
    return str(pdes).strip()


def row_to_modern(fields: list[str], values: list) -> dict | None:
    rec = dict(zip(fields, values))
    a = _num(rec.get("a"))
    e = _num(rec.get("e"))
    if a is None or e is None:
        return None
    arc = _num(rec.get("data_arc"))
    first = rec.get("first_obs")
    last = rec.get("last_obs")
    sigma = {}
    mapping = (
        ("sigma_a", "a_au"),
        ("sigma_e", "e"),
        ("sigma_i", "i_deg"),
        ("sigma_om", "om_deg"),
        ("sigma_w", "w_deg"),
        ("sigma_ma", "ma_deg"),
    )
    for src, dst in mapping:
        v = _num(rec.get(src))
        if v is not None and np.isfinite(v):
            sigma[dst] = v
    pdes = str(rec.get("pdes") or "").strip()
    full_name = str(rec.get("full_name") or pdes).strip()
    return {
        "pdes": pdes,
        "desig": designation(pdes, full_name),
        "full_name": full_name,
        "spkid": str(rec.get("spkid") or ""),
        "a_au": a,
        "e": e,
        "i_deg": _num(rec.get("i")),
        "om_deg": _num(rec.get("om")),
        "w_deg": _num(rec.get("w")),
        "ma_deg": _num(rec.get("ma")),
        "q_au": _num(rec.get("q")) if rec.get("q") not in (None, "") else a * (1.0 - e),
        "condition_code": _int(rec.get("condition_code")),
        "data_arc_days": arc,
        "n_obs": _num(rec.get("n_obs_used")),
        "first_obs": first,
        "last_obs": last,
        "source": rec.get("source"),
        "n_oppositions": n_oppositions_est(first, last, arc),
        "sigma": sigma,
    }


def fetch_cen_catalog(force: bool = False) -> dict:
    if CATALOG_PATH.exists() and not force:
        return json.loads(CATALOG_PATH.read_text(encoding="utf-8"))
    q = urllib.parse.urlencode(
        {"sb-class": "CEN", "full-prec": "true", "fields": FIELDS}
    )
    url = f"{QUERY_API}?{q}"
    with urllib.request.urlopen(url, timeout=120) as resp:
        raw = json.loads(resp.read().decode("utf-8"))
    fields = raw["fields"]
    objects = []
    for values in raw["data"]:
        rec = row_to_modern(fields, values)
        if rec is None:
            continue
        objects.append(rec)
    blob = {
        "source": "JPL SBDB Query API sb-class=CEN",
        "n": len(objects),
        "n_oppositions_note": (
            "SBDB Query has no n_opps field; n_oppositions is the inclusive "
            "calendar-year span of first_obs/last_obs (fallback: 1+floor(arc/365.25))."
        ),
        "objects": objects,
    }
    CATALOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    CATALOG_PATH.write_text(json.dumps(blob, indent=2), encoding="utf-8")
    return blob


def _arc_bin(arc: float | None, edges: np.ndarray) -> str:
    if arc is None or not np.isfinite(arc) or arc <= 0:
        return "arc_missing"
    # edges are log10 cut points for tertiles among valid arcs
    loga = math.log10(arc)
    if loga < edges[0]:
        return "arc_short"
    if loga < edges[1]:
        return "arc_medium"
    return "arc_long"


def _cc_bin(cc: int | None) -> str:
    if cc is None:
        return "U_missing"
    if cc <= 2:
        return "U_0-2"
    if cc <= 5:
        return "U_3-5"
    return "U_6-9"


def _opp_bin(n_opp: int | None) -> str:
    if n_opp is None:
        return "opp_missing"
    if n_opp <= 1:
        return "opp_1"
    if n_opp <= 3:
        return "opp_2-3"
    return "opp_4+"


def _eligible(obj: dict) -> bool:
    """Need angles + a,e to integrate; sigmas optional (floors exist)."""
    if obj.get("a_au") is None or obj.get("e") is None:
        return False
    if obj.get("i_deg") is None or obj.get("om_deg") is None:
        return False
    if obj.get("w_deg") is None or obj.get("ma_deg") is None:
        return False
    if not (5.0 <= float(obj["a_au"]) <= 35.0):
        return False
    if not (0.0 <= float(obj["e"]) < 0.95):
        return False
    if obj.get("data_arc_days") is None or float(obj["data_arc_days"]) <= 0:
        return False
    return True


def build_stratified_sample(catalog: dict, n_target: int = 175, seed: int = SEED) -> dict:
    rng = np.random.default_rng(seed)
    eligible = [o for o in catalog["objects"] if _eligible(o)]
    arcs = np.array([math.log10(float(o["data_arc_days"])) for o in eligible])
    edges = np.quantile(arcs, [1.0 / 3.0, 2.0 / 3.0])
    buckets: dict[str, list[dict]] = defaultdict(list)
    for o in eligible:
        key = "|".join(
            (
                _arc_bin(o.get("data_arc_days"), edges),
                _cc_bin(o.get("condition_code")),
                _opp_bin(o.get("n_oppositions")),
            )
        )
        rec = dict(o)
        rec["stratum"] = key
        rec["arc_bin"] = _arc_bin(o.get("data_arc_days"), edges)
        rec["cc_bin"] = _cc_bin(o.get("condition_code"))
        rec["opp_bin"] = _opp_bin(o.get("n_oppositions"))
        buckets[key].append(rec)

    n_target = min(n_target, len(eligible))
    keys = sorted(buckets)
    # At least one object per non-empty stratum, then proportional fill.
    chosen: list[dict] = []
    leftover_slots = n_target
    per_stratum_pick: dict[str, list[dict]] = {}
    for k in keys:
        pool = list(buckets[k])
        rng.shuffle(pool)
        take = 1 if leftover_slots > 0 else 0
        per_stratum_pick[k] = pool[:take]
        leftover_slots -= take
    remaining_pools = {k: buckets[k][1:] for k in keys}
    weights = np.array([len(buckets[k]) for k in keys], dtype=float)
    weights = weights / weights.sum()
    extra = leftover_slots
    # Largest remainder
    raw = extra * weights
    base = np.floor(raw).astype(int)
    rem = extra - int(base.sum())
    order = np.argsort(-(raw - base))
    alloc = {keys[i]: int(base[i]) for i in range(len(keys))}
    for j in range(rem):
        alloc[keys[int(order[j])]] += 1
    for k in keys:
        pool = remaining_pools[k]
        rng.shuffle(pool)
        n_take = min(alloc[k], len(pool))
        per_stratum_pick[k].extend(pool[:n_take])
    for k in keys:
        chosen.extend(per_stratum_pick[k])
    # If still short (tiny strata), fill from unused eligible
    if len(chosen) < n_target:
        used = {o["desig"] for o in chosen}
        rest = [o for o in eligible if o["desig"] not in used]
        rng.shuffle(rest)
        for o in rest:
            if len(chosen) >= n_target:
                break
            rec = dict(o)
            rec["stratum"] = "|".join(
                (
                    _arc_bin(o.get("data_arc_days"), edges),
                    _cc_bin(o.get("condition_code")),
                    _opp_bin(o.get("n_oppositions")),
                )
            )
            rec["arc_bin"] = _arc_bin(o.get("data_arc_days"), edges)
            rec["cc_bin"] = _cc_bin(o.get("condition_code"))
            rec["opp_bin"] = _opp_bin(o.get("n_oppositions"))
            chosen.append(rec)
    chosen = chosen[:n_target]
    chosen.sort(key=lambda o: (o["arc_bin"], o.get("data_arc_days") or 0, o["desig"]))
    blob = {
        "seed": seed,
        "n_target": n_target,
        "n_catalog": catalog["n"],
        "n_eligible": len(eligible),
        "log10_arc_tertile_edges": [float(x) for x in edges],
        "arc_tertile_days": [float(10 ** float(x)) for x in edges],
        "n_strata": len(keys),
        "stratum_counts_eligible": {k: len(buckets[k]) for k in keys},
        "stratum_counts_sample": dict(Counter(o["stratum"] for o in chosen)),
        "objects": chosen,
    }
    SAMPLE_PATH.write_text(json.dumps(blob, indent=2), encoding="utf-8")
    return blob


def pick_pilot_objects(sample: dict, n: int = 5) -> list[dict]:
    objs = sorted(sample["objects"], key=lambda o: float(o["data_arc_days"]))
    if len(objs) < n:
        return objs
    idx = np.linspace(0, len(objs) - 1, n).round().astype(int)
    picked = []
    seen = set()
    for i in idx:
        o = objs[int(i)]
        if o["desig"] in seen:
            continue
        seen.add(o["desig"])
        picked.append(o)
    # Ensure we have n distinct
    for o in objs:
        if len(picked) >= n:
            break
        if o["desig"] not in seen:
            picked.append(o)
            seen.add(o["desig"])
    return picked[:n]


def object_summary(blob: dict) -> dict:
    clones = blob.get("clones") or []
    n = len(clones)
    classes = [c.get("class") for c in clones]
    counts = {k: classes.count(k) for k in ("D", "R", "Q")}
    modal = max(counts, key=counts.get) if n else None
    lives = [float(c["lifetime_myr"]) for c in clones if c.get("lifetime_myr") is not None]
    n_cens = sum(1 for c in clones if c.get("survived_full"))
    n_esc = sum(1 for c in clones if c.get("escaped"))
    n_long = sum(1 for c in clones if c.get("long_lived") or float(c.get("lifetime_myr") or 0) >= 22.0)
    nom = next((c for c in clones if int(c.get("clone", -1)) == 0), None)
    nom_class = nom.get("class") if nom else None
    minority_rq = (counts["R"] + counts["Q"]) / n if n else None
    disagree_nominal = (
        float(nom_class != modal) if nom_class is not None and modal is not None else None
    )
    modern = blob.get("modern_nominal") or {}
    return {
        "desig": blob.get("desig"),
        "campaign": blob.get("campaign"),
        "clone_scheme": blob.get("clone_scheme", "diagonal"),
        "n_clones": n,
        "t_max_myr": blob.get("t_max_myr"),
        "modal_class": modal,
        "counts_D": counts["D"],
        "counts_R": counts["R"],
        "counts_Q": counts["Q"],
        "frac_modal": counts[modal] / n if n and modal else None,
        "frac_ge22": n_long / n if n else None,
        "minority_rq_frac": minority_rq,
        "median_lifetime_myr": float(np.median(lives)) if lives else None,
        "censoring_count": n_cens,
        "censoring_rate": n_cens / n if n else None,
        "n_escaped": n_esc,
        "nominal_class": nom_class,
        "disagree_vs_nominal": disagree_nominal,
        "data_arc_days": modern.get("data_arc_days"),
        "condition_code": modern.get("condition_code"),
        "n_oppositions": modern.get("n_oppositions"),
        "n_obs": modern.get("n_obs"),
        "arc_bin": modern.get("arc_bin"),
        "cc_bin": modern.get("cc_bin"),
        "opp_bin": modern.get("opp_bin"),
        "stratum": modern.get("stratum"),
        "a_au": modern.get("a_au"),
        "e": modern.get("e"),
        "elapsed_s": blob.get("elapsed_s"),
    }


def spearman_table(rows: list[dict]) -> list[dict]:
    from scipy.stats import spearmanr

    predictors = {
        "log10_arc_days": lambda r: math.log10(r["data_arc_days"])
        if r.get("data_arc_days") and r["data_arc_days"] > 0
        else None,
        "condition_code": lambda r: r.get("condition_code"),
        "n_oppositions": lambda r: r.get("n_oppositions"),
    }
    outcomes = {
        "minority_rq_frac": lambda r: r.get("minority_rq_frac"),
        "disagree_vs_nominal": lambda r: r.get("disagree_vs_nominal"),
        "censoring_rate": lambda r: r.get("censoring_rate"),
        "frac_ge22": lambda r: r.get("frac_ge22"),
    }
    out = []
    for pk, pf in predictors.items():
        for ok, of in outcomes.items():
            xs, ys = [], []
            for r in rows:
                x, y = pf(r), of(r)
                if x is None or y is None:
                    continue
                try:
                    xf, yf = float(x), float(y)
                except (TypeError, ValueError):
                    continue
                if not (np.isfinite(xf) and np.isfinite(yf)):
                    continue
                xs.append(xf)
                ys.append(yf)
            n = len(xs)
            if n < 5:
                out.append(
                    {"predictor": pk, "outcome": ok, "n": n, "rho": None, "p": None, "note": "n<5"}
                )
                continue
            rho, p = spearmanr(xs, ys)
            out.append(
                {
                    "predictor": pk,
                    "outcome": ok,
                    "n": n,
                    "rho": float(rho) if rho is not None and np.isfinite(rho) else None,
                    "p": float(p) if p is not None and np.isfinite(p) else None,
                }
            )
    return out


def write_summary_csv(rows: list[dict], path: Path) -> None:
    if not rows:
        return
    keys = list(rows[0].keys())
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def make_plots(phase1_rows: list[dict], cov_rows: list[dict] | None, figdir: Path) -> list[Path]:
    import matplotlib.pyplot as plt

    figdir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 12,
        }
    )
    written: list[Path] = []

    # (a) minority R/Q vs arc
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    y = np.array(
        [
            r["minority_rq_frac"]
            for r in phase1_rows
            if r.get("data_arc_days") is not None and r.get("minority_rq_frac") is not None
        ]
    )
    x = np.array(
        [
            r["data_arc_days"]
            for r in phase1_rows
            if r.get("data_arc_days") is not None and r.get("minority_rq_frac") is not None
        ]
    )
    ax.scatter(x, y, s=22, c="#1f4e79", alpha=0.75, edgecolors="none")
    if len(x) >= 5:
        lx = np.log10(x)
        coeff = np.polyfit(lx, y, 1)
        xs = np.logspace(np.log10(x.min()), np.log10(x.max()), 80)
        ax.plot(xs, np.polyval(coeff, np.log10(xs)), color="#c45c26", lw=1.6, label="linear fit in log10(arc)")
        ax.legend(frameon=False)
    ax.set_xscale("log")
    ax.set_xlabel("Arc length (days)")
    ax.set_ylabel("R/Q clone fraction")
    ax.set_title("Phase 1: fragility vs arc length")
    ax.set_ylim(-0.05, 1.05)
    p = figdir / "fig_fragility_vs_arc.png"
    fig.tight_layout()
    fig.savefig(p, dpi=160)
    plt.close(fig)
    written.append(p)

    # (b) diagonal vs covariance
    if cov_rows:
        fig, ax = plt.subplots(figsize=(5.6, 5.2))
        dmap = {r["desig"]: r for r in phase1_rows}
        xs, ys, labels = [], [], []
        for c in cov_rows:
            d = dmap.get(c["desig"])
            if not d:
                continue
            xs.append(d["minority_rq_frac"])
            ys.append(c["minority_rq_frac"])
            labels.append(c["desig"])
        ax.plot([0, 1], [0, 1], color="0.7", lw=1)
        ax.scatter(xs, ys, s=36, c="#c45c26", zorder=3)
        ax.set_xlabel("Diagonal R/Q fraction")
        ax.set_ylabel("Covariance R/Q fraction")
        ax.set_title("Short-arc subsample: diagonal vs covariance")
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)
        ax.set_aspect("equal")
        p = figdir / "fig_diag_vs_cov.png"
        fig.tight_layout()
        fig.savefig(p, dpi=160)
        plt.close(fig)
        written.append(p)

    # (c) modal agreement with nominal
    fig, ax = plt.subplots(figsize=(5.8, 4.0))
    agree = [
        1.0 - float(r["disagree_vs_nominal"])
        for r in phase1_rows
        if r.get("disagree_vs_nominal") is not None
    ]
    ax.hist(agree, bins=[-0.05, 0.5, 1.05], color="#1f4e79", edgecolor="white")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["nominal ≠ mode", "nominal = mode"])
    ax.set_ylabel("Objects")
    ax.set_title("Phase 1: modal class vs clone-0 (nominal)")
    p = figdir / "fig_modal_agreement.png"
    fig.tight_layout()
    fig.savefig(p, dpi=160)
    plt.close(fig)
    written.append(p)
    return written

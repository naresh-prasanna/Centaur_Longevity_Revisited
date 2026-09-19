"""Assemble full-scale results into professional matrices + status markdown."""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np

from classify_core import reclassify_from_diagnostics

ROOT = Path(__file__).resolve().parents[1]
FS = ROOT / "results" / "fullscale"
RES = ROOT / "results"


def flat_of(rec: dict) -> bool:
    if rec.get("a_mean") is not None and rec.get("a_std") is not None:
        return float(rec["a_std"]) / (float(rec["a_mean"]) + 1e-9) < 0.03
    return bool(rec.get("flat"))


def label(rec: dict, rule: str = "v2") -> str:
    r = dict(rec)
    r["flat"] = flat_of(rec)
    return reclassify_from_diagnostics(r, rule)


def load_campaign(name: str) -> list[dict]:
    d = FS / name
    if not d.exists():
        return []
    out = []
    for f in sorted(d.glob("*.json")):
        out.append(json.loads(f.read_text(encoding="utf-8")))
    return out


def summarize_objs(objs: list[dict]) -> dict:
    if not objs:
        return {"n": 0, "flip_rate": None, "objects": []}
    rows = []
    for o in objs:
        # ensure v2 labels
        clones = o.get("clones") or []
        v2s = [label(c, "v2") for c in clones] if clones else [o.get("modal_class")]
        counts = dict(Counter(v2s))
        modal = Counter(v2s).most_common(1)[0][0]
        b = o.get("bailey_class") or o.get("bailey", {}).get("bailey_class")
        rows.append(
            {
                "desig": o["desig"],
                "bailey": b,
                "modal_v2": modal,
                "counts_v2": counts,
                "flip_v2": modal != b,
                "n_clones": len(clones) or o.get("n_clones"),
                "t_max_myr": o.get("t_max_myr"),
                "frac_long_v2": sum(1 for c in clones if c.get("lifetime_myr", 0) >= 22) / max(1, len(clones)),
                "median_life": float(np.median([c.get("lifetime_myr", 0) for c in clones])) if clones else None,
            }
        )
    return {
        "n": len(rows),
        "flip_rate": sum(1 for r in rows if r["flip_v2"]) / len(rows),
        "objects": rows,
    }


def main():
    qr = summarize_objs(load_campaign("qr20"))
    d5 = summarize_objs(load_campaign("d5"))
    th = summarize_objs(load_campaign("th173_40"))
    surv = summarize_objs(load_campaign("surv100"))

    # sensitivity summaries if present
    sens = {}
    for name in ("sens_q15", "sens_nohill"):
        sens[name] = summarize_objs(load_campaign(name))

    expcz = {}
    p = FS / "summary_exp_c_angles_zero.json"
    if p.exists():
        expcz = json.loads(p.read_text(encoding="utf-8"))

    # pilot / professional for comparison
    prof_path = RES / "professional_matrices.json"
    prof = json.loads(prof_path.read_text(encoding="utf-8")) if prof_path.exists() else {}

    out = {
        "default_rule": "v2",
        "qr20": qr,
        "d5_census": d5,
        "th173_40": th,
        "surv100": surv,
        "sensitivity": sens,
        "exp_c_angles_zero": {
            "n": expcz.get("n"),
            "n_orbit_drives": expcz.get("n_orbit_drives"),
        },
        "pilot_comparison": {
            "qr_flip_pilot_v2": prof.get("qr_pilot", {}).get("flip_rate_v2"),
            "qr_flip_full_v2": qr.get("flip_rate"),
            "d_stay_pilot_v2": prof.get("short_d_control", {}).get("stay_d_rate_v2"),
            "d_flip_census_v2": d5.get("flip_rate"),
        },
        "completeness": {
            "qr20_n": qr.get("n"),
            "qr20_target": 10,
            "d5_n": d5.get("n"),
            "d5_target": 51,
            "th173_40_n": th.get("n"),
            "surv100_n": surv.get("n"),
        },
    }
    dest = FS / "FULLSCALE_MATRICES.json"
    dest.write_text(json.dumps(out, indent=2), encoding="utf-8")

    lines = [
        "# Full-scale campaign status",
        "",
        f"Q/R 20-clone (40 Myr): **{qr.get('n')}/10** objects · flip_v2 = "
        f"**{(qr.get('flip_rate') or 0)*100:.0f}%**" if qr.get("n") else "Q/R 20-clone: not started",
        "",
        f"D census 5-clone (40 Myr): **{d5.get('n')}/51** · flip_v2 = "
        f"**{(d5.get('flip_rate') or 0)*100:.1f}%**" if d5.get("n") else "D census: not started",
        "",
        f"TH173 40 Myr: {th.get('n')} · surv100: {surv.get('n')}",
        "",
        "## Q/R hardened",
        "",
        "| Desig | BM09 | modal v2 | counts | frac≥22 Myr |",
        "|-------|------|----------|--------|-------------|",
    ]
    for r in qr.get("objects") or []:
        lines.append(
            f"| {r['desig']} | {r['bailey']} | {r['modal_v2']} | `{r['counts_v2']}` | {r['frac_long_v2']*100:.0f}% |"
        )
    lines += ["", "## D census flips (if any)", ""]
    flips = [r for r in (d5.get("objects") or []) if r["flip_v2"]]
    if not flips:
        lines.append("_none yet / none_")
    else:
        for r in flips:
            lines.append(f"- **{r['desig']}** D→{r['modal_v2']} `{r['counts_v2']}`")

    lines += [
        "",
        "## Completeness",
        "",
        "```json",
        json.dumps(out["completeness"], indent=2),
        "```",
        "",
        "Rebuild paper when qr20≥10 and d5≥40: `python src/build_fullscale_paper.py`",
    ]
    (FS / "FULLSCALE_STATUS.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(out["completeness"], indent=2))
    print("flip qr20", qr.get("flip_rate"), "d5", d5.get("flip_rate"))
    print("wrote", dest)


if __name__ == "__main__":
    main()

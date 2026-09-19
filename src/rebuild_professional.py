"""Rebuild professional result matrices under classifier v1 (legacy) and v2 (default).

Uses stored integration diagnostics only — no re-integration required.
Writes:
  results/professional_matrices.json
  results/PROFESSIONAL_STATUS.md
  results/paper_summary.json  (updated with dual-rule honesty)
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from classify_core import reclassify_from_diagnostics

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
DATA = ROOT / "data"


def flat_of(rec: dict) -> bool:
    if rec.get("a_mean") is not None and rec.get("a_std") is not None:
        return float(rec["a_std"]) / (float(rec["a_mean"]) + 1e-9) < 0.03
    return bool(rec.get("flat"))


def label(rec: dict, rule: str) -> str:
    # Recompute flat with current definition when a moments exist
    r = dict(rec)
    r["flat"] = flat_of(rec)
    return reclassify_from_diagnostics(r, rule=rule)


def modal(classes: list[str]) -> str | None:
    if not classes:
        return None
    return Counter(classes).most_common(1)[0][0]


def load_table():
    return json.loads((DATA / "bailey_table2.json").read_text(encoding="utf-8"))["objects"]


def load_modern():
    return json.loads((DATA / "modern_elements.json").read_text(encoding="utf-8"))


def batch_object(path: Path) -> dict | None:
    d = json.loads(path.read_text(encoding="utf-8"))
    if "clones" not in d:
        return None
    des = d.get("desig")
    if not des:
        stem = path.stem.replace("_", " ")
        des = stem
    clones = []
    for c in d["clones"]:
        if "lifetime_myr" not in c:
            continue
        clones.append(
            {
                "clone": c.get("clone"),
                "life": c["lifetime_myr"],
                "flat": flat_of(c),
                "nonlinear": c.get("nonlinear"),
                "H": c.get("H"),
                "r2": c.get("r2"),
                "stored_class": c.get("class"),
                "v1": label(c, "v1"),
                "v2": label(c, "v2"),
            }
        )
    return {"desig": des, "clones": clones, "file": path.name}


def confusion(rows: list[dict], key: str) -> dict:
    mat: dict[str, dict[str, int]] = {}
    for r in rows:
        b = r["bailey"]
        m = r[key]
        mat.setdefault(b, {})
        mat[b][m] = mat[b].get(m, 0) + 1
    return mat


def main():
    table = {o["desig"]: o for o in load_table()}
    modern = load_modern()

    # --- batch pilots ---
    batch_rows = []
    for path in sorted((RES / "batch").glob("*.json")):
        obj = batch_object(path)
        if not obj or not obj["clones"]:
            continue
        des = obj["desig"]
        # map filename stem if needed
        if des not in table:
            for cand in table:
                if cand.replace(" ", "_") == path.stem:
                    des = cand
                    break
        if des not in table:
            continue
        bclass = table[des]["bailey_class"]
        v1s = [c["v1"] for c in obj["clones"]]
        v2s = [c["v2"] for c in obj["clones"]]
        m1, m2 = modal(v1s), modal(v2s)
        batch_rows.append(
            {
                "desig": des,
                "bailey": bclass,
                "n_clones": len(obj["clones"]),
                "counts_v1": dict(Counter(v1s)),
                "counts_v2": dict(Counter(v2s)),
                "modal_v1": m1,
                "modal_v2": m2,
                "flip_v1": m1 != bclass,
                "flip_v2": m2 != bclass,
                "arc": table[des].get("arc"),
                "modern_a": modern.get(des, {}).get("a_au"),
                "modern_e": modern.get(des, {}).get("e"),
                "bailey_a": table[des]["a_au"],
                "bailey_e": table[des]["e"],
                "U": modern.get(des, {}).get("U"),
                "arc_days": modern.get(des, {}).get("arc_days"),
            }
        )

    qr = [r for r in batch_rows if r["bailey"] in ("Q", "R")]
    sd = [r for r in batch_rows if r["bailey"] == "D"]

    # --- Exp C ---
    expc_path = RES / "exp_c_summary_qr40.json"
    expc = []
    if expc_path.exists():
        raw = json.loads(expc_path.read_text(encoding="utf-8"))
        for r in raw["results"]:
            b_ic, m_ic = r["bailey_ic"], r["modern_ic"]
            expc.append(
                {
                    "desig": r["desig"],
                    "bailey_label": r["bailey_class"],
                    "life_bailey": b_ic["lifetime_myr"],
                    "life_modern": m_ic["lifetime_myr"],
                    "v1_bailey": label(b_ic, "v1"),
                    "v1_modern": label(m_ic, "v1"),
                    "v2_bailey": label(b_ic, "v2"),
                    "v2_modern": label(m_ic, "v2"),
                    "orbit_drives_v1": label(b_ic, "v1") != label(m_ic, "v1"),
                    "orbit_drives_v2": label(b_ic, "v2") != label(m_ic, "v2"),
                    "a_std_bailey": b_ic.get("a_std"),
                    "a_std_modern": m_ic.get("a_std"),
                    "r2_bailey": b_ic.get("r2"),
                    "r2_modern": m_ic.get("r2"),
                }
            )

    # --- TH173 production (10 Myr: Q impossible; R also needs 22 Myr) ---
    th = json.loads((RES / "th173_production.json").read_text(encoding="utf-8"))
    th_clones = []
    for c in th.get("clones", []):
        th_clones.append(
            {
                "clone": c.get("clone"),
                "life": c["lifetime_myr"],
                "stored": c.get("class"),
                "v1": label(c, "v1"),
                "v2": label(c, "v2"),
            }
        )

    # validation
    val_path = RES / "classifier_validation.json"
    validation = json.loads(val_path.read_text(encoding="utf-8")) if val_path.exists() else {}

    out = {
        "default_rule": "v2",
        "classifier_note": (
            "v2 = BM09 lifetime-faithful (long+flat→Q; long→R; short→D). "
            "v1 = legacy Hurst gate (long+r2≥0.85→D) — deprecated for headlines."
        ),
        "validation": {
            "v2_ok": validation.get("v2_ok"),
            "pass_rate_v2": validation.get("pass_rate_v2"),
            "pass_rate_v1": validation.get("pass_rate_v1"),
            "n": validation.get("n"),
        },
        "qr_pilot": {
            "n": len(qr),
            "flip_rate_v1": sum(1 for r in qr if r["flip_v1"]) / max(1, len(qr)),
            "flip_rate_v2": sum(1 for r in qr if r["flip_v2"]) / max(1, len(qr)),
            "matrix_v1": confusion(qr, "modal_v1"),
            "matrix_v2": confusion(qr, "modal_v2"),
            "objects": qr,
        },
        "short_d_control": {
            "n": len(sd),
            "stay_d_rate_v1": sum(1 for r in sd if r["modal_v1"] == "D") / max(1, len(sd)),
            "stay_d_rate_v2": sum(1 for r in sd if r["modal_v2"] == "D") / max(1, len(sd)),
            "flip_rate_v1": sum(1 for r in sd if r["flip_v1"]) / max(1, len(sd)),
            "flip_rate_v2": sum(1 for r in sd if r["flip_v2"]) / max(1, len(sd)),
            "objects": sd,
        },
        "exp_c": {
            "n": len(expc),
            "n_orbit_drives_v1": sum(1 for r in expc if r["orbit_drives_v1"]),
            "n_orbit_drives_v2": sum(1 for r in expc if r["orbit_drives_v2"]),
            "n_v2_bailey_R": sum(1 for r in expc if r["v2_bailey"] == "R"),
            "n_v2_modern_R": sum(1 for r in expc if r["v2_modern"] == "R"),
            "n_v2_bailey_Q": sum(1 for r in expc if r["v2_bailey"] == "Q"),
            "n_v2_modern_Q": sum(1 for r in expc if r["v2_modern"] == "Q"),
            "objects": expc,
        },
        "th173_production": {
            "n": len(th_clones),
            "t_max_myr": th.get("t_max_myr"),
            "counts_v1": dict(Counter(c["v1"] for c in th_clones)),
            "counts_v2": dict(Counter(c["v2"] for c in th_clones)),
            "note": "T_max=10 Myr cannot assign Q or R (both need ≥22 Myr). Result is D-rejection of Q only.",
        },
        "claims_allowed": [
            "TH173 is not Q under modern elements (100 clones / 10 Myr; none reach Q criteria).",
            "Classifier v1 produced a false 100% Q/R→D headline; do not use as discovery.",
            "Under v2, Q/R pilot flip rate is reduced; several long survivors recover as R.",
            "Exp C under v1 was all-D (instrument); under v2, QP112 is R/R and 2 objects show orbit-driven class change.",
            "Short-arc D control mostly stays D under both rules; SN55 is the notable D→Q exception.",
        ],
        "claims_withdrawn": [
            "100% BM09 Q/R → D flip rate as a dynamical discovery under the legacy classifier.",
            "Experiment C all-D as evidence that orbit revision alone does not matter (v1 was degenerate for R).",
        ],
    }

    dest = RES / "professional_matrices.json"
    dest.write_text(json.dumps(out, indent=2), encoding="utf-8")

    # paper_summary dual-rule
    paper = {
        "default_rule": "v2",
        "n_objects": len(batch_rows),
        "qr_flip_rate_v1": out["qr_pilot"]["flip_rate_v1"],
        "qr_flip_rate_v2": out["qr_pilot"]["flip_rate_v2"],
        "qr_matrix_v2": out["qr_pilot"]["matrix_v2"],
        "short_d_stay_v2": out["short_d_control"]["stay_d_rate_v2"],
        "short_d_flip_rate_v2": out["short_d_control"]["flip_rate_v2"],
        "exp_c_orbit_drives_v2": out["exp_c"]["n_orbit_drives_v2"],
        "th173_production": out["th173_production"],
        "validation_v2_ok": out["validation"].get("v2_ok"),
        "objects": batch_rows,
        "flips_v2": [r for r in batch_rows if r["flip_v2"]],
        "claims_allowed": out["claims_allowed"],
        "claims_withdrawn": out["claims_withdrawn"],
    }
    (RES / "paper_summary.json").write_text(json.dumps(paper, indent=2), encoding="utf-8")

    # Markdown status
    lines = [
        "# Professional status — classifier locked",
        "",
        f"**Default rule: v2.** Synthetic battery v2_ok = `{out['validation'].get('v2_ok')}` "
        f"(pass rate {out['validation'].get('pass_rate_v2')}).",
        "",
        "## Rules",
        "",
        "- **v2 (production):** life≥22 Myr & flat → Q; life≥22 → R; else D.",
        "- **v1 (legacy/sensitivity):** long + (nonlinear | r²<0.85) → R; else long + r²≥0.85 → D. "
        "**Withdrawn for headlines** — forces long linear survivors to D.",
        "",
        "## Q/R pilot (modern ICs)",
        "",
        f"| Metric | v1 | v2 |",
        f"|--------|---:|---:|",
        f"| n | {out['qr_pilot']['n']} | {out['qr_pilot']['n']} |",
        f"| modal flip vs BM09 | {out['qr_pilot']['flip_rate_v1']*100:.0f}% | {out['qr_pilot']['flip_rate_v2']*100:.0f}% |",
        "",
        "| Desig | BM09 | modal v1 | modal v2 | counts v2 |",
        "|-------|------|----------|----------|-----------|",
    ]
    for r in qr:
        lines.append(
            f"| {r['desig']} | {r['bailey']} | {r['modal_v1']} | {r['modal_v2']} | `{r['counts_v2']}` |"
        )

    lines += [
        "",
        "## Short-arc D control",
        "",
        f"Stay-D rate: v1 {out['short_d_control']['stay_d_rate_v1']*100:.0f}% / "
        f"v2 {out['short_d_control']['stay_d_rate_v2']*100:.0f}% "
        f"(n={out['short_d_control']['n']}).",
        "",
        "## Experiment C (BM09 aei vs modern aei, same code)",
        "",
        f"orbit_drives_flip: v1 **{out['exp_c']['n_orbit_drives_v1']}/10** · "
        f"v2 **{out['exp_c']['n_orbit_drives_v2']}/10**",
        "",
        "| Desig | BM09 | v1 B/M | v2 B/M | life B/M |",
        "|-------|------|--------|--------|----------|",
    ]
    for r in expc:
        lines.append(
            f"| {r['desig']} | {r['bailey_label']} | {r['v1_bailey']}/{r['v1_modern']} | "
            f"{r['v2_bailey']}/{r['v2_modern']} | "
            f"{r['life_bailey']:.2f}/{r['life_modern']:.2f} |"
        )

    lines += [
        "",
        "## TH173 production",
        "",
        f"n={out['th173_production']['n']}, T_max={out['th173_production']['t_max_myr']} Myr. "
        f"counts v2 = `{out['th173_production']['counts_v2']}`. "
        f"{out['th173_production']['note']}",
        "",
        "## Allowed claims",
        "",
    ]
    for c in out["claims_allowed"]:
        lines.append(f"- {c}")
    lines += ["", "## Withdrawn claims", ""]
    for c in out["claims_withdrawn"]:
        lines.append(f"- ~~{c}~~")

    lines += [
        "",
        "## Still needed for a longer / RP-grade paper",
        "",
        "1. 20–100 clones on full Q/R set under v2",
        "2. 100 Myr integrations for objects that survive 40 Myr",
        "3. Escape-criterion and angle-completion sensitivity",
        "4. Full Table 2 D census under v2",
        "5. Honest journal rewrite using these matrices (not the legacy 100% flip headline)",
        "",
        f"Machine-readable: `{dest.name}`, `paper_summary.json`, `classifier_validation.json`.",
    ]
    status = RES / "PROFESSIONAL_STATUS.md"
    status.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"wrote {dest}")
    print(f"wrote {status}")
    print(
        f"QR flip v1/v2 = {out['qr_pilot']['flip_rate_v1']:.2f}/{out['qr_pilot']['flip_rate_v2']:.2f} "
        f"| ExpC orbit_drives v1/v2 = {out['exp_c']['n_orbit_drives_v1']}/{out['exp_c']['n_orbit_drives_v2']} "
        f"| val_v2_ok = {out['validation'].get('v2_ok')}"
    )


if __name__ == "__main__":
    main()

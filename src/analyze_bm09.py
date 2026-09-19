"""Analysis and figures for the BM09 Table 2 reclassification.

Produces the confusion matrix, the control/test decomposition with a Fisher exact
test, a per-object table in CSV and LaTeX, and three diagnostic figures.

The control/test split is the inferential core. Objects BM09 already had at multiple
oppositions have essentially unchanged orbits today, so their label changes measure
how much this pipeline differs from theirs. Short-arc objects carry that same
pipeline term plus the orbit revision. The difference in change rate between the two
groups is the orbit-driven component; the Fisher test asks whether it is real.

Usage:
  python src/analyze_bm09.py
  python src/analyze_bm09.py --campaign clones_censor_resolved
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "bm09_reclass"
FIGDIR = OUT / "figures"
CLASSES = ("D", "R", "Q")


def load_rows(campaign: str) -> list[dict]:
    rows = []
    for path in sorted((OUT / campaign).glob("*.json")):
        if path.name.endswith(".partial.json") or path.name.startswith("_"):
            continue
        blob = json.loads(path.read_text(encoding="utf-8"))
        if "clones" not in blob:
            continue
        m = blob["modern_nominal"]
        clones = blob["clones"]
        lives = [float(c["lifetime_myr"]) for c in clones if c.get("lifetime_myr") is not None]
        counts = blob["counts"]
        n = len(clones)
        shift = m.get("orbit_shift") or {}
        rows.append({
            "desig": blob["desig"],
            "bailey_class": m.get("bailey_class"),
            "modern_class": blob["modal_class"],
            "changed": int(m.get("bailey_class") != blob["modal_class"]),
            "frac_modal": blob["frac_modal"],
            "n_D": counts.get("D", 0), "n_R": counts.get("R", 0), "n_Q": counts.get("Q", 0),
            "n_clones": n,
            "median_lifetime_myr": float(np.median(lives)) if lives else None,
            "min_lifetime_myr": min(lives) if lives else None,
            "max_lifetime_myr": max(lives) if lives else None,
            "frac_ge22": sum(1 for x in lives if x >= 22.0) / n if n else None,
            "minority_rq_frac": (counts.get("R", 0) + counts.get("Q", 0)) / n if n else None,
            "n_censored": sum(1 for c in clones if c.get("survived_full")),
            "censoring_rate": sum(1 for c in clones if c.get("survived_full")) / n if n else None,
            "orbit_2007_quality": m.get("orbit_2007_quality"),
            "arc_2007_raw": m.get("arc_2007_raw"),
            "arc_2007_days": m.get("arc_2007_days"),
            "arc_2007_opp": m.get("arc_2007_opp"),
            "arc_days_now": m.get("data_arc_days"),
            "condition_code_now": m.get("condition_code"),
            "n_obs_now": m.get("n_obs"),
            "a_au": m.get("a_au"), "e": m.get("e"), "q_au": m.get("q_au"),
            "bailey_a_au": (m.get("bailey_elements") or {}).get("a_au"),
            "bailey_e": (m.get("bailey_elements") or {}).get("e"),
            "bailey_q_au": (m.get("bailey_elements") or {}).get("q_au"),
            "de": shift.get("de"), "frac_dq": shift.get("frac_dq"),
            "shift_index": shift.get("shift_index"),
        })
    return rows


def fisher_exact_2x2(a: int, b: int, c: int, d: int) -> tuple[float, float]:
    """Odds ratio and two-sided p for [[a,b],[c,d]]. Falls back to scipy if present."""
    try:
        from scipy.stats import fisher_exact

        orr, p = fisher_exact([[a, b], [c, d]])
        return float(orr), float(p)
    except Exception:
        pass
    n = a + b + c + d

    def hyp(x):
        return (math.comb(a + b, x) * math.comb(c + d, a + c - x)) / math.comb(n, a + c)

    lo = max(0, a + c - (c + d))
    hi = min(a + b, a + c)
    p_obs = hyp(a)
    p = sum(hyp(x) for x in range(lo, hi + 1) if hyp(x) <= p_obs + 1e-12)
    orr = float("inf") if b * c == 0 else (a * d) / (b * c)
    return orr, min(1.0, p)


def analyse(rows: list[dict]) -> dict:
    conf = {bc: {mc: 0 for mc in CLASSES} for bc in CLASSES}
    for r in rows:
        if r["bailey_class"] in conf and r["modern_class"] in CLASSES:
            conf[r["bailey_class"]][r["modern_class"]] += 1

    ctrl = [r for r in rows if r["orbit_2007_quality"] == "multi_opposition"]
    test = [r for r in rows if r["orbit_2007_quality"] == "short_arc"]
    ca = sum(r["changed"] for r in ctrl)
    ta = sum(r["changed"] for r in test)
    orr, p = fisher_exact_2x2(ta, len(test) - ta, ca, len(ctrl) - ca) if ctrl and test else (None, None)

    per_class = {}
    for bc in CLASSES:
        sub = [r for r in rows if r["bailey_class"] == bc]
        if sub:
            per_class[bc] = {
                "n": len(sub),
                "n_retained": sum(1 for r in sub if not r["changed"]),
                "retention_rate": sum(1 for r in sub if not r["changed"]) / len(sub),
                "reassigned_to": {mc: sum(1 for r in sub if r["modern_class"] == mc)
                                  for mc in CLASSES},
            }

    # Restricted comparison. BM09's R means resonance hopping; the v2 rule can only
    # emit R for clones living >=22 Myr, so the two R's are different quantities and
    # the R objects cannot enter a like-for-like test. D and Q are comparable.
    dq = [r for r in rows if r["bailey_class"] in ("D", "Q") and r["modern_class"] in ("D", "Q")]
    dq_ctrl = [r for r in dq if r["orbit_2007_quality"] == "multi_opposition"]
    dq_test = [r for r in dq if r["orbit_2007_quality"] == "short_arc"]
    dqa, dqb = sum(r["changed"] for r in dq_test), sum(r["changed"] for r in dq_ctrl)
    dq_or, dq_p = (fisher_exact_2x2(dqa, len(dq_test) - dqa, dqb, len(dq_ctrl) - dqb)
                   if dq_ctrl and dq_test else (None, None))

    out = {
        "n_objects": len(rows),
        "confusion_bailey_to_modern": conf,
        "restricted_DQ_axis": {
            "rationale": (
                "BM09's R (resonance hopping) and the v2 R (survives >=22 Myr but a(t) not "
                "flat) are different quantities, so R objects are excluded from the "
                "like-for-like comparison. D and Q are directly comparable."
            ),
            "n": len(dq),
            "n_changed": sum(r["changed"] for r in dq),
            "changed": [{"desig": r["desig"], "from": r["bailey_class"], "to": r["modern_class"],
                         "arc_2007": r["arc_2007_raw"], "arc_days_now": r["arc_days_now"],
                         "median_lifetime_myr": r["median_lifetime_myr"],
                         "frac_modal": r["frac_modal"]}
                        for r in dq if r["changed"]],
            "multi_opposition": {"n": len(dq_ctrl), "n_changed": dqb,
                                 "frac_changed": dqb / len(dq_ctrl) if dq_ctrl else None},
            "short_arc": {"n": len(dq_test), "n_changed": dqa,
                          "frac_changed": dqa / len(dq_test) if dq_test else None},
            "fisher": {"odds_ratio": dq_or, "p_two_sided": dq_p},
        },
        "per_bailey_class": per_class,
        "n_changed": sum(r["changed"] for r in rows),
        "frac_changed": sum(r["changed"] for r in rows) / len(rows) if rows else None,
        "control_multi_opposition": {
            "n": len(ctrl), "n_changed": ca,
            "frac_changed": ca / len(ctrl) if ctrl else None,
            "desigs_changed": [r["desig"] for r in ctrl if r["changed"]],
        },
        "test_short_arc": {
            "n": len(test), "n_changed": ta,
            "frac_changed": ta / len(test) if test else None,
            "desigs_changed": [r["desig"] for r in test if r["changed"]],
        },
        "fisher_short_arc_vs_control": {"odds_ratio": orr, "p_two_sided": p},
        "n_fully_censored": sum(1 for r in rows if r["n_censored"] == r["n_clones"]),
        "fully_censored_desigs": [r["desig"] for r in rows if r["n_censored"] == r["n_clones"]],
    }
    if ctrl and test:
        out["excess_change_attributable_to_orbit"] = (ta / len(test)) - (ca / len(ctrl))
    return out


def write_tables(rows: list[dict], res: dict, tag: str) -> list[Path]:
    written = []
    OUT.mkdir(parents=True, exist_ok=True)

    csv_path = OUT / f"bm09_per_object{tag}.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    written.append(csv_path)

    conf = res["confusion_bailey_to_modern"]
    lines = [
        "% Auto-generated by src/analyze_bm09.py -- do not hand-edit.",
        r"\begin{tabular}{lrrrr}",
        r"\hline",
        r"BM09 class & \multicolumn{3}{c}{This work (modal of 10 clones)} & Retained \\",
        r" & D & R & Q & \\",
        r"\hline",
    ]
    for bc in CLASSES:
        if bc not in res["per_bailey_class"]:
            continue
        pc = res["per_bailey_class"][bc]
        lines.append(
            f"{bc} (n={pc['n']}) & {conf[bc]['D']} & {conf[bc]['R']} & {conf[bc]['Q']} "
            f"& {pc['retention_rate']*100:.0f}\\% \\\\"
        )
    lines += [r"\hline", r"\end{tabular}"]
    tex_path = OUT / f"bm09_confusion{tag}.tex"
    tex_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    written.append(tex_path)

    (OUT / f"bm09_analysis{tag}.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
    written.append(OUT / f"bm09_analysis{tag}.json")
    return written


def make_figures(rows: list[dict], res: dict, tag: str) -> list[Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIGDIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "serif", "font.size": 10,
                         "axes.labelsize": 11, "axes.titlesize": 11})
    written = []

    # (a) change rate, control vs short-arc
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    ctrl, test = res["control_multi_opposition"], res["test_short_arc"]
    labels = [f"Multi-opposition\n(n={ctrl['n']})", f"Short arc\n(n={test['n']})"]
    vals = [(ctrl["frac_changed"] or 0) * 100, (test["frac_changed"] or 0) * 100]
    bars = ax.bar(labels, vals, color=["#4C72B0", "#C44E52"], width=0.55)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 1.5, f"{v:.0f}%", ha="center", fontsize=10)
    ax.set_ylabel("Objects whose class changed (\\%)")
    ax.set_title("2007 orbit quality vs.\\ classification change")
    ax.set_ylim(0, max(vals + [10]) * 1.35)
    p = res["fisher_short_arc_vs_control"]["p_two_sided"]
    if p is not None:
        ax.text(0.5, 0.94, f"Fisher exact $p$ = {p:.3f}", transform=ax.transAxes,
                ha="center", fontsize=9, style="italic")
    fig.tight_layout()
    f1 = FIGDIR / f"bm09_change_by_arc_quality{tag}.pdf"
    fig.savefig(f1, dpi=200); plt.close(fig); written.append(f1)

    # (b) clone lifetime spread per object, ordered by median
    sub = [r for r in rows if r["median_lifetime_myr"] is not None]
    sub.sort(key=lambda r: r["median_lifetime_myr"])
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    x = np.arange(len(sub))
    med = [r["median_lifetime_myr"] for r in sub]
    lo = [r["median_lifetime_myr"] - r["min_lifetime_myr"] for r in sub]
    hi = [r["max_lifetime_myr"] - r["median_lifetime_myr"] for r in sub]
    colors = {"D": "#4C72B0", "R": "#DD8452", "Q": "#55A868"}
    ax.errorbar(x, med, yerr=[lo, hi], fmt="none", ecolor="0.75", elinewidth=0.9, capsize=0)
    for cls in CLASSES:
        idx = [i for i, r in enumerate(sub) if r["modern_class"] == cls]
        if idx:
            ax.scatter(idx, [med[i] for i in idx], s=22, color=colors[cls],
                       label=f"{cls} (this work)", zorder=3)
    ax.axhline(22.0, ls="--", lw=1.0, color="0.35")
    ax.text(0.01, 22.6, "Q threshold, 22 Myr", fontsize=8.5, color="0.35")
    ax.set_xlabel("BM09 objects, ordered by median clone lifetime")
    ax.set_ylabel("Clone lifetime (Myr)")
    ax.set_title("Clone lifetime range per object (bars span min to max)")
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    fig.tight_layout()
    f2 = FIGDIR / f"bm09_lifetime_spread{tag}.pdf"
    fig.savefig(f2, dpi=200); plt.close(fig); written.append(f2)

    # (c) how far the orbit moved vs whether the label moved
    pts = [r for r in rows if r.get("shift_index") is not None]
    if pts:
        fig, ax = plt.subplots(figsize=(5.2, 3.8))
        for changed, color, lab in ((0, "#4C72B0", "class retained"), (1, "#C44E52", "class changed")):
            g = [r for r in pts if r["changed"] == changed]
            if g:
                ax.scatter([max(r["shift_index"], 1e-5) for r in g],
                           [r["minority_rq_frac"] for r in g],
                           s=26, color=color, alpha=0.8, label=lab)
        ax.set_xscale("log")
        ax.set_xlabel(r"Orbit shift since BM09,  $|\Delta e| + |\Delta q / q|$")
        ax.set_ylabel("Minority R/Q clone fraction")
        ax.set_title("Orbit revision vs.\\ classification instability")
        ax.legend(frameon=False, fontsize=9)
        fig.tight_layout()
        f3 = FIGDIR / f"bm09_shift_vs_fragility{tag}.pdf"
        fig.savefig(f3, dpi=200); plt.close(fig); written.append(f3)

    return written


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--campaign", default="clones")
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args()

    rows = load_rows(args.campaign)
    if not rows:
        print(f"no completed objects in {OUT / args.campaign}")
        return
    tag = "" if args.campaign == "clones" else f"_{args.campaign}"
    res = analyse(rows)
    written = write_tables(rows, res, tag)
    if not args.no_figures:
        written += make_figures(rows, res, tag)

    print(json.dumps(res, indent=2))
    print("\nwrote:")
    for p in written:
        print(f"  {p.relative_to(ROOT).as_posix()}")


if __name__ == "__main__":
    main()

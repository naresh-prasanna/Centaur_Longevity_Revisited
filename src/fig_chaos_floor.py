"""Headline figure: outcome dispersion is set by dynamics, not by measurement.

Left  - censoring-aware dispersion against orbit quality, coloured by
        inclination, the one covariate that does predict it.
Right - the same after conditioning on (a, e, i, q), with the envelope of
        slopes the data exclude at 95%.
"""
from __future__ import annotations

import glob
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
from chaos_floor import load_objects, pool  # noqa: E402
from chaos_floor_final import ARC_MIN, CONTROLS, is_centaur  # noqa: E402
from chaos_floor_km import FLOOR_LIFE, TCAP, km_dispersion  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
LABEL_BOX = dict(boxstyle="round,pad=0.35", fc="white", ec="0.8", alpha=0.92)
FIG = ROOT / "paper" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "STIXGeneral", "Times New Roman", "Times", "serif"],
    "mathtext.fontset": "stix",
    "font.size": 11, "axes.labelsize": 12, "axes.titlesize": 12,
    "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 9,
})


def build():
    clones, caps = {}, {}
    for f in (glob.glob("results/pop_fragility/phase1/*.json")
              + glob.glob("results/pop_fragility/phase3_merged/*.json")
              + glob.glob("results/bm09_reclass/clones/*.json")
              + glob.glob("results/bm09_reclass/clones_censor_resolved/*.json")):
        if "partial" in f:
            continue
        b = json.load(open(f, encoding="utf-8"))
        if "clones" not in b:
            continue
        cap = float(b.get("t_max_myr") or TCAP)
        if cap >= caps.get(b["desig"], -1):
            clones[b["desig"]] = b["clones"]
            caps[b["desig"]] = cap
    rows = pool(load_objects(["results/pop_fragility/phase1/*.json",
                              "results/pop_fragility/phase3_merged/*.json"], "survey"),
                load_objects(["results/bm09_reclass/clones/*.json",
                              "results/bm09_reclass/clones_censor_resolved/*.json"], "bm09"))
    out = []
    for r in rows:
        if not (all(r.get(k) is not None for k in CONTROLS) and is_centaur(r)
                and (r["arc_days"] or 0) >= ARC_MIN and r["desig"] in clones):
            continue
        cl = clones[r["desig"]]
        t = [min(float(c["lifetime_myr"]), caps.get(r["desig"], TCAP)) for c in cl]
        e = [0 if c.get("survived_full") else 1 for c in cl]
        if min(t) <= FLOOR_LIFE:
            continue
        d = km_dispersion(t, e)
        if d is None:
            continue
        r = dict(r); r["km"] = d
        out.append(r)
    return out


def main() -> None:
    res = json.loads((ROOT / "results" / "chaos_floor" / "chaos_floor_km.json").read_text())
    s = build()
    x = np.array([r["log10_sigma_a_rel"] for r in s])
    y = np.array([r["km"] for r in s])
    inc = np.array([r["i_deg"] for r in s])

    floor = float(y.mean())
    boot = np.array([np.mean(np.random.default_rng(i).choice(y, len(y)))
                     for i in range(4000)])
    flo, fhi = float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))
    sg = res["regression"]["terms"]["log10_sigma_a_rel"]
    lim = max(abs(sg["ci95"][0]), abs(sg["ci95"][1]))

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(10.0, 4.3))

    sc = axL.scatter(x, y, c=inc, cmap="cividis", s=30, alpha=0.9,
                     edgecolor="k", linewidth=0.3, zorder=3)
    axL.axhspan(flo, fhi, color="0.55", alpha=0.22, zorder=1)
    axL.axhline(floor, color="0.2", lw=1.6, zorder=2)
    cb = fig.colorbar(sc, ax=axL, pad=0.02)
    cb.set_label("inclination (deg)", fontsize=10)
    axL.set_xlabel(r"$\log_{10}\,(\sigma_a/a)$")
    axL.set_ylabel(r"dispersion of $\log_{10}$ lifetime  (dex)")
    axL.set_title("(a) observed", loc="left")
    axL.annotate("", xy=(0.06, 0.90), xytext=(0.40, 0.90), xycoords="axes fraction",
                 arrowprops=dict(arrowstyle="->", color="0.35", lw=1.1))
    axL.text(0.42, 0.90, "better orbits", transform=axL.transAxes,
             va="center", ha="left", fontsize=9, color="0.3")
    axL.text(0.97, 0.05, f"floor {floor:.2f} dex\n(factor {10**floor:.1f} in lifetime)",
             transform=axL.transAxes, va="bottom", ha="right", fontsize=9, color="0.2",
             bbox=LABEL_BOX, zorder=5)

    C = np.column_stack([[r[c] for r in s] for c in CONTROLS] + [np.ones(len(s))])

    def resid(v):
        b, *_ = np.linalg.lstsq(C, v, rcond=None)
        return v - C @ b

    xr, yr = resid(x), resid(y) + floor
    axR.scatter(xr, yr, s=30, alpha=0.85, color="#2b6cb0",
                edgecolor="k", linewidth=0.3, zorder=3)
    axR.axhspan(flo, fhi, color="0.55", alpha=0.22, zorder=1)
    axR.axhline(floor, color="0.2", lw=1.6, zorder=2)
    xs = np.linspace(xr.min(), xr.max(), 50)
    axR.plot(xs, floor + lim * xs, ls="--", lw=1.2, color="#b03030", zorder=2)
    axR.plot(xs, floor - lim * xs, ls="--", lw=1.2, color="#b03030", zorder=2)
    axR.set_xlabel(r"$\log_{10}\,(\sigma_a/a)$, residual after $(a,e,i,q)$")
    axR.set_ylabel("dispersion, residual  (dex)")
    axR.set_title("(b) conditioned on dynamical state", loc="left")
    axR.text(0.03, 0.96,
             f"slope ${sg['coef']:+.3f}$ dex/decade\n"
             f"95% CI $[{sg['ci95'][0]:+.3f}, {sg['ci95'][1]:+.3f}]$\n"
             f"$p = {sg['p']:.2f}$",
             transform=axR.transAxes, va="top", ha="left", fontsize=9,
             bbox=LABEL_BOX, zorder=5)
    axR.legend(loc="lower right", framealpha=0.92, edgecolor="0.8", handles=[
        Line2D([], [], color="0.2", lw=1.6, label="measured floor"),
        Line2D([], [], color="#b03030", ls="--", lw=1.2, label="slopes excluded at 95%"),
    ]).set_zorder(5)

    # A shared vertical scale keeps the two floors visually comparable, with
    # headroom at the top so the annotation never lands on a point.
    lo = min(y.min(), yr.min()) - 0.06
    hi = max(y.max(), yr.max()) + 0.30
    for ax in (axL, axR):
        ax.set_ylim(lo, hi)
        ax.grid(alpha=0.22, lw=0.6)
        ax.set_axisbelow(True)

    fig.suptitle("Classification uncertainty is set by dynamics, not by measurement "
                 f"($n = {len(s)}$ Centaurs)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out = FIG / "fig07_chaos_floor.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"wrote {out.relative_to(ROOT).as_posix()}  (n={len(s)}, floor={floor:.3f} dex, "
          f"sigma span {x.min():.1f} to {x.max():.1f})")


if __name__ == "__main__":
    main()

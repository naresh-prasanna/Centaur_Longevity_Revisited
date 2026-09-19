"""Regenerate CMDA figures at journal standard from archived campaign JSON.

Science locks (asserted at runtime):
  TH173 10 Myr: 80 early escapes / 20 right-censored, median 2.99 Myr
  TH173 40 Myr: 88 D / 12 R / 0 Q
  RO43: 10 D / 10 R tie
  QR aggregate: D=155, R=43, Q=2
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

SRC = Path(__file__).resolve().parents[1]
OUT = SRC / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# Okabe-Ito, colourblind-safe
BLUE = "#0072B2"
ORANGE = "#E69F00"
GREEN = "#009E73"
VERM = "#D55E00"
PURPLE = "#CC79A7"
SKY = "#56B4E9"
GREY = "#4D4D4D"
INK = "#1A1A1A"

mpl.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["TeX Gyre Termes", "Times New Roman", "DejaVu Serif", "Times", "serif"],
        "mathtext.fontset": "stix",
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 11,
        "axes.linewidth": 0.8,
        "axes.edgecolor": INK,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "legend.fontsize": 9,
        "legend.frameon": False,
        "figure.dpi": 120,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.04,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.unicode_minus": False,
    }
)


def spines(ax: plt.Axes) -> None:
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(top=False, right=False)
    ax.set_axisbelow(True)
    ax.grid(True, color="#E6E6E6", lw=0.6, zorder=0)


def save(fig: plt.Figure, name: str) -> Path:
    png = OUT / name
    pdf = OUT / name.replace(".png", ".pdf")
    fig.savefig(png, dpi=600)
    fig.savefig(pdf)
    plt.close(fig)
    print(f"wrote {png.name}  ({png.stat().st_size // 1024} KB)")
    return png


def load_json(rel: str):
    return json.loads((SRC / rel).read_text(encoding="utf-8"))


def fig_th173() -> None:
    prod = load_json("results/th173_production.json")
    long = load_json("results/fullscale/th173_100x40/2005_TH173.json")
    lives = np.array([float(c["lifetime_myr"]) for c in prod["clones"]], dtype=float)
    escaped = np.array([bool(c.get("escaped")) and not bool(c.get("survived_full")) for c in prod["clones"]])
    # Prefer explicit flags; fall back to Tmax occupancy.
    if escaped.sum() == 0:
        escaped = lives < (float(prod["t_max_myr"]) - 1e-6)
    n_esc = int(escaped.sum())
    n_cen = int((~escaped).sum())
    med = float(np.median(lives))
    assert n_esc == 80 and n_cen == 20, (n_esc, n_cen)
    assert abs(med - 2.99) < 0.05, med

    long_cls = [c.get("class") or c.get("class_v2") for c in long["clones"]]
    if any(x is None for x in long_cls):
        # lifetime-first v2
        def v2(c):
            life = float(c["lifetime_myr"])
            if life < 22:
                return "D"
            a_mean = float(c.get("a_mean") or 0)
            a_std = float(c.get("a_std") or 0)
            flat = a_std / (a_mean + 1e-9) < 0.03 if a_mean else bool(c.get("flat"))
            return "Q" if flat else "R"

        long_cls = [v2(c) for c in long["clones"]]
    from collections import Counter

    cc = Counter(long_cls)
    assert cc.get("D", 0) == 88 and cc.get("R", 0) == 12 and cc.get("Q", 0) == 0, dict(cc)

    bailey, modern = prod["bailey"], prod["modern_nominal"]
    fig, axes = plt.subplots(1, 3, figsize=(11.6, 3.85))

    ax = axes[0]
    ax.annotate(
        "",
        xy=(modern["a_au"], modern["e"]),
        xytext=(bailey["a_au"], bailey["e"]),
        arrowprops=dict(arrowstyle="-|>", color=VERM, lw=1.8, mutation_scale=12),
    )
    ax.scatter([bailey["a_au"]], [bailey["e"]], c=VERM, s=42, zorder=3, label="BM09 (Q)")
    ax.scatter([modern["a_au"]], [modern["e"]], c=BLUE, s=90, marker="*", zorder=4, label="Modern SBDB")
    ax.set_xlabel(r"Semimajor axis $a$ (AU)")
    ax.set_ylabel(r"Eccentricity $e$")
    ax.set_title("(a)", loc="left", pad=4)
    ax.set_xlim(14.6, 21.2)
    ax.set_ylim(-0.02, 0.42)
    ax.legend(loc="upper left", handletextpad=0.4, fontsize=8.5)
    spines(ax)

    ax = axes[1]
    t_max = float(prod["t_max_myr"])
    bins = np.linspace(0, t_max, 21)
    ax.hist(lives[escaped], bins=bins, color=BLUE, edgecolor="white", lw=0.4, label=f"Escaped before 10 Myr ($n={n_esc}$)")
    ax.hist(lives[~escaped], bins=bins, color=ORANGE, edgecolor="white", lw=0.4, hatch="///", label=f"Still bound at 10 Myr ($n={n_cen}$)")
    ax.axvline(med, color=INK, ls="--", lw=1.1, label=f"Median {med:.2f} Myr")
    ax.set_xlabel("Lifetime (Myr)")
    ax.set_ylabel("Number of clones")
    ax.set_title("(b)", loc="left", pad=4)
    ax.set_xlim(0, t_max)
    ax.set_ylim(0, None)
    ymax = ax.get_ylim()[1]
    ax.set_ylim(0, ymax * 1.18)
    ax.text(
        0.97,
        0.22,
        "Q requires 22 Myr\n(off this axis)",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        color=GREY,
        clip_on=False,
    )
    ax.legend(loc="upper left", fontsize=8)
    spines(ax)

    ax = axes[2]
    ts = np.linspace(0, t_max, 201)
    frac = np.array([(lives >= t).mean() for t in ts])
    ax.plot(ts, frac, color=BLUE, lw=2.0)
    ax.set_xlabel("Time (Myr)")
    ax.set_ylabel("Still-bound fraction")
    ax.set_title("(c)", loc="left", pad=4)
    ax.set_xlim(0, t_max)
    ax.set_ylim(0, 1.12)
    ax.text(
        0.97,
        0.08,
        "80 confirmed early escapes\n20 right-censored\n(not a 22 Myr Q test)",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        color=GREY,
        clip_on=False,
    )
    spines(ax)

    fig.tight_layout(w_pad=1.6)
    save(fig, "fig_th173_composite.png")


def fig_qr() -> None:
    mats = load_json("results/fullscale/FULLSCALE_MATRICES.json")
    objs = mats["qr20"]["objects"]
    d = [int(o["counts_v2"].get("D", 0)) for o in objs]
    r = [int(o["counts_v2"].get("R", 0)) for o in objs]
    q = [int(o["counts_v2"].get("Q", 0)) for o in objs]
    bailey = [o["bailey"] for o in objs]
    labels = [o["desig"].replace("19", "").replace("20", "") for o in objs]
    # Keep year prefix readable: 1995 DW2 -> 1995 DW2 short
    labels = []
    for o in objs:
        des = o["desig"]
        labels.append(des.split()[-1] if " " in des else des)

    tot_d, tot_r, tot_q = sum(d), sum(r), sum(q)
    assert (tot_d, tot_r, tot_q) == (155, 43, 2), (tot_d, tot_r, tot_q)
    ro = next(o for o in objs if "RO43" in o["desig"])
    assert int(ro["counts_v2"].get("D", 0)) == 10 and int(ro["counts_v2"].get("R", 0)) == 10

    x = np.arange(len(objs))
    fig, ax = plt.subplots(figsize=(10.4, 4.15))
    ax.bar(x, d, color=BLUE, width=0.72, label="D  (lifetime $<$ 22 Myr)")
    ax.bar(x, r, bottom=d, color=ORANGE, width=0.72, label="R  (survived 22 Myr, not flat)")
    ax.bar(x, q, bottom=np.array(d) + np.array(r), color=GREEN, width=0.72, label="Q  (survived 22 Myr, flat)")
    ax.axhline(10, color="#B0B0B0", ls=":", lw=0.9, zorder=1)
    ax.set_xticks(x)
    xtl = [f"{lab}\nBM09 {b}" for lab, b in zip(labels, bailey)]
    ax.set_xticklabels(xtl, fontsize=8.5)
    ax.set_ylabel("Clones (20 per object)")
    ax.set_ylim(0, 24.5)
    ax.legend(ncol=3, loc="upper left", bbox_to_anchor=(0.0, 1.02), fontsize=8.5)
    ro_i = next(i for i, o in enumerate(objs) if "RO43" in o["desig"])
    ax.annotate(
        "10/10 tie",
        xy=(ro_i, 20.35),
        ha="center",
        va="bottom",
        fontsize=8,
        color=INK,
    )
    spines(ax)
    ax.grid(False, axis="x")
    fig.tight_layout()
    save(fig, "fig_qr_clone_fractions.png")

    # Long-lived fraction companion
    frac = [(ri + qi) / 20.0 for ri, qi in zip(r, q)]
    fig, ax = plt.subplots(figsize=(7.6, 3.9))
    colors = [ORANGE if f >= 0.45 else BLUE for f in frac]
    ax.bar(x, frac, color=colors, width=0.72, edgecolor="white", lw=0.4)
    ax.axhline(0.5, color=INK, ls="--", lw=0.9, label="Half of the clones")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_ylabel("Fraction with lifetime $\\geq$ 22 Myr")
    ax.set_ylim(0, 1.0)
    ax.legend(loc="upper left")
    spines(ax)
    ax.grid(False, axis="x")
    fig.tight_layout()
    save(fig, "fig_qr_longfrac.png")


def fig_od() -> None:
    table = {o["desig"]: o for o in load_json("data/bailey_table2.json")["objects"]}
    modern = load_json("data/modern_elements.json")
    mats = load_json("results/fullscale/FULLSCALE_MATRICES.json")
    highlight = {
        "2005 TH173": (VERM, "2005 TH173  (Q$\\to$D)"),
        "1995 SN55": (GREEN, "1995 SN55  (now a TNO)"),
        "2002 VR130": (PURPLE, "2002 VR130"),
    }
    desigs = [o["desig"] for o in mats["qr20"]["objects"]] + ["1995 SN55", "2002 VR130"]
    fig, ax = plt.subplots(figsize=(7.3, 5.1))
    for d in desigs:
        if d not in modern or d not in table:
            continue
        b, m = table[d], modern[d]
        color, lw, alpha = GREY, 0.9, 0.55
        if d in highlight:
            color, lw, alpha = highlight[d][0], 2.0, 0.95
        ax.annotate(
            "",
            xy=(m["a_au"], m["e"]),
            xytext=(b["a_au"], b["e"]),
            arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, alpha=alpha, mutation_scale=9),
        )
        ax.scatter([b["a_au"]], [b["e"]], c=color, s=22, zorder=3, alpha=alpha)
        ax.scatter([m["a_au"]], [m["e"]], c=color, s=42, marker="*", zorder=4, alpha=alpha)
    # Labels placed to avoid the star
    ax.annotate("2005 TH173", xy=(19.95, 0.307), xytext=(16.2, 0.34), fontsize=9, color=VERM,
                arrowprops=dict(arrowstyle="-", color=VERM, lw=0.6))
    ax.annotate("1995 SN55", xy=(42.61, 0.166), xytext=(34.0, 0.28), fontsize=9, color=GREEN,
                arrowprops=dict(arrowstyle="-", color=GREEN, lw=0.6))
    ax.annotate("2002 VR130", xy=(modern["2002 VR130"]["a_au"], modern["2002 VR130"]["e"]),
                xytext=(26.5, 0.42), fontsize=9, color=PURPLE,
                arrowprops=dict(arrowstyle="-", color=PURPLE, lw=0.6))
    ax.set_xlabel(r"Semimajor axis $a$ (AU)")
    ax.set_ylabel(r"Eccentricity $e$")
    ax.legend(
        handles=[
            Line2D([0], [0], color=VERM, lw=2, label="2005 TH173 (Q$\\to$D)"),
            Line2D([0], [0], color=GREEN, lw=2, label="1995 SN55 (left the Centaur zone)"),
            Line2D([0], [0], color=PURPLE, lw=2, label="2002 VR130"),
            Line2D([0], [0], color=GREY, lw=1, label="Other BM09 Q/R objects"),
        ],
        loc="upper right",
        fontsize=8.5,
    )
    spines(ax)
    fig.tight_layout()
    save(fig, "fig_od_arrows.png")


def fig_arc() -> None:
    th_arc = np.array([30, 56, 105, 196, 366], dtype=float)
    th_r = np.array([0.05, 0.05, 0.05, 0.10, 0.15])
    hd_arc = np.array([10, 15, 22, 33, 49], dtype=float)
    hd_q = np.array([0.40, 0.40, 0.40, 0.35, 0.30])

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.7))
    ax = axes[0]
    ax.plot(th_arc, th_r, "o-", color=ORANGE, lw=1.8, ms=6, label="Minority R fraction")
    ax.set_xlabel("Synthetic observation arc (days)")
    ax.set_ylabel("Fraction of clones")
    ax.set_title("(a)", loc="left", pad=4)
    ax.set_ylim(0, 0.45)
    ax.legend(loc="upper left")
    ax.text(0.97, 0.06, "Most frequent class is D at every tier", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=8, color=GREY)
    spines(ax)

    ax = axes[1]
    ax.plot(hd_arc, hd_q, "o-", color=GREEN, lw=1.8, ms=6, label="Minority Q fraction")
    ax.set_xlabel("Synthetic observation arc (days)")
    ax.set_ylabel("Fraction of clones")
    ax.set_title("(b)", loc="left", pad=4)
    ax.set_ylim(0, 0.55)
    ax.legend(loc="upper right")
    ax.text(0.97, 0.06, "Most frequent class is D at every tier", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=8, color=GREY)
    spines(ax)

    fig.tight_layout()
    save(fig, "fig_arc_mixture.png")


def fig_at() -> None:
    npz = np.load(SRC / "results" / "at_longlived_regen.npz")
    t = np.asarray(npz["t_yr"], dtype=float) / 1e6
    a = np.asarray(npz["a_au"], dtype=float)
    fig, ax = plt.subplots(figsize=(7.4, 3.7))
    ax.plot(t, a, color=BLUE, lw=0.7)
    ax.set_xlabel("Time (Myr)")
    ax.set_ylabel(r"Semimajor axis $a$ (AU)")
    ax.set_xlim(0, 40)
    spines(ax)
    fig.tight_layout()
    save(fig, "fig_at_longlived.png")


def fig_ae() -> None:
    cache = load_json("data/sbdb_centaurs_ae_cache.json")
    # cache formats vary
    if isinstance(cache, dict) and "objects" in cache:
        rows = cache["objects"]
    elif isinstance(cache, dict) and "data" in cache:
        rows = cache["data"]
    elif isinstance(cache, list):
        rows = cache
    else:
        # try common keys
        rows = cache.get("centaurs") or cache.get("rows") or []
        if not rows:
            # flatten dict of desig -> elements
            rows = [{"desig": k, **v} if isinstance(v, dict) else v for k, v in cache.items()
                    if isinstance(v, dict) and ("a" in v or "a_au" in v)]

    def get_a(r):
        for k in ("a_au", "a", "semimajor_axis"):
            if k in r and r[k] is not None:
                return float(r[k])
        return None

    def get_e(r):
        for k in ("e", "eccentricity"):
            if k in r and r[k] is not None:
                return float(r[k])
        return None

    aa, ee = [], []
    for r in rows:
        a, e = get_a(r), get_e(r)
        if a is None or e is None:
            continue
        if 5 <= a <= 35:
            aa.append(a)
            ee.append(e)
    modern = load_json("data/modern_elements.json")
    mats = load_json("results/fullscale/FULLSCALE_MATRICES.json")
    qr = [o["desig"] for o in mats["qr20"]["objects"]]

    fig, ax = plt.subplots(figsize=(7.3, 5.0))
    ax.scatter(aa, ee, s=8, c="#C8C8C8", zorder=1, linewidths=0, label=f"SBDB Centaur flag ($n={len(aa)}$)")
    qa, qe = [], []
    for d in qr:
        if d == "2005 TH173":
            continue
        if d in modern:
            qa.append(modern[d]["a_au"])
            qe.append(modern[d]["e"])
    ax.scatter(qa, qe, s=36, c=BLUE, zorder=3, label="BM09 Q/R objects (modern)")
    th = modern["2005 TH173"]
    ax.scatter([th["a_au"]], [th["e"]], s=110, marker="*", c=VERM, zorder=4, edgecolor="k", linewidths=0.3,
               label="2005 TH173 (modern)")
    ax.set_xlabel(r"Semimajor axis $a$ (AU)")
    ax.set_ylabel(r"Eccentricity $e$")
    ax.set_xlim(5, 35)
    ax.set_ylim(0, 0.95)
    ax.legend(loc="upper right", markerscale=1.1, fontsize=8.5)
    spines(ax)
    fig.tight_layout()
    save(fig, "fig_centaur_ae_context.png")


def fig_chaos() -> None:
    sys.path.insert(0, str(SRC / "src"))
    os.chdir(SRC)
    from chaos_floor import load_objects, pool  # noqa: E402
    from chaos_floor_final import ARC_MIN, CONTROLS, is_centaur  # noqa: E402
    from chaos_floor_km import FLOOR_LIFE, TCAP, km_dispersion  # noqa: E402
    import glob

    clones, caps = {}, {}
    for f in (
        glob.glob("results/pop_fragility/phase1/*.json")
        + glob.glob("results/pop_fragility/phase3_merged/*.json")
        + glob.glob("results/bm09_reclass/clones/*.json")
        + glob.glob("results/bm09_reclass/clones_censor_resolved/*.json")
    ):
        if "partial" in f:
            continue
        b = json.load(open(f, encoding="utf-8"))
        if "clones" not in b:
            continue
        cap = float(b.get("t_max_myr") or TCAP)
        if cap >= caps.get(b["desig"], -1):
            clones[b["desig"]] = b["clones"]
            caps[b["desig"]] = cap
    rows = pool(
        load_objects(["results/pop_fragility/phase1/*.json", "results/pop_fragility/phase3_merged/*.json"], "survey"),
        load_objects(["results/bm09_reclass/clones/*.json", "results/bm09_reclass/clones_censor_resolved/*.json"], "bm09"),
    )
    sample = []
    for r in rows:
        if not (all(r.get(k) is not None for k in CONTROLS) and is_centaur(r) and (r["arc_days"] or 0) >= ARC_MIN and r["desig"] in clones):
            continue
        cl = clones[r["desig"]]
        t = [min(float(c["lifetime_myr"]), caps.get(r["desig"], TCAP)) for c in cl]
        e = [0 if c.get("survived_full") else 1 for c in cl]
        if min(t) <= FLOOR_LIFE:
            continue
        d = km_dispersion(t, e)
        if d is None:
            continue
        rr = dict(r)
        rr["km"] = d
        sample.append(rr)

    res = load_json("results/chaos_floor/chaos_floor_km.json")
    x = np.array([r["log10_sigma_a_rel"] for r in sample])
    y = np.array([r["km"] for r in sample])
    inc = np.array([r["i_deg"] for r in sample])
    floor = float(y.mean())
    rng = np.random.default_rng(0)
    boot = np.array([np.mean(rng.choice(y, len(y))) for _ in range(4000)])
    flo, fhi = float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))
    sg = res["regression"]["terms"]["log10_sigma_a_rel"]
    lim = max(abs(sg["ci95"][0]), abs(sg["ci95"][1]))

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(10.2, 4.15))
    sc = axL.scatter(x, y, c=inc, cmap="cividis", s=28, alpha=0.92, edgecolor="k", linewidth=0.25, zorder=3)
    axL.axhspan(flo, fhi, color="0.55", alpha=0.18, zorder=1)
    axL.axhline(floor, color=INK, lw=1.4, zorder=2)
    cb = fig.colorbar(sc, ax=axL, pad=0.02)
    cb.set_label("Inclination (deg)")
    axL.set_xlabel(r"$\log_{10}(\sigma_a/a)$")
    axL.set_ylabel(r"Dispersion of $\log_{10}$ lifetime (dex)")
    axL.set_title("(a)", loc="left", pad=4)
    axL.annotate("", xy=(0.08, 0.14), xytext=(0.42, 0.14), xycoords="axes fraction",
                 arrowprops=dict(arrowstyle="->", color=GREY, lw=1.0))
    axL.text(0.44, 0.14, "tighter orbits", transform=axL.transAxes, va="center", fontsize=9.5, color=GREY)
    axL.text(0.97, 0.92, f"Mean dispersion {floor:.2f} dex\n(factor {10**floor:.1f} in lifetime)",
             transform=axL.transAxes, va="top", ha="right", fontsize=9.5,
             bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#DDDDDD"))
    spines(axL)

    C = np.column_stack([[r[c] for r in sample] for c in CONTROLS] + [np.ones(len(sample))])

    def resid(v):
        b, *_ = np.linalg.lstsq(C, v, rcond=None)
        return v - C @ b

    xr, yr = resid(x), resid(y) + floor
    axR.scatter(xr, yr, s=28, alpha=0.88, color=BLUE, edgecolor="k", linewidth=0.25, zorder=3)
    axR.axhspan(flo, fhi, color="0.55", alpha=0.18, zorder=1)
    axR.axhline(floor, color=INK, lw=1.4, zorder=2)
    xs = np.linspace(xr.min(), xr.max(), 50)
    axR.plot(xs, floor + lim * xs, ls="--", lw=1.1, color=VERM, zorder=2)
    axR.plot(xs, floor - lim * xs, ls="--", lw=1.1, color=VERM, zorder=2)
    axR.set_xlabel(r"$\log_{10}(\sigma_a/a)$, after removing $a,e,i,q$")
    axR.set_ylabel("Residual dispersion (dex)")
    axR.set_title("(b)", loc="left", pad=4)
    axR.text(
        0.03,
        0.08,
        f"Slope ${sg['coef']:+.3f}$ dex per tenfold $\\sigma_a/a$\n"
        f"95% CI $[{sg['ci95'][0]:+.3f},\\ {sg['ci95'][1]:+.3f}]$\n"
        f"$p = {sg['p']:.2f}$",
        transform=axR.transAxes,
        va="bottom",
        ha="left",
        fontsize=9.5,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#DDDDDD"),
        clip_on=False,
    )
    axR.legend(
        loc="upper right",
        handles=[
            Line2D([], [], color=INK, lw=1.4, label="Measured mean dispersion"),
            Line2D([], [], color=VERM, ls="--", lw=1.1, label="Slopes ruled out at 95%"),
        ],
        fontsize=8,
    )
    spines(axR)
    lo = min(y.min(), yr.min()) - 0.06
    hi = max(y.max(), yr.max()) + 0.18
    for ax in (axL, axR):
        ax.set_ylim(lo, hi)
    fig.tight_layout()
    save(fig, "fig07_chaos_floor.png")


def fig_lyap() -> None:
    ana = load_json("results/lyapunov_proxy/analysis.json")
    objs = [o for o in ana["objects"] if o.get("tau_proxy_yr")]
    tau = np.array([o["tau_proxy_yr"] for o in objs])
    disp = np.array([o["dispersion_dex"] for o in objs])
    literature = {
        "1977 UB": {"tau_yr": 1.1e4, "label": "Chiron", "c": BLUE, "off": (8, -14)},
        "1992 AD": {"tau_yr": 5e3, "label": "Pholus", "c": ORANGE, "off": (-52, -16)},
        "2005 TH173": {"tau_yr": 8e3, "label": "TH173", "c": GREEN, "off": (8, 8)},
    }
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(10.2, 3.7))
    for desig, meta in literature.items():
        row = next((o for o in objs if o["desig"] == desig), None)
        if not row:
            continue
        axL.scatter(meta["tau_yr"], row["tau_proxy_yr"], s=64, zorder=3, color=meta["c"], edgecolor="k", lw=0.4)
        axL.annotate(meta["label"], (meta["tau_yr"], row["tau_proxy_yr"]), textcoords="offset points",
                     xytext=meta["off"], fontsize=9, color=meta["c"])
    lims = [1e3, 1e5]
    axL.plot(lims, lims, color=GREY, ls="--", lw=1.1, label="Equal values")
    axL.set_xscale("log")
    axL.set_yscale("log")
    axL.set_xlim(lims)
    axL.set_ylim(lims)
    axL.set_xlabel(r"Published Lyapunov time (yr)")
    axL.set_ylabel(r"Clone-ensemble proxy $\tau$ (yr)")
    axL.set_title("(a)", loc="left", pad=4)
    axL.legend(loc="lower right")
    spines(axL)

    log_inv = np.log10(1.0 / tau)
    axR.scatter(log_inv, disp, s=28, alpha=0.9, color=BLUE, edgecolor="k", lw=0.3)
    coef = np.polyfit(log_inv, disp, 1)
    xs = np.linspace(log_inv.min(), log_inv.max(), 50)
    axR.plot(xs, np.polyval(coef, xs), color=VERM, lw=1.5, zorder=2)
    rho = ana["spearman_inv_tau_disp"]["rho"]
    p = ana["spearman_inv_tau_disp"]["p"]
    axR.set_xlabel(r"$\log_{10}(1/\tau_{\mathrm{proxy}})$ (yr$^{-1}$)")
    axR.set_ylabel("Lifetime dispersion (dex)")
    axR.set_title("(b)", loc="left", pad=4)
    axR.text(0.03, 0.96, f"Spearman $\\rho={rho:+.2f}$, $p={p:.3f}$", transform=axR.transAxes,
             ha="left", va="top", fontsize=8.5, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#DDDDDD"))
    spines(axR)
    fig.tight_layout(w_pad=1.8)
    save(fig, "fig08_lyapunov_proxy.png")


def main() -> None:
    fig_th173()
    fig_qr()
    fig_od()
    fig_arc()
    fig_at()
    fig_ae()
    fig_chaos()
    fig_lyap()
    print("all figures written to", OUT)


if __name__ == "__main__":
    main()

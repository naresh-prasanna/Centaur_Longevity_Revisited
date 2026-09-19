"""Generate paper figures from pilot + Exp C JSON (matplotlib)."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "paper" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

# Match Times-like Icarus body: serif + larger tick/axis fonts for print.
plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["DejaVu Serif", "STIXGeneral", "Times New Roman", "Times", "serif"],
        "mathtext.fontset": "stix",
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.titlesize": 12,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 10,
        "figure.titlesize": 12,
    }
)


def load():
    table = {o["desig"]: o for o in json.loads((ROOT / "data" / "bailey_table2.json").read_text())["objects"]}
    modern = json.loads((ROOT / "data" / "modern_elements.json").read_text())
    qr = json.loads((ROOT / "results" / "batch_summary_qr_long40.json").read_text())
    sd = json.loads((ROOT / "results" / "batch_summary_short_d40.json").read_text())
    return table, modern, qr, sd


def fig_orbit_revisions(table, modern, qr, sd):
    fig, ax = plt.subplots(figsize=(6.5, 5))
    for r in qr["all"] + sd["all"]:
        d = r["desig"]
        if d not in modern:
            continue
        b, m = table[d], modern[d]
        color = {"Q": "C3", "R": "C1", "D": "C0"}[r["bailey"]]
        ax.annotate(
            "",
            xy=(m["a_au"], m["e"]),
            xytext=(b["a_au"], b["e"]),
            arrowprops=dict(arrowstyle="->", color=color, lw=1.0, alpha=0.7),
        )
        ax.scatter([b["a_au"]], [b["e"]], c=color, s=20, zorder=3)
        ax.scatter([m["a_au"]], [m["e"]], c=color, s=36, marker="*", zorder=4)
    ax.set_xlabel("Semimajor axis a (AU)")
    ax.set_ylabel("Eccentricity e")
    ax.set_title("Bailey 2009 → modern SBDB (arrows)")
    ax.legend(
        handles=[
            plt.Line2D([0], [0], color="C3", marker="o", label="Bailey Q"),
            plt.Line2D([0], [0], color="C1", marker="o", label="Bailey R"),
            plt.Line2D([0], [0], color="C0", marker="o", label="Bailey D (short-arc)"),
        ],
        frameon=False,
    )
    fig.tight_layout()
    fig.savefig(FIG / "fig02_orbit_revisions.png", dpi=160)
    plt.close(fig)


def fig_flip_matrix(qr, sd):
    labels = ["D", "R", "Q"]
    # rows = bailey, cols = modern modal
    M = np.zeros((3, 3), dtype=int)
    idx = {k: i for i, k in enumerate(labels)}
    for r in qr["all"] + sd["all"]:
        M[idx[r["bailey"]], idx[r["modern"]]] += 1
    fig, ax = plt.subplots(figsize=(4.5, 4))
    im = ax.imshow(M, cmap="Blues")
    ax.set_xticks(range(3), labels)
    ax.set_yticks(range(3), labels)
    ax.set_xlabel("Modern modal class")
    ax.set_ylabel("Bailey 2009 class")
    ax.set_title("Classification flip matrix (pilots)")
    for i in range(3):
        for j in range(3):
            ax.text(j, i, str(M[i, j]), ha="center", va="center", color="black")
    fig.colorbar(im, ax=ax, fraction=0.046)
    fig.tight_layout()
    fig.savefig(FIG / "fig04_flip_matrix.png", dpi=160)
    plt.close(fig)


def fig_lifetimes(qr, sd):
    def lives_from_summary(summary):
        # need per-object batch files
        out = []
        for r in summary["all"]:
            slug = r["desig"].replace(" ", "_")
            p = ROOT / "results" / "batch" / f"{slug}.json"
            if not p.exists():
                continue
            d = json.loads(p.read_text())
            out.extend([c["lifetime_myr"] for c in d["clones"]])
        return out

    qr_l = lives_from_summary(qr)
    sd_l = lives_from_summary(sd)
    fig, ax = plt.subplots(figsize=(6, 4))
    bins = np.linspace(0, 40, 21)
    ax.hist(qr_l, bins=bins, alpha=0.65, label=f"Bailey Q/R clones (n={len(qr_l)})", color="C1")
    ax.hist(sd_l, bins=bins, alpha=0.55, label=f"Short-arc D clones (n={len(sd_l)})", color="C0")
    ax.axvline(22, color="k", ls="--", lw=1, label="Long-lived cut (22 Myr)")
    ax.set_xlabel("Lifetime (Myr)")
    ax.set_ylabel("Count")
    ax.set_title("Clone lifetimes (40 Myr pilots)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "fig05_lifetimes.png", dpi=160)
    plt.close(fig)


def fig_survival(qr, sd):
    def surv_curve(summary):
        lives = []
        for r in summary["all"]:
            slug = r["desig"].replace(" ", "_")
            p = ROOT / "results" / "batch" / f"{slug}.json"
            if not p.exists():
                continue
            d = json.loads(p.read_text())
            lives.extend([c["lifetime_myr"] for c in d["clones"]])
        lives = np.asarray(lives, dtype=float)
        ts = np.linspace(0, 40, 81)
        frac = np.array([(lives >= t).mean() if len(lives) else 0 for t in ts])
        return ts, frac

    fig, ax = plt.subplots(figsize=(6, 4))
    t, f = surv_curve(qr)
    ax.plot(t, f, label="Bailey Q/R", color="C1", lw=2)
    t, f = surv_curve(sd)
    ax.plot(t, f, label="Short-arc D", color="C0", lw=2)
    ax.axvline(22, color="k", ls="--", lw=1)
    ax.set_xlabel("Time (Myr)")
    ax.set_ylabel("Surviving fraction")
    ax.set_ylim(0, 1.05)
    ax.set_title("Survival curves (clone ensembles)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "fig06_survival.png", dpi=160)
    plt.close(fig)


def fig_th173_composite():
    """RNAAS single figure: BM09→modern (a,e) + lifetime hist + survival."""
    prod = json.loads((ROOT / "results" / "th173_production.json").read_text())
    bailey = prod["bailey"]
    modern = prod["modern_nominal"]
    lives = np.asarray([c["lifetime_myr"] for c in prod["clones"]], dtype=float)
    t_max = float(prod["t_max_myr"])

    fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.5))
    fig.subplots_adjust(wspace=0.38)

    ax = axes[0]
    ax.annotate(
        "",
        xy=(modern["a_au"], modern["e"]),
        xytext=(bailey["a_au"], bailey["e"]),
        arrowprops=dict(arrowstyle="->", color="C3", lw=1.8),
    )
    ax.scatter([bailey["a_au"]], [bailey["e"]], c="C3", s=55, zorder=3, label="BM09 (Q)")
    ax.scatter(
        [modern["a_au"]],
        [modern["e"]],
        c="C0",
        s=70,
        marker="*",
        zorder=4,
        label="modern SBDB",
    )
    ax.set_xlabel(r"$a$ (AU)")
    ax.set_ylabel(r"$e$")
    ax.set_title("(a) Orbit revision", pad=8)
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    ax.set_xlim(14.5, 21.5)
    ax.set_ylim(-0.02, 0.38)

    ax = axes[1]
    bins = np.linspace(0, t_max, 21)
    ax.hist(lives, bins=bins, color="C0", alpha=0.85, edgecolor="white", lw=0.4)
    med = float(np.median(lives))
    ax.axvline(med, color="0.35", ls="--", lw=1.2, label=f"median {med:.1f} Myr")
    ax.set_xlabel("Lifetime (Myr)")
    ax.set_ylabel("Clones")
    ax.set_title("(b) 100-clone lifetimes", pad=8)
    ax.set_xlim(0, t_max)
    ax.text(
        0.97, 0.97,
        "Q cut = 22 Myr\n(off scale)",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        color="0.25",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.85", alpha=0.92),
    )
    ax.legend(frameon=False, fontsize=8, loc="upper left")

    ax = axes[2]
    ts = np.linspace(0, t_max, 101)
    frac = np.array([(lives >= t).mean() for t in ts])
    ax.plot(ts, frac, color="C0", lw=2)
    ax.axhline(0, color="0.7", lw=0.5)
    ax.set_xlabel("Time (Myr)")
    ax.set_ylabel("Surviving fraction")
    ax.set_title("(c) Survival", pad=8)
    ax.set_ylim(0, 1.05)
    ax.set_xlim(0, t_max)
    ax.text(
        0.98,
        0.92,
        f"modal D {prod['counts']['D']}/100\n"
        f"t_max={t_max:.0f} Myr (rejects Q)",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9,
    )

    fig.tight_layout(pad=1.0)
    out = FIG / "fig_th173_composite.png"
    fig.savefig(out, dpi=200, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    return out


def fig_qr_clone_fractions():
    """Stacked D/R/Q bars for 10 BM09 Q/R objects (full-scale 20 clones, v2)."""
    mats = json.loads((ROOT / "results" / "fullscale" / "FULLSCALE_MATRICES.json").read_text())
    objs = mats["qr20"]["objects"]
    labels = [o["desig"] for o in objs]
    d_counts = [int(o["counts_v2"].get("D", 0)) for o in objs]
    r_counts = [int(o["counts_v2"].get("R", 0)) for o in objs]
    q_counts = [int(o["counts_v2"].get("Q", 0)) for o in objs]
    bailey = [o["bailey"] for o in objs]

    x = np.arange(len(objs))
    fig, ax = plt.subplots(figsize=(9.5, 4.2))
    ax.bar(x, d_counts, color="#4C78A8", label="D", width=0.72)
    ax.bar(x, r_counts, bottom=d_counts, color="#F58518", label="R", width=0.72)
    bottom_rq = [d + r for d, r in zip(d_counts, r_counts)]
    ax.bar(x, q_counts, bottom=bottom_rq, color="#54A24B", label="Q", width=0.72)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"{lab}\n(BM09 {b})" for lab, b in zip(labels, bailey)],
        fontsize=9,
        rotation=0,
    )
    ax.set_ylabel("Clone count (n=20)")
    ax.set_ylim(0, 22)
    ax.set_title("BM09 Q/R sample: classifier v2 clone fractions (40 Myr)")
    ax.axhline(10, color="0.75", ls=":", lw=0.8)
    ax.legend(frameon=False, ncol=3, loc="upper right")
    ax.text(
        0.01,
        0.98,
        "Modal class = D for 10/10; aggregate D=155, R=43, Q=2",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        color="0.25",
    )
    fig.tight_layout()
    out = FIG / "fig_qr_clone_fractions.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out


def fig_od_arrows():
    """BM09 → modern (a,e) arrows for Q/R set + SN55 (and VR130)."""
    table = {o["desig"]: o for o in json.loads((ROOT / "data" / "bailey_table2.json").read_text())["objects"]}
    modern = json.loads((ROOT / "data" / "modern_elements.json").read_text())
    mats = json.loads((ROOT / "results" / "fullscale" / "FULLSCALE_MATRICES.json").read_text())
    highlight = {
        "2005 TH173": "C3",
        "1995 SN55": "C2",
        "2002 VR130": "C4",
    }
    desigs = [o["desig"] for o in mats["qr20"]["objects"]] + ["1995 SN55", "2002 VR130"]

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    for d in desigs:
        if d not in modern or d not in table:
            continue
        b, m = table[d], modern[d]
        color = highlight.get(d, "0.55")
        lw = 2.0 if d in highlight else 0.9
        alpha = 0.95 if d in highlight else 0.55
        ax.annotate(
            "",
            xy=(m["a_au"], m["e"]),
            xytext=(b["a_au"], b["e"]),
            arrowprops=dict(arrowstyle="->", color=color, lw=lw, alpha=alpha),
        )
        ax.scatter([b["a_au"]], [b["e"]], c=color, s=28, zorder=3, alpha=alpha)
        ax.scatter([m["a_au"]], [m["e"]], c=color, s=48, marker="*", zorder=4, alpha=alpha)
        if d in highlight:
            ax.annotate(
                d.replace("20", ""),
                xy=(m["a_au"], m["e"]),
                xytext=(6, 6),
                textcoords="offset points",
                fontsize=10,
                color=color,
            )
    ax.set_xlabel(r"Semimajor axis $a$ (AU)")
    ax.set_ylabel(r"Eccentricity $e$")
    ax.set_title(r"Orbit revision: BM09 Table 2 $\rightarrow$ modern SBDB")
    ax.legend(
        handles=[
            plt.Line2D([0], [0], color="C3", lw=2, label="2005 TH173 (Q→D)"),
            plt.Line2D([0], [0], color="C2", lw=2, label="1995 SN55 (D→Q)"),
            plt.Line2D([0], [0], color="C4", lw=2, label="2002 VR130 (D→R)"),
            plt.Line2D([0], [0], color="0.55", lw=1, label="BM09 Q/R others"),
        ],
        frameon=False,
        fontsize=10,
        loc="upper right",
    )
    fig.tight_layout()
    out = FIG / "fig_od_arrows.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out


def fig_arc_mixture():
    """Option A: TH173 minority-R / frac_long vs arc + HD12 minority-Q vs arc."""
    # Locked numbers from STRETCH_E_VERDICT.md / ARC_CLASS_PRODUCTION.md
    th_arc = np.array([30, 56, 105, 196, 366], dtype=float)
    th_r = np.array([0.05, 0.05, 0.05, 0.10, 0.15])
    th_fl = np.array([0.05, 0.05, 0.05, 0.10, 0.15])
    hd_arc = np.array([10, 15, 22, 33, 49], dtype=float)
    hd_q = np.array([0.40, 0.40, 0.40, 0.35, 0.30])
    hd_fl = np.array([0.40, 0.40, 0.40, 0.35, 0.30])

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8), sharey=False)

    ax = axes[0]
    ax.plot(th_arc, th_r, "o-", color="#F58518", lw=2, label=r"minority R frac")
    ax.plot(th_arc, th_fl, "s--", color="#4C78A8", lw=1.5, label=r"frac$\geq$22 Myr")
    ax.set_xlabel("Synthetic arc (d)")
    ax.set_ylabel("Fraction of clones")
    ax.set_title(r"(a) 2005 TH173 — rising mixture ($\rho=+0.89$)")
    ax.set_ylim(0, 0.45)
    ax.legend(frameon=False, fontsize=10)
    ax.text(
        0.98,
        0.05,
        "modal D at every tier",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        color="0.35",
    )

    ax = axes[1]
    ax.plot(hd_arc, hd_q, "o-", color="#54A24B", lw=2, label=r"minority Q frac")
    ax.plot(hd_arc, hd_fl, "s--", color="#4C78A8", lw=1.5, label=r"frac$\geq$22 Myr")
    ax.set_xlabel("Synthetic arc (d)")
    ax.set_ylabel("Fraction of clones")
    ax.set_title(r"(b) 1999 HD12 — falling Q minority ($\rho=-0.89$)")
    ax.set_ylim(0, 0.55)
    ax.legend(frameon=False, fontsize=10)
    ax.text(
        0.98,
        0.05,
        "modal D at every tier",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        color="0.35",
    )

    fig.suptitle(
        r"Option A arc$\rightarrow$class-mixture (n=20, $T=25$ Myr, v2; not a modal-flip curve)",
        fontsize=12,
        y=1.02,
    )
    fig.tight_layout()
    out = FIG / "fig_arc_mixture.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out


def main():
    # Legacy pilot figures (best-effort; skip if batch summaries missing)
    try:
        table, modern, qr, sd = load()
        fig_orbit_revisions(table, modern, qr, sd)
        fig_flip_matrix(qr, sd)
        fig_lifetimes(qr, sd)
        fig_survival(qr, sd)
    except Exception as exc:  # noqa: BLE001
        print(f"skip legacy pilot figures: {exc}")

    th173 = fig_th173_composite()
    qr_frac = fig_qr_clone_fractions()
    od = fig_od_arrows()
    arc = fig_arc_mixture()
    print(f"wrote figures to {FIG}")
    print(f"TH173 composite: {th173}")
    print(f"Q/R fractions:   {qr_frac}")
    print(f"OD arrows:       {od}")
    print(f"Arc mixture:     {arc}")


if __name__ == "__main__":
    main()

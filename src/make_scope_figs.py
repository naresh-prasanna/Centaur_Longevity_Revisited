"""Build scope-expansion figures (a(t), longfrac, SBDB a-e context)."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "paper" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

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
        "legend.fontsize": 9,
    }
)


def fig_centaur_ae_context() -> None:
    rows = json.loads((ROOT / "data" / "sbdb_centaurs_ae_cache.json").read_text())["objects"]
    modern = json.loads((ROOT / "data" / "modern_elements.json").read_text())
    qr = [
        "1995 DW2",
        "1998 QM107",
        "1998 TF35",
        "2000 FZ53",
        "2003 QP112",
        "2003 UW292",
        "2005 RL43",
        "2005 RO43",
        "2005 TH173",
        "2006 SX368",
    ]
    fig, ax = plt.subplots(figsize=(6.8, 5.0))
    ax.scatter(
        [o["a"] for o in rows],
        [o["e"] for o in rows],
        s=8,
        c="0.75",
        alpha=0.55,
        label=f"SBDB Centaurs (n={len(rows)})",
        zorder=1,
    )
    qa = [modern[d]["a_au"] for d in qr if d in modern]
    qe = [modern[d]["e"] for d in qr if d in modern]
    ax.scatter(
        qa,
        qe,
        s=42,
        c="#1f4e79",
        marker="o",
        edgecolors="white",
        linewidths=0.4,
        label="BM09 Q/R sample (modern)",
        zorder=3,
    )
    th = modern["2005 TH173"]
    ax.scatter(
        [th["a_au"]],
        [th["e"]],
        s=110,
        c="#c45c26",
        marker="*",
        edgecolors="k",
        linewidths=0.4,
        label="2005 TH173 (modern)",
        zorder=4,
    )
    ax.set_xlabel(r"Semimajor axis $a$ (AU)")
    ax.set_ylabel(r"Eccentricity $e$")
    ax.set_xlim(5, 35)
    ax.set_ylim(0, 0.95)
    ax.set_title(r"Modern TH173 and BM09 Q/R sample in SBDB Centaur $(a,e)$")
    ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    outp = FIG / "fig_centaur_ae_context.png"
    fig.savefig(outp, dpi=200)
    plt.close(fig)
    print("wrote", outp)


def fig_at_from_npz() -> None:
    p = ROOT / "results" / "at_longlived_regen.npz"
    if not p.exists():
        print("skip a(t): missing", p)
        return
    z = np.load(p, allow_pickle=True)
    desig = str(z["desig"])
    ci = int(z["clone"])
    cls = str(z["class_v2"])
    t = z["t_yr"] / 1e6
    a = z["a_au"]
    fig, ax = plt.subplots(figsize=(6.8, 3.6))
    ax.plot(t, a, color="#1f4e79", lw=0.65)
    ax.set_xlabel("Time (Myr)")
    ax.set_ylabel(r"Semimajor axis $a$ (AU)")
    ax.set_xlim(0, min(40.0, float(t[-1]) + 0.5))
    ax.set_title(f"{desig} clone {ci}: $a(t)$ over 40 Myr (v2 class {cls})")
    fig.tight_layout()
    outp = FIG / "fig_at_longlived.png"
    fig.savefig(outp, dpi=200)
    plt.close(fig)
    print("wrote", outp)


def fig_qr_longfrac() -> None:
    desigs = []
    fracs = []
    for p in sorted((ROOT / "results" / "fullscale" / "qr20").glob("*.json")):
        d = json.loads(p.read_text())
        desigs.append(d["desig"])
        n = len(d["clones"])
        fracs.append(sum(1 for c in d["clones"] if c.get("lifetime_myr", 0) >= 22) / n)
    labels = [d.split()[-1] for d in desigs]
    fig, ax = plt.subplots(figsize=(6.8, 3.8))
    colors = ["#c45c26" if f >= 0.4 else "#4c6a8a" for f in fracs]
    ax.bar(range(len(fracs)), fracs, color=colors, edgecolor="0.2", linewidth=0.4)
    ax.axhline(0.5, color="0.5", ls="--", lw=0.8)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel(r"Fraction with lifetime $\geq$ 22 Myr")
    ax.set_ylim(0, 1.05)
    ax.set_title("BM09 Q/R sample: long-lived clone fraction (20 clones, 40 Myr, v2)")
    fig.tight_layout()
    outp = FIG / "fig_qr_longfrac.png"
    fig.savefig(outp, dpi=200)
    plt.close(fig)
    print("wrote", outp, list(zip(desigs, fracs)))


if __name__ == "__main__":
    fig_at_from_npz()
    fig_qr_longfrac()
    fig_centaur_ae_context()

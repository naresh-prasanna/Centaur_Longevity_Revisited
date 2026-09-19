"""Figure: ensemble Lyapunov proxy vs literature and floor correlation."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "lyapunov_proxy"
FIG = ROOT / "paper" / "figures"

LITERATURE = {
    "1977 UB": {"tau_yr": 1.1e4, "label": "Chiron"},
    "1992 AD": {"tau_yr": 5e3, "label": "Pholus"},
    "2005 TH173": {"tau_yr": 8e3, "label": "TH173"},
}

COLORS = {"Chiron": "#2b6cb0", "Pholus": "#c05621", "TH173": "#2f855a"}


def _clean_spines(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main() -> None:
    ana = json.loads((OUT / "analysis.json").read_text(encoding="utf-8"))
    objs = [o for o in ana["objects"] if o.get("tau_proxy_yr")]
    tau = np.array([o["tau_proxy_yr"] for o in objs])
    disp = np.array([o["dispersion_dex"] for o in objs])

    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["DejaVu Serif", "STIXGeneral", "Times New Roman", "Times"],
        "mathtext.fontset": "stix",
        "font.size": 10,
    })

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(10.2, 3.6))

    # --- (a) validation ---
    for desig, meta in LITERATURE.items():
        row = next((o for o in objs if o["desig"] == desig), None)
        if not row:
            continue
        c = COLORS.get(meta["label"], "C0")
        axL.scatter(meta["tau_yr"], row["tau_proxy_yr"], s=70, zorder=3, color=c, edgecolor="k", linewidth=0.4)
        axL.annotate(
            meta["label"],
            (meta["tau_yr"], row["tau_proxy_yr"]),
            textcoords="offset points",
            xytext=(7, 5),
            fontsize=9,
            color=c,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.75", alpha=0.95),
        )
    lims = [1e3, 1e5]
    axL.plot(lims, lims, color="0.45", ls="--", lw=1.2, label="1:1")
    axL.set_xscale("log")
    axL.set_yscale("log")
    axL.set_xlim(lims)
    axL.set_ylim(lims)
    axL.set_xlabel("literature $\\tau_L$ (yr)")
    axL.set_ylabel("ensemble proxy $\\tau$ (yr)")
    axL.set_title("(a) validation objects ($n=3$)", loc="left", fontsize=10, pad=10)
    axL.legend(loc="lower right", frameon=False, fontsize=8)
    axL.grid(alpha=0.22, which="major")
    _clean_spines(axL)

    # --- (b) floor vs chaos rate (log axis for readable ticks) ---
    log_inv = np.log10(1.0 / tau)  # same rank order as 1/tau for Spearman
    axR.scatter(log_inv, disp, s=32, alpha=0.88, color="#2b6cb0", edgecolor="k", linewidth=0.35)
    if len(log_inv) >= 3:
        coef = np.polyfit(log_inv, disp, 1)
        xs = np.linspace(log_inv.min(), log_inv.max(), 50)
        axR.plot(xs, np.polyval(coef, xs), color="#b03030", lw=1.5, zorder=2)
    rho = ana["spearman_inv_tau_disp"]["rho"]
    p = ana["spearman_inv_tau_disp"]["p"]
    axR.set_xlabel("$\\log_{10}(1/\\tau_{\\mathrm{proxy}})$ (yr$^{-1}$)")
    axR.set_ylabel("lifetime dispersion (dex)")
    axR.set_title("(b) faster chaos $\\rightarrow$ higher floor ($n=39$)", loc="left", fontsize=10, pad=10)
    axR.xaxis.set_major_locator(MaxNLocator(nbins=5, prune="both"))
    axR.yaxis.set_major_locator(MaxNLocator(nbins=5))
    axR.grid(alpha=0.22, which="major")
    _clean_spines(axR)
    axR.text(
        0.03, 0.97,
        f"Spearman $\\rho={rho:+.2f}$, $p={p:.3f}$",
        transform=axR.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.75", alpha=0.95),
    )

    fig.subplots_adjust(left=0.08, right=0.97, bottom=0.20, top=0.88, wspace=0.38)

    FIG.mkdir(parents=True, exist_ok=True)
    path = FIG / "fig08_lyapunov_proxy.png"
    fig.savefig(path, dpi=220, bbox_inches="tight", pad_inches=0.10)
    plt.close(fig)
    print(f"wrote {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

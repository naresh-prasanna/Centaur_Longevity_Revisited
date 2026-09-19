"""Study roadmap figure for Introduction."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "paper" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "STIXGeneral", "Times New Roman"],
    "mathtext.fontset": "stix",
    "font.size": 9,
})


def box(ax, x, y, w, h, text, fc="#f7f7f7", ec="#333333", fontsize=8.5):
    p = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.02",
        linewidth=1.0,
        edgecolor=ec,
        facecolor=fc,
        transform=ax.transAxes,
    )
    ax.add_patch(p)
    ax.text(
        x + w / 2, y + h / 2, text,
        ha="center", va="center", fontsize=fontsize,
        transform=ax.transAxes, wrap=True,
    )


def arrow(ax, x0, y0, x1, y1, color="#444444"):
    a = FancyArrowPatch(
        (x0, y0), (x1, y1),
        arrowstyle="-|>",
        mutation_scale=10,
        linewidth=1.0,
        color=color,
        transform=ax.transAxes,
        shrinkA=2,
        shrinkB=2,
    )
    ax.add_patch(a)


def main() -> None:
    fig, ax = plt.subplots(figsize=(10.5, 5.2))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Title strip
    box(
        ax, 0.04, 0.90, 0.92, 0.07,
        "Two linked questions: (i) chaos-floor bound on orbit quality  "
        "(ii) BM09 Table-2 labels under a stated modern pipeline",
        fc="#e8eef5", ec="#2b6cb0", fontsize=9,
    )

    # Left column — chaos floor
    box(ax, 0.04, 0.72, 0.44, 0.14,
        "Chaos-floor track\n219 clone ensembles ($n=65$--96)  "
        "$\\rightarrow$ lifetime dispersion vs.\\ $\\sigma_a/a$",
        fc="#eef6ff", ec="#2b6cb0")
    arrow(ax, 0.26, 0.72, 0.26, 0.64)
    box(ax, 0.04, 0.54, 0.44, 0.09,
        "Flat floor once $(a,e,i,q)$ fixed; $a$, $i$ set height",
        fc="#eef6ff", ec="#2b6cb0", fontsize=8)
    arrow(ax, 0.26, 0.54, 0.26, 0.46)
    box(ax, 0.04, 0.36, 0.44, 0.09,
        "Lyapunov proxy pilot ($n=39$): faster chaos $\\rightarrow$ higher floor",
        fc="#eef6ff", ec="#2b6cb0", fontsize=8)

  # Right column — BM09
    box(ax, 0.52, 0.72, 0.44, 0.14,
        "BM09 reclassification track\n61 census $+$ 10 Q/R deep sample  "
        "(orbit $+$ integrator $+$ escape $+$ classifier change together)",
        fc="#faf0e6", ec="#c05621")
    arrow(ax, 0.74, 0.72, 0.74, 0.64)
    box(ax, 0.52, 0.54, 0.44, 0.09,
        "Experiment~C: Bailey IC vs.\\ modern IC (pipeline fixed)",
        fc="#faf0e6", ec="#c05621", fontsize=8)
    arrow(ax, 0.74, 0.54, 0.74, 0.46)
    box(ax, 0.52, 0.36, 0.44, 0.09,
        "Label change $\\neq$ ``astrometry alone'' without isolating branches",
        fc="#faf0e6", ec="#c05621", fontsize=8)

    # Bottom — shared limits
    box(
        ax, 0.04, 0.14, 0.92, 0.18,
        "Shared interpretational limits (Limitations section):\n"
        "• v2 R = lifetime ≥22 Myr, not BM09 resonance-hopping R "
        "(nine BM09 R objects excluded from like-for-like census)\n"
        "• Classifier v1→v2 after pilot stop-check; not preregistered; "
        "TH173 Q rejection invariant\n"
        "• Headline BM09 claim: labels do not generally reproduce under our "
        "operational rule—not “BM09 was wrong”",
        fc="#f5f5f5", ec="#555555", fontsize=8.2,
    )

    path = FIG / "fig_study_roadmap.png"
    fig.savefig(path, dpi=220, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    print(f"wrote {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

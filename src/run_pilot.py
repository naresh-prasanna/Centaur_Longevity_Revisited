"""Pilot / production: agreement vs arc for Centaurs (Option A)."""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from arc_core import (
    SEED,
    agreement_vs_full,
    arc_tiers_for,
    load_bailey,
    load_modern,
    run_tier,
)
from planet_ephem import DEFAULT_KERNEL_DIR, open_giant_ephem

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results" / "synthetic_arc"
PAPER = ROOT / "paper"
CACHE = ROOT / "data" / "j2000_giants_horizons_de441.json"

# Production set: leads + short-arc D + short R + long R + long D
PRODUCTION_DESIGS = [
    "2005 TH173",  # BM09 Q, arc 366 d
    "1995 DW2",  # BM09 R, arc 5576 d
    "1999 HD12",  # BM09 D, arc 49 d (≤50 d control)
    "1996 RX33",  # BM09 D, arc 13 d
    "2003 QP112",  # BM09 R, arc 54 d
    "1998 QM107",  # BM09 R, long arc
    "1977 UB",  # Chiron, BM09 D, very long arc
]


def partial_path(tag: str, desig: str) -> Path:
    return RES / f"arc_{tag}_partial_{desig.replace(' ', '_')}.json"


def write_progress(payload: dict) -> None:
    """Atomic-ish progress for progress_bar.py / PROGRESS.html."""
    path = RES / "progress.json"
    payload = {**payload, "ts": time.time()}
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def run_object(
    desig: str,
    n_clones: int,
    t_max_myr: float,
    n_tiers: int,
    kernel_dir: Path,
    tag: str = "prod",
    resume_partial: bool = True,
) -> dict:
    modern = load_modern(desig)
    bailey = load_bailey(desig)
    tiers = arc_tiers_for(modern, n_tiers=n_tiers)
    print(f"{desig}: arc_full={modern.get('data_arc_days')} tiers={tiers} bailey={bailey}", flush=True)
    tier_rows: list[dict] = []
    t0 = time.time()
    start_i = 0
    pp = partial_path(tag, desig)
    if resume_partial and pp.exists():
        try:
            prev = json.loads(pp.read_text(encoding="utf-8"))
            if prev.get("desig") == desig and prev.get("n_clones") == n_clones:
                tier_rows = list(prev.get("tiers") or [])
                start_i = len(tier_rows)
                print(f"  RESUME partial {pp.name}: {start_i}/{len(tiers)} tiers", flush=True)
        except Exception as exc:
            print(f"  partial load failed ({exc}); starting fresh", flush=True)
            tier_rows = []
            start_i = 0

    ephem_name = None
    with open_giant_ephem(kernel_dir, CACHE) as ephem:
        ephem_name = type(ephem).__name__
        print(f"  ephem={ephem_name}", flush=True)
        for i, arc in enumerate(tiers):
            if i < start_i:
                continue
            write_progress(
                {
                    "desig": desig,
                    "tier_index": i,
                    "n_tiers": len(tiers),
                    "arc_days": float(arc),
                    "status": "running",
                    "n_clones": n_clones,
                }
            )
            print(f"  arc={arc:.1f} d ...", flush=True)
            tr = run_tier(modern, arc, n_clones, t_max_myr * 1e6, ephem, seed=SEED)
            tier_rows.append(tr)
            # Persist after every tier so Cursor abort / crash does not redo work
            pp.write_text(
                json.dumps(
                    {
                        "desig": desig,
                        "bailey_class": bailey,
                        "n_clones": n_clones,
                        "t_max_myr": t_max_myr,
                        "tiers": tier_rows,
                    },
                    indent=2,
                ),
                encoding="utf-8",
            )
            write_progress(
                {
                    "desig": desig,
                    "tier_index": i + 1,
                    "n_tiers": len(tiers),
                    "arc_days": float(arc),
                    "status": "tier_done",
                    "modal": tr["modal"],
                    "n_clones": n_clones,
                }
            )
            print(
                f"    modal={tr['modal']} counts={tr['counts']} med_life={tr['median_life_myr']:.2f} Myr",
                flush=True,
            )
    curve = agreement_vs_full(tier_rows)
    if pp.exists():
        pp.unlink(missing_ok=True)
    return {
        "desig": desig,
        "bailey_class": bailey,
        "data_arc_days": modern.get("data_arc_days"),
        "n_clones": n_clones,
        "t_max_myr": t_max_myr,
        "seed": SEED,
        "classifier": "v2",
        "ephem": ephem_name,
        "elapsed_s": time.time() - t0,
        "curve": curve,
        "tiers": tier_rows,
    }


def write_results_md(summaries: list[dict], path: Path, title: str) -> None:
    lines = [
        f"# {title}",
        "",
        r"Classifier **v2**. Arc model: \(\sigma\propto\sqrt{A_{\mathrm{full}}/\mathrm{arc}}\).",
        "",
    ]
    for s in summaries:
        lines.append(
            f"## {s['desig']} (BM09 {s.get('bailey_class')}, full arc {s.get('data_arc_days')} d, "
            f"n={s.get('n_clones')}, T={s.get('t_max_myr')} Myr)"
        )
        lines.append("")
        lines.append("| arc (d) | σ scale | modal | agree w/ full | med life (Myr) | frac≥22 Myr |")
        lines.append("|--------:|--------:|:-----:|-------------:|---------------:|------------:|")
        for c in s["curve"]:
            lines.append(
                f"| {c['arc_days']:.0f} | {c['sigma_scale']:.2f} | {c['modal']} | "
                f"{c['clone_agreement']:.2f} | {c['median_life_myr']:.2f} | {c['frac_long']:.2f} |"
            )
        lines.append("")
    lines.append("See `CONTRACT.md`. Option I production remains paused.")
    path.write_text("\n".join(lines), encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--desigs", nargs="+", default=None)
    ap.add_argument("--production", action="store_true", help="20 clones, production object set")
    ap.add_argument("--n-clones", type=int, default=None)
    ap.add_argument("--t-max-myr", type=float, default=25.0)
    ap.add_argument("--n-tiers", type=int, default=5)
    ap.add_argument("--kernel-dir", default=str(DEFAULT_KERNEL_DIR))
    ap.add_argument("--resume", action="store_true", help="Skip objects with existing result JSON")
    args = ap.parse_args()

    if args.production:
        desigs = args.desigs or PRODUCTION_DESIGS
        n_clones = args.n_clones or 20
        tag = "prod"
        md_name = "ARC_CLASS_PRODUCTION.md"
        summary_name = "arc_prod_summary.json"
    else:
        desigs = args.desigs or ["2005 TH173", "1995 DW2"]
        n_clones = args.n_clones or 5
        tag = "pilot"
        md_name = "ARC_CLASS_PILOT.md"
        summary_name = "arc_pilot_summary.json"

    RES.mkdir(parents=True, exist_ok=True)
    PAPER.mkdir(parents=True, exist_ok=True)
    kdir = Path(args.kernel_dir)

    summaries = []
    for des in desigs:
        out = RES / f"arc_{tag}_{des.replace(' ', '_')}.json"
        if args.resume and out.exists():
            print(f"SKIP {des} (exists {out.name})", flush=True)
            summaries.append(json.loads(out.read_text(encoding="utf-8")))
            continue
        try:
            s = run_object(
                des,
                n_clones,
                args.t_max_myr,
                args.n_tiers,
                kdir,
                tag=tag,
                resume_partial=args.resume,
            )
        except Exception as exc:
            print(f"FAIL {des}: {exc}", flush=True)
            continue
        out.write_text(json.dumps(s, indent=2), encoding="utf-8")
        print(f"wrote {out}", flush=True)
        summaries.append(s)
        # Incremental summary so a crash still leaves a usable curve set
        (RES / summary_name).write_text(json.dumps(summaries, indent=2), encoding="utf-8")
        write_results_md(
            summaries,
            PAPER / md_name,
            "Option A — Arc→class stability production"
            if args.production
            else "Option A — Arc→class stability pilot",
        )

    if summaries:
        (RES / summary_name).write_text(json.dumps(summaries, indent=2), encoding="utf-8")
        write_results_md(
            summaries,
            PAPER / md_name,
            "Option A — Arc→class stability production"
            if args.production
            else "Option A — Arc→class stability pilot",
        )
        print("wrote", PAPER / md_name, flush=True)


if __name__ == "__main__":
    main()

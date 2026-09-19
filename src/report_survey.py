"""Population fragility survey: Spearman table and stratum breakdown.

Predictors are orbit-determination quality metrics; outcomes are measures of how
unstable an object's classification is across its clone ensemble.
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from pop_survey import spearman_table, write_summary_csv  # noqa: E402
from run_pop_fragility import OUT, load_campaign_rows  # noqa: E402

rows = load_campaign_rows("phase1")
print(f"Phase 1 objects: {len(rows)}\n")

sp = spearman_table(rows)
print(f"{'predictor':<18}{'outcome':<24}{'n':>5}{'rho':>9}{'p':>12}   sig")
print("-" * 72)
for r in sp:
    rho = "  --" if r.get("rho") is None else f"{r['rho']:+.3f}"
    p = r.get("p")
    ps = "  --" if p is None else (f"{p:.2e}" if p < 1e-3 else f"{p:.4f}")
    star = "***" if p is not None and p < 0.001 else "**" if p is not None and p < 0.01 \
        else "*" if p is not None and p < 0.05 else ""
    print(f"{r['predictor']:<18}{r['outcome']:<24}{r['n']:>5}{rho:>9}{ps:>12}   {star}")

(OUT / "spearman_phase1.json").write_text(json.dumps(sp, indent=2), encoding="utf-8")
write_summary_csv(rows, OUT / "summary_phase1.csv")

print("\nModal class distribution:")
dist = defaultdict(int)
for r in rows:
    dist[r["modal_class"]] += 1
for k in ("D", "R", "Q"):
    print(f"   {k}: {dist.get(k, 0)}")

print("\nBy arc-length stratum:")
by = defaultdict(list)
for r in rows:
    by[r.get("arc_bin")].append(r)
print(f"{'stratum':<12}{'n':>5}{'mean R/Q frac':>15}{'mean censor':>13}{'median life':>13}")
for k in sorted(by, key=lambda x: str(x)):
    g = by[k]
    def mean(key):
        v = [x[key] for x in g if x.get(key) is not None]
        return sum(v) / len(v) if v else float("nan")
    print(f"{str(k):<12}{len(g):>5}{mean('minority_rq_frac'):>15.3f}"
          f"{mean('censoring_rate'):>13.3f}{mean('median_lifetime_myr'):>13.2f}")

print(f"\nwrote {(OUT / 'spearman_phase1.json').as_posix()}")
print(f"wrote {(OUT / 'summary_phase1.csv').as_posix()}")

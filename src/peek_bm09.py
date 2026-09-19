"""Quick progress peek at the BM09 reclassification while it runs."""
import glob
import json
from collections import Counter

rows = []
for f in sorted(glob.glob("results/bm09_reclass/clones/*.json")):
    if "partial" in f:
        continue
    b = json.load(open(f, encoding="utf-8"))
    m = b["modern_nominal"]
    lf = sorted(c["lifetime_myr"] for c in b["clones"])
    rows.append({
        "desig": b["desig"],
        "bailey": m["bailey_class"],
        "modern": b["modal_class"],
        "frac_modal": b["frac_modal"],
        "med_life": lf[len(lf) // 2],
        "min_life": lf[0],
        "max_life": lf[-1],
        "arc": m["arc_2007_raw"],
        "quality": m["orbit_2007_quality"],
    })

print(f"{'desig':<12} {'BM09':<5} {'now':<5} {'frac':<5} {'med':>7} {'min':>7} {'max':>7}  {'arc2007':<9} flag")
for r in rows:
    flag = "CHANGED" if r["bailey"] != r["modern"] else ""
    print(f"{r['desig']:<12} {r['bailey']:<5} {r['modern']:<5} {r['frac_modal']:<5.1f} "
          f"{r['med_life']:>7.2f} {r['min_life']:>7.2f} {r['max_life']:>7.2f}  {r['arc']:<9} {flag}")

n = len(rows)
ch = sum(1 for r in rows if r["bailey"] != r["modern"])
print(f"\nn={n}  changed={ch} ({ch/n:.0%})" if n else "no rows yet")
for q in ("multi_opposition", "short_arc"):
    sub = [r for r in rows if r["quality"] == q]
    if sub:
        c = sum(1 for r in sub if r["bailey"] != r["modern"])
        print(f"  {q:<18} n={len(sub):<3} changed={c} ({c/len(sub):.0%})")
print("BM09 label mix so far:", dict(Counter(r["bailey"] for r in rows)))
print("modern label mix:     ", dict(Counter(r["modern"] for r in rows)))

"""Does the R (resonance-hopping) class reproduce under v1 vs v2 of the classifier?

Every clone record stores class_v1 and class_v2, so the BM09 sample can be re-scored
under either rule without re-integrating anything.
"""
import glob
import json
from collections import Counter

rows = []
for f in sorted(glob.glob("results/bm09_reclass/clones/*.json")):
    if "partial" in f:
        continue
    b = json.load(open(f, encoding="utf-8"))
    m = b["modern_nominal"]
    clones = b["clones"]

    def modal(key):
        c = Counter(cl.get(key) for cl in clones if cl.get(key))
        return c.most_common(1)[0][0] if c else None

    rows.append({
        "desig": b["desig"],
        "bailey": m["bailey_class"],
        "v1": modal("class_v1"),
        "v2": modal("class_v2"),
        "n_v1_R": sum(1 for c in clones if c.get("class_v1") == "R"),
        "n_v2_R": sum(1 for c in clones if c.get("class_v2") == "R"),
        "quality": m.get("orbit_2007_quality"),
    })

for rule in ("v1", "v2"):
    agree = sum(1 for r in rows if r["bailey"] == r[rule])
    print(f"{rule}: overall agreement with BM09 = {agree}/{len(rows)} ({agree/len(rows):.0%})")
    for bc in ("D", "R", "Q"):
        sub = [r for r in rows if r["bailey"] == bc]
        if sub:
            a = sum(1 for r in sub if r[rule] == bc)
            print(f"   BM09 {bc} (n={len(sub):>2}): reproduced {a:>2} ({a/len(sub):>4.0%})"
                  f"   -> {dict(Counter(r[rule] for r in sub))}")
    print()

print("BM09 R-class objects in detail:")
print(f"{'desig':<12} {'BM09':<5} {'v1':<4} {'v2':<4} {'#clones R (v1)':>15} {'#clones R (v2)':>15}")
for r in rows:
    if r["bailey"] == "R":
        print(f"{r['desig']:<12} {r['bailey']:<5} {str(r['v1']):<4} {str(r['v2']):<4} "
              f"{r['n_v1_R']:>15} {r['n_v2_R']:>15}")

tot_v1 = sum(r["n_v1_R"] for r in rows)
tot_v2 = sum(r["n_v2_R"] for r in rows)
n_clones = sum(10 for _ in rows)
print(f"\nR-labelled clones across the whole sample: v1={tot_v1}/{n_clones}, v2={tot_v2}/{n_clones}")

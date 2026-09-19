"""Feasibility probes for a predictive at-risk model, using only archived clones.

Three questions:
  1. How far is the population from the 22 Myr cut? If nothing is near it, no
     OD-quality covariate can move a class, and the null is structural.
  2. How much does the class split depend on where the cut is placed?
  3. How much within-object dispersion is there to model at all?
"""
import glob
import json
import math
import statistics as st

CUTS = [2, 5, 10, 15, 22, 30, 40]


def load(pat):
    objs = []
    for f in glob.glob(pat):
        if "partial" in f:
            continue
        b = json.load(open(f, encoding="utf-8"))
        if "clones" not in b:
            continue
        objs.append(b)
    return objs


for label, pat in (("SURVEY", "results/pop_fragility/phase1/*.json"),
                   ("BM09", "results/bm09_reclass/clones/*.json")):
    objs = load(pat)
    lives = [c["lifetime_myr"] for b in objs for c in b["clones"]]
    med_per_obj = [st.median([c["lifetime_myr"] for c in b["clones"]]) for b in objs]
    print(f"=== {label}: {len(objs)} objects, {len(lives)} clones ===")
    print(f"  clone lifetime: median {st.median(lives):.2f}  p90 {sorted(lives)[int(.9*len(lives))]:.2f}  max {max(lives):.2f}")
    print(f"  object median lifetime: max {max(med_per_obj):.2f} Myr  (cut is 22 Myr)")
    near = sum(1 for m in med_per_obj if 11 <= m <= 44)
    print(f"  objects with median within a factor 2 of the cut: {near}/{len(objs)}")

    print("  class split vs where the cut is placed (flat rule aside, D = life < cut):")
    for cut in CUTS:
        frac_long = sum(1 for x in lives if x >= cut) / len(lives)
        n_modal_long = sum(1 for b in objs
                           if sum(1 for c in b["clones"] if c["lifetime_myr"] >= cut) > len(b["clones"]) / 2)
        print(f"     cut {cut:>2} Myr -> {frac_long*100:5.1f}% of clones long-lived, "
              f"{n_modal_long:>3}/{len(objs)} objects modal non-D")

    # Within-object dispersion: the quantity OD quality should actually drive.
    disp = []
    for b in objs:
        v = [max(c["lifetime_myr"], 1e-4) for c in b["clones"]]
        lg = [math.log10(x) for x in v]
        if len(lg) > 1:
            disp.append(st.pstdev(lg))
    print(f"  within-object sd of log10(lifetime): median {st.median(disp):.2f} dex, "
          f"range {min(disp):.2f}-{max(disp):.2f} dex")
    print(f"  objects with dispersion > 0.5 dex: {sum(1 for d in disp if d > 0.5)}/{len(disp)}\n")

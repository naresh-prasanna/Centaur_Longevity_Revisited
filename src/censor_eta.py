"""How much censoring work remains, and which sampled objects produced no file."""
import glob
import json
from pathlib import Path

for label, pat in (("BM09", "results/bm09_reclass/clones/*.json"),
                   ("SURVEY", "results/pop_fragility/phase1/*.json")):
    n_obj = 0
    n_cens = 0
    per = []
    for f in glob.glob(pat):
        if "partial" in f:
            continue
        b = json.load(open(f, encoding="utf-8"))
        c = sum(1 for x in b.get("clones", []) if x.get("survived_full"))
        if c:
            n_obj += 1
            n_cens += c
            per.append((b["desig"], c))
    per.sort(key=lambda t: -t[1])
    print(f"{label}: {n_obj} objects carry {n_cens} censored clones")
    print(f"   worst: {per[:6]}")
    # A censored clone survived Tmax, so extending to 40 Myr costs the full window.
    print(f"   est. cost at ~215 s/clone on 8 workers: {n_cens*215/8/60:.0f} min\n")

sample = json.load(open("data/pop_fragility_sample.json", encoding="utf-8"))
want = {o["desig"] for o in sample["objects"]}
have = set()
for f in glob.glob("results/pop_fragility/phase1/*.json"):
    if "partial" in f:
        continue
    have.add(json.load(open(f, encoding="utf-8"))["desig"])
missing = sorted(want - have)
print(f"Phase 1 sampled {len(want)}, wrote {len(have)}, missing {len(missing)}")
for d in missing:
    p = Path("results/pop_fragility/phase1") / (d.replace(" ", "_") + ".partial.json")
    print(f"   {d}  partial_exists={p.exists()}")

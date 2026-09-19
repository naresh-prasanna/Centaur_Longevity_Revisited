"""Are all BM09 Table 2 objects still Centaurs under modern orbits?

A dynamical classification only means something if the object is still in the
population being classified. BM09 selected on 2007 orbits, some from arcs of days.
This applies the standard Centaur zone (perihelion beyond Jupiter, semimajor axis
interior to Neptune) to both the 2007 and the modern elements.
"""
import glob
import json

A_NEPTUNE = 30.07
Q_JUPITER = 5.20

rows = []
for f in sorted(glob.glob("results/bm09_reclass/clones/*.json")):
    if "partial" in f:
        continue
    b = json.load(open(f, encoding="utf-8"))
    m = b["modern_nominal"]
    old = m.get("bailey_elements") or {}
    a_new, q_new = m.get("a_au"), m.get("q_au")
    a_old, q_old = old.get("a_au"), old.get("q_au")
    if None in (a_new, q_new, a_old, q_old):
        continue

    def is_cen(a, q):
        return (q > Q_JUPITER) and (a < A_NEPTUNE)

    rows.append({
        "desig": b["desig"],
        "bailey_class": m.get("bailey_class"),
        "modern_class": b["modal_class"],
        "a_old": a_old, "q_old": q_old, "cen_old": is_cen(a_old, q_old),
        "a_new": a_new, "q_new": q_new, "cen_new": is_cen(a_new, q_new),
        "arc_2007": m.get("arc_2007_raw"),
        "quality": m.get("orbit_2007_quality"),
        "arc_now": m.get("data_arc_days"),
        "full_name": m.get("full_name"),
    })

lost = [r for r in rows if r["cen_old"] and not r["cen_new"]]
print(f"BM09 objects checked: {len(rows)}")
print(f"Centaur under 2007 elements: {sum(1 for r in rows if r['cen_old'])}")
print(f"Centaur under modern elements: {sum(1 for r in rows if r['cen_new'])}")
print(f"No longer Centaurs: {len(lost)}\n")

if lost:
    print(f"{'desig':<12} {'name':<28} {'a_2007':>7} {'q_2007':>7} {'a_now':>7} {'q_now':>7}  {'arc2007':<8} {'arc_now':>8}")
    for r in lost:
        print(f"{r['desig']:<12} {(r['full_name'] or '')[:28]:<28} "
              f"{r['a_old']:>7.2f} {r['q_old']:>7.2f} {r['a_new']:>7.2f} {r['q_new']:>7.2f}  "
              f"{r['arc_2007']:<8} {r['arc_now']:>8.0f}")
    print()
    by_q = {}
    for r in lost:
        by_q.setdefault(r["quality"], []).append(r["desig"])
    for k, v in by_q.items():
        print(f"  {k}: {len(v)} -> {v}")

json.dump({"n_checked": len(rows), "n_lost": len(lost), "lost": lost, "rows": rows},
          open("results/bm09_reclass/_MEMBERSHIP.json", "w", encoding="utf-8"), indent=2)
print("\nwrote results/bm09_reclass/_MEMBERSHIP.json")

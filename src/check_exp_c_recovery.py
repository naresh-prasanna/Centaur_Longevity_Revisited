"""Audit Experiment C label recovery against BM09's published classes.

Answers two questions the manuscript conflated:
  (a) how many objects have a branch label equal to their BM09 published class,
      under the integration-era rule (v1) and under v2;
  (b) how many objects change label between the Bailey-IC and modern-IC branches.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
LONG = 22.0


def classify_v2(rec: dict) -> str:
    life = float(rec.get("lifetime_myr") or 0.0)
    flat = bool(rec.get("flat"))
    if life >= LONG and flat:
        return "Q"
    if life >= LONG:
        return "R"
    return "D"


def main() -> None:
    src = RES / "exp_c_summary_qr40.json"
    d = json.loads(src.read_text())
    rows = d["results"]

    print(f"source        : {src.name}")
    print(f"label         : {d.get('label')}")
    print(f"T_max (Myr)   : {d.get('t_max_myr')}   angles: {d.get('angles')}   n: {d.get('n')}")
    print()
    hdr = (f"{'object':14} {'BM09':4} | {'v1 bail':7} {'v1 mod':6} {'v1 flip':7} "
           f"| {'v2 bail':7} {'v2 mod':6} {'v2 flip':7} | {'life bail':9} {'life mod':9}")
    print(hdr)
    print("-" * len(hdr))

    n = 0
    rec_v1_b = rec_v1_m = rec_v2_b = rec_v2_m = 0
    rec_v1_either = rec_v2_either = 0
    flip_v1 = flip_v2 = 0

    for r in rows:
        n += 1
        bm = r["bailey_class"]
        b, m = r["bailey_ic"], r["modern_ic"]
        v1b, v1m = b["class"], m["class"]
        v2b, v2m = classify_v2(b), classify_v2(m)

        rec_v1_b += v1b == bm
        rec_v1_m += v1m == bm
        rec_v2_b += v2b == bm
        rec_v2_m += v2m == bm
        rec_v1_either += (v1b == bm) or (v1m == bm)
        rec_v2_either += (v2b == bm) or (v2m == bm)
        flip_v1 += v1b != v1m
        flip_v2 += v2b != v2m

        print(f"{r['desig']:14} {bm:4} | {v1b:7} {v1m:6} {str(v1b != v1m):7} "
              f"| {v2b:7} {v2m:6} {str(v2b != v2m):7} "
              f"| {b['lifetime_myr']:9.2f} {m['lifetime_myr']:9.2f}")

    print()
    print(f"n objects                                  : {n}")
    print(f"RECOVERY of BM09 class, Bailey-IC branch   : v1 {rec_v1_b}/{n}   v2 {rec_v2_b}/{n}")
    print(f"RECOVERY of BM09 class, modern-IC branch   : v1 {rec_v1_m}/{n}   v2 {rec_v2_m}/{n}")
    print(f"RECOVERY on at least one branch            : v1 {rec_v1_either}/{n}   v2 {rec_v2_either}/{n}")
    print(f"BRANCH FLIP (Bailey-IC vs modern-IC label)  : v1 {flip_v1}/{n}   v2 {flip_v2}/{n}")
    print()
    print("objects recovering BM09 class under v2:")
    for r in rows:
        bm = r["bailey_class"]
        v2b, v2m = classify_v2(r["bailey_ic"]), classify_v2(r["modern_ic"])
        if v2b == bm or v2m == bm:
            which = ("Bailey-IC" if v2b == bm else "") + ("/modern-IC" if v2m == bm else "")
            print(f"  {r['desig']:14} BM09 {bm}  v2 {v2b}/{v2m}   match on {which.strip('/')}")
    print("objects changing label between branches under v2:")
    for r in rows:
        v2b, v2m = classify_v2(r["bailey_ic"]), classify_v2(r["modern_ic"])
        if v2b != v2m:
            print(f"  {r['desig']:14} {v2b} -> {v2m}")


if __name__ == "__main__":
    main()

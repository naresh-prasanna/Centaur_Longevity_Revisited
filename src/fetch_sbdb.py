"""Fetch modern JPL SBDB elements for Bailey Table 2 objects."""
from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / "data" / "bailey_table2.json"
OUT = ROOT / "data" / "modern_elements.json"
API = "https://ssd-api.jpl.nasa.gov/sbdb.api"


def fetch_one(desig: str) -> dict:
    q = urllib.parse.urlencode({"sstr": desig, "full-prec": "true", "phys-par": "false"})
    url = f"{API}?{q}"
    with urllib.request.urlopen(url, timeout=60) as resp:
        raw = json.loads(resp.read().decode("utf-8"))
    if "object" not in raw or "orbit" not in raw:
        raise RuntimeError(f"SBDB miss for {desig}: keys={list(raw.keys())}")
    obj = raw["object"]
    orb = raw["orbit"]
    elems = {e["name"]: e for e in orb["elements"]}
    # SBDB sigma often under orbit.elements[*].sigma or orbit.covariance
    sigma = {}
    for key, name in (
        ("a", "a_au"),
        ("e", "e"),
        ("i", "i_deg"),
        ("om", "om_deg"),
        ("w", "w_deg"),
        ("ma", "ma_deg"),
    ):
        el = elems.get(key, {})
        if "sigma" in el and el["sigma"] not in (None, ""):
            sigma[name] = float(el["sigma"])

    def val(key: str) -> float:
        return float(elems[key]["value"])

    a = val("a")
    e = val("e")
    cc = orb.get("condition_code")
    try:
        cc_out = int(cc) if cc is not None and str(cc).strip() != "" else None
    except (TypeError, ValueError):
        cc_out = cc

    return {
        "query": desig,
        "full_name": obj.get("fullname") or desig,
        "spkid": str(obj.get("spkid", "")),
        "pdes": obj.get("pdes", desig),
        "epoch_jd": float(orb["epoch"]),
        "a_au": a,
        "e": e,
        "i_deg": val("i"),
        "om_deg": val("om"),
        "w_deg": val("w"),
        "ma_deg": val("ma"),
        "q_au": a * (1.0 - e),
        "condition_code": cc_out,
        "data_arc_days": _num(orb.get("data_arc")),
        "n_obs": _num(orb.get("n_obs_used")),
        "sigma": sigma,
        "pe_used": orb.get("pe_used"),
        "source": "JPL SBDB",
    }


def _num(x):
    try:
        return float(x) if x is not None else None
    except (TypeError, ValueError):
        return None


def main():
    table = json.loads(TABLE.read_text(encoding="utf-8"))
    modern = {}
    if OUT.exists():
        modern = json.loads(OUT.read_text(encoding="utf-8"))

    errors = {}
    for row in table["objects"]:
        des = row["desig"]
        if des in modern and "a_au" in modern[des]:
            print(f"skip cached {des}", flush=True)
            continue
        try:
            modern[des] = fetch_one(des)
            print(
                f"ok {des}: a={modern[des]['a_au']:.3f} e={modern[des]['e']:.3f} "
                f"arc={modern[des].get('data_arc_days')} U={modern[des].get('condition_code')}",
                flush=True,
            )
            OUT.write_text(json.dumps(modern, indent=2), encoding="utf-8")
            time.sleep(0.35)
        except Exception as exc:  # noqa: BLE001
            errors[des] = str(exc)
            print(f"FAIL {des}: {exc}", flush=True)
            time.sleep(0.5)

    OUT.write_text(json.dumps(modern, indent=2), encoding="utf-8")
    summary = {
        "n_ok": sum(1 for v in modern.values() if "a_au" in v),
        "n_fail": len(errors),
        "errors": errors,
    }
    (ROOT / "data" / "modern_elements_fetch_log.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()

"""Emit every number the manuscript quotes for the new census, in one JSON.

Nothing here computes physics; it only aggregates archived per-object results so
the text and the data files cannot drift apart.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from analyze_bm09 import analyse, load_rows  # noqa: E402
from pop_survey import spearman_table  # noqa: E402
from run_pop_fragility import load_campaign_rows  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
A_NEPTUNE, Q_JUPITER = 30.07, 5.20


def is_centaur(a, q):
    return a is not None and q is not None and q > Q_JUPITER and a < A_NEPTUNE


bm = load_rows("clones")
res = analyse(bm)
membership = json.loads((ROOT / "results" / "bm09_reclass" / "_MEMBERSHIP.json")
                        .read_text(encoding="utf-8"))
lost = {r["desig"] for r in membership["lost"]}

# Like-for-like set: still a Centaur today, and on the comparable D/Q axis.
comparable = [r for r in bm
              if r["desig"] not in lost
              and r["bailey_class"] in ("D", "Q")
              and r["modern_class"] in ("D", "Q")]
comp_ctrl = [r for r in comparable if r["orbit_2007_quality"] == "multi_opposition"]
comp_test = [r for r in comparable if r["orbit_2007_quality"] == "short_arc"]

survey = load_campaign_rows("phase1")
sp = spearman_table(survey)
strata = {}
for r in survey:
    strata.setdefault(r.get("arc_bin"), []).append(r)


def mean(g, k):
    v = [x[k] for x in g if x.get(k) is not None]
    return sum(v) / len(v) if v else None


out = {
    "census": {
        "n_table2": 63,
        "n_excluded_neptune_trojan": 1,
        "n_excluded_no_modern_solution": 1,
        "n_integrated": res["n_objects"],
        "n_clones_per_object": 10,
        "t_max_myr": 25.0,
        "censor_extension_myr": 40.0,
        "seed": 87,
        "confusion": res["confusion_bailey_to_modern"],
        "per_class": res["per_bailey_class"],
        "n_changed_all_classes": res["n_changed"],
        "control_multi_opposition": res["control_multi_opposition"],
        "test_short_arc": res["test_short_arc"],
        "fisher_all_classes": res["fisher_short_arc_vs_correlation"]
        if "fisher_short_arc_vs_correlation" in res else res["fisher_short_arc_vs_control"],
    },
    "membership": {
        "n_checked": membership["n_checked"],
        "n_centaur_2007": sum(1 for r in membership["rows"] if r["cen_old"]),
        "n_centaur_modern": sum(1 for r in membership["rows"] if r["cen_new"]),
        "n_lost": membership["n_lost"],
        "lost": [{k: r[k] for k in ("desig", "full_name", "a_old", "q_old", "a_new",
                                    "q_new", "arc_2007", "arc_now", "quality")}
                 for r in membership["lost"]],
        "criterion": f"q > {Q_JUPITER} AU and a < {A_NEPTUNE} AU",
    },
    "like_for_like_DQ": {
        "rationale": ("Excludes the nine BM09 R objects (BM09's R is resonance hopping; "
                      "the v2 R requires survival to 22 Myr, so the two are not the same "
                      "quantity) and the two objects that are no longer Centaurs."),
        "n": len(comparable),
        "n_changed": sum(r["changed"] for r in comparable),
        "changed": [{"desig": r["desig"], "from": r["bailey_class"], "to": r["modern_class"],
                     "arc_2007": r["arc_2007_raw"], "frac_modal": r["frac_modal"],
                     "median_lifetime_myr": r["median_lifetime_myr"]}
                    for r in comparable if r["changed"]],
        "multi_opposition": {"n": len(comp_ctrl),
                             "n_changed": sum(r["changed"] for r in comp_ctrl)},
        "short_arc": {"n": len(comp_test),
                      "n_changed": sum(r["changed"] for r in comp_test)},
    },
    "R_class_non_equivalence": {
        "n_bm09_R": res["per_bailey_class"].get("R", {}).get("n"),
        "n_reproduced": res["per_bailey_class"].get("R", {}).get("n_retained"),
        "cause": ("assign_class_v2 returns R only for clones with lifetime >= 22 Myr; "
                  "BM09's R denotes resonance hopping irrespective of lifetime."),
    },
    "population_survey": {
        "n_objects": len(survey),
        "n_clones_per_object": 10,
        "t_max_myr": 25.0,
        "modal_class_distribution": {
            k: sum(1 for r in survey if r["modal_class"] == k) for k in ("D", "R", "Q")},
        "spearman": sp,
        "min_p": min([r["p"] for r in sp if r.get("p") is not None], default=None),
        "result": "null; no predictor reaches p < 0.14",
        "cause": ("every object is modal D because the v2 rule needs 22 Myr survival and "
                  "median clone lifetimes are 0.19-2.67 Myr at Tmax = 25 Myr, so the "
                  "outcome variable has no variance"),
        "by_arc_stratum": {
            str(k): {"n": len(g),
                     "mean_minority_rq_frac": mean(g, "minority_rq_frac"),
                     "mean_censoring_rate": mean(g, "censoring_rate"),
                     "mean_median_lifetime_myr": mean(g, "median_lifetime_myr")}
            for k, g in sorted(strata.items(), key=lambda t: str(t[0]))},
    },
}

cf = json.loads((ROOT / "results" / "chaos_floor" / "chaos_floor_final.json")
                .read_text(encoding="utf-8"))
ck = json.loads((ROOT / "results" / "chaos_floor" / "chaos_floor_km.json")
                .read_text(encoding="utf-8"))
out["chaos_floor"] = {
    "n_integrated": cf["n_integrated"],
    "n_complete_case": cf["n_sample"],
    "n_failing_membership_integrated": cf["n_failing_membership"],
    "floor_dex": cf["floor_dex"],
    "floor_ci95": cf["floor_ci95"],
    "floor_median_dex": cf["floor_median_dex"],
    "lifetime_factor": cf["lifetime_factor"],
    "sigma_range_log10": cf["sigma_range_log10"],
    "sigma_coef_dex_per_decade": cf["sigma_coef_dex_per_decade"],
    "sigma_coef_ci95": cf["sigma_coef_ci95"],
    "sigma_coef_p": cf["sigma_coef_p"],
    "excluded_effect_dex_per_decade": cf["excluded_effect_dex_per_decade"],
    "max_reducible_dispersion_dex": cf["max_reducible_dispersion_dex"],
    "groups": cf["groups"],
    "catalogue": cf["catalogue"],
    "km_n": ck["n_sample"],
    "km_disp_mean_dex": ck["km_disp_mean"],
    "km_disp_median_dex": ck["km_disp_median"],
    "km_sigma_coef": ck["regression"]["terms"]["log10_sigma_a_rel"]["coef"],
    "km_sigma_p": ck["regression"]["terms"]["log10_sigma_a_rel"]["p"],
    "km_i_deg_p": ck["regression"]["terms"]["i_deg"]["p"],
    "km_interaction_p": ck["interaction"]["terms"]["sigma_x_life"]["p"],
}

lp = json.loads((ROOT / "results" / "lyapunov_proxy" / "analysis.json")
                .read_text(encoding="utf-8"))
out["lyapunov_proxy"] = {
    "n_objects": lp["n_objects"],
    "n_with_tau": lp["n_with_tau"],
    "spearman_inv_tau_disp_rho": lp["spearman_inv_tau_disp"]["rho"],
    "spearman_inv_tau_disp_p": lp["spearman_inv_tau_disp"]["p"],
    "validation": lp["validation"],
}

path = ROOT / "results" / "PAPER_NUMBERS.json"
path.write_text(json.dumps(out, indent=2), encoding="utf-8")
print(json.dumps(out, indent=2))
print(f"\nwrote {path.relative_to(ROOT).as_posix()}")

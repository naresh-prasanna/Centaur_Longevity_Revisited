# Population fragility survey — run report

Honest record of what ran. Numbers are from archived JSON, not invented.

## Machine
{
  "cpu_count": 16,
  "workers": 8,
  "ram_gb": 7.34,
  "ram_avail_gb": 0.39
}

## Pilot
{
  "phase": "pilot",
  "objects": [
    "2012 TT72",
    "2017 QJ173",
    "2020 KH7",
    "2013 UU17",
    "2001 PT13"
  ],
  "n_objects": 5,
  "n_clones": 5,
  "t_max_myr": 25.0,
  "workers": 8,
  "wall_s": 0.007387876510620117,
  "n_clone_jobs": 25,
  "mean_s_per_clone": 10.021741161346435,
  "median_s_per_clone": 0.753000020980835,
  "mean_s_per_clone_myr": 3.8077349156996507
}

## Projection after pilot
{
  "phase1": {
    "n_obj": 175,
    "n_clones": 10,
    "t_max_myr": 25.0,
    "workers": 8,
    "s_per_clone_assumed": 10.021741161346435,
    "n_jobs": 1750,
    "wall_h": 0.6089599664012592,
    "wall_s": 2192.255879044533,
    "note": "as planned"
  },
  "phase2_h": 0.06263588225841522,
  "phase3_h_if_20pct_censored": 0.1437145520706971,
  "phases_1_3_h": 0.8153104007303715,
  "pilot_wall_h": 2.052187919616699e-06
}

## What actually completed
- Phase 1 objects: 172 (requested 175)
- Phase 2 covariance objects: 13
- Phase 3 objects with 40 Myr extensions: 36
- Phase 1 wall_h: 1.819568225675159
- Phase 2 wall_h: 0.014654081397586399
- Phase 3 wall_h: 1.9690331157048544
- Deviations: ['covariance fetch/run failures: {\'2012 TT72\': "\'NoneType\' object is not subscriptable", \'2014 OF331\': "\'NoneType\' object is not subscriptable", \'2015 KS161\': "cannot pickle \'_io.BufferedReader\' object", \'2015 PK79\': "\'elements\'", \'2015 TP10\': "\'NoneType\' object is not subscriptable", \'2016 LR88\': "\'elements\'", \'2018 PQ18\': "\'NoneType\' object is not subscriptable", \'2014 EB132\': "\'elements\'", \'2014 EV344\': "\'elements\'", \'2016 UP173\': "\'NoneType\' object is not subscriptable", \'2013 AE174\': "\'elements\'", \'2014 EU131\': "\'elements\'", \'2015 HW414\': "\'elements\'", \'2015 OA196\': "\'elements\'", \'2026 EZ79\': "cannot pickle \'_io.BufferedReader\' object", \'2014 EJ323\': "\'elements\'", \'2014 EM342\': "cannot pickle \'_io.BufferedReader\' object", \'2014 EP322\': "\'elements\'"}']

## Spearman (Phase 1)
- log10_arc_days vs minority_rq_frac: rho=0.0816224065905856 p=0.2871297331632323 n=172
- log10_arc_days vs disagree_vs_nominal: rho=0.04356267291611452 p=0.5704271066290859 n=172
- log10_arc_days vs censoring_rate: rho=0.07264380790498275 p=0.3436288665915382 n=172
- log10_arc_days vs frac_ge22: rho=0.0816224065905856 p=0.2871297331632323 n=172
- condition_code vs minority_rq_frac: rho=0.05602997074793203 p=0.4653647059743084 n=172
- condition_code vs disagree_vs_nominal: rho=-0.009286908356029359 p=0.9037614164622496 n=172
- condition_code vs censoring_rate: rho=0.06869155679351324 p=0.3705891264335212 n=172
- condition_code vs frac_ge22: rho=0.05602997074793203 p=0.4653647059743084 n=172
- n_oppositions vs minority_rq_frac: rho=0.11218183327437352 p=0.14287952990038727 n=172
- n_oppositions vs disagree_vs_nominal: rho=0.0754428490512747 p=0.32530655234646905 n=172
- n_oppositions vs censoring_rate: rho=0.10336100134871971 p=0.17724092514711307 n=172
- n_oppositions vs frac_ge22: rho=0.11218183327437352 p=0.14287952990038727 n=172

## Diagonal vs covariance (short-arc)
- n=13; mean (cov − diag) R/Q fraction = 0.0000
- Positive mean: covariance yields more R/Q clones than diagonal.
- Negative mean: diagonal overestimates the R/Q minority.

## Figures
- `C:/Users/Hp/Downloads/Horizon87-v1.0-horizon87-v1-clean/option_j_th173_class/results/pop_fragility/figures/fig_fragility_vs_arc.png`
- `C:/Users/Hp/Downloads/Horizon87-v1.0-horizon87-v1-clean/option_j_th173_class/results/pop_fragility/figures/fig_diag_vs_cov.png`
- `C:/Users/Hp/Downloads/Horizon87-v1.0-horizon87-v1-clean/option_j_th173_class/results/pop_fragility/figures/fig_modal_agreement.png`

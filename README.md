# Centaur Longevity Revisited — Reproducibility Archive

[![DOI](https://zenodo.org/badge/1376715397.svg)](https://doi.org/10.5281/zenodo.22840759)

REBOUND-based N-body integration code and archived production data
reclassifying the Bailey & Malhotra (2009) Centaur sample on modern JPL SBDB
orbits.

## Paper

Naresh Prasanna and Chrisphin Karthick. "Centaur Longevity Revisited: Lifetime
Dispersion and the Bailey & Malhotra Classifications under Modern Orbital
Solutions." Submitted to *Celestial Mechanics and Dynamical Astronomy*.

## Contents

- `src/` — REBOUND integration driver, clone generator,
  classifier (production rule v2 and the v1 Hurst-gated comparison), campaign
  scheduler, statistical aggregation, and figure-generation scripts.
- `data/` — cached BM09 Table 2 and modern SBDB element files (snapshot
  2026-07-26), plus the J2000 giant-planet cache used to initialize WHFast.
- `results/` — archived campaign JSON: TH173 (10 Myr and 40 Myr
  runs), the 20-clone BM09 Q/R matrix, the 61-object BM09 Table 2 census, the
  TH173 SBDB-covariance Q-test, the orbit-update control, the 219-ensemble
  pooled dispersion analysis, the 39-object Lyapunov proxy pilot, the
  Appendix B synthetic-arc-length scan, and the Appendix C 172-object
  population survey. See Table 4 of the paper for what each campaign covers.

Campaign folders under `results/`:

| Paper campaign | Archive path |
|---|---|
| TH173, 100-clone, 10 Myr | `results/th173_production.json` |
| TH173, 100-clone, 40 Myr | `results/fullscale/th173_100x40/` |
| Q/R full-scale, 20-clone, 40 Myr | `results/fullscale/qr20/` |
| BM09 Table 2 census + 40 Myr extensions | `results/bm09_reclass/` |
| TH173 SBDB-covariance Q-test, 20-clone, 10 Myr | `results/fullscale/cov_th173_q10/` |
| Orbit-update control (single-nominal branches) | `results/exp_c/` |
| 219-ensemble dispersion analysis | `results/chaos_floor/` |
| Lyapunov proxy, 39 objects | `results/lyapunov_proxy/` |
| Appendix B synthetic-arc scan | `results/synthetic_arc/` |
| Appendix C population survey | `results/pop_fragility/phase1/` |

## Reproducing the results

All production integrations use random seed 87 and a JPL SBDB snapshot dated
2026-07-26.

Regenerate every tabulated number from the archived per-object JSON (no
re-integration):

```
python src/paper_numbers.py
python src/validate_classifier.py
```

Rebuild Figs. 1–9 from the archived JSON:

```
python src/rebuild_figures.py
```

To re-run integrations, DE441 split kernels must be available locally (they
are not in this archive). Pass the kernel directory to the campaign scripts,
for example:

```
python src/run_th173.py --n-clones 100 --t-max-myr 10 --label production --kernel-dir /path/to/de441
python src/run_fullscale.py --phase A --workers 2 --kernel-dir /path/to/de441
python src/run_bm09_reclass.py --kernel-dir /path/to/de441
python src/run_covariance_campaign.py --campaign cov_th173_q10 --kernel-dir /path/to/de441
python src/run_exp_c.py --kernel-dir /path/to/de441
python src/run_limitation_suite.py --kernel-dir /path/to/de441
```

`environment.txt` records the production Python / REBOUND / seed settings.

## Requirements

See `requirements.txt`. REBOUND is the core dependency (WHFast and Mercurius
integrators).

## Citation

If you use this code or data, please cite the paper above. A Zenodo DOI for
this archive is forthcoming and will be added here once minted.

## License

MIT — see `LICENSE`.

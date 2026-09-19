# TH173 integrator check (WHFast vs Mercurius)

**Call:** gap closed

## Setup

- Object: 2005 TH173 (modern SBDB)
- Clones: first 20 from `sample_clones_from_modern` (seed 87)
- t_max: 10.0 Myr (same as production)
- Integrators: whfast, mercurius
- Ephem: `CachedJ2000Ephem`
- Escape cuts / sampling / forces: unchanged from `classify_core.integrate_clone`

## Key numbers

- Median lifetime WHFast: **3.857 Myr**
- Median lifetime Mercurius: **5.562 Myr**
- Clones with life < 22 Myr under **both**: **20** / 20
- Clones Q under **both**: **0** / 20
- v2 class agreement: **20** / 20 (1.00)
- Spearman ρ(lifetime WHFast, Mercurius): **-0.436**
- Risk clones (WHFast→D but Mercurius long-lived/Q/R): **0**

## Interpretation

Per-clone lifetimes are **not** reproducible across integrators (\(\rho<0\)): that is expected Centaur chaos (Lyapunov times \(\sim10^4\)–\(10^5\) yr), not a bug. The closed gap is **ensemble/class-level** only — both integrators give 20/20 D and 0 Q; “0 risk clones” does **not** mean bit-for-bit trajectory agreement.

If Mercurius confirms short lives and no Q, the WHFast Q-rejection gap is closed for this referee concern. If Mercurius yields long-lived/Q survivors where WHFast gave D, the production headline is at risk.

**Verdict:** gap closed (class-level)

Elapsed: 1481.8 s

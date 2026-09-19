# 9-hour sequential compute status

**Window:** 9 hours, **one Python job at a time**, `--workers 1`.  
**Started:** 2026-08-14 evening.

## What happened earlier

Parallel launch of 6 heavy jobs caused Windows paging-file / CLR spawn failures. DE441 `.bsp` kernels are missing, so jobs use `CachedJ2000Ephem`. Partial progress: d20 finished **2/51** objects; TH173 100×40, covariance, Exp C clones, and Mercurius 100 did **not** finish.

## Queue order (fits what it can in 9h)

| Priority | Job | Expected | Notes |
|----------|-----|----------|--------|
| 1 | TH173 100×40 Myr | ~3–5 h | Highest scientific value |
| 2 | Mercurius 100×10 Myr | ~2 h | Resumes from 1/100 partial |
| 3 | Covariance TH173 100×40 | ~3–5 h | Starts if time remains |
| 4 | Covariance QP112 20×40 | ~0.5–1 h | |
| 5 | Exp C 20-clone ensembles | ~6–10 h | Likely truncated by deadline |
| 6 | D census 20×40 remainder | ~16 h full | 2/51 already done; resume if time |
| last | v1/v2 re-label | minutes | if anything finished |

Arc 10-tier and full cov_qr20 are **not** in this window (too large after the TH173 jobs).

Monitor: `results/limitation_suite/nine_hour_status.json` and `results/limitation_suite/*.log`.

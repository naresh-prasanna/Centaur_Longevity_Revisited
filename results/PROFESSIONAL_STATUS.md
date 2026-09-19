# Professional status — classifier locked

**Default rule: v2.** Synthetic battery v2_ok = `True` (pass rate 1.0).

## Rules

- **v2 (production):** life≥22 Myr & flat → Q; life≥22 → R; else D.
- **v1 (legacy/sensitivity):** long + (nonlinear | r²<0.85) → R; else long + r²≥0.85 → D. **Withdrawn for headlines** — forces long linear survivors to D.

## Q/R pilot (modern ICs)

| Metric | v1 | v2 |
|--------|---:|---:|
| n | 10 | 10 |
| modal flip vs BM09 | 100% | 70% |

| Desig | BM09 | modal v1 | modal v2 | counts v2 |
|-------|------|----------|----------|-----------|
| 1995 DW2 | R | D | D | `{'D': 18, 'R': 2}` |
| 1998 QM107 | R | D | D | `{'D': 18, 'R': 2}` |
| 1998 TF35 | R | D | D | `{'D': 3}` |
| 2000 FZ53 | R | D | R | `{'R': 2, 'D': 1}` |
| 2003 QP112 | R | D | R | `{'R': 2, 'D': 1}` |
| 2003 UW292 | R | D | D | `{'D': 3}` |
| 2005 RL43 | R | D | D | `{'D': 3}` |
| 2005 RO43 | R | D | R | `{'D': 1, 'R': 2}` |
| 2005 TH173 | Q | D | D | `{'D': 3}` |
| 2006 SX368 | R | D | D | `{'R': 1, 'D': 2}` |

## Short-arc D control

Stay-D rate: v1 95% / v2 95% (n=21).

## Experiment C (BM09 aei vs modern aei, same code)

orbit_drives_flip: v1 **0/10** · v2 **2/10**

| Desig | BM09 | v1 B/M | v2 B/M | life B/M |
|-------|------|--------|--------|----------|
| 1995 DW2 | R | D/D | D/D | 1.21/0.30 |
| 1998 QM107 | R | D/D | D/D | 10.22/0.43 |
| 1998 TF35 | R | D/D | D/D | 1.14/3.82 |
| 2000 FZ53 | R | D/D | D/R | 7.65/22.69 |
| 2003 QP112 | R | D/D | R/R | 40.00/40.00 |
| 2003 UW292 | R | D/D | D/D | 2.06/2.42 |
| 2005 RL43 | R | D/D | D/D | 7.28/6.77 |
| 2005 RO43 | R | D/D | D/D | 17.70/11.08 |
| 2005 TH173 | Q | D/D | D/D | 8.74/21.38 |
| 2006 SX368 | R | D/D | D/R | 6.49/26.04 |

## TH173 production

n=100, T_max=10.0 Myr. counts v2 = `{'D': 100}`. T_max=10 Myr cannot assign Q or R (both need ≥22 Myr). Result is D-rejection of Q only.

## Allowed claims

- TH173 is not Q under modern elements (100 clones / 10 Myr; none reach Q criteria).
- Classifier v1 produced a false 100% Q/R→D headline; do not use as discovery.
- Under v2, Q/R pilot flip rate is reduced; several long survivors recover as R.
- Exp C under v1 was all-D (instrument); under v2, QP112 is R/R and 2 objects show orbit-driven class change.
- Short-arc D control mostly stays D under both rules; SN55 is the notable D→Q exception.

## Withdrawn claims

- ~~100% BM09 Q/R → D flip rate as a dynamical discovery under the legacy classifier.~~
- ~~Experiment C all-D as evidence that orbit revision alone does not matter (v1 was degenerate for R).~~

## Still needed for a longer / RP-grade paper

1. 20–100 clones on full Q/R set under v2
2. 100 Myr integrations for objects that survive 40 Myr
3. Escape-criterion and angle-completion sensitivity
4. Full Table 2 D census under v2
5. Honest journal rewrite using these matrices (not the legacy 100% flip headline)

Machine-readable: `professional_matrices.json`, `paper_summary.json`, `classifier_validation.json`.

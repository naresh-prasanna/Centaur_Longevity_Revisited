# Classifier integrity stop-check

## Verdict

You were right to stop. The classifier is **not** literally stuck on D — R/Q fire on a handful of real clones — but it is **structurally D-biased** for the objects that matter (long-lived BM09-R survivors). Experiment C’s all-D table is therefore **part instrument failure**, not a clean dynamical finding.

**Paper claims that treat “100% Q/R → D” as discovery are on hold.**

---

## Direct answer: has R or Q ever appeared?

**Yes, but almost never — and almost never on the BM09 Q/R targets.**

Global walk (`CLASSIFIER_AUDIT.json`): D:286 / R:3 / Q:5 across 294 labeled fields. Strip metadata false-Qs (BM09 label dicts with `life=null`) → **5 real clone hits**:

| Object | Clone | Class | Life | Context |
|--------|------:|-------|-----:|---------|
| 1995 SN55 | 0 | Q | 40 | **D-control** (not a BM09 R/Q target) |
| 1995 SN55 | 1 | Q | 40 | D-control |
| 1995 SN55 | 2 | R | 40 | D-control |
| 2003 QC112 | 1 | R | 40 | short-arc D control |
| 2005 RO43 | 1 | R | 40 | only BM09-R target with any R under v1 |

Among **lifetime ≥ 22 Myr** records: **25/30 still labeled D** under current rules.

So: not a null instrument in the absolute sense — SN55 proves Q and R can fire — but for Experiment C / Q/R pilot, the R channel is effectively dead. That **is** the “null instrument for R/Q under this setup” signature you flagged.

---

## Why Exp C is all-D (two mechanisms)

### 1. Early escape (legitimate D under BM09 lifetime cut)

Most Exp C objects die well before 22 Myr on **both** ICs (e.g. DW2 ~1 Myr, TF35 ~1–4 Myr, UW292 ~2 Myr, TH173 bailey ~9 Myr). Those **should** be D. That part is dynamics / pipeline lifetime, not Hurst gating.

### 2. Long survivors forced to D (instrument bug relative to BM09)

Decision tree in `classify_bailey` (v1):

```text
Q  if flat and life ≥ 22
R  if life ≥ 22 and (nonlinear or r² < 0.85)
D  if r² ≥ 0.85          ← fires on long-lived linear wanderers
D  if life < 22
```

`nonlinear` requires `SSR_lin > 1.5 SSR_quad` **and** `r² < 0.92`. Survivors with excellent linear Hurst fits never get R.

Smoking gun — **2003 QP112** (BM09 **R**), Exp C, both ICs survive 40 Myr:

| IC | life | H | r² | nonlinear | flat | **v1** |
|----|-----:|--:|----:|:---------:|:----:|:------:|
| BM09 aei | 40 | 0.41 | 0.974 | no | no | **D** |
| modern | 40 | 0.45 | 0.994 | no | no | **D** |

Modern `a_std ≈ 167 AU` — not “diffusive escape,” just excellent R² → forced D. Same label under BM09’s own (a,e,i). **That cannot be an orbit-flip result.**

BM09’s own text correlated long life ↔ R and short life ↔ D. Our rule 3 does the opposite for long-lived *linear* series.

---

## Offline fix test: classifier v2

BM09-faithful lifetime rule (Hurst nonlinearity kept as diagnostic only):

1. life ≥ 22 and flat → **Q**
2. life ≥ 22 → **R**
3. else → **D**

`src/reclass_v2.py` → `results/CLASSIFIER_V2_RECLASS.json`

### Experiment C under v1 vs v2

| Object | BM09 | v1 B/M | v2 B/M | life B/M (Myr) |
|--------|------|--------|--------|----------------|
| 1995 DW2 | R | D/D | D/D | 1.2 / 0.3 |
| 1998 QM107 | R | D/D | D/D | 10.2 / 0.4 |
| 1998 TF35 | R | D/D | D/D | 1.1 / 3.8 |
| 2000 FZ53 | R | D/D | **D/R** | 7.6 / 22.7 |
| 2003 QP112 | R | D/D | **R/R** | 40 / 40 |
| 2003 UW292 | R | D/D | D/D | 2.1 / 2.4 |
| 2005 RL43 | R | D/D | D/D | 7.3 / 6.8 |
| 2005 RO43 | R | D/D | D/D | 17.7 / 11.1 |
| 2005 TH173 | Q | D/D | D/D | 8.7 / 21.4 |
| 2006 SX368 | R | D/D | **D/R** | 6.5 / 26.0 |

| Metric | v1 | v2 |
|--------|---:|---:|
| orbit_drives_flip | **0/10** | **2/10** |
| Bailey-IC → R | 0 | 1 (QP112) |
| Modern-IC → R | 0 | 3 |

**QP112 becoming R/R under v2** is the decisive check: the all-D Exp C table was the classifier refusing R for long survivors, not “both orbit sets are D.”

### Q/R pilot modal flip rate

| | v1 | v2 |
|--|---:|---:|
| Modal flip vs BM09 | **10/10 (100%)** | **7/10 (70%)** |
| Objects recovering modal R under v2 | 0 | FZ53, QP112, RO43 |

Headline “100% flip” **does not survive** a BM09-aligned R definition. Remaining ~70% flips are mostly early escapes (still need pipeline / IC / angle scrutiny — not automatically “fragile BM09 labels”).

### What still holds

- **TH173 Q rejection** (100 clones, 10 Myr): still valid — Q requires ≥22 Myr; modern orbit is not Bailey’s nearly circular a≈15.7.
- **SN55 D→Q** in the short-arc control: still the clearest positive that the code can emit Q.
- Early-escape majority in Exp C: real under *this* integrator/escape cut — but that is a different claim than “R→D reclassification.”

---

## On page count / Research Papers

A **7-page compressed note is not RP-ready** if the core classifier has not been stress-tested. More pages only help if they are **real analyses**, not appendices:

Must-have before more prose:
1. Adopt or sensitivity-test classifier (v1 vs v2 vs BM09 nonlinear definition); publish both tables
2. Synthetic Q/R/D toys that pass/fail each rule set
3. Per-clone class fractions (not only modal), with R/Q hit rates
4. 20–100 clones on Q/R set; 100 Myr subset for survivors
5. Full Table 2 D census under the chosen rule
6. Escape-criterion and angle-completion sensitivity
7. Honest rewrite: TH173 orbit revision = strong; “100% R flip” = withdrawn pending calibration

Until (1)–(3) are done, stop writing Results text that depends on v1 R/Q labels.

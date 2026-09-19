# Full-scale campaign status

Q/R 20-clone (40 Myr): **10/10** objects · flip_v2 = **100%**

D census 5-clone (40 Myr): **51/51** · flip_v2 = **3.9%**

TH173 40 Myr: 1 · surv100: 15

## Q/R hardened

| Desig | BM09 | modal v2 | counts | frac≥22 Myr |
|-------|------|----------|--------|-------------|
| 1995 DW2 | R | D | `{'D': 18, 'R': 2}` | 10% |
| 1998 QM107 | R | D | `{'D': 18, 'R': 2}` | 10% |
| 1998 TF35 | R | D | `{'D': 16, 'R': 4}` | 20% |
| 2000 FZ53 | R | D | `{'R': 4, 'D': 16}` | 20% |
| 2003 QP112 | R | D | `{'R': 7, 'D': 11, 'Q': 2}` | 45% |
| 2003 UW292 | R | D | `{'D': 17, 'R': 3}` | 15% |
| 2005 RL43 | R | D | `{'D': 16, 'R': 4}` | 20% |
| 2005 RO43 | R | D | `{'D': 10, 'R': 10}` | 50% |
| 2005 TH173 | Q | D | `{'D': 17, 'R': 3}` | 15% |
| 2006 SX368 | R | D | `{'R': 4, 'D': 16}` | 20% |

## D census flips (if any)

- **1995 SN55** D→Q `{'Q': 4, 'R': 1}`
- **2002 VR130** D→R `{'D': 1, 'R': 4}`

## Completeness

```json
{
  "qr20_n": 10,
  "qr20_target": 10,
  "d5_n": 51,
  "d5_target": 51,
  "th173_40_n": 1,
  "surv100_n": 15
}
```

Rebuild paper when qr20≥10 and d5≥40: `python src/build_fullscale_paper.py`

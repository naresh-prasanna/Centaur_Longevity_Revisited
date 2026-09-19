"""Integrate TH173 clones; Bailey-style Hurst / lifetime classification."""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import rebound

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(ROOT / "src"))

from elements import BAILEY, MODERN, sample_clones, write_snapshot  # noqa: E402
from classify_core import classify_bailey  # noqa: E402
from planet_ephem import (  # noqa: E402
    DAYS_PER_YEAR,
    GIANT_NAMES,
    JD_J2000,
    DE441Ephem,
    gm_sun_au3_yr2,
    heliocentric_state_au_day as giant_state,
    masses_msun,
)

SEED = 87
G_AU_YR = gm_sun_au3_yr2()


def hill_radius(a: float, m: float) -> float:
    return a * (m / 3.0) ** (1.0 / 3.0)


def integrate_clone(
    elems: dict,
    ephem: DE441Ephem,
    t_max_yr: float,
    dt: float = 0.25,
    sample_yr: float = 300.0,
) -> dict:
    """Return a-series (yr), lifetime_yr, escaped flag."""
    m = masses_msun()
    state = giant_state(ephem, JD_J2000)

    sim = rebound.Simulation()
    sim.G = G_AU_YR
    sim.integrator = "whfast"
    sim.dt = dt
    sim.add(m=1.0)
    a_p = []
    for name in GIANT_NAMES:
        x, v_day = state[name]
        v = v_day * DAYS_PER_YEAR
        sim.add(m=m[name], x=x[0], y=x[1], z=x[2], vx=v[0], vy=v[1], vz=v[2])
        a_p.append(float(np.linalg.norm(x)))
    sim.add(
        m=0.0,
        a=elems["a"],
        e=elems["e"],
        inc=elems["inc"],
        Omega=elems["Omega"],
        omega=elems["omega"],
        M=elems["M"],
    )
    sim.N_active = 5
    sim.move_to_com()

    a_p = np.asarray(a_p, dtype=np.float64)
    m_p = np.array([m[n] for n in GIANT_NAMES], dtype=np.float64)
    tp_i = 5

    times: list[float] = []
    a_vals: list[float] = []
    escaped = False
    t_samp = 0.0
    try:
        while sim.t < t_max_yr:
            t_target = min(t_samp, t_max_yr)
            if t_target > sim.t:
                sim.integrate(t_target)

            primary = sim.particles[0]
            p = sim.particles[tp_i]
            orb = p.orbit(primary=primary)
            sun = np.array([primary.x, primary.y, primary.z])
            r = float(np.linalg.norm(np.array([p.x, p.y, p.z]) - sun))
            q = orb.a * (1.0 - orb.e) if np.isfinite(orb.a) else 0.0
            planet_h = (
                np.array([[sim.particles[k].x, sim.particles[k].y, sim.particles[k].z] for k in range(1, 5)])
                - sun
            )
            tp_h = np.array([p.x, p.y, p.z]) - sun
            hit = any(
                np.linalg.norm(tp_h - planet_h[k]) < hill_radius(float(a_p[k]), float(m_p[k]))
                for k in range(4)
            )
            if (not np.isfinite(orb.a)) or r > 1.0e4 or q < 2.5 or hit:
                escaped = True
                break

            times.append(float(sim.t))
            a_vals.append(float(orb.a))
            t_samp += sample_yr
            if t_samp > t_max_yr and sim.t >= t_max_yr - 0.5 * dt:
                break
    except Exception:
        escaped = True

    lifetime = float(sim.t)
    return {
        "t_yr": np.asarray(times, dtype=np.float64),
        "a_au": np.asarray(a_vals, dtype=np.float64),
        "lifetime_yr": lifetime,
        "escaped": bool(escaped),
        "survived_full": (not escaped) and lifetime >= t_max_yr * 0.999,
    }




def run_ensemble(n_clones: int, t_max_yr: float, kernel_dir: Path, label: str):
    write_snapshot(ROOT / "data" / "th173_elements.json")
    (ROOT / "data").mkdir(exist_ok=True)
    (ROOT / "results").mkdir(exist_ok=True)

    clones = sample_clones(n_clones, seed=SEED)
    # also include exact modern nominal as clone 0 replacement
    clones[0] = {
        "a": MODERN["a_au"],
        "e": MODERN["e"],
        "inc": np.deg2rad(MODERN["i_deg"]),
        "Omega": np.deg2rad(MODERN["om_deg"]),
        "omega": np.deg2rad(MODERN["w_deg"]),
        "M": np.deg2rad(MODERN["ma_deg"]),
    }

    rows = []
    t0 = time.time()
    with DE441Ephem(kernel_dir) as ephem:
        for i, el in enumerate(clones):
            r = integrate_clone(el, ephem, t_max_yr)
            c = classify_bailey(r["a_au"], r["t_yr"], r["lifetime_yr"])
            row = {
                "clone": i,
                "a0": el["a"],
                "e0": el["e"],
                "escaped": r["escaped"],
                "survived_full": r["survived_full"],
                **c,
            }
            rows.append(row)
            print(
                f"  clone {i+1}/{n_clones} class={c['class']} "
                f"life={c['lifetime_myr']:.2f} Myr H={c['H']} r2={c['r2']}",
                flush=True,
            )

    classes = [r["class"] for r in rows]
    counts = {k: classes.count(k) for k in ("D", "R", "Q")}
    modal = max(counts, key=counts.get)
    frac_modal = counts[modal] / len(rows)
    frac_Q = counts["Q"] / len(rows)
    flip = frac_Q < (2.0 / 3.0)  # contract: ≥2/3 Q → reject flip
    out = {
        "label": label,
        "seed": SEED,
        "n_clones": n_clones,
        "t_max_myr": t_max_yr / 1.0e6,
        "bailey": BAILEY,
        "modern_nominal": {k: MODERN[k] for k in ("a_au", "e", "i_deg", "q_au", "condition_code")},
        "counts": counts,
        "modal_class": modal,
        "frac_modal": frac_modal,
        "frac_Q": frac_Q,
        "H1_flip": bool(flip),
        "call": (
            f"FLIP vs Bailey Q — modal={modal}, frac_Q={frac_Q:.2f}"
            if flip
            else f"NO FLIP — frac_Q={frac_Q:.2f} ≥ 2/3"
        ),
        "elapsed_s": time.time() - t0,
        "clones": rows,
    }
    path = ROOT / "results" / f"th173_{label}.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({k: out[k] for k in out if k != "clones"}, indent=2), flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-clones", type=int, default=20)
    ap.add_argument("--t-max-myr", type=float, default=10.0)
    ap.add_argument("--label", default="pilot")
    ap.add_argument("--kernel-dir", default=r"C:\Users\Hp\Downloads\Zenith")
    args = ap.parse_args()
    run_ensemble(args.n_clones, args.t_max_myr * 1.0e6, Path(args.kernel_dir), args.label)


if __name__ == "__main__":
    main()

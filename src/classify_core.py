"""Shared Bailey-style integrator + Hurst classifier (multi-object)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import rebound

sys.path.insert(0, str(Path(__file__).resolve().parent))

from planet_ephem import (  # noqa: E402
    DAYS_PER_YEAR,
    GIANT_NAMES,
    JD_J2000,
    CachedJ2000Ephem,
    DE441Ephem,
    gm_sun_au3_yr2,
    heliocentric_state_au_day as giant_state,
    masses_msun,
    open_giant_ephem,
)

SUPPORTED_INTEGRATORS = ("whfast", "ias15", "mercurius")

G_AU_YR = gm_sun_au3_yr2()
LONG_LIVED_MYR = 22.0


def hill_radius(a: float, m: float) -> float:
    return a * (m / 3.0) ** (1.0 / 3.0)


def elems_from_modern(modern: dict) -> dict:
    return {
        "a": float(modern["a_au"]),
        "e": float(modern["e"]),
        "inc": float(np.deg2rad(modern["i_deg"])),
        "Omega": float(np.deg2rad(modern["om_deg"])),
        "omega": float(np.deg2rad(modern["w_deg"])),
        "M": float(np.deg2rad(modern["ma_deg"])),
    }


def sample_clones_from_modern(modern: dict, n: int, seed: int = 87) -> list[dict]:
    """Diagonal-sigma clones; clone 0 = exact nominal."""
    rng = np.random.default_rng(seed)
    s = modern.get("sigma") or {}
    defaults = {
        "a_au": 0.05,
        "e": 0.01,
        "i_deg": 0.05,
        "om_deg": 0.05,
        "w_deg": 0.5,
        "ma_deg": 0.5,
    }
    # If SBDB omitted sigmas, use modest floors so n>1 still explores
    sa = float(s.get("a_au", defaults["a_au"]))
    se = float(s.get("e", defaults["e"]))
    si = float(s.get("i_deg", defaults["i_deg"]))
    so = float(s.get("om_deg", defaults["om_deg"]))
    sw = float(s.get("w_deg", defaults["w_deg"]))
    sm = float(s.get("ma_deg", defaults["ma_deg"]))

    out = [elems_from_modern(modern)]
    for _ in range(max(0, n - 1)):
        e = float(np.clip(rng.normal(modern["e"], se), 0.0, 0.95))
        a = float(max(5.0, rng.normal(modern["a_au"], sa)))
        out.append(
            {
                "a": a,
                "e": e,
                "inc": float(np.deg2rad(rng.normal(modern["i_deg"], si))),
                "Omega": float(np.deg2rad(rng.normal(modern["om_deg"], so))),
                "omega": float(np.deg2rad(rng.normal(modern["w_deg"], sw))),
                "M": float(np.deg2rad(rng.normal(modern["ma_deg"], sm))),
            }
        )
    return out


def _configure_integrator(sim: rebound.Simulation, integrator: str, dt: float) -> None:
    """Configure REBOUND 5.x integrator object (sim.integrator.<field>)."""
    name = integrator.lower().strip()
    if name not in SUPPORTED_INTEGRATORS:
        raise ValueError(f"unsupported integrator={integrator!r}; expected one of {SUPPORTED_INTEGRATORS}")
    sim.integrator = name
    sim.dt = float(dt)
    # REBOUND 5: settings live on sim.integrator (not ri_*).
    if name == "mercurius":
        # Hybrid WHFast↔IAS15; default r_crit_hill=3. Keep safe_mode for clean sampling.
        sim.integrator.r_crit_hill = 3.0
        sim.integrator.safe_mode = 1
    elif name == "ias15":
        # Bound adaptive cost on deep encounters while keeping high accuracy.
        sim.integrator.epsilon = 1e-6
        sim.integrator.min_dt = 1e-4 * float(dt)
    elif name == "whfast":
        # Match production: plain WHFast at fixed dt (safe_mode default is fine).
        pass


def integrate_clone(
    elems: dict,
    ephem: DE441Ephem | CachedJ2000Ephem,
    t_max_yr: float,
    dt: float = 0.25,
    sample_yr: float = 300.0,
    r_escape_au: float = 1.0e4,
    q_min_au: float = 2.5,
    use_hill: bool = True,
    integrator: str = "whfast",
) -> dict:
    m = masses_msun()
    state = giant_state(ephem, JD_J2000)

    sim = rebound.Simulation()
    sim.G = G_AU_YR
    _configure_integrator(sim, integrator, dt)
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
            hit = False
            if use_hill:
                planet_h = (
                    np.array(
                        [
                            [sim.particles[k].x, sim.particles[k].y, sim.particles[k].z]
                            for k in range(1, 5)
                        ]
                    )
                    - sun
                )
                tp_h = np.array([p.x, p.y, p.z]) - sun
                hit = any(
                    np.linalg.norm(tp_h - planet_h[k]) < hill_radius(float(a_p[k]), float(m_p[k]))
                    for k in range(4)
                )
            if (not np.isfinite(orb.a)) or r > r_escape_au or q < q_min_au or hit:
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
        "integrator": integrator.lower().strip(),
        "escape_cfg": {
            "r_escape_au": r_escape_au,
            "q_min_au": q_min_au,
            "use_hill": use_hill,
        },
    }


def integrate_clone_snapshots(
    elems: dict,
    ephem: DE441Ephem | CachedJ2000Ephem,
    t_max_yr: float,
    snapshot_myr: list[float],
    dt: float = 0.25,
    sample_yr: float = 300.0,
    r_escape_au: float = 1.0e4,
    q_min_au: float = 2.5,
    use_hill: bool = True,
    integrator: str = "whfast",
) -> dict:
    """Integrate one clone; record $a$ at fixed times (Myr) while still bound."""
    del sample_yr  # snapshot-only path; no 300 yr series stored
    snaps_myr = sorted({float(x) for x in snapshot_myr if 0 < x <= t_max_yr / 1e6})
    snap_yr = [t * 1e6 for t in snaps_myr]

    m = masses_msun()
    state = giant_state(ephem, JD_J2000)
    sim = rebound.Simulation()
    sim.G = G_AU_YR
    _configure_integrator(sim, integrator, dt)
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

    captured: list[tuple[float, float]] = []
    escaped = False

    def _escaped_now() -> bool:
        primary = sim.particles[0]
        p = sim.particles[tp_i]
        orb = p.orbit(primary=primary)
        sun = np.array([primary.x, primary.y, primary.z])
        r = float(np.linalg.norm(np.array([p.x, p.y, p.z]) - sun))
        q = orb.a * (1.0 - orb.e) if np.isfinite(orb.a) else 0.0
        if (not np.isfinite(orb.a)) or r > r_escape_au or q < q_min_au:
            return True
        if use_hill:
            planet_h = (
                np.array(
                    [
                        [sim.particles[k].x, sim.particles[k].y, sim.particles[k].z]
                        for k in range(1, 5)
                    ]
                )
                - sun
            )
            tp_h = np.array([p.x, p.y, p.z]) - sun
            if any(
                np.linalg.norm(tp_h - planet_h[k]) < hill_radius(float(a_p[k]), float(m_p[k]))
                for k in range(4)
            ):
                return True
        return False

    try:
        for t_myr, t_goal in zip(snaps_myr, snap_yr):
            if escaped:
                break
            if t_goal > sim.t:
                sim.integrate(t_goal)
            if _escaped_now():
                escaped = True
                break
            orb = sim.particles[tp_i].orbit(primary=sim.particles[0])
            if np.isfinite(orb.a):
                captured.append((t_myr, float(orb.a)))
    except Exception:
        escaped = True

    lifetime = float(sim.t)
    return {
        "t_yr": np.asarray([], dtype=np.float64),
        "a_au": np.asarray([], dtype=np.float64),
        "lifetime_yr": lifetime,
        "escaped": bool(escaped),
        "survived_full": (not escaped) and lifetime >= t_max_yr * 0.999,
        "integrator": integrator.lower().strip(),
        "snapshots_myr": captured,
        "escape_cfg": {
            "r_escape_au": r_escape_au,
            "q_min_au": q_min_au,
            "use_hill": use_hill,
        },
    }


def hurst_curve(a: np.ndarray, t: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if len(a) < 20:
        return np.array([]), np.array([])
    n = len(a)
    ws = np.unique(np.geomspace(10, max(10, n // 4), num=12).astype(int))
    log_w = []
    log_s = []
    for w in ws:
        if w < 5 or w >= n:
            continue
        sigs = []
        for i0 in range(0, n - w, max(1, w // 4)):
            sigs.append(np.std(a[i0 : i0 + w]))
        if not sigs:
            continue
        sbar = float(np.mean(sigs))
        if sbar <= 0:
            continue
        log_w.append(np.log10(w * 300.0))
        log_s.append(np.log10(sbar))
    return np.asarray(log_w), np.asarray(log_s)


DEFAULT_RULE = "v2"  # BM09 lifetime-faithful; v1 kept for sensitivity only


def hurst_diagnostics(a: np.ndarray, t: np.ndarray) -> dict:
    """Hurst / flatness diagnostics shared by all rule sets."""
    log_w, log_s = hurst_curve(a, t)
    H = np.nan
    r2 = np.nan
    nonlinear = False
    rel_std = float(np.std(a) / (np.mean(a) + 1e-9)) if len(a) else np.nan
    # Flatness is an a(t) property. Do NOT require |H|≪1: near-constant series
    # make log-σ ill-conditioned and can yield spurious H.
    flat = bool(np.isfinite(rel_std) and rel_std < 0.03 and len(a) >= 20)

    if len(log_w) >= 4:
        coeff = np.polyfit(log_w, log_s, 1)
        H = float(coeff[0])
        pred = np.polyval(coeff, log_w)
        ss_res = float(np.sum((log_s - pred) ** 2))
        ss_tot = float(np.sum((log_s - log_s.mean()) ** 2)) + 1e-30
        r2 = 1.0 - ss_res / ss_tot
        coeff2 = np.polyfit(log_w, log_s, 2)
        pred2 = np.polyval(coeff2, log_w)
        ss2 = float(np.sum((log_s - pred2) ** 2))
        nonlinear = (ss_res > 1.5 * ss2 + 1e-6) and (r2 < 0.92)

    return {
        "H": H,
        "r2": r2,
        "nonlinear": bool(nonlinear),
        "flat": flat,
        "rel_std": rel_std,
    }


def assign_class_v1(life_myr: float, flat: bool, nonlinear: bool, r2: float, H: float) -> str:
    """Legacy Hurst-gated rules (pre-stopcheck). Biases long-lived linear survivors to D."""
    if flat and life_myr >= LONG_LIVED_MYR:
        return "Q"
    if life_myr >= LONG_LIVED_MYR and (nonlinear or (np.isfinite(r2) and r2 < 0.85)):
        return "R"
    if np.isfinite(r2) and r2 >= 0.85 and np.isfinite(H):
        return "D"
    if life_myr < LONG_LIVED_MYR:
        return "D"
    return "R" if nonlinear else "D"


def assign_class_v2(life_myr: float, flat: bool) -> str:
    """BM09-faithful lifetime correlation: long+flat→Q; long→R; short→D.

    Matches Bailey & Malhotra's empirical link (short life ↔ D; long life ↔ R),
    with flat a(t) reserved for Q. Hurst nonlinearity is diagnostic, not a gate.
    """
    if life_myr >= LONG_LIVED_MYR and flat:
        return "Q"
    if life_myr >= LONG_LIVED_MYR:
        return "R"
    return "D"


def classify_bailey(
    a: np.ndarray,
    t: np.ndarray,
    lifetime_yr: float,
    rule: str = DEFAULT_RULE,
) -> dict:
    """Classify a clone. Default rule is v2. Always returns class_v1 and class_v2."""
    life_myr = lifetime_yr / 1.0e6
    diag = hurst_diagnostics(a, t)
    class_v1 = assign_class_v1(
        life_myr, diag["flat"], diag["nonlinear"], diag["r2"], diag["H"]
    )
    class_v2 = assign_class_v2(life_myr, diag["flat"])
    if rule not in ("v1", "v2"):
        raise ValueError(f"unknown classifier rule: {rule}")
    cls = class_v1 if rule == "v1" else class_v2
    return {
        "class": cls,
        "class_v1": class_v1,
        "class_v2": class_v2,
        "rule": rule,
        "H": diag["H"],
        "r2": diag["r2"],
        "nonlinear": diag["nonlinear"],
        "flat": diag["flat"],
        "rel_std": diag["rel_std"],
        "lifetime_myr": life_myr,
        "long_lived": life_myr >= LONG_LIVED_MYR,
        "a_mean": float(np.mean(a)) if len(a) else np.nan,
        "a_std": float(np.std(a)) if len(a) else np.nan,
    }


def reclassify_from_diagnostics(rec: dict, rule: str = DEFAULT_RULE) -> str:
    """Re-label a stored clone record that already has lifetime/flat/r2/H fields."""
    life = float(rec.get("lifetime_myr") or 0.0)
    flat = bool(rec.get("flat"))
    # Prefer stored flat; if absent, fall back to rel_std / a_std
    if "flat" not in rec:
        a_mean = rec.get("a_mean")
        a_std = rec.get("a_std")
        if a_mean is not None and a_std is not None:
            flat = float(a_std) / (float(a_mean) + 1e-9) < 0.03
    if rule == "v2":
        return assign_class_v2(life, flat)
    return assign_class_v1(
        life,
        flat,
        bool(rec.get("nonlinear")),
        float(rec["r2"]) if rec.get("r2") is not None else np.nan,
        float(rec["H"]) if rec.get("H") is not None else np.nan,
    )

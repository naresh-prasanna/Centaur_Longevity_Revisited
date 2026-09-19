"""DE441 heliocentric giant-planet positions (AU) for Option I."""
from __future__ import annotations

from pathlib import Path

import numpy as np
from jplephem.spk import SPK

# Zenith split kernels (same handoff as option_b_recurrence)
DEFAULT_KERNEL_DIR = Path(r"C:\Users\Hp\Downloads\Zenith")
SPLIT_JD = 2440416.5
JD_J2000 = 2451545.0
DAYS_PER_YEAR = 365.25

# AU in km (IAU)
KM_PER_AU = 149597870.700

# DE441 part-2 sun/planet coverage ends ~year 17191 ≈ +15.191 kyr from J2000
DE441_T_MAX_YR = 15190.0

# J2000 mean obliquity — rotate ICRF/equatorial SPK vectors into ecliptic
EPS_J2000 = np.deg2rad(23.4392911)
_CE, _SE = np.cos(EPS_J2000), np.sin(EPS_J2000)


def eq_to_ecl(x: np.ndarray, v: np.ndarray | None = None):
    """Rotate equatorial (ICRF) position/velocity into J2000 ecliptic."""
    x = np.asarray(x, dtype=np.float64)
    xe = np.array([x[0], x[1] * _CE + x[2] * _SE, -x[1] * _SE + x[2] * _CE])
    if v is None:
        return xe
    v = np.asarray(v, dtype=np.float64)
    ve = np.array([v[0], v[1] * _CE + v[2] * _SE, -v[1] * _SE + v[2] * _CE])
    return xe, ve


GIANT_NAIF = {
    "jupiter": 5,
    "saturn": 6,
    "uranus": 7,
    "neptune": 8,
}
GIANT_NAMES = ("jupiter", "saturn", "uranus", "neptune")

# GM in AU^3 / day^2 (DE440/441 style solar-system values)
GM_SUN_AU3_DAY2 = 0.2959122082855911e-3
GM_PLANET_AU3_DAY2 = {
    "jupiter": 0.282534590952422e-6,
    "saturn": 0.845970607324503e-7,
    "uranus": 0.129202482578296e-7,
    "neptune": 0.152435734788511e-7,
}


def gm_sun_au3_yr2() -> float:
    return GM_SUN_AU3_DAY2 * (DAYS_PER_YEAR**2)


def gm_planets_au3_yr2() -> dict[str, float]:
    s = DAYS_PER_YEAR**2
    return {k: v * s for k, v in GM_PLANET_AU3_DAY2.items()}


def masses_msun() -> dict[str, float]:
    """Planet masses in solar masses (GM_p / GM_sun)."""
    gms = gm_sun_au3_yr2()
    return {k: gm_planets_au3_yr2()[k] / gms for k in GIANT_NAMES}


class DE441Ephem:
    """Open split DE441 kernels; heliocentric giant XYZ in AU (ecliptic)."""

    def __init__(self, kernel_dir: Path | None = None):
        self.kernel_dir = Path(kernel_dir) if kernel_dir else DEFAULT_KERNEL_DIR
        self.s1 = SPK.open(str(self.kernel_dir / "de441_part-1.bsp"))
        self.s2 = SPK.open(str(self.kernel_dir / "de441_part-2.bsp"))

    def close(self):
        self.s1.close()
        self.s2.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def _spk_for_jd(self, jd: float) -> SPK:
        return self.s1 if jd < SPLIT_JD else self.s2

    def _bary_km(self, spk: SPK, body: int, jd: float) -> np.ndarray:
        return np.asarray(spk[0, body].compute(jd), dtype=np.float64)

    def _bary_km_state(self, spk: SPK, body: int, jd: float):
        x, v = spk[0, body].compute_and_differentiate(jd)
        return np.asarray(x, dtype=np.float64), np.asarray(v, dtype=np.float64)

    def heliocentric_au(self, jd: float) -> dict[str, np.ndarray]:
        spk = self._spk_for_jd(jd)
        sun = self._bary_km(spk, 10, jd)
        out = {}
        for name, nid in GIANT_NAIF.items():
            out[name] = eq_to_ecl((self._bary_km(spk, nid, jd) - sun) / KM_PER_AU)
        return out

    def heliocentric_state_au_day(self, jd: float) -> dict[str, tuple[np.ndarray, np.ndarray]]:
        """Heliocentric ecliptic (x AU, v AU/day)."""
        spk = self._spk_for_jd(jd)
        sun_x, sun_v = self._bary_km_state(spk, 10, jd)
        out = {}
        for name, nid in GIANT_NAIF.items():
            x, v = self._bary_km_state(spk, nid, jd)
            x = (x - sun_x) / KM_PER_AU
            v = (v - sun_v) / KM_PER_AU
            out[name] = eq_to_ecl(x, v)
        return out

    def heliocentric_au_array(self, jd: np.ndarray) -> dict[str, np.ndarray]:
        """Return dict name -> (3, N) AU."""
        jd = np.asarray(jd, dtype=np.float64)
        out = {n: np.empty((3, len(jd)), dtype=np.float64) for n in GIANT_NAMES}
        for i, t in enumerate(jd):
            pos = self.heliocentric_au(float(t))
            for n in GIANT_NAMES:
                out[n][:, i] = pos[n]
        return out

    def mean_semi_major_au(self, year0: float, year1: float, n_samp: int = 512) -> dict[str, float]:
        """Mean heliocentric distance over [year0, year1] as circular-control a."""
        year1 = min(year1, 2000.0 + DE441_T_MAX_YR - 1.0)
        year0 = min(year0, year1 - 1.0)
        jd = JD_J2000 + np.linspace(year0 - 2000.0, year1 - 2000.0, n_samp) * DAYS_PER_YEAR
        pos = self.heliocentric_au_array(jd)
        means = {}
        for n in GIANT_NAMES:
            r = np.sqrt((pos[n] ** 2).sum(axis=0))
            means[n] = float(r.mean())
        return means


def heliocentric_state_au_day(
    ephem: "DE441Ephem | CachedJ2000Ephem", jd: float
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Module-level wrapper for Option J: giant_state(ephem, jd)."""
    return ephem.heliocentric_state_au_day(jd)


class CachedJ2000Ephem:
    """Kernel-free giant ICs at J2000 from a Horizons/DE441 JSON cache.

    Only supports heliocentric_state_au_day at JD_J2000 (Option J WHFast init).
    """

    def __init__(self, cache_path: Path | None = None):
        if cache_path is None:
            cache_path = (
                Path(__file__).resolve().parents[1]
                / "data"
                / "j2000_giants_horizons_de441.json"
            )
        self.cache_path = Path(cache_path)
        import json

        raw = json.loads(self.cache_path.read_text(encoding="utf-8"))
        self.jd = float(raw["jd_tdb"])
        self._state: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        for name in GIANT_NAMES:
            b = raw["bodies"][name]
            x = np.asarray(b["x_au"], dtype=np.float64)
            v = np.asarray(b["v_au_day"], dtype=np.float64)
            self._state[name] = (x, v)

    def close(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def heliocentric_state_au_day(self, jd: float) -> dict[str, tuple[np.ndarray, np.ndarray]]:
        if abs(float(jd) - self.jd) > 1e-6:
            raise ValueError(f"CachedJ2000Ephem only supports JD={self.jd}, got {jd}")
        return {k: (v[0].copy(), v[1].copy()) for k, v in self._state.items()}


def open_giant_ephem(kernel_dir: Path | None = None, cache_path: Path | None = None):
    """Prefer local DE441 kernels; fall back to J2000 JSON cache."""
    kdir = Path(kernel_dir) if kernel_dir else DEFAULT_KERNEL_DIR
    if (kdir / "de441_part-1.bsp").exists() and (kdir / "de441_part-2.bsp").exists():
        return DE441Ephem(kdir)
    return CachedJ2000Ephem(cache_path)


def _rv_to_elements(r: np.ndarray, v: np.ndarray, mu: float) -> dict:
    """Osculating Keplerian elements from heliocentric r (AU), v (AU/yr)."""
    rr = np.linalg.norm(r)
    vv = np.linalg.norm(v)
    h = np.cross(r, v)
    hh = np.linalg.norm(h)
    n_vec = np.array([-h[1], h[0], 0.0])
    nn = np.linalg.norm(n_vec)
    e_vec = np.cross(v, h) / mu - r / rr
    e = float(np.linalg.norm(e_vec))
    a = 1.0 / (2.0 / rr - vv * vv / mu)
    inc = float(np.arccos(np.clip(h[2] / hh, -1.0, 1.0)))
    if nn > 1e-14:
        Omega = float(np.arctan2(n_vec[1], n_vec[0]))
        omega = float(np.arctan2(np.dot(np.cross(n_vec, e_vec), h) / hh, np.dot(n_vec, e_vec)))
    else:
        Omega = 0.0
        omega = float(np.arctan2(e_vec[1], e_vec[0]))
    if e > 1e-12:
        cos_nu = np.clip(np.dot(e_vec, r) / (e * rr), -1.0, 1.0)
        nu = float(np.arccos(cos_nu))
        if np.dot(r, v) < 0:
            nu = 2 * np.pi - nu
        E = 2 * np.arctan2(np.sqrt(1 - e) * np.sin(nu / 2), np.sqrt(1 + e) * np.cos(nu / 2))
        # better: eccentric anomaly from true anomaly
        E = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(nu / 2))
        M = E - e * np.sin(E)
    else:
        M = float(np.arctan2(r[1], r[0]) - Omega - omega)
    return {"a": float(a), "e": e, "inc": inc, "Omega": Omega, "omega": omega, "M0": float(M)}


def _elements_to_r(elem: dict, t_yr: float, t0_yr: float, mu: float) -> np.ndarray:
    """Position at t_yr from elements referred to epoch t0_yr (M = M0 + n(t-t0))."""
    a, e = elem["a"], elem["e"]
    n = np.sqrt(mu / a**3)
    M = elem["M0"] + n * (t_yr - t0_yr)
    # Solve Kepler
    E = M
    for _ in range(12):
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    cosE, sinE = np.cos(E), np.sin(E)
    x_orb = a * (cosE - e)
    y_orb = a * np.sqrt(1 - e * e) * sinE
    cosO, sinO = np.cos(elem["Omega"]), np.sin(elem["Omega"])
    cosw, sinw = np.cos(elem["omega"]), np.sin(elem["omega"])
    cosi, sini = np.cos(elem["inc"]), np.sin(elem["inc"])
    R11 = cosO * cosw - sinO * sinw * cosi
    R12 = -cosO * sinw - sinO * cosw * cosi
    R21 = sinO * cosw + cosO * sinw * cosi
    R22 = -sinO * sinw + cosO * cosw * cosi
    R31 = sinw * sini
    R32 = cosw * sini
    return np.array([R11 * x_orb + R12 * y_orb, R21 * x_orb + R22 * y_orb, R31 * x_orb + R32 * y_orb])


class EphemTable:
    """DE441 table with two-body Keplerian continuation beyond kernel coverage."""

    def __init__(self, ephem: DE441Ephem, t_years: np.ndarray):
        self.t = np.asarray(t_years, dtype=np.float64)
        jd = JD_J2000 + self.t * DAYS_PER_YEAR
        raw = ephem.heliocentric_au_array(jd)
        self.xyz = np.stack([raw[n] for n in GIANT_NAMES], axis=0)  # (4,3,N)
        self.t0_cont = float(self.t[-1])
        mu = gm_sun_au3_yr2()
        state = ephem.heliocentric_state_au_day(JD_J2000 + self.t0_cont * DAYS_PER_YEAR)
        self.cont_elems = []
        for name in GIANT_NAMES:
            x, v_day = state[name]
            v_yr = v_day * DAYS_PER_YEAR
            self.cont_elems.append(_rv_to_elements(x, v_yr, mu))
        self.mu = mu

    def positions(self, t_year: float) -> np.ndarray:
        """Return (4, 3) AU at time t_year from J2000."""
        if t_year <= self.t[0]:
            return self.xyz[:, :, 0].copy()
        if t_year <= self.t[-1]:
            i = int(np.searchsorted(self.t, t_year) - 1)
            i = max(0, min(i, len(self.t) - 2))
            t0, t1 = self.t[i], self.t[i + 1]
            u = (t_year - t0) / (t1 - t0)
            return (1.0 - u) * self.xyz[:, :, i] + u * self.xyz[:, :, i + 1]
        # Beyond DE441: Keplerian continuation from last tabulated epoch
        out = np.zeros((4, 3), dtype=np.float64)
        for k, elem in enumerate(self.cont_elems):
            out[k] = _elements_to_r(elem, t_year, self.t0_cont, self.mu)
        return out


def build_ephem_table(ephem: DE441Ephem, t_end: float, sample_yr: float = 0.25) -> EphemTable:
    """Tabulate DE441 giants up to min(t_end, DE441_T_MAX_YR); continue analytically after."""
    t_tab = min(float(t_end), DE441_T_MAX_YR)
    n = int(np.ceil(t_tab / sample_yr)) + 1
    t = np.linspace(0.0, t_tab, n)
    return EphemTable(ephem, t)


def save_ephem_table(table: EphemTable, path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        t=table.t,
        xyz=table.xyz,
        t0_cont=np.array([table.t0_cont]),
        mu=np.array([table.mu]),
        cont_a=np.array([e["a"] for e in table.cont_elems]),
        cont_e=np.array([e["e"] for e in table.cont_elems]),
        cont_inc=np.array([e["inc"] for e in table.cont_elems]),
        cont_Omega=np.array([e["Omega"] for e in table.cont_elems]),
        cont_omega=np.array([e["omega"] for e in table.cont_elems]),
        cont_M0=np.array([e["M0"] for e in table.cont_elems]),
    )


def load_ephem_table(path: Path) -> EphemTable:
    z = np.load(path)
    table = object.__new__(EphemTable)
    table.t = z["t"]
    table.xyz = z["xyz"]
    table.t0_cont = float(z["t0_cont"][0])
    table.mu = float(z["mu"][0])
    table.cont_elems = [
        {
            "a": float(z["cont_a"][k]),
            "e": float(z["cont_e"][k]),
            "inc": float(z["cont_inc"][k]),
            "Omega": float(z["cont_Omega"][k]),
            "omega": float(z["cont_omega"][k]),
            "M0": float(z["cont_M0"][k]),
        }
        for k in range(4)
    ]
    return table


def circular_n_yr(a_au: float) -> float:
    """Mean motion (rad/yr) for circular Keplerian orbit about the Sun."""
    return float(np.sqrt(gm_sun_au3_yr2() / a_au**3))


def circular_positions(t_year: float, a_giants: dict[str, float], phases: np.ndarray) -> np.ndarray:
    """(4, 3) AU for circular coplanar giants at time t_year."""
    out = np.zeros((4, 3), dtype=np.float64)
    for k, name in enumerate(GIANT_NAMES):
        a = a_giants[name]
        n = circular_n_yr(a)
        M = float(phases[k]) + n * t_year
        out[k, 0] = a * np.cos(M)
        out[k, 1] = a * np.sin(M)
        out[k, 2] = 0.0
    return out

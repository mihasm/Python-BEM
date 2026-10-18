# Python BEM - Blade Element Momentum Theory Software.

# Copyright (C) 2022 M. Smrekar

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
Loads the experimental validation cases in validation/cases/ and runs the BEM solver on them.

A case is a folder with case.json, geometry.csv, airfoils/ and experiment/ (see validation/README.md).
"""

import contextlib
import io
import json
import math
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

# BEM_SOLVER_ROOT runs the cases against another copy of the solver, e.g. an export of an older commit
# (git archive <rev> | tar -x -C /tmp/bem-old), to compare solver versions on the same data.
ROOT = Path(os.environ.get("BEM_SOLVER_ROOT") or Path(__file__).resolve().parents[1]).resolve()
CASES_DIR = Path(__file__).resolve().parent / "cases"
(ROOT / "foils").mkdir(exist_ok=True)  # the solver writes airfoil .dat files there
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from calculation import Calculator  # noqa: E402
from montgomerie import Montgomerie  # noqa: E402

# Solver configurations compared by the runner. Keys are the run_array arguments they change.
# method: 0 classic, 1 Spera, 2 Buhl (AeroDyn), 3 Buhl (QBlade), 4 Wilson-Walker, 5 modified ABS
# tip_loss_mode: 0 none, 1 Prandtl, 2 Shen, 3 Adkins; hub_loss_mode: 0 none, 1 Prandtl, 2 Shen
# rotational_augmentation_correction_method: 0 Snel, 1 Du-Selig, 2 Chaviaropoulos-Hansen, 3 Lindenburg,
#   4 Dumitrescu-Cardos, 5 Snel (propeller), 6 Gur-Rosen (propeller)
VARIANTS = {
    "buhl_prandtl": {"method": 2, "tip_loss_mode": 1, "hub_loss_mode": 1},
    "qblade_prandtl": {"method": 3, "tip_loss_mode": 1, "hub_loss_mode": 1},
    "classic_prandtl": {"method": 0, "tip_loss_mode": 1, "hub_loss_mode": 1},
    "spera_prandtl": {"method": 1, "tip_loss_mode": 1, "hub_loss_mode": 1},
    "wilson_walker_prandtl": {"method": 4, "tip_loss_mode": 1, "hub_loss_mode": 1},
    "abs_prandtl": {"method": 5, "tip_loss_mode": 1, "hub_loss_mode": 1},
    "buhl_shen": {"method": 2, "tip_loss_mode": 2, "hub_loss_mode": 2},
    "buhl_noloss": {"method": 2, "tip_loss_mode": 0, "hub_loss_mode": 0},
    "buhl_prandtl_snel": {"method": 2, "tip_loss_mode": 1, "hub_loss_mode": 1,
                          "rotational_augmentation_correction": True,
                          "rotational_augmentation_correction_method": 0},
    "buhl_prandtl_dusel": {"method": 2, "tip_loss_mode": 1, "hub_loss_mode": 1,
                           "rotational_augmentation_correction": True,
                           "rotational_augmentation_correction_method": 1},
    "buhl_prandtl_chav": {"method": 2, "tip_loss_mode": 1, "hub_loss_mode": 1,
                          "rotational_augmentation_correction": True,
                          "rotational_augmentation_correction_method": 2},
    "buhl_prandtl_lind": {"method": 2, "tip_loss_mode": 1, "hub_loss_mode": 1,
                          "rotational_augmentation_correction": True,
                          "rotational_augmentation_correction_method": 3},
    # propellers: the solver accepts only the classic momentum balance (method 0)
    "prop_prandtl": {"method": 0, "tip_loss_mode": 1, "hub_loss_mode": 1},
    "prop_adkins": {"method": 0, "tip_loss_mode": 3, "hub_loss_mode": 1},
    "prop_noloss": {"method": 0, "tip_loss_mode": 0, "hub_loss_mode": 0},
    "prop_prandtl_mach": {"method": 0, "tip_loss_mode": 1, "hub_loss_mode": 1, "mach_number_correction": True},
    "prop_prandtl_snel": {"method": 0, "tip_loss_mode": 1, "hub_loss_mode": 1,
                          "rotational_augmentation_correction": True,
                          "rotational_augmentation_correction_method": 5},
    "prop_prandtl_gur": {"method": 0, "tip_loss_mode": 1, "hub_loss_mode": 1,
                         "rotational_augmentation_correction": True,
                         "rotational_augmentation_correction_method": 6},
}
PROPELLER_VARIANTS = [k for k in VARIANTS if k.startswith("prop_")]
TURBINE_VARIANTS = [k for k in VARIANTS if not k.startswith("prop_")]
DEFAULT_VARIANTS = ["buhl_prandtl", "classic_prandtl", "buhl_shen", "prop_prandtl", "prop_adkins"]


def variants_for(case, variants):
    """The requested variants that apply to the case's rotor type."""
    allowed = PROPELLER_VARIANTS if case.is_propeller else TURBINE_VARIANTS
    return [v for v in variants if v in allowed]

SOLVER_DEFAULTS = {
    "method": 2,
    "tip_loss_mode": 1,
    "hub_loss_mode": 1,
    "mach_number_correction": False,
    "cascade_correction": False,
    "rotational_augmentation_correction": False,
    "rotational_augmentation_correction_method": 0,
    "skewed_wake_correction": False,
    "max_iterations": 300,
    "convergence_limit": 1e-5,
    "relaxation_factor": 0.3,
    "a_initial": 0.3,
    "aprime_initial": 0.01,
    "fix_reynolds": False,
    "reynolds": 1e6,
    "use_minimization_solver": False,
    "invert_alpha": False,
    "print_out": False,
    "print_all": False,
    "blade_design": 0,  # solid cross-section
    "blade_thickness": 0.001,
    "mass_density": 1000.0,
}

POLAR_NCRIT = 0.0  # every polar is stored under one Ncrit key; the solver filters by it


@dataclass
class Experiment:
    file: str
    kind: str
    data: pd.DataFrame
    meta: dict
    defaults: dict = field(default_factory=dict)

    @property
    def name(self):
        return Path(self.file).stem


@dataclass
class Case:
    path: Path
    meta: dict
    geometry: pd.DataFrame
    airfoils: dict
    experiments: list

    @property
    def id(self):
        return self.meta.get("id", self.path.name)

    @property
    def is_propeller(self):
        return int(self.meta.get("turbine_type", 0)) == 1


def read_csv(path):
    """Reads a case CSV; header lines starting with # are comments."""
    df = pd.read_csv(path, comment="#", skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    return df


def read_csv_comments(path):
    out = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                out.append(line[1:].strip())
            elif line.strip():
                break
    return out


def read_selig_dat(path):
    """Reads Selig-format airfoil coordinates (Lednicer format is converted)."""
    rows = []
    with open(path, encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()
    for line in lines[1:]:
        parts = line.replace(",", " ").split()
        if len(parts) < 2:
            continue
        try:
            rows.append((float(parts[0]), float(parts[1])))
        except ValueError:
            continue
    xy = np.array(rows)
    if len(xy) and xy[:, 0].max() > 1.5:
        # Lednicer: first row holds the point counts of upper and lower side
        n_up, n_low = int(xy[0, 0]), int(xy[0, 1])
        up, low = xy[1:1 + n_up], xy[1 + n_up:1 + n_up + n_low]
        xy = np.vstack([up[::-1], low[1:]])
    return xy[:, 0], xy[:, 1]


def load_case(path):
    path = Path(path)
    with open(path / "case.json", encoding="utf-8") as f:
        meta = json.load(f)
    geometry_file = path / (meta.get("geometry_file") or "geometry.csv")
    geometry = read_csv(geometry_file) if geometry_file.exists() else None
    airfoils = {}
    for name, a in meta.get("airfoils", {}).items():
        entry = {"meta": a}
        if a.get("coordinates"):
            entry["x"], entry["y"] = read_selig_dat(path / a["coordinates"])
        if a.get("polars"):
            entry["polars"] = read_csv(path / a["polars"])
        airfoils[name] = entry
    experiments = []
    for e in meta.get("experiments", []):
        f = path / e["file"]
        if not f.exists():
            continue
        data = read_csv(f) if f.suffix == ".csv" else None
        experiments.append(Experiment(file=e["file"], kind=e.get("kind", "performance"), data=data, meta=e,
                                      defaults=e.get("defaults", {})))
    return Case(path=path, meta=meta, geometry=geometry, airfoils=airfoils, experiments=experiments)


def list_cases(names=None, cases_dir=CASES_DIR):
    out = []
    for p in sorted(Path(cases_dir).iterdir()):
        if (p / "case.json").exists() and (names is None or p.name in names):
            out.append(p)
    return out


# ---------------------------------------------------------------- airfoils

def _circle_coordinates(n=60):
    t = np.linspace(0, 2 * np.pi, n)
    return 0.5 + 0.5 * np.cos(t), 0.5 * np.sin(t)


def _extrapolate_polar(x, y, alpha, cl, cd, re, cd90=1.5):
    """Resamples a measured polar to 1 deg and extends it to +-180 deg with the program's Montgomerie method."""
    order = np.argsort(alpha)
    alpha, cl, cd = np.asarray(alpha)[order], np.asarray(cl)[order], np.asarray(cd)[order]
    alpha, idx = np.unique(alpha, return_index=True)
    cl, cd = cl[idx], cd[idx]
    grid = np.arange(math.ceil(alpha.min()), math.floor(alpha.max()) + 1, 1.0)
    cl_g, cd_g = np.interp(grid, alpha, cl), np.interp(grid, alpha, cd)
    if len(grid) < 8 or not (grid.min() < 0 < grid.max()):
        raise ValueError("polar at Re=%g covers too little of the attached range (%s..%s deg)" %
                         (re, alpha.min(), alpha.max()))
    m = Montgomerie(x=list(x), y=list(y), alpha=grid, Cl=cl_g, Cd=cd_g, Re=re, A=-5, Am=8, B=5, Bm=5,
                    m_CD90=cd90, slope=0.106)
    with contextlib.redirect_stdout(io.StringIO()):
        a_ex, cl_ex, cd_ex = m.calculate_extrapolation()
    a_ex, cl_ex, cd_ex = np.asarray(a_ex, float), np.asarray(cl_ex, float), np.asarray(cd_ex, float)
    inside = (a_ex >= grid.min()) & (a_ex <= grid.max())
    a_out = np.concatenate([a_ex[~inside], grid])
    cl_out = np.concatenate([cl_ex[~inside], cl_g])
    cd_out = np.concatenate([cd_ex[~inside], cd_g])
    o = np.argsort(a_out)
    a_out, cl_out, cd_out = a_out[o], cl_out[o], cd_out[o]
    a_out, idx = np.unique(a_out, return_index=True)
    return a_out, cl_out[idx], cd_out[idx]


def build_airfoils(case, extrapolate=True, cd90=1.5):
    """Returns the airfoils dict in the format the solver (utils.get_curves_functions) expects."""
    out = {}
    for name, a in case.airfoils.items():
        polars = a.get("polars")
        is_circle = polars is None or name.lower().startswith(("circle", "cylinder"))
        if "x" in a:
            x, y = np.asarray(a["x"]), np.asarray(a["y"])
        else:
            x, y = _circle_coordinates()
        gathered, re_list, smin, smax = [], [], [], []
        if polars is None:
            cd_cyl = float(a["meta"].get("cd", 1.0))
            for al in (-180.0, 0.0, 180.0):
                gathered.append([1e6, POLAR_NCRIT, al, 0.0, cd_cyl])
            re_list, smin, smax = [1e6], [-5.0], [10.0]
        else:
            col = {c.lower(): c for c in polars.columns}
            re_col = col.get("re")
            for re, grp in (polars.groupby(re_col) if re_col else [(1e6, polars)]):
                alpha = grp[col["alpha_deg"]].to_numpy(float)
                cl = grp[col["cl"]].to_numpy(float)
                cd = grp[col["cd"]].to_numpy(float)
                if is_circle or np.allclose(cl, 0):
                    al, cl_, cd_ = np.array([-180.0, 0.0, 180.0]), np.zeros(3), np.full(3, np.mean(cd))
                elif extrapolate:
                    al, cl_, cd_ = _extrapolate_polar(x, y, alpha, cl, cd, float(re), cd90=cd90)
                else:
                    al, cl_, cd_ = alpha, cl, cd
                for i in range(len(al)):
                    gathered.append([float(re), POLAR_NCRIT, float(al[i]), float(cl_[i]), float(cd_[i])])
                re_list.append(float(re))
                attached = (alpha > -15) & (alpha < 25)
                smax.append(float(alpha[attached][np.argmax(cl[attached])]) if attached.any() else 10.0)
                smin.append(float(alpha[attached][np.argmin(cl[attached])]) if attached.any() else -5.0)
        thickness = float(np.max(y) - np.min(y)) if "x" in a else 1.0
        out[name] = {
            "x": list(map(float, x)), "y": list(map(float, y)),
            "max_thickness": thickness,
            "link": a["meta"].get("source", ""),
            "interp_function_cl": None, "interp_function_cd": None,
            "curves": [], "gathered_curves": gathered,
            "centroid_x": 0.25, "centroid_y": 0.0,
            "ncrit_selected": POLAR_NCRIT, "extrapolation_bool": 2,
            "stall_angles": [re_list, smin, smax],
        }
    return out


# ---------------------------------------------------------------- geometry

def blade_sections(case, n_sections=12, spacing="cosine"):
    """
    Divides the blade from its first tabulated station to the tip into annuli and interpolates the geometry at
    their midpoints. Returns r, dr, c, theta, foils.
    """
    g = case.geometry
    r_tab = g["r_m"].to_numpy(float)
    c_tab = g["chord_m"].to_numpy(float)
    t_tab = g["twist_deg"].to_numpy(float)
    f_tab = list(g["airfoil"].astype(str).str.strip())
    R = float(case.meta["R"])
    r0 = max(float(r_tab[0]), float(case.meta.get("Rhub", 0.0)))
    if spacing == "cosine":
        # half-cosine: fine at the tip (where loads change fastest), coarse at the root
        s = np.sin(np.linspace(0, np.pi / 2, n_sections + 1))
    else:
        s = np.linspace(0, 1, n_sections + 1)
    edges = r0 + (R - r0) * s
    r = 0.5 * (edges[1:] + edges[:-1])
    dr = np.diff(edges)
    c = np.interp(r, r_tab, c_tab)
    theta = np.interp(r, r_tab, t_tab)
    foils = []
    for _r in r:
        i = int(np.argmin(np.abs(r_tab - _r)))
        foils.append(f_tab[i])
    # a 'transition' row needs named airfoils on both sides; keep the ends named
    if foils[0] == "transition":
        foils[0] = next(f for f in f_tab if f != "transition")
    if foils[-1] == "transition":
        foils[-1] = next(f for f in reversed(f_tab) if f != "transition")
    return r, dr, c, theta, foils


# ---------------------------------------------------------------- running

class CaseRunner:
    def __init__(self, case, variant="buhl_prandtl", n_sections=12, solver_overrides=None):
        self.case = case
        self.variant = variant
        self.solver = {**SOLVER_DEFAULTS, **VARIANTS[variant], **(solver_overrides or {})}
        self.r, self.dr, self.c, self.theta, self.foils = blade_sections(case, n_sections)
        self.airfoils = build_airfoils(case)
        fluid = case.meta.get("fluid", {})
        self.rho = float(fluid.get("rho", 1.225))
        self.nu = float(fluid.get("kin_viscosity", 1.4607e-5))
        water = str(fluid.get("medium", "air")).lower() in ("water", "sea water", "seawater", "fresh water")
        self.speed_of_sound = float(fluid.get("speed_of_sound", 1480.0 if water else 343.0))
        self.calculator = None

    def _args(self, v, rpm, pitch, yaw, rho, nu):
        m = self.case.meta
        return {
            **self.solver,
            "theta": self.theta, "c": self.c, "r": self.r, "dr": self.dr, "foils": self.foils,
            "B": int(m["B"]), "R": float(m["R"]), "Rhub": float(m.get("Rhub", 0.0)),
            "rpm": float(rpm), "v": float(v), "pitch": float(pitch), "yaw_angle": float(yaw),
            "turbine_type": int(m.get("turbine_type", 0)),
            "rho": rho, "kin_viscosity": nu, "speed_of_sound": self.speed_of_sound,
            "airfoils": self.airfoils, "return_print": [], "return_results": [],
        }

    def run_point(self, v, rpm, pitch=0.0, yaw=0.0, rho=None, nu=None):
        """One operating point. Returns the solver's results dict, or None when it gave no result."""
        args = self._args(v, rpm, pitch, yaw, rho or self.rho, nu or self.nu)
        if self.calculator is None:
            self.calculator = Calculator(args)
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                res = self.calculator.run_array(**args)
        except Exception as e:  # noqa: BLE001 - a diverging point is reported, not fatal
            return {"error": "%s: %s" % (type(e).__name__, e)}
        return res

    def conditions(self, exp, row):
        """Operating point of one experiment row: (v, rpm, pitch, yaw, rho, nu)."""
        d = {**self.case.meta.get("defaults", {}), **exp.defaults}
        for k in ("v_ms", "rpm", "pitch_deg", "yaw_deg", "rho", "kin_viscosity"):
            if k in row and pd.notna(row[k]):
                d[k] = float(row[k])
        R = float(self.case.meta["R"])
        v, rpm = d.get("v_ms"), d.get("rpm")
        if self.case.is_propeller:
            J = row.get("J", row.get("j"))
            if v is None and J is not None and rpm is not None:
                v = float(J) * rpm / 60 * 2 * R
            if rpm is None and J is not None and v is not None and J > 0:
                rpm = v / (float(J) * 2 * R) * 60
        else:
            tsr = row.get("tsr")
            if rpm is None and tsr is not None and v is not None:
                rpm = float(tsr) * v / R * 60 / (2 * math.pi)
            if v is None and tsr is not None and rpm is not None:
                v = rpm * 2 * math.pi / 60 * R / float(tsr)
        if v is None or rpm is None:
            raise ValueError("%s: cannot determine v and rpm for row %s" % (exp.file, dict(row)))
        return (v, rpm, d.get("pitch_deg", 0.0), d.get("yaw_deg", 0.0), d.get("rho", self.rho),
                d.get("kin_viscosity", self.nu))


def bending_moment_about(res, r_ref):
    """Flap (out-of-plane) and edge (in-plane) moment of one blade about radius r_ref, from the section forces."""
    r = np.asarray(res["r"])
    fn = np.asarray(res["dFn"])
    ft = np.asarray(res["dFt"])
    out = r > r_ref
    return float(np.sum(fn[out] * (r[out] - r_ref))), float(np.sum(ft[out] * (r[out] - r_ref)))


def summarize(res, case, r_ref=None):
    """Scalar results of one point, in the units the experiment files use."""
    if res is None or "error" in res:
        return {"converged": False, "error": None if res is None else res["error"]}
    out = {"converged": True, "cp": res["cp"], "ct": res["ct"], "cq": res["cq"], "power_w": res["power"],
           "thrust_n": res["thrust"], "torque_nm": res["Q"], "tsr_calc": res["TSR"], "J_calc": res["J"],
           "eta": res["eff"] if case.is_propeller else np.nan,
           "stall_fraction": res["blade_stall_percentage"]}
    if r_ref is not None:
        out["root_flap_moment_nm"], out["root_edge_moment_nm"] = bending_moment_about(res, r_ref)
    return out

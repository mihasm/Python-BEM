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

"""Checks that every validation case in validation/cases is complete and that the solver runs on it."""

import os
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "validation"))

import bem_validation as bv  # noqa: E402

CASES = bv.list_cases()
AERO_KINDS = ("performance", "power_curve", "spanwise", "loads")


def aero_cases():
    out = []
    for p in CASES:
        try:
            c = bv.load_case(p)
        except Exception:  # noqa: BLE001 - a case that does not load fails test_case_is_consistent
            continue
        if any(e.kind in AERO_KINDS for e in c.experiments):
            out.append(p)
    return out


@pytest.mark.parametrize("path", CASES, ids=[p.name for p in CASES])
def test_case_is_consistent(path):
    case = bv.load_case(path)
    for key in ("B", "R", "turbine_type"):
        assert key in case.meta, "case.json lacks %s" % key
    assert case.experiments, "no experiment file found"
    for e in case.meta["experiments"]:
        assert (path / e["file"]).exists(), "missing %s" % e["file"]
    if not any(e.kind in AERO_KINDS for e in case.experiments):
        return
    g = case.geometry
    R = float(case.meta["R"])
    assert np.all(np.diff(g["r_m"]) > 0), "radius must increase"
    assert g["r_m"].iloc[-1] <= R * 1.001
    assert np.all(g["chord_m"] > 0) and g["chord_m"].max() < R
    for f in set(g["airfoil"]) - {"transition"}:
        assert f in case.airfoils, "geometry names airfoil %s which case.json does not define" % f
    airfoils = bv.build_airfoils(case)
    for name, a in airfoils.items():
        assert len(a["gathered_curves"]) >= 3, name


@pytest.mark.parametrize("path", aero_cases(), ids=[p.name for p in aero_cases()])
def test_solver_runs_on_first_point(path):
    os.chdir(ROOT)
    case = bv.load_case(path)
    runner = bv.CaseRunner(case, "prop_prandtl" if case.is_propeller else "buhl_prandtl", n_sections=12)
    exp = next(e for e in case.experiments if e.kind in AERO_KINDS)
    # a well-loaded point: the solver drops propeller points with negative thrust by design
    row = exp.data.loc[exp.data["ct"].idxmax()] if "ct" in exp.data else exp.data.iloc[len(exp.data) // 2]
    v, rpm, pitch, yaw, rho, nu = runner.conditions(exp, row)
    res = runner.run_point(v, rpm, pitch, yaw, rho, nu)
    assert res is not None and "error" not in res, res
    assert np.isfinite(res["cp"]) and np.isfinite(res["ct"])

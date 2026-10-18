"""Runs the program's Calculator.get_crossection_data on every exact row of experiment/section_properties.csv (theta 0)
and prints the relative errors of A, Ix, Iy (audit aid; the comparison harness proper is elsewhere).
Run from anywhere: python3 validation/cases/structural_analytical/extraction/check_program_sections.py"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.dirname(HERE)
sys.path.insert(0, os.path.abspath(os.path.join(CASE, "..", "..", "..")))
import calculation  # noqa: E402

ref = pd.read_csv(os.path.join(CASE, "experiment", "section_properties.csv"), comment="#")
ref = ref[~ref.method.str.startswith("approx") & ~ref.method.str.startswith("tabulated") & (ref.theta_deg == 0)]
calc = calculation.Calculator.__new__(calculation.Calculator)
for _, r in ref.iterrows():
    d = np.loadtxt(os.path.join(CASE, r.airfoil_file), skiprows=1)
    o = calc.get_crossection_data(r.c_m, 0.0, d[:, 0], d[:, 1], 1 if r.construction == "hollow" else 0, r.wall_m)
    print("%-34s A %+.4f  Ix %+.4f  Iy %+.4f  (tolerance %.0e)" % (
        r.case_id, o["A"] / r.A_m2 - 1, o["Ix"] / r.Ix_m4 - 1, o["Iy"] / r.Iy_m4 - 1, r.tolerance_rel))

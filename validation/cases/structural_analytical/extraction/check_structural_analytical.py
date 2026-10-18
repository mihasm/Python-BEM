"""Sanity plot of the structural_analytical reference data -> extraction/check_structural_analytical.png"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.dirname(HERE)
sec = pd.read_csv(os.path.join(CASE, "experiment", "section_properties.csv"), comment="#")
beam = pd.read_csv(os.path.join(CASE, "experiment", "beam_closed_form.csv"), comment="#")
cen = pd.read_csv(os.path.join(CASE, "experiment", "centrifugal_closed_form.csv"), comment="#")

fig, ax = plt.subplots(2, 3, figsize=(15, 8.5))
sym = sec[sec.case_id.str.match(r"naca00\d\d_solid$")]
cam = sec[sec.case_id.str.match(r"naca[2-6]4\d\d_solid$")]
a = ax[0, 0]
a.plot(sym.t_c, sym.A_m2, "o-", label="NACA 00xx exact A/c^2")
a.plot(sym.t_c, sym.xfoil_A_m2, "x", label="XFOIL BEND")
a.plot(sym.t_c, 0.60 * sym.t_c, "--", label="Drela 0.60 t/c")
a.plot(cam.t_c, cam.A_m2, "s", mfc="none", label="cambered (2412..6412)")
a.set_xlabel("t/c")
a.set_ylabel("A / c^2")
a.legend(fontsize=8)
a.set_title("solid NACA 4-digit: area")
a = ax[0, 1]
a.plot(sym.t_c, sym.Ix_m4 / sym.t_c ** 3, "o-", label="exact Ix/(c^4 tau^3)")
a.plot(sym.t_c, sym.xfoil_Ix_m4 / sym.t_c ** 3, "x", label="XFOIL")
a.axhline(0.036, ls="--", color="C2", label="Drela K_I = 0.036")
a.plot(cam.t_c, cam.K_I, "s", mfc="none", label="cambered K_I")
a.set_ylim(0.03, 0.045)
a.set_xlabel("t/c")
a.legend(fontsize=8)
a.set_title("flapwise inertia coefficient")
a = ax[0, 2]
for code, m in (("0012", "o"), ("4412", "s"), ("0018", "^")):
    h = sec[sec.case_id.str.match(r"naca%s_hollow_w[0-9.]+$" % code)]
    a.plot(h.wall_m * 1000, h.Ix_m4 / sec[sec.case_id == "naca%s_solid" % code].Ix_m4.iloc[0], m + "-",
           label="NACA %s exact" % code)
    a.plot(h.wall_m * 1000, h.thinwall_Ix_m4 / sec[sec.case_id == "naca%s_solid" % code].Ix_m4.iloc[0], m + ":",
           mfc="none", label="thin-wall approx")
a.set_xlabel("wall [mm] (chord 1 m)")
a.set_ylabel("Ix_hollow / Ix_solid")
a.legend(fontsize=7)
a.set_title("hollow sections")
a = ax[1, 0]
for k, g in beam.groupby("case_id"):
    a.plot(g.r_m, g.moment_flap_nm, label=k)
a.set_xlabel("r [m]")
a.set_ylabel("M [N m]")
a.legend(fontsize=8)
a.set_title("cantilever bending moment")
a = ax[1, 1]
for k, g in beam.groupby("case_id"):
    a.plot(g.r_m, g.deflection_flap_m * 1000, label=k)
a.set_xlabel("r [m]")
a.set_ylabel("w [mm]")
a.legend(fontsize=8)
a.set_title("deflection (E=70 GPa, NACA0012 c=0.1)")
a = ax[1, 2]
for k, g in cen.groupby("case_id"):
    a.plot(g.r_m, g.stress_centrifugal_pa / 1e6, label=k)
r = np.linspace(0.1, 1, 50)
a.plot(r, 2700 * (600 * 2 * np.pi / 60) ** 2 * (1 - r ** 2) / 2 / 1e6, "k:", label="rho w^2 (R^2-r^2)/2")
a.set_xlabel("r [m]")
a.set_ylabel("sigma_c [MPa]")
a.legend(fontsize=8)
a.set_title("centrifugal stress, 600 rpm")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "check_structural_analytical.png"), dpi=100)

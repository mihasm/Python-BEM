"""Sanity checks + overview plot of the nrel_phase_vi case -> check_nrel_phase_vi.png"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.dirname(HERE)
c = json.load(open(os.path.join(CASE, "case.json")))
R = c["R"]
g = pd.read_csv(os.path.join(CASE, c["geometry_file"]), comment="#")
pc = pd.read_csv(os.path.join(CASE, "experiment", "power_curve.csv"), comment="#")
ld = pd.read_csv(os.path.join(CASE, "experiment", "loads.csv"), comment="#")
sp = pd.read_csv(os.path.join(CASE, "experiment", "spanwise_cn_ct.csv"), comment="#")
dut = pd.read_csv(os.path.join(CASE, "airfoils", "s809_polars.csv"), comment="#")
osu = pd.read_csv(os.path.join(CASE, "airfoils", "s809_osu_polars.csv"), comment="#")
xy = np.loadtxt(os.path.join(CASE, "airfoils", "s809.dat"), skiprows=1)

# --- checks
assert abs(g.r_m.iloc[-1] - R) < 1e-9, "radius must end at R"
assert g.r_m.is_monotonic_increasing
assert set(g.airfoil) <= set(c["airfoils"]) | {"transition"}
s = g[g.airfoil == "s809"]
assert (np.diff(s.twist_deg) <= 0).all(), "twist must decrease outboard on the S809 part"
assert abs(s.twist_deg.iloc[-1]) < 1e-9
assert (g.chord_m > 0.1).all() and (g.chord_m < 0.8).all()
assert (pc.cp < 16 / 27).all() and (pc.cp > 0).all()
assert np.allclose(pc.power_w, pc.torque_nm * 72 * np.pi / 30, rtol=1e-3)
for col, ref in (("cn_chord", "cn_simms2001"), ("ct_chord", "ct_simms2001")):
    d = (sp[col] - sp[ref]).abs().dropna()
    assert d.max() < 0.01, (col, d.max())
assert (abs(ld.flap_moment_gauge_axes_nm - ld.root_flap_lindenburg2003_nm) < 10).all()
for f in [e["file"] for e in c["experiments"]]:
    assert os.path.exists(os.path.join(CASE, f)), f
# thrust consistency: flap moment / lever arm ~ thrust per blade (rough, centre of pressure ~0.7 R)
T_est = 2 * ld.root_flap_moment_nm / (0.7 * R - 0.432)
ct_est = T_est / (0.5 * ld.rho * ld.v_ms ** 2 * np.pi * R ** 2)
print("rough rotor CT from flap moment (c.o.p. 0.7R):", dict(zip(ld.v_ms, ct_est.round(3))))
print("geometry rows:", len(g), " Cp max:", pc.cp.max().round(3), " t/c(s809 coords) =",
      round(np.max(np.interp(np.linspace(0, 1, 400), xy[:32][::-1, 0], xy[:32][::-1, 1])
                   - np.interp(np.linspace(0, 1, 400), xy[31:, 0], xy[31:, 1])), 4))
print("all checks passed")

# --- plot
fig, ax = plt.subplots(3, 3, figsize=(16, 13))
a = ax[0, 0]
a.plot(g.r_m, g.chord_m, "o-")
for _, r in g.iterrows():
    a.annotate(r.airfoil[0], (r.r_m, r.chord_m), fontsize=7, xytext=(2, 4), textcoords="offset points")
a.set(xlabel="r [m]", ylabel="chord [m]", title="Chord (c=circle, t=transition, s=S809)")
a = ax[0, 1]
a.plot(g.r_m, g.twist_deg, "o-", label="twist_deg (rel. to tip)")
a.plot(g.r_m, g.twist_deg + 3.0, "--", label="twist + pitch 3 deg (sequence S)")
a.axhline(0, color="k", lw=0.5)
a.set(xlabel="r [m]", ylabel="deg", title="Twist")
a.legend()
a = ax[0, 2]
a.plot(xy[:, 0], xy[:, 1], ".-")
a.set_aspect("equal")
a.set(title="S809 coordinates (Hand 2001 Table A-2)")
a = ax[1, 0]
for re_, d in dut.groupby("re"):
    a.plot(d.alpha_deg, d.cl, "o-", ms=3, label=f"DUT Re={re_/1e6:g}e6")
for re_, d in osu.groupby("re"):
    a.plot(d.alpha_deg, d.cl, "x--", ms=3, lw=0.7, label=f"OSU Re={re_/1e6:g}e6")
a.set(xlabel="alpha [deg]", ylabel="cl", title="S809 measured cl")
a.legend(fontsize=6)
a = ax[1, 1]
for re_, d in dut.groupby("re"):
    a.plot(d.cd, d.cl, "o-", ms=3, label=f"DUT {re_/1e6:g}e6")
for re_, d in osu[osu.alpha_deg < 16].groupby("re"):
    a.plot(d.cd, d.cl, "x--", ms=3, lw=0.7, label=f"OSU {re_/1e6:g}e6")
a.set(xlabel="cd", ylabel="cl", title="S809 drag polar", xlim=(0, 0.15))
a = ax[1, 2]
m = pc.campaign.notna()
a.errorbar(pc.v_ms[m], pc.torque_nm[m], yerr=pc.torque_nm_std[m], fmt="ks", capsize=3, label="Simms 2001 Fig.8 (+-1 std, Sorensen)")
a.plot(pc.v_ms[~m], pc.torque_nm[~m], "r^", label="Lindenburg 2003 Fig.6.2 (extra speeds)")
a.set(xlabel="v [m/s]", ylabel="LSS torque [N m]", title="Sequence S torque")
a2 = a.twinx()
a2.plot(pc.v_ms, pc.cp, "g:", label="Cp")
a2.set_ylabel("Cp", color="g")
a.legend(fontsize=7)
a = ax[2, 0]
a.errorbar(ld.v_ms, ld.flap_moment_gauge_axes_nm, yerr=ld.flap_moment_gauge_axes_nm_std, fmt="bo-", capsize=3, label="root flap (r=0.432 m)")
a.errorbar(ld.v_ms, ld.edge_moment_gauge_axes_nm, yerr=ld.edge_moment_gauge_axes_nm_std, fmt="ms-", capsize=3, label="root edge (+-std incl. gravity)")
a.set(xlabel="v [m/s]", ylabel="moment [N m]", title="Root bending moments (tip-chord axes)")
a.legend(fontsize=7)
a = ax[2, 1]
for case, d in sp.groupby("campaign"):
    a.errorbar(d.r_R, d.cn_chord, yerr=d.cn_chord_std, fmt="o-", capsize=2, label=case)
a.set(xlabel="r/R", ylabel="cn", title="Normal force coefficient")
a.legend(fontsize=7)
a = ax[2, 2]
for case, d in sp.groupby("campaign"):
    a.errorbar(d.r_R, d.ct_chord, yerr=d.ct_chord_std, fmt="o-", capsize=2, label=case)
a.axhline(0, color="k", lw=0.5)
a.set(xlabel="r/R", ylabel="ct (chordwise, +LE)", title="Tangential force coefficient")
fig.suptitle("nrel_phase_vi - check of case data")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "check_nrel_phase_vi.png"), dpi=90)
print("plot written")

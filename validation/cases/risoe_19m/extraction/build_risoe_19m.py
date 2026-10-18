"""Build the Risoe 19 m / 100 kW stall-regulated test turbine case (IEA Annex XIV/XVIII field database).

All numbers are TABULATED in the sources (nothing digitised):
  raw/PV_r            - measured power curve, electrical + mechanical (IEA XVIII database file 'PV_r' = Fig. G.5)
  raw/cna*_rot.r,
  raw/cta*_rot.r      - rotating-blade cn/ct vs alpha at 37/68/98 % span (database files, Figs. G.6/G.7)
  Table G.2 (planform), G.3 (2D lift), G.4 (2D drag) of ECN-C--02-016 (Final report IEA Annex XVIII), typed below
  raw/naca634221-il.dat - NACA 63(4)-221 coordinates (airfoiltools.com) used as the shape of the 63-2xx family
Run: python3 build_risoe_19m.py
"""
import os
import subprocess
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.dirname(HERE)
RAW = os.path.join(HERE, "raw")
XF_DIR = os.path.join(HERE, "xfoil")
XFOIL = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(CASE))), "xfoil_executables", "xfoil")


def xfoil_negative(name, re_=1.5e6, ncrit=9):
    """XFOIL polar 0 .. -12 deg of the generated section (negative side is not in the report tables).
    Cached in extraction/xfoil/<name>.pol; the log is kept next to it."""
    os.makedirs(XF_DIR, exist_ok=True)
    pol = os.path.join(XF_DIR, name + ".pol")
    if not os.path.exists(pol):
        # XFOIL truncates long file names -> run inside XF_DIR with short relative names
        import shutil
        shutil.copy(os.path.join(CASE, "airfoils", name + ".dat"), os.path.join(XF_DIR, name + ".dat"))
        cmds = f"PLOP\nG F\n\nLOAD {name}.dat\nPANE\nOPER\nVPAR\nN {ncrit}\n\nVISC {re_:.0f}\nITER 200\nPACC\n{name}.pol\n\nASEQ 0 -12 -1\n\nQUIT\n"
        with open(os.path.join(XF_DIR, name + ".log"), "w") as log:
            subprocess.run(["perl", "-e", "alarm 150; exec @ARGV", XFOIL], input=cmds.encode(), stdout=log, stderr=log, cwd=XF_DIR)
    rows = []
    for line in open(pol).read().splitlines()[12:]:
        v = line.split()
        if len(v) >= 3:
            rows.append([float(v[0]), float(v[1]), float(v[2])])
    a = np.array(sorted({r[0]: r for r in rows}.values()))
    return a
REP = ("Schepers J.G. et al. (2002) Final report of IEA Annex XVIII 'Enhanced Field Rotor Aerodynamics Database', "
       "ECN-C--02-016, Appendix G (updated Aug 2003); data: Zenodo doi:10.5281/zenodo.10619556 (CC-BY-4.0)")

R = 9.5
RPM = 47.5
PITCH = 1.7           # mean of the measured blade tip angles 1.8, 1.8, 1.5 deg (App. G.1)
RHO = 1.25            # see notes: not stated for PV_r; Jan-1993 campaigns had 1.266-1.279

# ---------------------------------------------------------------- Table G.2 (typed from the report, p. 294/297)
tabG2 = np.array([
    [2.700, 1.090, 15.000, 24.600],
    [3.550, 1.005, 9.500, 20.700],
    [4.400, 0.925, 6.100, 18.700],
    [5.250, 0.845, 3.900, 17.600],
    [6.100, 0.765, 2.400, 16.600],
    [6.950, 0.685, 1.500, 15.600],
    [7.800, 0.605, 0.900, 14.600],
    [8.650, 0.525, 0.400, 13.600],
    [9.500, 0.445, 0.000, 12.600],
])
# ---------------------------------------------------------------- Tables G.3 / G.4 (typed), t/c = 12, 15, 18, 21, 25 %
TC = np.array([12.0, 15.0, 18.0, 21.0, 25.0])
alpha = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50, 60, 70, 80, 90], float)
CL = np.array([
    [.16, .16, .16, .16, .16],
    [.38, .38, .38, .38, .38],
    [.60, .60, .60, .60, .60],
    [.82, .82, .82, .82, .81],
    [1.04, 1.03, 1.01, 1.00, .98],
    [1.26, 1.23, 1.19, 1.14, 1.08],
    [1.43, 1.36, 1.29, 1.20, 1.12],
    [1.50, 1.42, 1.32, 1.21, 1.12],
    [1.34, 1.30, 1.25, 1.18, 1.09],
    [.80, 1.05, 1.10, 1.12, 1.06],
    [.70, .85, 1.00, 1.05, 1.04],
    [.92, .95, .97, .99, 1.01],
    [.99, 1.00, 1.01, 1.02, 1.03],
    [.99, 1.00, 1.01, 1.02, 1.03],
    [.92, .92, .92, .92, .92],
    [.75, .75, .75, .75, .75],
    [.51, .51, .51, .51, .51],
    [.24, .24, .24, .24, .24],
    [-.02, -.02, -.02, -.02, -.02],
])
CD = np.array([
    [.0048, .0052, .0058, .0064, .0073],
    [.0048, .0052, .0058, .0063, .0072],
    [.0056, .0058, .0061, .0065, .0073],
    [.0080, .0080, .0080, .0079, .0078],
    [.0110, .0110, .0108, .0105, .0102],
    [.0158, .0156, .0154, .0152, .0150],
    [.0240, .0230, .0225, .0218, .0210],
    [.0400, .0380, .0360, .0340, .0320],
    [.0700, .0700, .0700, .0700, .0700],
    [.1500, .1500, .1500, .1500, .1500],
    [.2150, .2200, .2250, .2300, .2350],
    [.3700, .3600, .3500, .3400, .3300],
    [.4700, .4700, .4600, .4500, .4500],
    [.6600, .6600, .6600, .6600, .6600],
    [.8500, .8500, .8500, .8500, .8500],
    [1.0100, 1.0100, 1.0100, 1.0100, 1.0100],
    [1.1500, 1.1500, 1.1500, 1.1500, 1.1500],
    [1.2500, 1.2500, 1.2500, 1.2500, 1.2500],
    [1.3000, 1.3000, 1.3000, 1.3000, 1.3000],
])
# keep the typed tables in the audit trail
tab = pd.DataFrame({"alpha_deg": np.repeat(alpha, 5), "tc_pct": np.tile(TC, len(alpha)),
                    "cl": CL.ravel(), "cd": CD.ravel()})
with open(os.path.join(HERE, "tables_G3_G4_typed.csv"), "w") as f:
    f.write(f"# Tables G.3 (2D lift) and G.4 (2D drag) for the Risoe test blade, typed from {REP}, p. 294-295\n")
    f.write("# 'basically 2D data from Abbott & von Doenhoff but with corrections for the actual Reynolds number' (App. G.2.2)\n")
    tab.to_csv(f, index=False)

def name_of(r):
    return f"risoe_r{int(round(r * 1000)):04d}"

# ---------------------------------------------------------------- geometry
geo = pd.DataFrame({"r_m": tabG2[:, 0], "chord_m": tabG2[:, 1], "twist_deg": tabG2[:, 2],
                    "airfoil": [name_of(r) for r in tabG2[:, 0]], "thickness_ratio": np.round(tabG2[:, 3] / 100, 4)})
with open(os.path.join(CASE, "geometry.csv"), "w") as f:
    f.write("# Risoe 19 m test turbine (LM 8.2 m blades + 1 m root extensions, 3 blades, NACA 63-2xx) - planform\n")
    f.write(f"# Source: Table G.2 'Planform data for the 19 m rotor with LM 8.2 m blades', {REP}, p. 294 (tabulated).\n")
    f.write("# Twist is given relative to the tip station (twist 0 at r = 9.5 m), positive towards feather = solver convention;\n")
    f.write("#   the measured tip angles (1.8, 1.8, 1.5 deg, App. G.1) are the global pitch -> pitch_deg = 1.7 in experiment files.\n")
    f.write("# The profiled blade starts at r = 2.7 m (blade 8.2 m incl. 1.4 m non-profiled root + 1 m extension); inner part not given.\n")
    f.write("# Cross-check: Table E.1 of the same report (same LM 8.2 m blade on the IC/RAL 17 m rotor, distance from blade root)\n")
    f.write("#   agrees except chord 1.000 vs 1.005 m (2nd station) and 0.770 vs 0.765 m (5th station).\n")
    f.write("# airfoil: one polar set per station, linearly interpolated in t/c from Tables G.3/G.4 (see airfoils/*_polars.csv).\n")
    geo.to_csv(f, index=False)

# ---------------------------------------------------------------- airfoils: coordinates (63-2xx shape scaled in thickness)
def read_xy(fn):
    pts = []
    for line in open(fn):
        p = line.split()
        if len(p) == 2:
            try:
                pts.append([float(p[0]), float(p[1])])
            except ValueError:
                pass
    return np.array(pts)

base = read_xy(os.path.join(RAW, "naca634221-il.dat"))
ile = int(np.argmin(base[:, 0]))
up = base[:ile + 1][::-1]      # LE -> TE
lo = base[ile:]                # LE -> TE
x = np.unique(np.concatenate([up[:, 0], lo[:, 0]]))
yu = np.interp(x, up[:, 0], up[:, 1]); yl = np.interp(x, lo[:, 0], lo[:, 1])
camber = 0.5 * (yu + yl); thick = yu - yl
t_base = thick.max()

def shape_for(tc):
    s = tc / t_base
    yu2 = camber + 0.5 * thick * s
    yl2 = camber - 0.5 * thick * s
    return np.vstack([np.column_stack([x[::-1], yu2[::-1]]), np.column_stack([x[1:], yl2[1:]])])

xf_notes = {}
for r, c, tw, tpct in tabG2:
    nm = name_of(r)
    xy = shape_for(tpct / 100)
    with open(os.path.join(CASE, "airfoils", nm + ".dat"), "w") as f:
        f.write(f"{nm} NACA 63-2xx t/c={tpct:.1f}% (GENERATED: NACA 63(4)-221 camber + thickness scaled; exact LM section not published)\n")
        for a, b in xy:
            f.write(f"{a:9.6f} {b:9.6f}\n")

    # polars interpolated in t/c, stored 0..30 deg (the published range with flow physics; 40-90 deg rows are a
    # flat-plate-like extrapolation the runner replaces with its own); no negative alpha published
    m = alpha <= 30
    cl = np.array([np.interp(tpct, TC, row) for row in CL[m]])
    cd = np.array([np.interp(tpct, TC, row) for row in CD[m]])
    p = pd.DataFrame({"re": 1.5e6, "alpha_deg": alpha[m], "cl": np.round(cl, 4), "cd": np.round(cd, 5),
                      "source": f"Risoe Tables G.3/G.4 (Abbott & von Doenhoff based, Re-corrected) interpolated to t/c={tpct}%"})
    # negative side from XFOIL on the generated section, made continuous at alpha = 0 with the table:
    # cl shifted by a constant (|shift| < 0.01), cd scaled by table/XFOIL at 0 deg
    xf = xfoil_negative(nm)
    x0 = xf[np.argmin(np.abs(xf[:, 0]))]
    assert abs(x0[0]) < 1e-6, nm
    dcl = cl[0] - x0[1]
    kcd = cd[0] / x0[2]
    neg = xf[xf[:, 0] < -0.5]
    pn = pd.DataFrame({"re": 1.5e6, "alpha_deg": neg[:, 0], "cl": np.round(neg[:, 1] + dcl, 4),
                       "cd": np.round(neg[:, 2] * kcd, 5),
                       "source": f"XFOIL 6.97, Ncrit=9, Re 1.5e6, generated t/c={tpct}% section; cl{dcl:+.4f}, cd*{kcd:.3f} for continuity at 0 deg"})
    p = pd.concat([pn, p]).sort_values("alpha_deg")
    xf_notes[nm] = (dcl, kcd)
    with open(os.path.join(CASE, "airfoils", nm + "_polars.csv"), "w") as f:
        f.write(f"# {nm}: 2D polars for t/c = {tpct} % - linear interpolation in t/c between the columns 12/15/18/21/25 %\n")
        f.write(f"#   of Tables G.3 (lift) and G.4 (drag), {REP}, p. 294-295 (tabulated; typed copy in extraction/tables_G3_G4_typed.csv).\n")
        f.write("# The tables are 'basically 2D data from Abbott & von Doenhoff with corrections for the actual Reynolds number'.\n")
        f.write("# Only alpha = 0..30 deg is published with flow physics (40-90 deg rows dropped); NO negative alpha published:\n")
        f.write(f"#   alpha -12..-1 deg are COMPUTED with XFOIL 6.97 (Ncrit 9, Re 1.5e6) on the generated section, cl shifted by\n")
        f.write(f"#   {xf_notes[nm][0]:+.4f} and cd scaled by {xf_notes[nm][1]:.3f} to join the table at 0 deg (see 'source' column;\n")
        f.write("#   raw XFOIL output in extraction/xfoil/). XFOIL gives cl(0) = 0.155-0.159 vs 0.16 in the table, i.e. consistent.\n")
        f.write("# re = 1.5e6 is NOMINAL: the report does not state the Re; at 47.5 rpm the blade runs at Re ~1.1e6 (root) .. 1.7e6.\n")
        p.to_csv(f, index=False)

# ---------------------------------------------------------------- experiment: power curve
pv = np.loadtxt(os.path.join(RAW, "PV_r"), skiprows=5)
omega = RPM * 2 * np.pi / 60
A = np.pi * R ** 2
df = pd.DataFrame({"v_ms": pv[:, 0], "power_w": pv[:, 2] * 1000, "power_elec_w": pv[:, 1] * 1000})
df["tsr"] = np.round(omega * R / df.v_ms, 4)
df["cp"] = np.round(df.power_w / (0.5 * RHO * df.v_ms ** 3 * A), 4)
df["torque_nm"] = np.round(df.power_w / omega, 1)
with open(os.path.join(CASE, "experiment", "power_curve.csv"), "w") as f:
    f.write("# Risoe 19 m test turbine - measured power curve at the 47.5 rpm (95 kW) setting, fixed pitch (stall regulated)\n")
    f.write(f"# Source: IEA Annex XVIII database file 'PV_r' (Risoe data, plotted as Fig. G.5 of {REP}). Tabulated, copied.\n")
    f.write("# power_w = 'P_mechanical' column (shaft power), power_elec_w = 'P_elec' column. Field measurement, binned;\n")
    f.write("#   wind speed from the met mast (hub height 29.3 m); number of samples per bin and density normalisation NOT stated.\n")
    f.write(f"# tsr, cp, torque_nm DERIVED here: tsr = Omega*R/v (R = 9.5 m, 47.5 rpm), cp = P_mech/(0.5*rho*v^3*pi*R^2) with rho = {RHO}\n")
    f.write("#   (assumption, see case.json notes), torque = P_mech/Omega.\n")
    f.write("# Extraction uncertainty: none (tabulated). Measurement uncertainty not published; field data (turbulence, tilt 5 deg).\n")
    df.to_csv(f, index=False)

# ---------------------------------------------------------------- experiment: rotating cn/ct vs alpha
def read_coef(fn):
    rows = []
    for line in open(fn):
        p = line.split()
        try:
            rows.append([float(v) for v in p])
        except ValueError:
            continue
    return rows

stations = {"37": 3.555, "68": 6.475, "98": 9.275}   # segment radii from NY_LOG (R1, R2, R3)
for st, r in stations.items():
    cn = read_coef(os.path.join(RAW, f"cna{ {'37': 1, '68': 2, '98': 3}[st]}_rot.r"))
    cnd = pd.DataFrame([row[:3] if len(row) >= 3 else row[:2] + [np.nan] for row in cn],
                       columns=["alpha_deg", "cn_chord", "cn_chord_std"])
    out = cnd.copy()
    if st != "37":
        ct = read_coef(os.path.join(RAW, f"cta{ {'68': 2, '98': 3}[st]}_rot.r"))
        ctd = pd.DataFrame([row[:3] for row in ct], columns=["alpha_deg", "ct_chord", "ct_chord_std"])
        out = cnd.merge(ctd, on="alpha_deg", how="outer").sort_values("alpha_deg")
    out.insert(0, "r_m", r)
    out.insert(1, "r_R", round(r / R, 3))
    fn = os.path.join(CASE, "experiment", f"rotating_coefficients_{st}pct.csv")
    with open(fn, "w") as f:
        f.write(f"# Risoe 19 m - rotating-blade chord-normal (cn_chord) and chordwise (ct_chord) force coefficients vs alpha, {st} % span (r = {r} m)\n")
        f.write(f"# Source: IEA Annex XVIII database files cna*_rot_r / cta*_rot_r (Figs. G.6/G.7 of {REP}). Tabulated, copied;\n")
        f.write("#   *_std = standard deviation within the alpha bin (as published). Forces measured directly on blade segments\n")
        f.write("#   (segment radii 3.555 / 6.475 / 9.275 m), not from pressures; 47.5 rpm, normal operation, January 1993.\n")
        if st != "37":
            f.write("# CAUTION: alpha is the five-hole-pitot inflow angle at r = 6.83 m minus a constant 4.7 deg chosen so that the measured\n")
            f.write("#   power curve vs alpha matched a BEM computation (NY_LOG: 'subtracted a constant 4.7 deg'); not an independent alpha.\n")
        else:
            f.write("# CAUTION: alpha from the pitot at r = 3.92 m; Risoe derives alpha with a constant BEM-matched offset (report sect. 6.2);\n")
            f.write("#   the offset used for this station is not stated.\n")
        if st == "37":
            f.write("# 37 %: updated Aug 2003 with the pitot at position 2 (r = 3.92 m); values 'read from figure' by Risoe; no ct, no std.\n")
        f.write("# cn_chord = force normal to the chord, ct_chord = force along the chord (positive towards the leading edge), both divided\n#   by the local dynamic pressure and chord; reference velocity = measured inflow velocity scaled to the segment radius.\n")
        out.to_csv(f, index=False)

# ---------------------------------------------------------------- check plot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(2, 3, figsize=(15, 8.5))
ax[0, 0].plot(geo.r_m, geo.chord_m, "o-"); ax[0, 0].set_ylabel("chord [m]")
a = ax[0, 0].twinx(); a.plot(geo.r_m, geo.twist_deg, "rs-"); a.set_ylabel("twist [deg] (red)")
ax[0, 0].set_xlabel("r [m]"); ax[0, 0].grid()
for r, c, tw, tpct in tabG2[::2]:
    xy = read_xy(os.path.join(CASE, "airfoils", name_of(r) + ".dat"))
    ax[0, 1].plot(xy[:, 0], xy[:, 1], label=f"{name_of(r)} t/c {tpct}%")
ax[0, 1].set_aspect("equal"); ax[0, 1].legend(fontsize=7)
for r, c, tw, tpct in tabG2[::2]:
    p = pd.read_csv(os.path.join(CASE, "airfoils", name_of(r) + "_polars.csv"), comment="#")
    ax[0, 2].plot(p.alpha_deg, p.cl, "o-", ms=3, label=f"t/c {tpct}%")
    ax[0, 2].plot(p.alpha_deg, p.cd * 5, "--", lw=0.8)
ax[0, 2].set_xlabel("alpha"); ax[0, 2].set_ylabel("cl (solid), 5*cd (dashed)"); ax[0, 2].legend(fontsize=7); ax[0, 2].grid()
ax[1, 0].plot(df.v_ms, df.power_w / 1e3, "o-", label="mechanical"); ax[1, 0].plot(df.v_ms, df.power_elec_w / 1e3, "s-", label="electrical")
ax[1, 0].set_xlabel("v [m/s]"); ax[1, 0].set_ylabel("P [kW]"); ax[1, 0].legend(); ax[1, 0].grid()
ax[1, 1].plot(df.tsr, df.cp, "o-"); ax[1, 1].axhline(16 / 27, color="gray", ls=":")
ax[1, 1].set_xlabel("tsr"); ax[1, 1].set_ylabel(f"cp (rho={RHO})"); ax[1, 1].grid()
for st in stations:
    d = pd.read_csv(os.path.join(CASE, "experiment", f"rotating_coefficients_{st}pct.csv"), comment="#")
    ax[1, 2].plot(d.alpha_deg, d.cn_chord, "o-", ms=3, label=f"cn {st}%")
p = pd.read_csv(os.path.join(CASE, "airfoils", name_of(6.95) + "_polars.csv"), comment="#")
ax[1, 2].plot(p.alpha_deg, p.cl, "k--", label="2D cl t/c 15.6%")
ax[1, 2].set_xlabel("alpha"); ax[1, 2].legend(fontsize=7); ax[1, 2].grid()
fig.suptitle("Risoe 19 m test turbine - check"); fig.tight_layout()
fig.savefig(os.path.join(HERE, "check_risoe_19m.png"), dpi=110)

assert (np.diff(geo.twist_deg) < 0).all() and geo.r_m.iloc[-1] == R
assert (df.cp < 16 / 27).all()
print(df.to_string())

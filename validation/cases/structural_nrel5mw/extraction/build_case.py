"""
Builds the case `structural_nrel5mw` from
  - Jonkman et al. 2009, NREL/TP-500-38060: Table 2-1 (distributed blade structural properties), Table 2-2 (integrated
    blade properties), Table 3-1 (blade aerodynamic properties) - text of the tables as extracted by pdftotext is kept in
    extraction/jonkman2009_table2-1_pdftotext.txt and ..._table3-1_pdftotext.txt
  - OpenFAST r-test (glue-codes/openfast/5MW_Baseline, github.com/OpenFAST/r-test, main branch, fetched 2026-10-18):
    ElastoDyn and BeamDyn blade files (second source for Table 2-1), airfoil coordinates and polars (files in
    sources/openfast_r-test/).
Run: python3 validation/cases/structural_nrel5mw/extraction/build_case.py  (numpy, pandas, matplotlib)
"""
import os
import re

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.dirname(HERE)
RT = os.path.join(CASE, "sources", "openfast_r-test")

# ---------------------------------------------------------------------------------------------------------------------
# Table 2-1 (primary source)
# ---------------------------------------------------------------------------------------------------------------------
cols21 = ["r_m", "blfract", "aerocent", "strctwist_deg", "mass_per_length_kg_m", "ei_flap_nm2", "ei_edge_nm2",
          "gj_nm2", "ea_n", "alpha", "flap_inertia_kg_m", "edge_inertia_kg_m", "precrv_m", "preswp_m", "flpcg_m",
          "edgcg_m", "flpea_m", "edgea_m"]
rows = []
with open(os.path.join(HERE, "jonkman2009_table2-1_pdftotext.txt")) as f:
    for line in f:
        tok = line.split()
        if len(tok) == 18 and re.match(r"^\d+\.\d+$", tok[0]):
            vals = [float(t.replace("E+6", "e6")) for t in tok]
            rows.append(vals)
t21 = pd.DataFrame(rows, columns=cols21)
assert len(t21) == 49, len(t21)
t21 = t21.drop(columns=["alpha", "precrv_m", "preswp_m", "flpcg_m", "flpea_m", "edgea_m"])  # all zero in the table

# second source: OpenFAST r-test ElastoDyn (BlFract, StrcTwst, BMassDen, FlpStff, EdgStff) and BeamDyn (6x6 matrices)
ed = []
with open(os.path.join(RT, "NRELOffshrBsline5MW_Blade.dat")) as f:
    lines = f.read().splitlines()
i0 = [k for k, l in enumerate(lines) if "DISTRIBUTED BLADE PROPERTIES" in l][0] + 3
for l in lines[i0:i0 + 49]:
    ed.append([float(v) for v in l.split()])
ed = np.array(ed)
bd_K, bd_M, bd_eta = [], [], []
with open(os.path.join(RT, "NRELOffshrBsline5MW_BeamDyn_Blade.dat")) as f:
    txt = f.read().splitlines()
k = [i for i, l in enumerate(txt) if "Distributed Properties" in l][0] + 1
while len(bd_eta) < 49:
    while not txt[k].strip():
        k += 1
    bd_eta.append(float(txt[k]))
    K = np.array([[float(v) for v in txt[k + 1 + j].split()] for j in range(6)])
    M = np.array([[float(v) for v in txt[k + 8 + j].split()] for j in range(6)])
    bd_K.append(K)
    bd_M.append(M)
    k += 14
bd_K, bd_M = np.array(bd_K), np.array(bd_M)
checks = {
    "blfract vs ElastoDyn": np.max(np.abs(t21.blfract - ed[:, 0])),
    "strctwist vs ElastoDyn": np.max(np.abs(t21.strctwist_deg - ed[:, 1])),
    "mass vs ElastoDyn (rel)": np.max(np.abs(t21.mass_per_length_kg_m / ed[:, 2] - 1)),
    "EI_flap vs ElastoDyn (rel)": np.max(np.abs(t21.ei_flap_nm2 / ed[:, 3] - 1)),
    "EI_edge vs ElastoDyn (rel)": np.max(np.abs(t21.ei_edge_nm2 / ed[:, 4] - 1)),
    "EA vs BeamDyn K33 (rel)": np.max(np.abs(t21.ea_n / bd_K[:, 2, 2] - 1)),
    "EI_edge vs BeamDyn K44 (rel)": np.max(np.abs(t21.ei_edge_nm2 / bd_K[:, 3, 3] - 1)),
    "EI_flap vs BeamDyn K55 (rel)": np.max(np.abs(t21.ei_flap_nm2 / bd_K[:, 4, 4] - 1)),
    "GJ vs BeamDyn K66 (rel)": np.max(np.abs(t21.gj_nm2 / bd_K[:, 5, 5] - 1)),
    "mass vs BeamDyn M11 (rel)": np.max(np.abs(t21.mass_per_length_kg_m / bd_M[:, 0, 0] - 1)),
    "edge inertia vs BeamDyn M44 (rel)": np.max(np.abs(t21.edge_inertia_kg_m / bd_M[:, 3, 3] - 1)),
    "flap inertia vs BeamDyn M55 (rel)": np.max(np.abs(t21.flap_inertia_kg_m / bd_M[:, 4, 4] - 1)),
}
for k_, v in checks.items():
    print("%-36s max diff %.3g" % (k_, v))

# integrated properties (Table 2-2) from the table with the 4.536 % mass scaling
s = (t21.r_m - 1.5).values
mden = t21.mass_per_length_kg_m.values * 1.04536
m_tot = np.trapezoid(mden, s)
m1 = np.trapezoid(mden * s, s)
m2 = np.trapezoid(mden * s ** 2, s)
print("integrated mass %.0f kg (Table 2-2: 17,740), first moment %.0f (363,231), second %.0f (11,776,047), CM %.3f m "
      "(20.475)" % (m_tot, m1, m2, m1 / m_tot))

hdr = """# NREL 5-MW reference blade - distributed structural properties (kind: structural; reference-model data, NOT a measurement)
# Source: Jonkman, Butterfield, Musial & Scott 2009, NREL/TP-500-38060, Table 2-1 (p. 5), tabulated (pdftotext of the
#   report, extraction/jonkman2009_table2-1_pdftotext.txt). Composite (GFRP) blade - the properties come from the DOWEC
#   study (Lindenburg/Kooijman 2003); no single E or rho applies.
# Cross-checked against OpenFAST r-test 5MW_Baseline ElastoDyn and BeamDyn blade files (sources/openfast_r-test): identical
#   (max relative differences printed by extraction/build_case.py; BeamDyn K33=EA, K44=EI_edge, K55=EI_flap, K66=GJ,
#   M44=edge inertia, M55=flap inertia).
# Columns: r_m from the rotor apex along the pitch axis (root at 1.5 m, tip at 63 m); blfract (r-1.5)/61.5;
#   aerocent (aerodynamic-centre chord fraction, FAST convention); strctwist_deg structural twist = aerodynamic twist;
#   mass_per_length_kg_m AS TABULATED (the model scales it by +4.536 % to the 17,740 kg blade mass - see
#   integrated_properties.csv); ei_flap_nm2 / ei_edge_nm2 about the principal structural axes oriented by strctwist;
#   gj_nm2 torsional stiffness; ea_n axial stiffness (estimated by Jonkman as 1e7 x average mass moment of inertia);
#   flap_inertia_kg_m / edge_inertia_kg_m mass moments of inertia per length; edgcg_m chordwise CM offset from the pitch
#   axis, positive toward the trailing edge. (Alpha, PrecrvRef, PreswpRef, FlpcgOf, FlpEAOf, EdgEAOf are zero in the table.)
"""
with open(os.path.join(CASE, "experiment", "distributed_properties.csv"), "w") as f:
    f.write(hdr)
    t21.to_csv(f, index=False, float_format="%.6g")

integ_rows = [
    dict(quantity="blade_length_m", value=61.5, trapezoid_table21=s[-1], note="w.r.t. root along preconed axis"),
    dict(quantity="mass_scaling_factor", value=1.04536, note="applied to mass_per_length_kg_m"),
    dict(quantity="blade_mass_kg", value=17740.0, trapezoid_table21=m_tot),
    dict(quantity="first_mass_moment_root_kgm", value=363231.0, trapezoid_table21=m1, note="w.r.t. blade root"),
    dict(quantity="second_mass_moment_root_kgm2", value=11776047.0, trapezoid_table21=m2, note="w.r.t. blade root"),
    dict(quantity="cm_from_root_m", value=20.475, trapezoid_table21=m1 / m_tot, note="along preconed axis"),
]
INTEG_HDR = ("# NREL 5-MW blade - integrated structural properties (kind: structural)\n"
             "# value: Jonkman et al. 2009, Table 2-2 (p. 6), tabulated.\n"
             "# trapezoid_table21: trapezoidal integral of Table 2-1 x 1.04536 (differs by -0.7 %).\n"
             "# midpoint_17_elements: Table 2-1 x 1.04536 linearly interpolated to the 17 element centres of Table 3-1\n"
             "#   (geometry.csv r_m) times the element lengths dr_m - the lumped-mass rule FAST uses AND the program uses\n"
             "#   (section mass at the element centre); it reproduces Table 2-2 to 5 digits. A program run with this mass\n"
             "#   distribution on these 17 sections must give these numbers (mass, and with omega the root centrifugal force\n"
             "#   F = omega^2 * sum(m_i dr_i r_i) about the rotor axis - row centrifugal_root_force_n_at_12.1rpm).\n")
# ---------------------------------------------------------------------------------------------------------------------
# geometry: Table 3-1
# ---------------------------------------------------------------------------------------------------------------------
afmap = {"Cylinder1.dat": "circle", "Cylinder2.dat": "circle_cd035", "DU40_A17.dat": "du40", "DU35_A17.dat": "du35",
         "DU30_A17.dat": "du30", "DU25_A17.dat": "du25", "DU21_A17.dat": "du21", "NACA64_A17.dat": "naca64_618"}
tc = {"circle": 1.0, "circle_cd035": 1.0, "du40": 0.405, "du35": 0.35, "du30": 0.30, "du25": 0.25, "du21": 0.21,
      "naca64_618": 0.18}
geo = []
with open(os.path.join(HERE, "jonkman2009_table3-1_pdftotext.txt")) as f:
    for line in f:
        tok = line.split()
        if len(tok) == 6 and tok[0].isdigit() and tok[5] in afmap:
            geo.append(dict(r_m=float(tok[1]), chord_m=float(tok[4]), twist_deg=float(tok[2]),
                            airfoil=afmap[tok[5]], thickness_ratio=tc[afmap[tok[5]]], dr_m=float(tok[3])))
geo = pd.DataFrame(geo)
assert len(geo) == 17 and abs(geo.dr_m.sum() - 61.5) < 1e-3
mm = np.interp(geo.r_m, t21.r_m, t21.mass_per_length_kg_m) * 1.04536
M_mid = (mm * geo.dr_m).sum()
m1_mid = (mm * geo.dr_m * (geo.r_m - 1.5)).sum()
m2_mid = (mm * geo.dr_m * (geo.r_m - 1.5) ** 2).sum()
mid = {"blade_length_m": geo.dr_m.sum(), "blade_mass_kg": M_mid, "first_mass_moment_root_kgm": m1_mid,
       "second_mass_moment_root_kgm2": m2_mid, "cm_from_root_m": m1_mid / M_mid}
for r_ in integ_rows:
    r_["midpoint_17_elements"] = mid.get(r_["quantity"], np.nan)
om = 12.1 * 2 * np.pi / 60
integ_rows.append(dict(quantity="centrifugal_root_force_n_at_12.1rpm", value=np.nan,
                       midpoint_17_elements=om ** 2 * (mm * geo.dr_m * geo.r_m).sum(),
                       trapezoid_table21=om ** 2 * np.trapezoid(mden * (s + 1.5), s),
                       note="derived here (not in the source): omega^2 sum(m dr r), rated 12.1 rpm, no precone"))
integ = pd.DataFrame(integ_rows)[["quantity", "value", "midpoint_17_elements", "trapezoid_table21", "note"]]
with open(os.path.join(CASE, "experiment", "integrated_properties.csv"), "w") as f:
    f.write(INTEG_HDR)
    integ.to_csv(f, index=False, float_format="%.7g")
print(integ.to_string())
with open(os.path.join(CASE, "geometry.csv"), "w") as f:
    f.write("# NREL 5-MW reference blade, Jonkman et al. 2009, Table 3-1 (p. 7), tabulated. r_m = RNodes (element centres,\n"
            "# from the rotor apex along the pitch axis; root 1.5 m, tip 63 m); twist_deg = AeroTwst (chord vs rotor plane at\n"
            "# zero pitch, positive towards feather - already the solver's wind-turbine convention, no conversion);\n"
            "# dr_m = DRNodes element length. thickness_ratio: DU 00-W2-401 0.405, DU 00-W2-350 0.35, DU 97-W-300 0.30,\n"
            "# DU 91-W2-250 0.25, DU 93-W-210 0.21, NACA 64-618 0.18 (airfoil names; Jonkman 2009 §3), cylinders 1.0.\n")
    geo.to_csv(f, index=False, float_format="%.6g")

# ---------------------------------------------------------------------------------------------------------------------
# airfoils: coordinates and polars from r-test
# ---------------------------------------------------------------------------------------------------------------------
pol_all = []
for src, name in afmap.items():
    base = src.replace(".dat", "")
    with open(os.path.join(RT, base + "_coords.txt")) as f:
        L = f.read().splitlines()
    xy = []
    started = False
    for l in L:
        if l.strip().startswith("!") or not l.strip():
            continue
        tok = l.split()
        try:
            a, b = float(tok[0]), float(tok[1])
        except (ValueError, IndexError):
            continue
        xy.append((a, b))
    xy = xy[1:]  # first numeric pair is the airfoil reference point
    with open(os.path.join(CASE, "airfoils", name + ".dat"), "w") as f:
        f.write("%s (OpenFAST r-test 5MW_Baseline/Airfoils/%s_coords.txt)\n" % (name, base))
        for a, b in xy:
            f.write("%.6f %.6f\n" % (a, b))
    with open(os.path.join(RT, src)) as f:
        L = f.read().splitlines()
    n_alf = [int(l.split()[0]) for l in L if "NumAlf" in l][0]
    re_m = [float(l.split()[0]) for l in L if re.search(r"\bRe\b", l) and "Reynolds" in l][0]
    i = [k for k, l in enumerate(L) if "NumAlf" in l][0] + 1
    tab = []
    while len(tab) < n_alf:
        l = L[i]
        i += 1
        if l.strip().startswith("!") or not l.strip():
            continue
        tab.append([float(v) for v in l.split()[:4]])
    tab = np.array(tab)
    if name.startswith("circle"):
        keep = tab  # constant cd, cl = 0 (3 rows: -180, 0, 180)
    else:
        keep = tab[(tab[:, 0] >= -20) & (tab[:, 0] <= 30)]
    pol = pd.DataFrame(dict(re=re_m * 1e6, alpha_deg=keep[:, 0], cl=keep[:, 1], cd=keep[:, 2], cm=keep[:, 3],
                            source="Jonkman 2009 airfoil table %s (DOWEC 2D data + AirfoilPrep 3D stall-delay and "
                                   "Viterna corrections), OpenFAST r-test copy; -20..30 deg kept" % src))
    with open(os.path.join(CASE, "airfoils", name + "_polars.csv"), "w") as f:
        f.write("# %s polar from OpenFAST r-test 5MW_Baseline/Airfoils/%s (the NREL 5-MW model's own table, Jonkman 2009\n"
                "# §3: 2D DOWEC data corrected for rotational stall delay and with Viterna drag; cylinders cl=0, cd 0.50\n"
                "# (Cylinder1) / 0.35 (Cylinder2)). Re as stated in the file (0.75e6, nominal). Full -180..180 table in\n"
                "# sources/openfast_r-test/.\n" % (name, src))
        pol.to_csv(f, index=False, float_format="%.6g")
    pol_all.append(pol.assign(name=name))

# ---------------------------------------------------------------------------------------------------------------------
# derived: equivalent homogeneous properties of the program's section model (information only)
# ---------------------------------------------------------------------------------------------------------------------
def poly_props(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    xj, yj = np.roll(x, -1), np.roll(y, -1)
    cr = x * yj - xj * y
    A = cr.sum() / 2
    cx, cy = ((x + xj) * cr).sum() / (6 * A), ((y + yj) * cr).sum() / (6 * A)
    Iyy = (cr * (y * y + y * yj + yj * yj)).sum() / 12 - A * cy ** 2
    Ixx = (cr * (x * x + x * xj + xj * xj)).sum() / 12 - A * cx ** 2
    per = np.hypot(xj - x, yj - y).sum()
    return abs(A), abs(Iyy), abs(Ixx), per


der = []
for _, g in geo.iterrows():
    d = np.loadtxt(os.path.join(CASE, "airfoils", g.airfoil + ".dat"), skiprows=1)
    A, Ixf, Iye, per = poly_props(d[:, 0], d[:, 1])
    c = g.chord_m
    m_ = np.interp(g.r_m, t21.r_m, t21.mass_per_length_kg_m)
    eif = np.interp(g.r_m, t21.r_m, t21.ei_flap_nm2)
    eie = np.interp(g.r_m, t21.r_m, t21.ei_edge_nm2)
    der.append(dict(r_m=g.r_m, airfoil=g.airfoil, chord_m=c, A_solid_m2=A * c ** 2, Ix_solid_m4=Ixf * c ** 4,
                    Iy_solid_m4=Iye * c ** 4, perimeter_m=per * c, mass_per_length_kg_m=m_, ei_flap_nm2=eif,
                    ei_edge_nm2=eie, rho_eq_solid_kg_m3=m_ / (A * c ** 2), E_eq_flap_solid_pa=eif / (Ixf * c ** 4),
                    E_eq_edge_solid_pa=eie / (Iye * c ** 4), wall_eq_m_rho1850=m_ / (1850.0 * per * c)))
der = pd.DataFrame(der)
with open(os.path.join(CASE, "experiment", "derived_equivalent_homogeneous.csv"), "w") as f:
    f.write("# DERIVED HERE, NOT FROM THE SOURCE (information only) - what a homogeneous section model (the program's) needs\n"
            "# to reproduce the NREL 5-MW blade. Section geometry: the airfoil outline of geometry.csv scaled by the chord, not\n"
            "# rotated (chord frame; Table 2-1 EI are about principal axes turned by the structural twist - small difference).\n"
            "# Mass/EI: Table 2-1 linearly interpolated to r_m (mass WITHOUT the 4.536 %% scaling).\n"
            "# rho_eq_solid = m'/A_solid; E_eq_*_solid = EI/I_solid; wall_eq_m_rho1850 = thin-shell wall thickness of a GFRP\n"
            "# shell (rho 1850 kg/m3, assumed) that carries m' (m' / (rho * perimeter)).\n")
    der.to_csv(f, index=False, float_format="%.5g")
print(der[["r_m", "airfoil", "rho_eq_solid_kg_m3", "E_eq_flap_solid_pa", "E_eq_edge_solid_pa", "wall_eq_m_rho1850"]].to_string())

# ---------------------------------------------------------------------------------------------------------------------
# check plot
# ---------------------------------------------------------------------------------------------------------------------
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(2, 2, figsize=(12, 8))
ax[0, 0].plot(geo.r_m, geo.chord_m, "o-", label="chord [m]")
ax[0, 0].plot(geo.r_m, geo.twist_deg, "s-", label="twist [deg]")
ax[0, 0].plot(t21.r_m, t21.strctwist_deg, "k:", label="StrcTwst (Table 2-1)")
ax[0, 0].set_xlabel("r [m]")
ax[0, 0].legend()
ax[0, 0].set_title("geometry (Table 3-1)")
ax[0, 1].semilogy(t21.r_m, t21.ei_flap_nm2, label="EI flap")
ax[0, 1].semilogy(t21.r_m, t21.ei_edge_nm2, label="EI edge")
ax[0, 1].semilogy(t21.r_m, t21.gj_nm2, label="GJ")
ax[0, 1].semilogy(t21.r_m, t21.ea_n, label="EA [N]")
ax[0, 1].semilogy(1.5 + np.array(bd_eta) * 61.5, bd_K[:, 4, 4], "k.", ms=3, label="BeamDyn K55")
ax[0, 1].legend()
ax[0, 1].set_title("stiffness (Table 2-1)")
ax[1, 0].plot(t21.r_m, t21.mass_per_length_kg_m, "o-", ms=3, label="Table 2-1")
ax[1, 0].plot(1.5 + ed[:, 0] * 61.5, ed[:, 2], "k.", ms=3, label="ElastoDyn")
ax[1, 0].set_ylabel("mass per length [kg/m]")
ax[1, 0].set_xlabel("r [m]")
ax[1, 0].legend()
for p in pol_all:
    if not p.name.iloc[0].startswith("circle"):
        ax[1, 1].plot(p.alpha_deg, p.cl, label=p.name.iloc[0])
ax[1, 1].set_xlabel("alpha [deg]")
ax[1, 1].set_ylabel("cl")
ax[1, 1].legend(fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "check_structural_nrel5mw.png"), dpi=110)
print(geo.to_string())

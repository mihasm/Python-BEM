"""Build the TUM G1 validation case from the raw TWEET-IE blind-test files in extraction/raw/.

Everything here is tabulated data (no digitising). Run from anywhere:
    python3 build_tum_g1.py
Writes geometry.csv, airfoils/*.dat, airfoils/*_polars.csv, experiment/*.csv and extraction/check_tum_g1.png.

Raw files (copied unchanged from Zenodo 10.5281/zenodo.13972361 'Data files.zip' and
10.5281/zenodo.15607777 'NTUA Dataset.zip', both CC-BY-4.0):
  raw/AeroDyn_blade.dat              - TUM AeroDyn15 blade table, 60 nodes, polar file per node
  raw/Blade_design_G1_RG14.xlsx      - TUM design table, 41 stations (chord, twist, t/c, pitch axis)
  raw/RG14 airfoil.txt               - RG14 coordinates
  raw/polars/section_1..10.dat       - TUM polars, 4 Re (0.047, 0.057, 0.069, 0.084 million)
  raw/case0/*.dat                    - measured Cp(/Ct) vs TSR, single turbines, TUM wind tunnel
  raw/ntua/*.dat                     - one measured operating point each, NTUA wind tunnel (Phase II)
"""
import os
import re
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.dirname(HERE)
RAW = os.path.join(HERE, "raw")

HUB_R = 0.0535      # ElastoDyn HubRad of the TUM OpenFAST model (m)
TIP_R = 0.5495      # ElastoDyn TipRad (m); docs give D = 1.1 m
R_DOC = 0.55        # radius used by the data providers for Cp, Ct, TSR (D = 1.1 m)
U_TUM = 5.74        # m/s, Pitot 4D upstream, uncorrected (Phase I documentation, Cases 1-2; Case 0 consistent)

SRC_ZEN1 = "TWEET-IE blind test Phase I, Zenodo doi:10.5281/zenodo.13972361 (CC-BY-4.0)"
SRC_ZEN2 = "TWEET-IE blind test Phase II, Zenodo doi:10.5281/zenodo.15607777 (CC-BY-4.0)"


# ---------------------------------------------------------------- geometry
def read_aerodyn_blade(fn):
    rows = []
    for line in open(fn):
        p = line.split()
        if len(p) == 7 and re.match(r"^-?\d", p[0]) and p[6].startswith("section_"):
            rows.append([float(p[0]), float(p[4]), float(p[5]), p[6]])
    return pd.DataFrame(rows, columns=["span", "twist", "chord", "polar"])


def read_design_xlsx(fn):
    x = pd.read_excel(fn, header=None)
    d = x.iloc[3:, :9].copy()
    d.columns = ["r", "span", "chord", "twist", "thick_pct", "pitch_axis_pct", "ac_pct", "af_id", "profile"]
    for c in ["r", "span", "chord", "twist", "thick_pct"]:
        d[c] = d[c].astype(float)
    return d.reset_index(drop=True)


ad = read_aerodyn_blade(os.path.join(RAW, "AeroDyn_blade.dat"))
des = read_design_xlsx(os.path.join(RAW, "Blade_design_G1_RG14.xlsx"))

AFNAME = {f"section_{k}": (f"g1_sec{k}" if 2 <= k <= 9 else None) for k in range(1, 11)}
AFNAME["section_1"] = "circle"
AFNAME["section_10"] = "rg14"

geo = pd.DataFrame({
    "r_m": np.round(HUB_R + ad.span.values, 4),
    "chord_m": ad.chord.values,
    "twist_deg": ad.twist.values,
    "airfoil": [AFNAME[p] for p in ad.polar],
    # t/c interpolated from the 41-station design table on the blade-length coordinate
    "thickness_ratio": np.round(np.interp(ad.span.values, des.span.values, des.thick_pct.values) / 100.0, 4),
})
hdr = [
    "# TUM G1 scaled wind turbine (D = 1.1 m, 3 blades, RG14 airfoil) - blade geometry",
    f"# Source: TUM AeroDyn15 blade input 'AeroDyn_blade.dat' (60 nodes) from {SRC_ZEN1},",
    "#   folder 'G1 Wind Turbine Model/Blade distributions'. Tabulated, copied unchanged (chord, twist).",
    f"# r_m = HubRad + BlSpn with HubRad = {HUB_R} m (TUM OpenFAST ElastoDyn); last node r = {TIP_R} m = TipRad.",
    "# twist_deg = AeroDyn BlTwist = local chord angle to the rotor plane WITHOUT pitch (OpenFAST convention,",
    "#   alpha = phi - (twist + pitch)); identical to the solver convention. Global pitch 0.4 deg is in the experiment files.",
    "# airfoil: TUM polar assignment of AeroDyn_blade.dat: section_1 = circle (cylinder root), section_2..9 = root",
    "#   transition (g1_sec2..g1_sec9), section_10 = RG14 (rg14).",
    "# thickness_ratio: interpolated from the 41-station design table 'Blade_design_G1_RG14.xlsx' (same source),",
    "#   chord/twist of the two tables agree to <0.001 m / <0.1 deg (see check_tum_g1.png).",
]
with open(os.path.join(CASE, "geometry.csv"), "w") as f:
    f.write("\n".join(hdr) + "\n")
    geo.to_csv(f, index=False)

# ---------------------------------------------------------------- airfoils
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


rg14 = read_xy(os.path.join(RAW, "RG14 airfoil.txt"))
ile = int(np.argmin(rg14[:, 0]))
# circle of diameter = chord on the same x distribution (upper surface y>0 from TE to LE, then lower)
xc = rg14[:, 0]
yc = 0.5 * np.sqrt(np.clip(1.0 - (2.0 * xc - 1.0) ** 2, 0.0, None))
yc[ile:] *= -1.0
yc[ile] = 0.0
circle = np.column_stack([xc, yc])
t_rg14 = float(np.max(np.interp(np.linspace(0, 1, 2001), rg14[:ile + 1][::-1, 0], rg14[:ile + 1][::-1, 1])
                      - np.interp(np.linspace(0, 1, 2001), rg14[ile:, 0], rg14[ile:, 1])))


def write_dat(name, xy, title):
    with open(os.path.join(CASE, "airfoils", name + ".dat"), "w") as f:
        f.write(title + "\n")
        for x, y in xy:
            f.write(f"{x:9.6f} {y:9.6f}\n")


write_dat("rg14", rg14, "RG14 (TUM G1 blade airfoil, coordinates from TWEET-IE dataset 'RG14 airfoil.txt')")
write_dat("circle", circle, "circle (cylindrical blade root, generated, D = chord)")

# transition sections: no coordinates are published; generated as a blend circle <-> RG14 whose
# maximum thickness equals the t/c interpolated from the design table at the node where TUM uses that polar
first_node = {}
for _, g in geo.iterrows():
    first_node.setdefault(g.airfoil, g)
blend_info = {}
for k in range(2, 10):
    name = f"g1_sec{k}"
    t = float(first_node[name].thickness_ratio)
    w = min(max((t - t_rg14) / (1.0 - t_rg14), 0.0), 1.0)
    xy = (1 - w) * rg14 + w * circle
    write_dat(name, xy, f"{name} G1 root transition, GENERATED blend circle/RG14 with t/c = {t:.3f} (no published coordinates)")
    blend_info[name] = (t, w)

# ---------------------------------------------------------------- polars
def read_aerodyn_polars(fn):
    lines = open(fn).read().splitlines()
    out = []
    i = 0
    while i < len(lines):
        if "Reynolds number in millions" in lines[i]:
            re_ = float(lines[i].split()[0]) * 1e6
            while "NumAlf" not in lines[i]:
                i += 1
            n = int(lines[i].split()[0])
            i += 1
            while lines[i].strip().startswith("!"):
                i += 1
            for k in range(n):
                a, cl, cd, cm = [float(v) for v in lines[i + k].split()[:4]]
                out.append((re_, a, cl, cd, cm))
            i += n
        else:
            i += 1
    return pd.DataFrame(out, columns=["re", "alpha_deg", "cl", "cd", "cm"])


for k in range(1, 11):
    name = AFNAME[f"section_{k}"]
    p = read_aerodyn_polars(os.path.join(RAW, "polars", f"section_{k}.dat"))
    # store only the finely tabulated range -20..+20 deg (0.2 deg steps); outside it TUM gives a coarse
    # 10-deg extrapolation which the runner replaces by its own (Montgomerie) extrapolation
    p = p[(p.alpha_deg >= -20.0) & (p.alpha_deg <= 20.0)].copy()
    # the 4th column of the TUM files duplicates cl exactly (cm == cl in every row) -> not a moment coefficient; dropped
    assert np.allclose(p.cm, p.cl) or k == 1
    p = p.drop(columns="cm")
    p["re"] = p.re.round().astype(int)
    p["source"] = f"TUM polars section_{k}.dat (TWEET-IE Phase I dataset); nominal G1 polars identified from rotor measurements"
    with open(os.path.join(CASE, "airfoils", name + "_polars.csv"), "w") as f:
        f.write(f"# {name}: TUM-supplied polars '{'Polars/section_%d.dat' % k}' from {SRC_ZEN1}\n")
        f.write("# Tabulated, copied unchanged for -20..+20 deg (0.2 deg steps) at Re = 47 000, 57 000, 69 000, 84 000.\n")
        f.write("# The dataset does not state their origin; Re values and shape match the 'nominal' G1 polars p* of\n")
        f.write("#   Wang, Campagnolo & Bottasso (2020) WES 5:1537, Fig. 5, which were identified from G1 rotor\n")
        f.write("#   measurements with the method of Bottasso et al. (2014a) - i.e. NOT independent 2D wind-tunnel data.\n")
        if k == 1:
            f.write("# section_1 = cylinder root: cl = 0, cd = 1.0 constant (TUM value).\n")
        elif 2 <= k <= 9:
            t, w = blend_info[name]
            f.write(f"# Root transition section (TUM 'transition' polar); nominal t/c at its first node = {t:.3f}.\n")
        f.write("# cm is not given: the 4th column of the TUM file is an exact copy of cl (file error), so it was dropped.\n")
        f.write("# Columns: re, alpha_deg, cl, cd, source\n")
        p.to_csv(f, index=False)

# ---------------------------------------------------------------- experiments
def read_case0(fn):
    rows = [l.split() for l in open(fn) if l.strip() and l.split()[0].isdigit()]
    return np.array(rows, float)


def tsr_from_rpm(rpm, U):
    return rpm * 2 * np.pi / 60.0 * R_DOC / U


exp_dir = os.path.join(CASE, "experiment")
common_hdr = [
    f"# Data: {SRC_ZEN1}, folder 'CASE 0 (Single Wind Turbines)'. Tabulated, copied unchanged",
    "#   (power, cp, ct and their combined standard uncertainties = the '_err' columns; rpm).",
    "# Test: TUM closed-loop boundary-layer wind tunnel, test section 2.7 m (W) x 1.8 m (H); single G1 turbine;",
    "#   constant tunnel speed, rotor speed varied. U_inf = 5.74 m/s from a Pitot tube 4D upstream, UNCORRECTED.",
    "#   Blockage ratio A_rotor/A_tunnel = 20 %; NO blockage correction applied (Phase I documentation p. 4).",
    "#   Cp = P/(0.5 rho U^3 pi R^2), Ct = T/(0.5 rho U^2 pi R^2) with R = 0.55 m (D = 1.1 m). Inflow is a",
    "#   low-turbulence boundary-layer profile (hub height 0.82 m); see raw inflow profiles in the Zenodo record.",
    "# tsr: recomputed = (rpm*2*pi/60)*0.55/5.74 (the source gives TSR rounded to 0.1; tsr_source keeps that value;",
    "#   it agrees within 0.015 with the set-points 3.0, 3.5, ... 5.25, 5.5, 5.75, ... 8.25, 8.75 ..., which the source",
    "#   prints rounded, e.g. 5.25 -> 5.3).",
    "# rho: not tabulated; rho = P/(cp*0.5*U^3*A) gives 1.13-1.17 (cp rounded to 3 decimals) -> defaults rho = 1.145.",
    "# pitch_deg: 0.4 deg (Phase I documentation: 'blade pitch ~0.4 deg'; OpenFAST BlPitch = 0.4; CAD 'pitch_0d4').",
]

files = {
    "performance_tum_upstream_turbine.csv": ("Single_Upstream turbine_upstream position.dat", True,
        "Upstream G1 turbine alone at the upstream position (paper symbol: red circle). Cp AND Ct."),
    "performance_tum_downstream_turbine_upstream_pos.csv": ("Single_Downstream turbine_upstream position.dat", False,
        "Second (nominally identical) G1 turbine alone at the upstream position - rotor-to-rotor repeatability. Cp only."),
    "performance_tum_downstream_turbine_downstream_pos.csv": ("Single_Downstream turbine_downstream position.dat", False,
        "Second G1 turbine alone at the downstream position (no turbine upstream; local inflow differs slightly). Cp only."),
}
exp_frames = {}
for out, (raw, has_ct, desc) in files.items():
    a = read_case0(os.path.join(RAW, "case0", raw))
    rpm = a[:, -1]
    df = pd.DataFrame({
        "tsr": np.round(tsr_from_rpm(rpm, U_TUM), 4),
        "cp": a[:, 5], "cp_err": a[:, 8],
    })
    if has_ct:
        df["ct"] = a[:, 9]
        df["ct_err"] = a[:, 12]
    df["power_w"] = a[:, 1]
    df["power_w_err"] = a[:, 4]
    df["rpm"] = rpm
    df["tsr_source"] = a[:, -2]
    assert (np.abs(df.tsr - df.tsr_source) < 0.07).all(), out
    with open(os.path.join(exp_dir, out), "w") as f:
        f.write(f"# TUM G1 - {desc}\n# Raw file: '{raw}'\n")
        f.write("\n".join(common_hdr) + "\n")
        if has_ct:
            f.write("# Ct: thrust from the tower-base load cell. Its type-A uncertainty grows strongly above tsr ~ 8.5\n")
            f.write("#   (ct_err = 0.6 at tsr 9.5 and 2.25 at tsr 10.0 in the source) - treat those Ct points as unreliable.\n")
        f.write("# Extraction uncertainty: none (tabulated). Measurement uncertainty: the *_err columns (combined, as published).\n")
        df.to_csv(f, index=False)
    exp_frames[out] = df

# NTUA single points (Phase II), larger tunnel, 11 % blockage
ntua_rows = []
for raw, ti in [("Reference_Low_Turb_5D_WT_1.dat", "low"), ("Reference_High_Turb_5D_WT_1.dat", "high")]:
    a = read_case0(os.path.join(RAW, "ntua", raw))[0]
    ntua_rows.append({"tsr": a[-2], "cp": a[5], "cp_err": a[8], "power_w": a[1], "power_w_err": a[4],
                      "rpm": a[-1], "v_ms": round(a[-1] * 2 * np.pi / 60 * R_DOC / a[-2], 3),
                      "turbulence": ti})
ntua = pd.DataFrame(ntua_rows)
with open(os.path.join(exp_dir, "performance_ntua_upstream_turbine_point.csv"), "w") as f:
    f.write("# TUM G1 - upstream turbine, baseline (no wake control) operating point in the NTUA wind tunnel\n")
    f.write(f"# Data: {SRC_ZEN2}, 'NTUA Dataset/Case 1 - Baseline/<Low|High> Turbulence/Reference_*_5D_WT_1.dat'.\n")
    f.write("#   Tabulated, copied unchanged (power, cp, combined uncertainty, tsr, rpm). A downstream G1 runs 5D behind.\n")
    f.write("# Test: NTUA closed-circuit tunnel, section 2.5 m (H) x 3.5 m (W); blockage 11 %, no blockage correction;\n")
    f.write("#   U from a Pitot tube 2D upstream. 'high' turbulence = screen + passive grid. Same G1 rotor, pitch ~0.4 deg.\n")
    f.write("# v_ms: NOT tabulated; derived here as (rpm*2*pi/60)*0.55/tsr (tsr is rounded to 0.1 -> v uncertain by ~0.6 %).\n")
    f.write("# Note the much lower Cp than in the TUM tunnel at the same tsr (0.39 vs 0.50): blockage 11 % vs 20 % and\n")
    f.write("#   different inflow; useful as a reminder of the size of tunnel effects on this rotor.\n")
    ntua.to_csv(f, index=False)

# ---------------------------------------------------------------- check plot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(2, 3, figsize=(15, 8.5))
a0 = ax[0, 0]
a0.plot(geo.r_m, geo.chord_m, "b.-", label="geometry.csv (AeroDyn 60 nodes)")
a0.plot(des.r, des.chord, "r+", label="design xlsx 41 stations")
a0.set_xlabel("r [m]"); a0.set_ylabel("chord [m]"); a0.legend(fontsize=8); a0.grid()
a1 = ax[0, 1]
a1.plot(geo.r_m, geo.twist_deg, "b.-", label="geometry.csv")
a1.plot(des.r, des.twist, "r+", label="design xlsx")
a1.set_xlabel("r [m]"); a1.set_ylabel("twist [deg]"); a1.legend(fontsize=8); a1.grid()
a2 = ax[0, 2]
for k in [1, 2, 4, 6, 8]:
    nm = f"g1_sec{k}" if 2 <= k <= 9 else "circle"
    xy = read_xy(os.path.join(CASE, "airfoils", nm + ".dat"))
    a2.plot(xy[:, 0], xy[:, 1], lw=0.8, label=nm)
a2.plot(rg14[:, 0], rg14[:, 1], "k", lw=1.5, label="rg14")
a2.set_aspect("equal"); a2.legend(fontsize=7); a2.set_title("airfoils")
a3 = ax[1, 0]
pr = pd.read_csv(os.path.join(CASE, "airfoils", "rg14_polars.csv"), comment="#")
for re_, g in pr.groupby("re"):
    a3.plot(g.alpha_deg, g.cl, label=f"Re {re_:.0f}")
a3b = a3.twinx()
for re_, g in pr.groupby("re"):
    a3b.plot(g.alpha_deg, g.cd, "--", lw=0.8)
a3.set_xlabel("alpha [deg]"); a3.set_ylabel("cl (solid)"); a3b.set_ylabel("cd (dashed)"); a3.legend(fontsize=7); a3.grid()
a4 = ax[1, 1]
for out, df in exp_frames.items():
    a4.errorbar(df.tsr, df.cp, yerr=df.cp_err, fmt="o-", ms=3, label=out.replace("performance_", "").replace(".csv", ""))
a4.errorbar(ntua.tsr, ntua.cp, yerr=ntua.cp_err, fmt="ks", label="NTUA point")
a4.axhline(16 / 27, color="gray", ls=":")
a4.set_xlabel("tsr"); a4.set_ylabel("cp"); a4.legend(fontsize=7); a4.grid()
a5 = ax[1, 2]
df = exp_frames["performance_tum_upstream_turbine.csv"]
a5.errorbar(df.tsr, df.ct, yerr=np.minimum(df.ct_err, 0.3), fmt="o-", ms=3)
a5.set_ylim(0, 1.2); a5.set_xlabel("tsr"); a5.set_ylabel("ct (err bars clipped at 0.3)"); a5.grid()
fig.suptitle("TUM G1 validation case - check")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "check_tum_g1.png"), dpi=110)

# ---------------------------------------------------------------- sanity
assert abs(geo.r_m.iloc[-1] - TIP_R) < 1e-4
assert (np.diff(geo.twist_deg) < 0).all()
assert (df.cp < 16 / 27).all()
print("t/c RG14 from coordinates:", round(t_rg14, 4))
print("blends:", {k: (round(v[0], 3), round(v[1], 3)) for k, v in blend_info.items()})
print(geo.head(12).to_string()); print("...", len(geo), "stations")
print("design vs aerodyn chord max diff:",
      np.abs(np.interp(des.span, ad.span, ad.chord) - des.chord)[3:].max(),
      " twist max diff:", np.abs(np.interp(des.span, ad.span, ad.twist) - des.twist)[3:].max())

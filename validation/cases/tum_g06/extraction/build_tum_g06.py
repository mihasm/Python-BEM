"""Build the TUM G06 case (Nanos et al. 2022, WES 7:1263) from the digitised Fig. 2 and Fig. 8.
Run digitise_fig2.py and digitise_fig8.py (run_method4) first; this script only converts and documents.
"""
import os
import shutil
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.dirname(HERE)
G1 = os.path.join(os.path.dirname(CASE), "tum_g1")
R = 0.3
SRC = "Nanos E.M., Bottasso C.L., Manolas D.I., Riziotis V.A. (2022), Wind Energ. Sci. 7, 1263-1287, doi:10.5194/wes-7-1263-2022"

# ---------------------------------------------------------------- geometry (Fig. 2, digitised)
g = pd.read_csv(os.path.join(HERE, "fig2_digitised.csv"))
geo = pd.DataFrame({"r_m": np.round(g.r_R * R, 5), "chord_m": np.round(g.c_R * R, 5),
                    "twist_deg": np.round(g.twist_deg, 3), "airfoil": "rg14", "thickness_ratio": 0.0848})
with open(os.path.join(CASE, "geometry.csv"), "w") as f:
    f.write("# TUM G06 scaled wind turbine (D = 0.6 m, 3 blades, RG-14 over the whole span) - blade planform\n")
    f.write(f"# Source: {SRC}, Fig. 2 'Chord and twist distributions along the blade' (c/R and twist vs r/R).\n")
    f.write("# DIGITISED from the embedded raster (2067x1400 px) with extraction/digitise_fig2.py; overlay: extraction/fig2_overlay.png.\n")
    f.write("# Extraction uncertainty: c/R +-0.0015 (+-0.45 mm), twist +-0.2 deg (line width ~8 px = 0.0017 c/R).\n")
    f.write("# r_m = r/R * 0.3; chord_m = c/R * 0.3. Twist as plotted ('Blade twist', positive = towards feather, i.e. the\n")
    f.write("#   solver convention; pitch beta of Fig. 8 is added on top). The blade starts at r/R = 0.06: inboard the figure shows\n")
    f.write("#   a cylindrical root (c/R = 0.067 for r/R < 0.04) and a jump to the maximum chord at r/R ~0.05, which is not modelled.\n")
    f.write("# airfoil: RG-14 everywhere (paper Table 1 / Sect. 2.3.1); t/c of RG-14 = 0.0848.\n")
    geo.to_csv(f, index=False)

# ---------------------------------------------------------------- airfoil (same RG-14 data as tum_g1)
shutil.copy(os.path.join(G1, "airfoils", "rg14.dat"), os.path.join(CASE, "airfoils", "rg14.dat"))
p = pd.read_csv(os.path.join(G1, "airfoils", "rg14_polars.csv"), comment="#")
p["source"] = "TUM G1 RG14 polars (TWEET-IE dataset section_10.dat) reused for G06"
with open(os.path.join(CASE, "airfoils", "rg14_polars.csv"), "w") as f:
    f.write("# RG-14 polars REUSED from the TUM G1 case (validation/cases/tum_g1/airfoils/rg14_polars.csv): TUM-supplied polars\n")
    f.write("#   'section_10.dat' of the TWEET-IE blind-test dataset, Zenodo doi:10.5281/zenodo.13972361, Re 47k-84k.\n")
    f.write("# The G06 paper uses polars tuned from G06 measurements (Sect. 3.2.2, Fig. 9) which are NOT published; the G1\n")
    f.write("#   polars are the closest public RG-14 data from the same group and Reynolds range (G06 Re ~4e4-8e4).\n")
    f.write("# Tabulated (copied), -20..+20 deg.\n")
    p.to_csv(f, index=False)

# ---------------------------------------------------------------- experiments (Fig. 8, digitised)
d = pd.read_csv(os.path.join(HERE, "fig8_markers_final.csv"))
BETA = {"c0": 0.0, "c1": 2.0, "c2": 3.0, "c3": 4.0, "c4": 6.0}
VEL = {"c0": 6.0, "c1": 7.0, "c2": 8.0, "c3": 9.0, "c4": 10.0}
METH = {2: "marker", 3: "marker", 4: "line"}


def merge(pcp, pct, key):
    rows = []
    for ck in BETA:
        a = d[(d.panel == pcp) & (d.curve == ck)].sort_values("tsr")
        b = d[(d.panel == pct) & (d.curve == ck)].sort_values("tsr")
        used = set()
        for _, ra in a.iterrows():
            row = {"curve": ck, "tsr": ra.tsr, "cp": ra.value, "ct": np.nan, "cp_read": METH[int(ra.method)], "ct_read": ""}
            if len(b):
                j = (b.tsr - ra.tsr).abs().idxmin()
                if abs(b.loc[j, "tsr"] - ra.tsr) < 0.12 and j not in used:
                    used.add(j)
                    row["ct"] = b.loc[j, "value"]; row["ct_read"] = METH[int(b.loc[j, "method"])]
                    row["tsr"] = round(0.5 * (ra.tsr + b.loc[j, "tsr"]), 4)
            rows.append(row)
        for j, rb in b.iterrows():
            if j not in used:
                rows.append({"curve": ck, "tsr": rb.tsr, "cp": np.nan, "ct": rb.value, "cp_read": "", "ct_read": METH[int(rb.method)]})
    df = pd.DataFrame(rows)
    if key == "pitch":
        df.insert(0, "pitch_deg", df.curve.map(BETA)); df["v_ms"] = 10.0
    else:
        df.insert(0, "v_ms", df.curve.map(VEL)); df["pitch_deg"] = 0.0
    df = df.drop(columns="curve").sort_values([df.columns[0], "tsr"])
    df["cp"] = df.cp.round(4); df["ct"] = df.ct.round(4)
    return df


hdr_common = [
    f"# Source: {SRC}, Fig. 8. DIGITISED from the embedded raster (PDF p. 12, 2067x1279 px) with extraction/digitise_fig8.py",
    "#   (gridline calibration, colour masks per curve, marker detection; overlay extraction/fig8_overlay_final.png).",
    "# cp_read/ct_read: 'marker' = centre of the plotted marker found automatically; 'line' = marker hidden under other curves,",
    "#   value read from the curve line at the 0.5-spaced tsr set-point (line interpolation between neighbouring markers).",
    "# Extraction uncertainty: tsr +-0.03, cp +-0.004, ct +-0.005 ('marker'); about twice that for 'line'.",
    "# Test: UTD BLAST boundary-layer wind tunnel, 2.8 m x 2.1 m, blockage 4.8 % (uncorrected), low turbulence (TI 0.15 %),",
    "#   uniform inflow (shear exponent 0) - paper Table 2. Torque from shaft strain gauges, thrust from the tower-base load cell",
    "#   (paper Sect. 2.5). Cp = P/(0.5 rho U^3 pi R^2), R = 0.3 m. Published measurement uncertainty ~1 % on loads, pitch +-0.3 deg.",
    "# Rows with cp or ct empty: the marker was not visible/detected in that panel.",
]
for fn, (pcp, pct, key, title) in {
    "performance_pitch_sweep_10ms.csv": ("a", "b", "pitch", "Fig. 8a,b: Cp, Ct vs tsr at U = 10 m/s for pitch beta = 0, 2, 3, 4, 6 deg"),
    "performance_speed_sweep_beta0.csv": ("d", "e", "speed", "Fig. 8d,e: Cp, Ct vs tsr at beta = 0 for U = 6, 7, 8, 9, 10 m/s (Reynolds-number effect)"),
}.items():
    df = merge(pcp, pct, key)
    with open(os.path.join(CASE, "experiment", fn), "w") as f:
        f.write(f"# TUM G06 - {title}\n")
        f.write("\n".join(hdr_common) + "\n")
        df.to_csv(f, index=False)
    print(fn, len(df))

# ---------------------------------------------------------------- check plot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(2, 3, figsize=(15, 8.5))
ax[0, 0].plot(geo.r_m, geo.chord_m, "o-", ms=3); ax[0, 0].set_ylabel("chord [m]"); ax[0, 0].set_xlabel("r [m]"); ax[0, 0].grid()
ax[0, 1].plot(geo.r_m, geo.twist_deg, "o-", ms=3); ax[0, 1].set_ylabel("twist [deg]"); ax[0, 1].set_xlabel("r [m]"); ax[0, 1].grid()
re_tip = 7.5 * 10 * geo.chord_m / 1.5e-5 * (geo.r_m / R)
ax[0, 2].plot(geo.r_m, 7.5 * 10 * (geo.r_m / R) * geo.chord_m / 1.5e-5, label="U=10, tsr 7.5")
ax[0, 2].plot(geo.r_m, 7.5 * 6 * (geo.r_m / R) * geo.chord_m / 1.5e-5, label="U=6, tsr 7.5")
ax[0, 2].set_ylabel("approx. chord Re (Omega r c / nu)"); ax[0, 2].legend(); ax[0, 2].grid()
for k, (fn, grp) in enumerate([("performance_pitch_sweep_10ms.csv", "pitch_deg"), ("performance_speed_sweep_beta0.csv", "v_ms")]):
    df = pd.read_csv(os.path.join(CASE, "experiment", fn), comment="#")
    for gv, gg in df.groupby(grp):
        ax[1, 0 + k].plot(gg.tsr, gg.cp, "o-", ms=3, label=f"{grp}={gv}")
        ax[1, 2].plot(gg.tsr, gg.ct, "s-" if k == 0 else "^--", ms=3, label=f"{grp}={gv}")
    ax[1, k].set_xlabel("tsr"); ax[1, k].set_ylabel("cp"); ax[1, k].legend(fontsize=7); ax[1, k].grid()
ax[1, 2].set_xlabel("tsr"); ax[1, 2].set_ylabel("ct"); ax[1, 2].legend(fontsize=6); ax[1, 2].grid()
fig.suptitle("TUM G06 validation case - check"); fig.tight_layout()
fig.savefig(os.path.join(HERE, "check_tum_g06.png"), dpi=110)
assert geo.r_m.max() <= R and (geo.chord_m > 0).all()

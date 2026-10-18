"""Builds the ntnu_bt1 case files (geometry, airfoil, polars, experiment CSVs) from the source texts and the
extraction outputs in this folder. Run from anywhere:  python3 build_case.py
Inputs (all in ../sources unless noted):
  bt5inv.txt            pdftotext -layout of Invitation_Blind_test_5.pdf (Saetran et al. 2017, zenodo 1218555)
  krogstad_lund_2012_we482.txt   pdftotext -layout of Krogstad & Lund (2012)
  sarlak18.txt          pdftotext -layout of Sarlak et al. (2018) Energies 11:1418
  wes-3-329-2018.txt    pdftotext -layout of Bartl et al. (2018a)
  ./krogstad_eriksen2013_fig3_raw.csv, ./krogstad_adaramola2012_fig4_raw.csv, ./krogstad_adaramola2012_fig5_raw.csv,
  ./bartl2018b_fig1_raw.csv, ./bartl2019_fig4_raw.csv   (made by the extract_/digitise_ scripts here)
"""
import os
import re
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.dirname(HERE)
SRC = os.path.join(CASE, "sources")
MINUS = "\u2212"


def txt(name):
    with open(os.path.join(SRC, name), encoding="utf-8") as f:
        return f.read()


def num(s):
    # pdftotext renders the minus sign of Krogstad & Lund (2012) as the control character \x02 (font MT2SYS,
    # MathTime symbol font, code 2 = minus; checked with PyMuPDF) (and the decimal
    # point after it as ':'), e.g. '\x020:71674' = -0.71674
    return float(s.replace(MINUS, "-").replace("\x02", "-").replace(":", "."))


def write_csv(path, header_lines, df, float_format=None):
    with open(path, "w", encoding="utf-8") as f:
        for h in header_lines:
            f.write("# " + h + "\n")
        df.to_csv(f, index=False, float_format=float_format, lineterminator="\n")
    print("wrote", os.path.relpath(path, CASE), len(df), "rows")


# ------------------------------------------------------------------ geometry (BT5 invitation Table C.1)
t = txt("bt5inv.txt")
blk = t[t.index("C.1. NTNU blade"):t.index("Table C.1:")]
geo = []
for line in blk.splitlines():
    p = line.split()
    if len(p) == 4 and p[0].isdigit():
        geo.append((int(p[0]), num(p[1]), num(p[2]), num(p[3])))
geo = pd.DataFrame(geo, columns=["element", "r_m", "chord_m", "twist_src"])
assert len(geo) == 30

# cross-check against Krogstad & Lund (2012) Table I (r/R, c/R, twist; R = 0.45 m)
kl = txt("krogstad_lund_2012_we482.txt")
blk = kl[kl.index("Table I. Blade definition."):kl.index("Wind Energ. 2012; 15:443", kl.index("Table I. Blade definition."))]
kl_rows = []
for line in blk.splitlines():
    p = line.split()
    if len(p) in (2, 3) and re.match(r"^[0-9.]+$", p[0]):
        tw = p[2] if len(p) == 3 else None
        kl_rows.append((num(p[0]) * 0.45, num(p[1]) * 0.45, num(tw) if tw not in (None, "–") else np.nan))
kl_df = pd.DataFrame(kl_rows, columns=["r", "c", "tw"])
assert len(kl_df) == 30, len(kl_df)
dr = np.abs(kl_df.r.values - geo.r_m.values).max()
dc = np.abs(kl_df.c.values - geo.chord_m.values).max()
dt = np.nanmax(np.abs(kl_df.tw.values[3:] - geo.twist_src.values[3:]))
print("KL2012 Table I vs BT5 Table C.1: max |dr| = %.2e m, |dc| = %.2e m, |dtwist| = %.2e deg" % (dr, dc, dt))

airfoil = ["circle"] * 3 + ["transition"] + ["s826"] * 26
twist = geo.twist_src.copy()
twist[:3] = 38.0      # 120 deg in the source is only a flag for "cylinder"; cl = 0 there, see header
thick = [1.0] * 3 + [np.nan] + [0.14] * 26
G = pd.DataFrame({"r_m": geo.r_m, "chord_m": geo.chord_m, "twist_deg": twist, "airfoil": airfoil,
                  "thickness_ratio": thick})
write_csv(os.path.join(CASE, "geometry.csv"), [
    "NTNU 'Blind test' rotor blade (turbine T2 of BT1-BT5, LARS1 of BT5; same blades on T1 of BT2-BT4).",
    "Source: Saetran, Muehle, Bartl, Schottler, Hoelling, Adaramola (2017) Invitation to the 'Blind test 5' workshop,",
    "  Appendix C, Table C.1 (tabulated), https://doi.org/10.5281/zenodo.1218555 ; identical to Krogstad & Lund (2012)",
    "  Wind Energy 15:443-457, Table I (given there as r/R, c/R with R = 0.45 m): max differences "
    "%.1e m in r, %.1e m in c, %.1e deg in twist." % (dr, dc, dt),
    "r_m: centre of 30 blade elements of 0.015 m from the rotor axis (outer edge of element 30 = 0.450 m).",
    "twist_deg: angle chord vs rotor plane (source: 'twist angle measured with respect to the rotor plane'),",
    "  blade pitch 0 in all tests; twist is about the quarter chord (K&L 2012 sec. 2.2). No conversion needed.",
    "Rows 1-3: circular cylinder d = 13.5 mm (chord column = diameter) used to fix the blade in the hub; the",
    "  source marks them with twist 120 deg (a flag). Here twist is set to 38.0 deg (= first transition row) so that",
    "  interpolation between rows does not produce spurious angles; irrelevant for the cylinder (cl = 0).",
    "Row 4 (r = 0.055 m): 'linear transition' from the cylinder (r <= 0.049 m) to the S826 (r >= 0.0675 m)",
    "  (K&L 2012: transition between r/R = 0.109 and 0.151). Rows 5-30: NREL S826, t/c = 0.14, sharp square tip.",
], G, float_format="%.6g")

# ------------------------------------------------------------------ S826 coordinates (BT5 invitation Table B.1)
blk = t[t.index("B.1 NTNU blade: NREL S826 airfoil"):t.index("Table B.1:")]
up, lo = [], []
for line in blk.splitlines():
    p = line.split()
    try:
        vals = [num(v) for v in p]
    except ValueError:
        continue
    if len(vals) == 4:
        up.append(vals[:2]); lo.append(vals[2:])
    elif len(vals) == 2:
        up.append(vals)
up, lo = np.array(up), np.array(lo)
assert len(up) == 34 and len(lo) == 32
sel = np.vstack([up[::-1], lo[1:]])
with open(os.path.join(CASE, "airfoils", "s826.dat"), "w") as f:
    f.write("NREL S826 (Somers 2005, NREL/SR-500-36344) as distributed in the NTNU Blind test 5 invitation, Table B.1\n")
    for x, y in sel:
        f.write("%9.5f %10.6f\n" % (x, y))
print("wrote airfoils/s826.dat", len(sel), "points; max t/c from coordinates:",
      round(float(np.interp(np.linspace(0, 1, 400), up[:, 0], up[:, 1]).max()
                  - np.interp(np.linspace(0, 1, 400), lo[:, 0], lo[:, 1]).min()), 4))

# ------------------------------------------------------------------ polars
# (1) DTU, Sarlak et al. 2018, Tables A2 (cl) and A3 (cd), Re 40k..400k  -> primary
s = txt("sarlak18.txt")
res = [40e3, 60e3, 80e3, 100e3, 120e3, 150e3, 200e3, 300e3, 400e3]


def table(start, end):
    b = s[s.index(start):s.index(end, s.index(start))]
    out = []
    for line in b.splitlines():
        line = re.sub(r"\([^)]*\)", "", line)          # drop the Re=40k downstroke values in parentheses
        p = line.split()
        if len(p) == 10:
            try:
                out.append([num(v) for v in p])
            except ValueError:
                pass
    return np.array(out)


A2 = table("Table A2. Lift coefficients", "Energies 2018, 11, 1418")
A3 = table("Table A3. Drag coefficients", "Energies 2018, 11, 1418")
assert A2.shape == (26, 10) and A3.shape == (26, 10) and np.allclose(A2[:, 0], A3[:, 0])
rows = []
for j, re_ in enumerate(res):
    for i in range(len(A2)):
        rows.append(dict(re=re_, alpha_deg=A2[i, 0], cl=A2[i, j + 1], cd=A3[i, j + 1],
                         source="exp DTU Sarlak 2018 Tab.A2/A3" + (" (upstroke)" if j == 0 else "")))
P = pd.DataFrame(rows)
write_csv(os.path.join(CASE, "airfoils", "s826_polars.csv"), [
    "ALTERNATIVE polar set (case.json alternative_polars): NREL S826 measured polars, DTU red wind tunnel (closed test section 0.5 x 0.5 m), 2D wing c = 0.100 m, span 0.499 m,",
    "  inflow TI < 0.2 %; lift from force balance, drag from wake rake (alpha <= ~15 deg reliable, beyond less so).",
    "Source: Sarlak, Frere, Mikkelsen, Sorensen (2018) Energies 11(6):1418, doi:10.3390/en11061418, Appendix",
    "  Tables A2 (cl) and A3 (cd), tabulated (copied exactly). Re = 40e3 ... 400e3. At Re = 40e3 the upstroke",
    "  (increasing alpha) values are used; the downstroke values (static hysteresis near stall) are omitted.",
    "Tunnel corrections: the paper states they are negligible (sec. 3.1.2); the tabulated values are as given.",
    "Note: this is the same DTU facility/model as Sarmast & Mikkelsen (2012/13) 'The experimental results of the NREL",
    "  S826 airfoil at low Reynolds numbers' (the report itself could not be downloaded); see case.json notes.",
], P, float_format="%.4g")

# (2) NTNU, Bartl et al. 2019 Fig. 4 (digitised), alternative
D = pd.read_csv(os.path.join(HERE, "bartl2019_fig4_raw.csv"))
D = D[(D.alpha_deg >= -10) & (D.alpha_deg <= 20)].copy()
D["source"] = "exp NTNU Bartl 2019 Fig.4 digitised; cd " + np.where(D.alpha_deg >= 13, "force balance", "wake rake")
D.loc[(D.cl_digitised == "occluded") | (D.cd_digitised == "occluded"), "source"] += " (marker hidden: mean of overlying curves)"
write_csv(os.path.join(CASE, "airfoils", "s826_polars_ntnu_bartl2019.csv"), [
    "DEFAULT polar set of the case: NREL S826 measured at NTNU (closed-return tunnel 2.7 x 1.8 m), 2D wing",
    "  c = 0.45 m between end dummies, TI(Re=70k) = 0.7 %. Re = 0.5e5 ... 6.0e5.",
    "Source: Bartl, Sagmo, Bracchi, Saetran (2019) Eur. J. Mech. B/Fluids 75:180-192, doi:10.1016/j.euromechflu.2018.10.002,",
    "  Fig. 4 (a) lift, (b) drag - DIGITISED from the raster figure (script extraction/digitise_bartl2019_fig4.py,",
    "  overlay extraction/overlay_bartl2019_fig4.png). cd: wake rake for alpha < 13 deg, force balance for alpha >= 13 deg",
    "  (as recommended by the paper p.16). Stored range -10..20 deg (curves are crowded beyond).",
    "Uncertainty from digitising: cl +-0.01 (visible markers), +-0.03 (hidden markers); cd +-0.002 below 12 deg,",
    "  +-0.01..0.03 between 13 and 20 deg where 16 lines overlap. Paper: lift uncertainty up to 3.5 %, drag up to 20 %.",
    "Re = 0.5e5 shows laminar-separation effects (lift collapse 2..8 deg).",
], D[["re", "alpha_deg", "cl", "cd", "source"]], float_format="%.4g")

# (3) organisers' XFOIL set, Table D.1 BT5 invitation (also distributed in BT4)
blk = t[t.index("D.1 NREL S826 airfoil"):t.index("Table D.1:")]
x = []
for line in blk.splitlines():
    p = line.split()
    if len(p) == 8:
        try:
            v = [num(q) for q in p]
        except ValueError:
            continue
        x.append(v[:3]); x.append(v[4:7])
X = pd.DataFrame(x, columns=["alpha_deg", "cl", "cd"]).sort_values("alpha_deg")
X.insert(0, "re", 1.0e5)
X["source"] = "XFOIL (organisers BT4/BT5 invitation Table D.1)"
write_csv(os.path.join(CASE, "airfoils", "s826_polars_bt5_xfoil_re1e5.csv"), [
    "ALTERNATIVE polar set: the 'standard' S826 polar distributed by the blind-test organisers, computed with XFOIL",
    "  at Re = 1.0e5 (Ncrit not stated in the document) for -40..50 deg.",
    "Source: Saetran et al. (2017) Invitation to the Blind test 5 workshop, Appendix D, Table D.1 (tabulated),",
    "  https://doi.org/10.5281/zenodo.1218555 . Values beyond ~20 deg are XFOIL post-stall and not reliable.",
], X, float_format="%.5g")

# ------------------------------------------------------------------ experiments
EXP = os.path.join(CASE, "experiment")
R_T2 = 0.447
R_KL = 0.45


def to_case_R(df, R_src=R_KL, R_case=R_T2):
    """Source normalised with R_src -> the case's fixed definitions with R_case (README: one column dictionary).
    tsr = Omega*R/v scales with R, cp and ct with 1/R^2. The source values are kept as tsr_src/cp_src/ct_src."""
    for c in ("tsr", "cp", "ct", "cp_err", "ct_err"):
        if c in df.columns:
            df[c + "_src"] = df[c]
    k = R_case / R_src
    df["tsr"] = (df["tsr_src"] * k).round(4)
    for c in ("cp", "ct", "cp_err", "ct_err"):
        if c in df.columns:
            df[c] = (df[c + "_src"] / k ** 2).round(5)
    return df


CONV_NOTE = ("CONVERTED to this case's R = 0.447 m (D = 0.894 m): tsr = tsr_src*0.447/0.45, cp = cp_src*(0.45/0.447)^2,"
             " ct = ct_src*(0.45/0.447)^2 (factor 1.0135); tsr_src/cp_src/ct_src = values as published (R = 0.45 m).")


def rpm_of(tsr, v, R):
    return np.round(tsr * v / R * 60 / (2 * np.pi), 1)




def merge_cp_ct(F, by, ct_scale=1.0):
    """Pairs the CP and CT markers of the same operating point (same series, tsr within 0.06).
    ct_scale: the CT panel's tsr is divided by it first (see the Fig. 5 note); the CP tsr is kept when paired."""
    cp = F[F.panel == "cp"].rename(columns={"value": "cp", "tsr": "tsr_cp"}).drop(columns="panel")
    ct = F[F.panel == "ct"].rename(columns={"value": "ct", "tsr": "tsr_ct"}).drop(columns="panel")
    ct["tsr_ct"] = (ct["tsr_ct"] / ct_scale).round(3)
    out = []
    for key, g in cp.groupby(by):
        h = ct[ct[by] == key]
        used = set()
        for _, r in g.iterrows():
            d = (h.tsr_ct - r.tsr_cp).abs()
            j = d.idxmin() if len(d) else None
            if j is not None and d[j] < 0.06 and j not in used:
                used.add(j)
                out.append({by: key, "tsr": round((r.tsr_cp + h.tsr_ct[j]) / 2, 3) if ct_scale == 1.0 else r.tsr_cp,
                            "cp": r.cp, "ct": h.ct[j]})
            else:
                out.append({by: key, "tsr": r.tsr_cp, "cp": r.cp, "ct": np.nan})
        for j, r in h.iterrows():
            if j not in used:
                out.append({by: key, "tsr": r.tsr_ct, "cp": np.nan, "ct": r.ct})
    return pd.DataFrame(out).sort_values([by, "tsr"]).reset_index(drop=True)

# (1) Krogstad & Lund 2012 Table II
blk = kl[kl.index("Table II. Measured power and thrust coefficients."):]
blk = blk[:blk.index("was, the case")]
tab = []
for line in blk.splitlines():
    p = line.split()
    if len(p) == 3:
        try:
            tab.append([num(p[0]), num(p[1]), num(p[2])])
        except ValueError:
            pass
T = pd.DataFrame(tab, columns=["tsr", "cp", "ct"])
assert len(T) == 22
T["v_ms"] = 10.0
T["rpm"] = rpm_of(T.tsr, 10.0, R_KL)
T["pitch_deg"] = 0.0
T["cp_err"] = (0.04 * T.cp.abs()).round(4)
T["ct_err"] = (0.02 * T.ct).round(4)
T = to_case_R(T)
write_csv(os.path.join(EXP, "performance_krogstad_lund2012_tableII.csv"), [
    "Measured CP and CT vs TSR, NTNU model turbine, pitch 0, yaw 0, uniform inflow TI ~0.3 %, U_inf ~10 m/s (nominal).",
    "Source: Krogstad & Lund (2012) Wind Energy 15:443-457, doi:10.1002/we.482, Table II (TABULATED, 22 rows copied exactly;",
    "  last row CP = -0.013258). The paper normalises with D = 0.9 m (R = 0.45 m); rpm = tsr_src*U/0.45 (measured speed).",
    CONV_NOTE,
    "CT: rotor thrust from 6-component balance with the drag of the model without blades (CD = 0.137) subtracted.",
    "NOT corrected for wind-tunnel blockage (closed test section 2.7 m x 2.0 m, rotor area 11.8 % of the section).",
    "Uncertainty: CT 0.5 % at design TSR rising to ~2 % at the lowest TSR (K&L 2012); CP +-4 % (as given for the same",
    "  rig in Krogstad, Saetran & Adaramola 2015, JFS 52:65-80). The *_err columns use 4 % CP and 2 % CT.",
], T, float_format="%.6g")

# (2) Krogstad & Eriksen 2013 Fig. 3 measured points (vector-exact)
F = pd.read_csv(os.path.join(HERE, "krogstad_eriksen2013_fig3_raw.csv"))
F["series"] = 0
M = merge_cp_ct(F, "series").drop(columns="series")
assert M.cp.notna().all() and M.ct.notna().all() and len(M) == 28
M["v_ms"] = 10.0
M["rpm"] = rpm_of(M.tsr, 10.0, R_T2)
M["pitch_deg"] = 0.0
M["cp_err"] = (0.04 * M.cp.abs()).round(4)
M["ct_err"] = (0.02 * M.ct).round(4)
write_csv(os.path.join(EXP, "performance_bt1_krogstad_eriksen2013_fig3.csv"), [
    "Measured CP and CT vs TSR of the BT1 turbine (NTNU T2, D = 0.894 m), pitch 0, yaw 0, U_ref = 10 m/s, TI 0.3 %.",
    "Source: Krogstad & Eriksen (2013) Renewable Energy 50:325-333, doi:10.1016/j.renene.2012.06.044, Fig. 3 (a),(b),",
    "  'Measured' series (large black dots). Extracted EXACTLY from the vector PDF drawing (PyMuPDF marker centres,",
    "  script extraction/extract_krogstad_eriksen2013_fig3.py, overlay extraction/overlay_krogstad_eriksen2013_fig3.png);",
    "  extraction error < 0.005 in tsr and < 0.001 in CP/CT. Normalised with D = 0.894 m; rpm = tsr*U/0.447.",
    "CT is rotor thrust (tower/nacelle drag subtracted). Not blockage-corrected (rotor+tower ~12 % of the section).",
    "Note: differs from Krogstad & Lund (2012) Table II by up to +0.05 in CT and +0.016 in CP (other test series / D).",
    "*_err: 4 % CP, 2 % CT (measurement uncertainty stated for this rig in Krogstad et al. 2015, JFS 52:65-80).",
], M, float_format="%.6g")

# (3) Reynolds-number dependence, Krogstad & Adaramola 2012 Fig. 4
F = pd.read_csv(os.path.join(HERE, "krogstad_adaramola2012_fig4_raw.csv"))
M = merge_cp_ct(F, "v_ms")
M["rpm"] = rpm_of(M.tsr, M.v_ms, R_KL)
M["pitch_deg"] = 0.0
M = to_case_R(M[["tsr", "cp", "ct", "v_ms", "rpm", "pitch_deg"]])
# CT of this figure lies 0.10-0.17 ABOVE the rotor CT of BT1 Fig. 3 and K&L Table II at every TSR, also at
# U = 10.3 m/s (CT = 0.28 vs 0.15 at TSR 1) - the offset matches the tower+nacelle drag (CD = 0.137 on the rotor
# area, K&L 2012). It is therefore kept only as ct_incl_tower (not the dictionary's rotor 'ct').
M = M.rename(columns={"ct": "ct_incl_tower"})
write_csv(os.path.join(EXP, "performance_re_dependence_krogstad_adaramola2012_fig4.csv"), [
    "Reynolds-number dependence: measured CP and CT vs TSR at U_inf = 7.22, 7.72, 8.24, 8.76, 9.33, 9.75, 10.30 m/s.",
    "NTNU model turbine (the paper states D = 0.90 m and normalises with it; rpm = tsr_src*U/0.45), pitch 0, yaw 0.",
    CONV_NOTE,
    "Source: Krogstad & Adaramola (2012) Wind Energy 15:743-756, doi:10.1002/we.502, Fig. 4 (a),(b). Extracted EXACTLY",
    "  from the vector PDF drawing (series told apart by marker shape; script extraction/extract_krogstad_adaramola2012_fig4.py,",
    "  overlay extraction/overlay_krogstad_adaramola2012_fig4.png); extraction error < 0.01 in tsr, < 0.001 in CP/CT.",
    "Empty cp: only at the highest TSR of each series, where CP < 0 lies below the plotted range (0..0.5).",
    "Paper: CP, CT independent of U for U > 9 m/s. Not blockage-corrected (12.4 %).",
    "!! ct_incl_tower: the published CT of this figure is 0.10-0.17 higher than the rotor CT of the same rotor in",
    "  Krogstad & Eriksen 2013 Fig. 3 and Krogstad & Lund 2012 Table II at all TSR (0.28 vs 0.15 at TSR 1, U = 10.3 m/s);",
    "  the offset equals the tower+nacelle drag (CD = 0.137, K&L 2012), so this CT very probably INCLUDES the tower and",
    "  nacelle drag although the paper's text says it was subtracted. Hence it is not stored as rotor 'ct'; use it only",
    "  for the Reynolds-number trend (differences between speeds). CP agrees with Fig. 3 within ~0.03.",
], M, float_format="%.6g")

# (4) yaw vs TSR, Krogstad & Adaramola 2012 Fig. 5
F = pd.read_csv(os.path.join(HERE, "krogstad_adaramola2012_fig5_raw.csv"))
# The CT panel (c) is plotted at tsr values 1.0240 times those of the CP and Cq panels (a),(b) for every series
# (median ratio of matched markers, scatter +-0.001) - same measurement points, different tsr scaling in the figure.
# The CP/Cq tsr is kept (it agrees with Fig. 4 and with the runaway TSR 11.4 quoted in the text).
CT_SCALE = 1.0240
M = merge_cp_ct(F, "yaw_deg", ct_scale=CT_SCALE)
M["v_ms"] = 10.0
M["rpm"] = rpm_of(M.tsr, 10.0, R_KL)
M["pitch_deg"] = 0.0
M = to_case_R(M[["tsr", "cp", "ct", "yaw_deg", "v_ms", "rpm", "pitch_deg"]])
write_csv(os.path.join(EXP, "performance_yaw_krogstad_adaramola2012_fig5.csv"), [
    "Yawed rotor: measured CP and CT vs TSR at yaw angles 0, 10, 20, 30, 40, 50 deg; pitch 0; U_inf ~10 m/s.",
    "NTNU model turbine, D = 0.90 m as stated in the paper (coefficients normalised with the full rotor area pi*D^2/4",
    "  and U_inf, not with the projected area); rpm = tsr_src*U/0.45.",
    CONV_NOTE,
    "Source: Krogstad & Adaramola (2012) Wind Energy 15:743-756, doi:10.1002/we.502, Fig. 5 panels (a) CP and (c) CT",
    "  (the caption's panel order is wrong; the axis titles are used). Extracted EXACTLY from the vector PDF drawing",
    "  (series by marker style; extraction/extract_krogstad_adaramola2012_fig5.py, overlay_krogstad_adaramola2012_fig5.png);",
    "  extraction error < 0.01 in tsr, < 0.001 in CP/CT. Empty cp/ct: marker missing in that panel.",
    "!! In the PDF the CT panel's markers sit at tsr values exactly 1.024x those of the CP and Cq panels (all six",
    "  series, scatter +-0.001), although they are the same operating points. tsr here = the CP-panel value (CT tsr",
    "  divided by 1.024); this choice agrees with the paper's Fig. 4. If the CT panel is right instead, tsr is 2.4 % low.",
    "Thrust: the paper defines CT = 2D/(rho pi R^2 U^2) with D the 'drag force' from the force balance, tower+nacelle",
    "  drag subtracted. Whether D is along the wind or along the (yawed) rotor axis is not stated - treat CT at large",
    "  yaw with care (factor cos(yaw) between the two). Not blockage-corrected (12.4 %).",
], M, float_format="%.6g")

# (5) yaw sweep at TSR 6, Bartl et al. 2018b Fig. 1 (inflow A) + Bartl et al. 2018a Table 2
F = pd.read_csv(os.path.join(HERE, "bartl2018b_fig1_raw.csv"))
F = F[F.inflow == "A"]
cp = F[F.panel == "cp"][["yaw_deg", "value"]].rename(columns={"value": "cp"})
ct = F[F.panel == "ct"][["yaw_deg", "value"]].rename(columns={"value": "ct"})
M = pd.merge(cp, ct, on="yaw_deg", how="outer").sort_values("yaw_deg")
M.insert(0, "tsr", 6.0)
M["v_ms"] = 10.0
M["rpm"] = rpm_of(6.0, 10.0, R_T2)
M["pitch_deg"] = 0.0
M["cp_err"] = 0.017
M["ct_err"] = 0.018
write_csv(os.path.join(EXP, "yaw_sweep_bartl2018_tsr6_fig1.csv"), [
    "Yaw sweep at TSR = 6.0 (x = yaw_deg): measured CP and CT of the NTNU rotor (D = 0.894 m, LARS1/T1 rig of",
    "  Bartl et al. 2018), pitch 0, U_ref = 10 m/s, inflow A = uniform, TI = 0.23 %. Positive yaw as defined in the paper.",
    "Source: Bartl, Muehle, Schottler, ... Saetran (2018) Wind Energ. Sci. 3:489-502, doi:10.5194/wes-3-489-2018,",
    "  Fig. 1 (a),(b), red series 'Inflow A'. DIGITISED from the raster figure (extraction/digitise_bartl2018b_fig1.py,",
    "  overlay extraction/overlay_bartl2018b_fig1.png); digitising error +-0.004. Empty value: marker hidden by",
    "  the inflow-B marker (cp at +40, ct at -40, -30, +10 deg); for those angles inflow A ~ inflow B within 0.005.",
    "*_err = measurement uncertainty stated in the paper (CP 0.011..0.017, CT 0.013..0.018; the larger values used).",
    "Cross-check Bartl et al. (2018a) Table 2 (file yaw_bartl2018a_table2_inflowA.csv): yaw 0: CP 0.468, CT 0.893 (other test",
    "  series and D = 0.90 m normalisation) vs 0.467 / 0.833 here.",
    "Thrust = force normal to the rotor plane, tower/nacelle drag subtracted. Not blockage-corrected (12.8 %).",
], M, float_format="%.6g")

w = txt("wes-3-329-2018.txt")
blk = w[w.index("Table 2. Turbine performance"):w.index("2.2   Measurement techniques")]
rows = []
for line in blk.splitlines():
    p = line.split()
    if len(p) == 7 and re.match(r"^[+\-−]?\d+$", p[0]):
        y = num(p[0].replace("+", ""))
        v = [num(q) for q in p[1:]]
        for lab, ti, (a, b) in (("A", 0.23, v[0:2]), ("B", 10.0, v[2:4]), ("C", 10.0, v[4:6])):
            rows.append(dict(tsr=6.0, yaw_deg=y, cp=a, ct=b, inflow=lab, ti_percent=ti))
M = pd.DataFrame(rows).sort_values(["inflow", "yaw_deg"])
M_all_note = "Values AS PUBLISHED (normalised with D = 0.90 m); this file is not listed in case.json."
M["v_ms"] = 10.0
M["rpm"] = rpm_of(6.0, 10.0, R_KL)
M["pitch_deg"] = 0.0
write_csv(os.path.join(EXP, "yaw_bartl2018a_table2.csv"), [
    "CP and CT at TSR = 6.0 and yaw 0, +30, -30 deg for three inflows: A uniform TI 0.23 % (compare with the solver),",
    "  B uniform grid turbulence TI 10 %, C sheared (power-law alpha = 0.11) TI 10 %. Pitch 0, U_ref = 10 m/s.",
    "Source: Bartl, Muehle, Schottler, Saetran, Peinke, Adaramola, Hoelling (2018) Wind Energ. Sci. 3:329-343,",
    "  doi:10.5194/wes-3-329-2018, Table 2 (TABULATED). Rotor D = 0.90 m as stated there (rpm = 6*U/0.45).",
    "Only the inflow-A rows correspond to the uniform low-turbulence inflow the solver assumes.",
    M_all_note,
], M[["tsr", "yaw_deg", "cp", "ct", "inflow", "ti_percent", "v_ms", "rpm", "pitch_deg"]], float_format="%.6g")
os.replace(os.path.join(EXP, "yaw_bartl2018a_table2.csv"), os.path.join(EXP, "yaw_bartl2018a_table2_all_inflows.csv"))
A = to_case_R(M[M.inflow == "A"][["tsr", "yaw_deg", "cp", "ct", "v_ms", "rpm", "pitch_deg"]].copy())
write_csv(os.path.join(EXP, "yaw_bartl2018a_table2_inflowA.csv"), [
    "CP and CT at TSR = 6.0, yaw 0 / +30 / -30 deg, uniform low-turbulence inflow A (TI 0.23 %), pitch 0, U_ref = 10 m/s.",
    "Source: Bartl et al. (2018) Wind Energ. Sci. 3:329-343, doi:10.5194/wes-3-329-2018, Table 2, inflow-A columns",
    "  (TABULATED). Rotor D = 0.90 m as stated there (rpm = 6*U/0.45). " + CONV_NOTE,
    "Uncertainty CP +-0.011..0.017, CT +-0.013..0.018",
    "  (same rig, Bartl et al. 2018, WES 3:489). Thrust = force normal to the rotor plane, tower/nacelle drag subtracted.",
    "The inflow B (TI 10 %) and C (shear) rows are in yaw_bartl2018a_table2_all_inflows.csv (not listed in case.json).",
], A, float_format="%.6g")

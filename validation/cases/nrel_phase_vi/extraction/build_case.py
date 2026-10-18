"""Build the NREL UAE Phase VI validation case files from the sources and the raw extractions.

Inputs (all in ../sources, git-ignored, or in this folder):
  29955.txt          pdftotext -layout of Hand et al. 2001, NREL/TP-500-29955 (Table A-1 blade geometry, Table A-2 S809)
  osu.txt            pdftotext -layout of Ramsay et al. 1995, NREL/TP-442-7817 (Appendix B, Tables B1-B4 clean)
  Selig_PROPID_wtprep_airfoil1_S809_DUT.html   S809 Delft (DUT) tabulated polars, Re = 1.0 ... 3.0e6
  ahlund.txt         pdftotext of Åhlund 2004 FOI-R--1243--SE (Table 27: DUT Re=1e6, used as cross-check)
  simms2001_measured_raw.csv, sorensen2002_measured_raw.csv, lindenburg2003_measured_raw.csv  (this folder;
                     produced by extract_simms2001.py, extract_sorensen2002.py, extract_lindenburg2003.py)
Outputs: ../geometry.csv, ../airfoils/*, ../experiment/*.csv
Run: python3 build_case.py
"""
import os
import re

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.dirname(HERE)
SRC = os.path.join(CASE, "sources")
R = 5.029
RPM = 72.0
PITCH = 3.0          # tip pitch of sequence S, deg
R_GAUGE = 0.432      # root strain gauges, m from the rotor axis
OMEGA = RPM * 2 * np.pi / 60

FLOAT = re.compile(r"-?\d+\.\d+|-?\d+")


def write_csv(path, header_lines, df, float_format=None):
    with open(path, "w") as f:
        for h in header_lines:
            f.write("# " + h + "\n")
        df.to_csv(f, index=False, float_format=float_format, lineterminator="\n")


# ----------------------------------------------------------------------------------------------- geometry
def hand_table_a1():
    txt = open(os.path.join(SRC, "29955.txt"), encoding="utf-8", errors="replace").read().splitlines()
    start = next(i for i, l in enumerate(txt) if "Table A-1. Blade chord and twist distributions" in l
                 and "...." not in l)
    rows = []
    for l in txt[start + 1:start + 60]:
        s = l.strip()
        if not s or not re.match(r"^\d+\.\d+\s", s):
            continue
        f = [float(x) for x in FLOAT.findall(s)]
        if len(f) < 6:
            continue
        r, chord, twist, thick, axis = f[0], f[3], f[4], f[5], f[6] if len(f) > 6 else np.nan
        rows.append(dict(r_m=r, chord_m=chord, twist_hand_deg=twist, thick=thick, twist_axis_pct=axis,
                         thickness_is_pct="%" in s))
        if r >= 5.531:
            break
    return pd.DataFrame(rows)


def build_geometry():
    t = hand_table_a1()
    t.to_csv(os.path.join(HERE, "hand2001_tableA1_parsed.csv"), index=False)
    # twist and chord at the standard-tip radius R = 5.029 m by linear interpolation between 5.000 and 5.305 m
    tw_tip = float(np.interp(R, t.r_m, t.twist_hand_deg))
    ch_tip = float(np.interp(R, t.r_m, t.chord_m))
    g = t[t.r_m < R].copy()
    g = pd.concat([g, pd.DataFrame([dict(r_m=R, chord_m=ch_tip, twist_hand_deg=tw_tip, thick=20.95,
                                         thickness_is_pct=True)])], ignore_index=True)
    # thickness ratio
    tc = []
    foil = []
    for _, row in g.iterrows():
        if row.thickness_is_pct:
            tc.append(0.2095)
            foil.append("s809")
        else:
            tc.append(row.thick / row.chord_m)
            if row.r_m <= 0.8836:
                foil.append("circle")
            elif row.r_m < 1.257:
                foil.append("transition")
            else:
                foil.append("s809")
    g["thickness_ratio"] = np.round(tc, 4)
    g["airfoil"] = foil
    g.loc[g.r_m >= 1.2574, "airfoil"] = "s809"
    g.loc[(g.r_m >= 1.2574), "thickness_ratio"] = 0.2095
    # solver convention: twist relative to the tip chord (pitch is measured at the tip, 0 deg at tip)
    g["twist_deg"] = np.round(g.twist_hand_deg - tw_tip, 4)
    out = g[["r_m", "chord_m", "twist_deg", "airfoil", "thickness_ratio"]].copy()
    out["chord_m"] = out.chord_m.round(4)
    hdr = [
        "NREL UAE Phase VI blade (standard tip, R = 5.029 m), sequence S configuration.",
        "Source: Hand et al. 2001, NREL/TP-500-29955, Appendix A, Table A-1 'Blade chord and twist distributions'",
        "(p. 64 of the report), TABULATED; rows parsed from the PDF text by extraction/build_case.py.",
        f"Last row r = {R} m (tip of the standard tip piece) interpolated linearly between the 5.000 and 5.305 m rows:",
        f"chord {ch_tip:.4f} m, Hand twist {tw_tip:.4f} deg (OpenFAST r-test UAE_VI uses -1.815 deg there).",
        "TWIST CONVERSION: Hand gives twist positive towards feather, zero at the 3.772 m station (75% span of a",
        "5.029-m blade); the test pitch angle is the angle of the TIP chord (Hand Fig. 19 'pitch angle convention,",
        f"0 deg at tip'). Here twist_deg = twist_Hand - twist_Hand(r=R) = twist_Hand + {-tw_tip:.4f}, i.e. relative to the",
        "tip chord, so that local blade angle = twist_deg + pitch_deg with pitch_deg = 3.0 for sequence S.",
        "(Equivalent: OpenFAST r-test keeps Hand's twist and sets BlPitch = 4.815 deg.) Sign: positive towards",
        "feather = the solver's turbine convention (alpha = phi - (twist + pitch)).",
        "Root: blade attaches to the hub at r = 0.508 m; cylinder 0.508-0.883 m (dia. 0.218 m tapering to 0.183 m;",
        "the twist values given there are irrelevant for a circle and are only shifted like the rest);",
        "cylinder-to-S809 transition 0.883-1.257 m (rows 1.0085, 1.0675, 1.1335 m are 'transition'; their thickness",
        "ratio = Hand's thickness [m] / chord). S809 (t/c = 20.95%) from 1.2575 m to the tip.",
        "Pitch/twist axis at 30% chord outboard of 1.257 m (35.9-31.9% in the transition, 50% on the cylinder).",
        "Columns: r_m [m], chord_m [m], twist_deg [deg], airfoil, thickness_ratio [-]",
    ]
    write_csv(os.path.join(CASE, "geometry.csv"), hdr, out)
    return out, tw_tip


# ----------------------------------------------------------------------------------------------- airfoils
def s809_coordinates():
    txt = open(os.path.join(SRC, "29955.txt"), encoding="utf-8", errors="replace").read().splitlines()
    start = next(i for i, l in enumerate(txt) if l.strip().startswith("Table A-2. Airfoil profile coordinates")
                 and "...." not in l)
    up, lo = [], []
    for l in txt[start + 1:start + 60]:
        f = FLOAT.findall(l)
        if len(f) == 4:
            a, b, c, d = map(float, f)
            up.append((a, b))
            lo.append((c, d))
        if "Source" in l:
            break
    # last line of Table A-2 is "1.00000 0.00000 | 0.00000 0.00000": upper TE and the leading-edge point
    le = lo.pop()           # (0, 0)
    assert le == (0.0, 0.0), le
    upper = sorted(up, key=lambda p: -p[0])        # TE -> LE
    lower = sorted(lo, key=lambda p: p[0])         # LE -> TE
    pts = upper + [le] + lower
    with open(os.path.join(CASE, "airfoils", "s809.dat"), "w") as f:
        f.write("NREL S809 (Hand et al. 2001 NREL/TP-500-29955 Table A-2 = Somers 1997 NREL/SR-440-6918 Table 2)\n")
        for x, y in pts:
            f.write(f"{x:9.5f} {y:9.5f}\n")
    # cross-check with Somers 1997 Table 2 (OCR text): compare a few values
    return np.array(pts)


def circle_files():
    th = np.linspace(0, 2 * np.pi, 73)
    with open(os.path.join(CASE, "airfoils", "circle.dat"), "w") as f:
        f.write("circle (root cylinder of the Phase VI blade, unit diameter)\n")
        for t in th:
            f.write(f"{0.5 + 0.5 * np.cos(t):9.5f} {0.5 * np.sin(t):9.5f}\n")
    rows = []
    for re_ in (1e5, 5e5):
        for a in np.arange(-180, 181, 10):
            rows.append(dict(re=re_, alpha_deg=float(a), cl=0.0, cd=1.0, cm=0.0,
                             source="assumed: circular cylinder, subcritical Re (~1e5-4e5), cd=1.0 as Lindenburg 2003 App. A.2.3"))
    hdr = ["Root cylinder of the Phase VI blade (r = 0.508-0.883 m, diameter 0.218-0.183 m).",
           "NOT measured: constant cd = 1.0, cl = 0, cm = 0 at all angles. Cylinder Reynolds number V_rel*d/nu ~ 1e5-4e5",
           "(subcritical; smooth-cylinder cd ~ 1.0-1.2). Lindenburg 2003 (ECN-C--03-025, App. A.2.3) used cd = 1.0",
           "(range 0.35 rough ... 1.2 smooth); the OpenFAST r-test cylinder.dat uses cd = 0.3. Choice here: 1.0."]
    write_csv(os.path.join(CASE, "airfoils", "circle_polars.csv"), hdr, pd.DataFrame(rows))


def dut_polars():
    txt = open(os.path.join(SRC, "Selig_PROPID_wtprep_airfoil1_S809_DUT.html")).read()
    lines = [l.strip() for l in re.sub(r"<[^>]+>", "", txt).splitlines() if l.strip()]
    i = 2
    nre = int(lines[i]); i += 1
    rows = []
    for _ in range(nre):
        re_ = float(lines[i]); i += 1
        n = int(lines[i]); i += 1
        tab = []
        for _ in range(n):
            cl, cd, a = map(float, lines[i].split()); i += 1
            tab.append((a, cl, cd))
        # The PROPID file pads every table with synthetic post-stall points (cl held constant at the last value,
        # cd = 0.1 / 0.175 / 0.275 / 0.363 at alpha = 16 / 20 / 25 / 27.5 deg). Those are not measured: drop them.
        pad = {(16.0, 0.1), (20.0, 0.175), (25.0, 0.275), (27.5, 0.363)}
        for a, cl, cd in tab:
            if (a, cd) in pad:
                continue
            rows.append(dict(re=re_, alpha_deg=a, cl=cl, cd=cd,
                             source="measured, TU Delft LST (Somers 1997 NREL/SR-440-6918); table from Selig PROPID/WTPREP airfoil1"))
    df = pd.DataFrame(rows)
    # cross-check Re = 1e6 with Åhlund 2004 Table 27
    t = open(os.path.join(SRC, "ahlund.txt"), encoding="utf-8", errors="replace").read().splitlines()
    s = next(k for k, l in enumerate(t) if "Table 27 - Wind tunnel profile coefficients from DUT" in l)
    ah = []
    for l in t[s + 2:s + 30]:
        f = FLOAT.findall(l)
        if len(f) == 3:
            ah.append(tuple(map(float, f)))
    ah = pd.DataFrame(ah, columns=["alpha_deg", "cl", "cd"])
    d1 = df[df.re == 1e6].reset_index(drop=True)
    m = []
    for _, r in ah.iterrows():
        k = (d1.alpha_deg - r.alpha_deg).abs().idxmin()
        if abs(d1.alpha_deg[k] - r.alpha_deg) < 0.05:
            m.append((r.alpha_deg, d1.alpha_deg[k], r.cl - d1.cl[k], r.cd - d1.cd[k]))
    m = np.array(m)
    print(f"DUT 1e6 vs Åhlund Table 27: {len(m)} matched points, max |dalpha| {np.max(np.abs(m[:, 0] - m[:, 1])):.3f},"
          f" max |dcl| {np.max(np.abs(m[:, 2])):.4f}, max |dcd| {np.max(np.abs(m[:, 3])):.4f}")
    hdr = [
        "S809 MEASURED 2-D polars, TU Delft Low Speed Laboratory low-turbulence tunnel (Somers 1997, NREL/SR-440-6918,",
        "tests 1986; the report itself only plots them, Fig. 14). TABULATED values taken from M. Selig's PROPID/WTPREP",
        "data set 'airfoil1: S809 Experimental data from Delft wind tunnel'",
        "(https://m-selig.ae.illinois.edu/propid/raj/wtprep/polardata-html/airfoil1.html), Re = 1.0, 1.5, 2.0, 2.5, 3.0e6.",
        "The synthetic post-stall padding points of that file (cl held, cd = 0.1/0.175/0.275/0.363 at 16/20/25/27.5 deg)",
        "were removed: only measured angles are kept (Re 1e6: -3.09 ... 17.23 deg).",
        "Cross-check: Re = 1e6 agrees with Åhlund 2004 (FOI-R--1243--SE) Table 27 (DUT, Re = 1e6): "
        f"{len(m)} points, |dalpha| <= {np.max(np.abs(m[:, 0] - m[:, 1])):.2f} deg, |dcl| <= {np.max(np.abs(m[:, 2])):.3f},"
        f" |dcd| <= {np.max(np.abs(m[:, 3])):.4f}.",
        "Åhlund also lists cl at 18.19-20.16 deg with cd frozen at 0.1853 (not measured) - not included.",
        "Clean (free transition). cd = wake-rake profile drag. No cm in this table.",
        "Rotor Re at the pressure stations ~0.7-1.1e6 (Lindenburg 2003 App. A); the lowest DUT Re is 1.0e6.",
    ]
    write_csv(os.path.join(CASE, "airfoils", "s809_polars.csv"), hdr, df)
    return df


def osu_polars():
    t = open(os.path.join(SRC, "osu.txt"), encoding="utf-8", errors="replace").read().splitlines()
    nominal = {"B1": 0.75e6, "B2": 1.0e6, "B3": 1.25e6, "B4": 1.5e6}
    rows = []
    for tab, re_nom in nominal.items():
        # data pages start after the list of tables; the table header contains "Table B1. S809, Clean"
        idx = [k for k, l in enumerate(t) if re.search(rf"Table {tab}\. S809, Clean", l) and "...." not in l]
        seen = set()
        for k0 in idx:
            for l in t[k0 + 1:k0 + 45]:
                if f"End of Table {tab}" in l or (re.search(r"Table B\d\.", l) and tab not in l):
                    break
                f = l.split()
                if len(f) == 7 and re.match(r"^\d+$", f[0]):
                    run = int(f[0])
                    if run in seen:
                        continue
                    seen.add(run)
                    a, cl, cdp, cm, rea = map(float, f[1:6])
                    cdw = float(f[6]) if f[6] != "--" else np.nan
                    rows.append(dict(re=re_nom, alpha_deg=a, cl=cl, cd=cdw if cdw == cdw else cdp, cm=cm,
                                     cdp=cdp, cdw=cdw, re_actual=rea * 1e6, run=run,
                                     source=f"measured, OSU 3x5 tunnel, Ramsay et al. 1995 NREL/TP-442-7817 Table {tab} (clean); "
                                            + ("cd = Cdw (wake rake)" if cdw == cdw else "cd = Cdp (pressure drag; no wake data)")))
    df = pd.DataFrame(rows).sort_values(["re", "alpha_deg"]).reset_index(drop=True)
    print("OSU points per Re:", df.groupby("re").size().to_dict())
    hdr = [
        "S809 MEASURED 2-D polars, Ohio State University AARL 3x5-ft tunnel, steady, CLEAN model:",
        "Ramsay, Hoffmann & Gregorek 1995 (rev. 12/1999), NREL/TP-442-7817, Appendix B, Tables B1-B4, TABULATED",
        "(parsed from the PDF text by extraction/build_case.py). Nominal Re 0.75, 1.0, 1.25, 1.5e6 in column 're';",
        "the actual Re of each run is in 're_actual'.",
        "cd rule: wake-rake drag Cdw where the report gives it (|alpha| up to ~10 deg), otherwise pressure drag Cdp",
        "(skin friction missing, so cd is somewhat low just past 10 deg; in deep stall Cdp ~ total drag). cdp and",
        "cdw columns keep both. cm = quarter-chord moment from pressures.",
        "Known issue: OSU vs Delft data differ by an angle-of-attack offset of roughly 0.5-1 deg (Lindenburg 2003",
        "shifted OSU by -0.53 deg; NREL's AeroDyn files of the OpenFAST r-test UAE_VI shift it by -1 deg). NOT shifted here.",
        "Alternative to the Delft set (s809_polars.csv); reaches -20...+40 deg (deep stall), useful because the",
        "Phase VI inboard sections run far beyond 20 deg at high wind speed.",
    ]
    write_csv(os.path.join(CASE, "airfoils", "s809_osu_polars.csv"), hdr, df)
    return df


# ----------------------------------------------------------------------------------------------- experiment
LINDENBURG_RHO = {7.0: 1.246, 10.0: 1.246, 13.1: 1.227, 15.1: 1.224, 20.1: 1.221, 25.1: 1.220}
NOMINAL = [7, 10, 13, 15, 20, 25]
LIND_V = {7: 7.0, 10: 10.0, 13: 13.1, 15: 15.1, 20: 20.1, 25: 25.1}


def experiments():
    sim = pd.read_csv(os.path.join(HERE, "simms2001_measured_raw.csv"))
    sor = pd.read_csv(os.path.join(HERE, "sorensen2002_measured_raw.csv"))
    lin = pd.read_csv(os.path.join(HERE, "lindenburg2003_measured_raw.csv"))
    lin = lin[~lin.legend]

    def simms(q):
        d = sim[sim.quantity == q].copy()
        d["vn"] = d.v_ms.round().astype(int)
        g = d.groupby("vn").agg(v=("v_ms", "mean"), val=("value", "mean"), spread=("value", lambda x: x.max() - x.min()),
                                scale=("value_per_pt", "first"))
        return g

    # ---- power curve: torque
    st = simms("lss_torque_nm")
    so = sor[sor.quantity == "lss_torque_nm"].set_index(sor[sor.quantity == "lss_torque_nm"].x.round().astype(int))
    lt = lin[lin.quantity == "lss_torque_knm"].copy()
    rows = []
    for vn in NOMINAL:
        v = round(float(st.loc[vn, "v"]), 2)
        tq = float(st.loc[vn, "val"])
        lv = lt.iloc[(lt.v_ms - v).abs().argmin()]
        rows.append(dict(v_ms=v, v_nominal=vn, campaign=f"S{vn:02d}00000", rpm=RPM, pitch_deg=PITCH, yaw_deg=0.0,
                         rho=LINDENBURG_RHO[LIND_V[vn]], rho_source="Lindenburg2003 p.62 table (per campaign)",
                         torque_nm=round(tq, 1),
                         torque_nm_err=round(max(3.0, abs(tq - float(so.loc[vn, "value"])), abs(tq - 1000 * lv.value)), 1),
                         torque_nm_std=round(float(so.loc[vn, "std"]), 1),
                         torque_source="Simms2001 Fig.8 (vector, 4 panels mean); checked vs Sorensen2002 Fig.3 and Lindenburg2003 Fig.6.2",
                         torque_sorensen2002_nm=round(float(so.loc[vn, "value"]), 1),
                         torque_lindenburg2003_nm=round(1000 * float(lv.value), 1)))
    main_v = [r["v_ms"] for r in rows]
    known = sorted(LINDENBURG_RHO.items())
    for _, lv in lt.iterrows():
        if min(abs(lv.v_ms - v) for v in main_v) < 0.2:
            continue
        v = round(float(lv.v_ms), 1)
        rho = float(np.interp(v, [k for k, _ in known], [r for _, r in known]))
        rows.append(dict(v_ms=v, v_nominal=int(round(v)), campaign="", rpm=RPM, pitch_deg=PITCH, yaw_deg=0.0,
                         rho=round(rho, 3), rho_source="ASSUMED: interpolated from Lindenburg2003 campaign table (not given)",
                         torque_nm=round(1000 * float(lv.value), 1), torque_nm_err=round(max(3.0, 0.5 * 1000 * float(lv.value_per_pt)), 1),
                         torque_nm_std=np.nan,
                         torque_source="Lindenburg2003 Fig.6.2 'Measured' (vector); only source for this wind speed",
                         torque_sorensen2002_nm=np.nan, torque_lindenburg2003_nm=round(1000 * float(lv.value), 1)))
    pc = pd.DataFrame(rows).sort_values("v_ms").reset_index(drop=True)
    pc["power_w"] = (pc.torque_nm * OMEGA).round(1)
    pc["power_w_err"] = (pc.torque_nm_err * OMEGA).round(1)
    pc["tsr"] = (OMEGA * R / pc.v_ms).round(4)
    pc["cp"] = (pc.power_w / (0.5 * pc.rho * pc.v_ms ** 3 * np.pi * R ** 2)).round(4)
    cols = ["v_ms", "v_nominal", "campaign", "rpm", "pitch_deg", "yaw_deg", "rho", "torque_nm", "torque_nm_err",
            "torque_nm_std", "power_w", "power_w_err", "tsr", "cp", "rho_source", "torque_source",
            "torque_sorensen2002_nm", "torque_lindenburg2003_nm"]
    hdr = [
        "NREL UAE Phase VI, NASA Ames 80x120 ft, sequence S: upwind, rigid hub, 0 deg cone, 0 deg yaw, tip pitch 3 deg,",
        "72 rpm, no five-hole probes (plugs taped). Mean LOW-SPEED SHAFT TORQUE (strain gauge LSSTQ, N m) vs tunnel speed.",
        "Sources (all measured data, read from VECTOR figures with PyMuPDF -> exact plotted coordinates, not raster):",
        "  Simms et al. 2001 NREL/TP-500-29494 Fig. 8 p.21 (6 points, value = mean of the 4 identical panel copies,",
        "    panel spread <= 2.1 N m) -> torque_nm for 7/10/13/15/20/25 m/s;",
        "  Sorensen, Michelsen & Schreck 2002 Wind Energy 5:151 Fig. 3 (same 6 points with +-1 std over a revolution)",
        "    -> torque_nm_std and the column torque_sorensen2002_nm (cross-check);",
        "  Lindenburg 2003 ECN-C--03-025 Fig. 6.2 p.63 (15 points, 6-25 m/s) -> torque_lindenburg2003_nm, and the ONLY",
        "    source for the extra speeds 6, 8, 9, 11, 12.1, 14.1, 16.1, 17.1, 18.1 m/s (assumed to be sequence S, 0 deg",
        "    yaw, as the chapter states; their campaign names and densities are not given).",
        "The three sources agree within 4 N m at the six common speeds.",
        "v_ms: tunnel speed as plotted by Simms (6 main points; Lindenburg's table rounds them to 7.0/10.0/13.1/15.1/20.1/25.1)",
        "or by Lindenburg (extra points). rho: Lindenburg 2003 p.62 campaign table (main points); interpolated for extras.",
        "torque_nm_err: +-(max of the three-source disagreement, 3 N m) = extraction/consistency uncertainty, NOT the",
        "measurement uncertainty (not stated in the sources). torque_nm_std: +-1 std of the 30-s record (Sorensen Fig. 3).",
        "power_w = torque_nm * Omega, Omega = 72 rpm (DERIVED, aerodynamic shaft power, not electrical).",
        "tsr = Omega R / v and cp = P / (0.5 rho v^3 pi R^2) with R = 5.029 m (DERIVED).",
        "Blockage: not corrected; Simms 2001 Fig. 2 estimates U_tunnel/U_free > 0.99 (negligible), Lindenburg <= 3%.",
        "rpm: nominal 72 (synchronous 71.63 rpm + slip; per-campaign rpm not available here).",
    ]
    write_csv(os.path.join(CASE, "experiment", "power_curve.csv"), hdr, pc[cols])

    # ---- loads: root flap / edge bending at r = 0.432 m
    sf = simms("root_flap_moment_nm")
    sof = sor[sor.quantity == "root_flap_moment_nm"].copy()
    sof.index = sof.x.round().astype(int)
    soe = sor[sor.quantity == "root_edge_moment_nm"].copy()
    soe.index = soe.x.round().astype(int)
    lf = lin[lin.quantity == "root_flap_moment_knm"]
    rows = []
    for vn in NOMINAL:
        v = round(float(sf.loc[vn, "v"]), 2)
        m = float(sf.loc[vn, "val"])
        lv = lf.iloc[(lf.v_ms - v).abs().argmin()]
        cands = [abs(m - 1000 * lv.value)]
        if vn in sof.index:
            cands.append(abs(m - sof.loc[vn, "value"]))
        rows.append(dict(v_ms=v, v_nominal=vn, campaign=f"S{vn:02d}00000", rpm=RPM, pitch_deg=PITCH, yaw_deg=0.0,
                         rho=LINDENBURG_RHO[LIND_V[vn]], gauge_r_m=0.432,
                         flap_moment_gauge_axes_nm=round(m, 1),
                         flap_moment_gauge_axes_nm_err=round(max(5.0, max(cands)), 1),
                         flap_moment_gauge_axes_nm_std=round(float(sof.loc[vn, "std"]), 1) if vn in sof.index else np.nan,
                         edge_moment_gauge_axes_nm=round(float(soe.loc[vn, "value"]), 1) if vn in soe.index else np.nan,
                         edge_moment_gauge_axes_nm_err=round(max(5.0, 0.5 * float(soe.loc[vn, "value_per_pt"]) + float(soe.loc[vn, "bar_asym"])), 1) if vn in soe.index else np.nan,
                         edge_moment_gauge_axes_nm_std=round(float(soe.loc[vn, "std"]), 1) if vn in soe.index else np.nan,
                         root_flap_sorensen2002_nm=round(float(sof.loc[vn, "value"]), 1) if vn in sof.index else np.nan,
                         root_flap_lindenburg2003_nm=round(1000 * float(lv.value), 1)))
    ld = pd.DataFrame(rows)
    # Canonical columns (validation/README.md): moments about the ROTOR-PLANE axes. The gauge axes are turned by
    # the tip pitch theta: M_oop = M_flap cos(theta) - M_edge sin(theta), M_ip = M_edge cos(theta) + M_flap sin(theta).
    # 13 m/s has no measured edge moment; for the axis rotation only, it is interpolated linearly in v between the
    # neighbouring speeds (its weight is sin(3 deg) = 0.052, so the effect on the flap value is ~1 N m).
    edge_for_rotation = ld.set_index("v_ms").edge_moment_gauge_axes_nm.interpolate(method="index").to_numpy()
    th = np.radians(ld.pitch_deg.to_numpy())
    ld["root_flap_moment_nm"] = (ld.flap_moment_gauge_axes_nm * np.cos(th) - edge_for_rotation * np.sin(th)).round(1)
    ld["root_flap_moment_nm_err"] = ld.flap_moment_gauge_axes_nm_err + np.where(
        ld.edge_moment_gauge_axes_nm.isna(), 5.0, 0.0)
    ld["root_edge_moment_nm"] = (ld.edge_moment_gauge_axes_nm * np.cos(th) + ld.flap_moment_gauge_axes_nm * np.sin(th)).round(1)
    ld["root_edge_moment_nm_err"] = ld.edge_moment_gauge_axes_nm_err
    ld = ld[["v_ms", "v_nominal", "campaign", "rpm", "pitch_deg", "yaw_deg", "rho", "gauge_r_m",
             "root_flap_moment_nm", "root_flap_moment_nm_err", "root_edge_moment_nm", "root_edge_moment_nm_err",
             "flap_moment_gauge_axes_nm", "flap_moment_gauge_axes_nm_err", "flap_moment_gauge_axes_nm_std",
             "edge_moment_gauge_axes_nm", "edge_moment_gauge_axes_nm_err", "edge_moment_gauge_axes_nm_std",
             "root_flap_sorensen2002_nm", "root_flap_lindenburg2003_nm"]]
    hdr = [
        "NREL UAE Phase VI, sequence S (upwind, rigid, 0 deg cone, 0 deg yaw, tip pitch 3 deg, 72 rpm).",
        "Mean BLADE ROOT BENDING MOMENTS [N m] from the blade-root strain gauges (channels B1RFB/B1REB/B3RFB/B3REB),",
        "located on the 80-mm pitch shaft at r = 0.432 m from the rotor axis (Hand 2001 p.34; Lindenburg 2003 p.62)",
        "-> column gauge_r_m = 0.432.",
        "COLUMNS: root_flap_moment_nm / root_edge_moment_nm are the canonical ROTOR-PLANE moments (out-of-plane,",
        "downwind positive / in-plane, direction of rotation positive) about r = gauge_r_m, CONVERTED from the measured",
        "gauge-axes values *_gauge_axes_nm (the figures' values, kept unchanged). At 13 m/s the edge moment needed for",
        "the conversion is interpolated between 10 and 15 m/s (effect ~1 N m; +5 N m added to root_flap_moment_nm_err).",
        "GAUGE AXES: the gauges rotate with blade pitch and are aligned with the TIP CHORD (Hand 2001 p.34), i.e. with",
        "pitch 3 deg the 'flap' axis is 3 deg off the rotor plane (Lindenburg 2003 p.62 says the same).",
        "With theta = tip pitch (3 deg, towards feather), M_oop = out-of-plane (thrust-direction) moment, M_ip = in-plane",
        "(rotation-direction) moment, both about r = 0.432 m:  M_flap = M_oop cos(theta) + M_ip sin(theta),",
        "M_edge = M_ip cos(theta) - M_oop sin(theta).  Signs (Hand 2001 p.34): flap positive for a DOWNWIND force,",
        "edge positive for a force in the DIRECTION OF ROTATION, for the upwind configuration.",
        "Gravity/centrifugal: 0 deg cone and tilt -> no mean gravity flap moment; gravity edge moment averages to zero",
        "over a revolution (it is the large +-std of the edge moment).",
        "Sources (measured, read from VECTOR figures with PyMuPDF):",
        "  flap_moment_gauge_axes_nm: Simms et al. 2001 NREL/TP-500-29494 Fig. 9 p.22 (mean of 4 identical panel copies);",
        "  cross-checks: Sorensen et al. 2002 Fig. 4 (no 13 m/s point) and Lindenburg 2003 Fig. 6.3 - agree within 4 N m;",
        "  flap_moment_gauge_axes_nm_std, edge_moment_gauge_axes_nm(_std): Sorensen et al. 2002 Figs. 4 and 5",
        "  (mean = '+' marker, std = half error bar). Root edge moment has NO 13 m/s point in any source found (empty).",
        "_err = extraction/consistency uncertainty (flap: max source disagreement, min 5 N m; edge: single source,",
        "half a point of plot resolution + marker/bar-centre offset, min 5 N m), NOT measurement uncertainty.",
        "_std = +-1 standard deviation over the record (cyclic variation incl. gravity for the edge moment).",
    ]
    write_csv(os.path.join(CASE, "experiment", "loads.csv"), hdr, ld)

    # ---- spanwise Cn / Ct at the five pressure stations
    chord_at = lambda r: float(np.interp(r, GEO.r_m, GEO.chord_m))
    stations = {0.30: 1.510, 0.47: 2.343, 0.63: 3.185, 0.80: 4.023, 0.95: 4.780}
    rr = {0.30: "030", 0.63: "063", 0.95: "095"}
    rows = []
    cases = ["S0700000", "S1000000", "S1300000", "S1500000", "S2000000", "S2500000"]
    for k, case in enumerate(cases):
        vn = NOMINAL[k]
        v = round(float(st.loc[vn, "v"]), 2)
        for rR, rm in stations.items():
            cn = sor[(sor.quantity == "cn") & (sor.case == case)]
            ct = sor[(sor.quantity == "ct") & (sor.case == case)]
            a = cn.iloc[(cn.x - rR).abs().argmin()]
            b = ct.iloc[(ct.x - rR).abs().argmin()]
            row = dict(campaign=case, v_ms=v, v_nominal=vn, rpm=RPM, pitch_deg=PITCH, yaw_deg=0.0,
                       rho=LINDENBURG_RHO[LIND_V[vn]], tsr=round(OMEGA * R / v, 4),
                       r_R=rR if rR not in (0.47, 0.63) else {0.47: 0.466, 0.63: 0.633}[rR], r_m=rm,
                       chord_m=round(chord_at(rm), 4),
                       cn_chord=round(float(a.value), 4), cn_chord_std=round(float(a["std"]), 4),
                       ct_chord=round(float(b.value), 4), ct_chord_std=round(float(b["std"]), 4))
            if rR in rr:
                s = simms("cn_" + rr[rR]).loc[vn]
                t2 = simms("ct_" + rr[rR]).loc[vn]
                cm = simms("cm_" + rr[rR]).loc[vn]
                q = simms("q_local_pa_" + rr[rR]).loc[vn]
                row.update(cm=round(float(cm.val), 4), q_local_pa=round(float(q.val), 1),
                           cn_simms2001=round(float(s.val), 4), ct_simms2001=round(float(t2.val), 4))
                row["fn_chord_n_per_m"] = round(row["cn_chord"] * row["q_local_pa"] * row["chord_m"], 1)
                row["ft_chord_n_per_m"] = round(row["ct_chord"] * row["q_local_pa"] * row["chord_m"], 2)
            rows.append(row)
    sp = pd.DataFrame(rows)
    sp["cn_chord_err"] = 0.005
    sp["ct_chord_err"] = 0.002
    dcn = (sp.cn_chord - sp.cn_simms2001).abs().max()
    dct = (sp.ct_chord - sp.ct_simms2001).abs().max()
    cols = ["campaign", "v_ms", "v_nominal", "rpm", "pitch_deg", "yaw_deg", "rho", "tsr", "r_R", "r_m", "chord_m",
            "cn_chord", "cn_chord_err", "cn_chord_std", "ct_chord", "ct_chord_err", "ct_chord_std", "cm", "q_local_pa",
            "fn_chord_n_per_m", "ft_chord_n_per_m", "cn_simms2001", "ct_simms2001"]
    hdr = [
        "NREL UAE Phase VI, sequence S (upwind, rigid, 0 deg yaw, tip pitch 3 deg, 72 rpm): time-averaged sectional",
        "force coefficients at the five pressure-tap stations r/R = 0.30, 0.466, 0.633, 0.80, 0.95 (r = 1.510, 2.343,",
        "3.185, 4.023, 4.780 m; Hand 2001 Table 7/Lindenburg 2003 App. A).",
        "DEFINITIONS (Hand 2001 pp.51-53, Eqs. 8-10, Fig. 25): cn_chord = normal force perpendicular to the CHORD,",
        "ct_chord = tangent force PARALLEL to the chord, positive towards the LEADING EDGE; both from integrated surface",
        "pressures (no skin friction) and normalised by the local STAGNATION dynamic pressure QNORM measured at the",
        "stagnation tap (q_local_pa) and the local chord: Fn' = cn_chord*q_local*c, Ft' = ct_chord*q_local*c.",
        "In BEM terms: cn_chord = cl cos(alpha) + cd sin(alpha), ct_chord = cl sin(alpha) - cd cos(alpha) (cd = pressure",
        "drag only). cm about the quarter chord (pressure), positive nose-up.",
        "Sources (measured, read from VECTOR figures with PyMuPDF):",
        "  cn_chord(_std), ct_chord(_std) (all 5 stations): Sorensen, Michelsen & Schreck 2002, Wind Energy 5:151,",
        "    Figs. 6 and 7 (circles = measured mean, dotted bars = +-1 std over the record).",
        "  cm, q_local_pa (0.30/0.63/0.95 only): Simms et al. 2001 NREL/TP-500-29494 Figs. 12/13, 16/17, 20/21.",
        "  cn_simms2001, ct_simms2001: Simms 2001 Figs. 10/11, 14/15, 18/19 - independent cross-check of the same",
        f"    data: max |dcn| = {dcn:.4f}, max |dct| = {dct:.4f}.",
        "  fn_chord_n_per_m, ft_chord_n_per_m = cn_chord*q_local*c, ct_chord*q_local*c (DERIVED, where q_local exists);",
        "    chord_m from geometry.csv.",
        "cn_chord_err / ct_chord_err: extraction uncertainty (read precision < 0.001; set to the Simms/Sorensen spread",
        "  rounded up) - NOT the measurement uncertainty. *_std = cyclic/turbulent variation over the record.",
        "Note: the 10 m/s case (S1000000) shows an anomalous ct at 0.466 R (0.007) - real, in both sources.",
        "Note (Simms 2001 p.17): the 0.30 R station is disturbed by the hub instrumentation boxes.",
    ]
    write_csv(os.path.join(CASE, "experiment", "spanwise_cn_ct.csv"), hdr, sp[cols])
    return pc, ld, sp


if __name__ == "__main__":
    os.makedirs(os.path.join(CASE, "airfoils"), exist_ok=True)
    os.makedirs(os.path.join(CASE, "experiment"), exist_ok=True)
    GEO, TW_TIP = build_geometry()
    print(GEO.to_string())
    print("tip twist (Hand) =", TW_TIP)
    pts = s809_coordinates()
    circle_files()
    dut_polars()
    osu_polars()
    pc, ld, sp = experiments()
    print(pc.to_string())
    print(ld.to_string())
    print(sp.to_string())

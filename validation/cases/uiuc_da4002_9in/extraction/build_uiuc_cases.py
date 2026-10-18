#!/usr/bin/env python3
"""
Builds the UIUC small-propeller validation cases (validation/cases/uiuc_*) from the raw downloads:

  * UIUC Propeller Database (Brandt, Deters, Ananda, Dantsker & Selig), whole archive
    https://m-selig.ae.illinois.edu/props/download/UIUC-propDB.zip  (md5 a41e484f1fd0fb6ff80b76e27410808b)
  * APC geometry files PE0-FILES_WEB-202602.zipx  (https://www.apcprop.com/technical-information/file-downloads/)
  * UIUC Low-Speed Airfoil Tests, Vol. 3 (Clark-Y measured polar)  https://m-selig.ae.illinois.edu/pd/pub/lsat/volume03.zip
  * airfoil coordinates: SDA1045/1075/1100 from Deters, Ananda & Selig (2014) Appendix A (pdftotext of the paper),
    E63 and Clark-Y from the UIUC airfoil coordinate database, NACA 4412 from the analytic 4-digit definition.

Polars are computed with XFOIL 6.97 (repo binary xfoil_executables/xfoil) at Ncrit = 9 (primary) and 5 (alternative).

Usage:  python3 build_uiuc_cases.py [--raw /tmp/uiuc] [--only case_id ...] [--skip-xfoil]
The raw folder layout is the one produced by the download step (see extraction/README_extraction.md).
"""
import argparse
import concurrent.futures as cf
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import splprep, splev

HERE = Path(__file__).resolve()
if HERE.parent.name == "extraction" and HERE.parents[2].name == "cases":
    CASES = HERE.parents[2]
else:
    CASES = Path("/Users/miha/Documents/dev/Python-BEM/validation/cases")
REPO = CASES.parents[1]
XFOIL = REPO / "xfoil_executables" / "xfoil"
IN2M = 0.0254
RE_LIST = [20e3, 30e3, 40e3, 60e3, 80e3, 100e3, 150e3, 200e3]
NCRITS = [9.0, 5.0]
XFOIL_WORKERS = 8
FLUID = {"rho": 1.17, "kin_viscosity": 1.55e-5, "medium": "air"}

# ------------------------------------------------------------------------------------------------ sources
SOURCES = {
    "uiucpropdb": {"key": "uiucpropdb",
                   "citation": "Brandt, J.B., Deters, R.W., Ananda, G.K., Dantsker, O.D., and Selig, M.S., UIUC Propeller "
                               "Database, Vols 1-4, University of Illinois at Urbana-Champaign, Department of Aerospace "
                               "Engineering (downloaded 2026-10-18, archive UIUC-propDB.zip md5 "
                               "a41e484f1fd0fb6ff80b76e27410808b)",
                   "url": "https://m-selig.ae.illinois.edu/props/propDB.html"},
    "deters2014": {"key": "deters2014",
                   "citation": "Deters, R.W., Ananda, G.K., and Selig, M.S., 'Reynolds Number Effects on the Performance of "
                               "Small-Scale Propellers', AIAA Aviation 2014, AIAA Paper 2014-2151, Atlanta, GA, June 2014",
                   "doi": "10.2514/6.2014-2151",
                   "url": "https://m-selig.ae.illinois.edu/pubs/DetersAnandaSelig-2014-AIAA-2014-2151.pdf"},
    "deters2014phd": {"key": "deters2014phd",
                      "citation": "Deters, R.W., 'Performance and Slipstream Characteristics of Small-Scale Propellers at "
                                  "Low Reynolds Numbers', PhD dissertation, University of Illinois at Urbana-Champaign, 2014",
                      "url": "https://m-selig.ae.illinois.edu/props/volume-2/propDB-volume-2.html"},
    "brandt2011": {"key": "brandt2011",
                   "citation": "Brandt, J.B., and Selig, M.S., 'Propeller Performance Data at Low Reynolds Numbers', 49th "
                               "AIAA Aerospace Sciences Meeting, AIAA Paper 2011-1255, Orlando, FL, January 2011",
                   "doi": "10.2514/6.2011-1255",
                   "url": "https://m-selig.ae.illinois.edu/pubs/BrandtSelig-2011-AIAA-2011-1255-LRN-Propellers.pdf"},
    "brandt2005": {"key": "brandt2005",
                   "citation": "Brandt, J.B., 'Small-Scale Propeller Performance at Low Speeds', MS thesis, University of "
                               "Illinois at Urbana-Champaign, 2005",
                   "url": "https://m-selig.ae.illinois.edu/props/volume-1/Brandt_2005_UIUC-MS-Thesis.pdf"},
    "selig2015corr": {"key": "selig2015corr",
                      "citation": "Selig, M.S., and Ananda, G.K., 'Low Reynolds Number Propeller Performance Data: Wind "
                                  "Tunnel Corrections for Motor Fixture Drag', UIUC report AE 15-01, January 22, 2015",
                      "url": "https://m-selig.ae.illinois.edu/props/volume-1/uiuc-props-wind-tunnel-correction-new.pdf"},
    "apcpe0": {"key": "apcpe0",
               "citation": "APC Propellers (Landing Products Inc.), propeller geometry files PE0-FILES_WEB-202602.zipx "
                           "(file version v2025-1001, simulation date 02/24/2026), downloaded 2026-10-18",
               "url": "https://www.apcprop.com/technical-information/file-downloads/"},
    "apceng": {"key": "apceng",
               "citation": "APC Propellers, 'Engineering' technical information page: the dominant basis of most APC "
                           "airfoils is similar to NACA 4412 and Clark-Y",
               "url": "https://www.apcprop.com/technical-information/engineering"},
    "mccrink2017": {"key": "mccrink2017",
                    "citation": "McCrink, M.H., and Gregory, J.W., 'Blade Element Momentum Modeling of Low-Reynolds "
                                "Electric Propulsion Systems', Journal of Aircraft 54(1), 2017, pp. 163-176",
                    "doi": "10.2514/1.C033622",
                    "url": "https://arc.aiaa.org/doi/10.2514/1.C033622"},
    "lsat3": {"key": "lsat3",
              "citation": "Lyon, C.A., Broeren, A.P., Giguere, P., Gopalarathnam, A., and Selig, M.S., Summary of "
                          "Low-Speed Airfoil Data, Vol. 3, SoarTech Publications, Virginia Beach, VA, 1998 (UIUC LSAT "
                          "data files volume03.zip, DRAG03.TXT 'Clark-Y (B)')",
              "url": "https://m-selig.ae.illinois.edu/pd.html"},
    "uiucads": {"key": "uiucads",
                "citation": "UIUC Airfoil Coordinates Database (Selig), files e63.dat and clarky.dat",
                "url": "https://m-selig.ae.illinois.edu/ads/coord_database.html"},
    "xfoil": {"key": "xfoil",
              "citation": "Drela, M., 'XFOIL: An Analysis and Design System for Low Reynolds Number Airfoils', in Low "
                          "Reynolds Number Aerodynamics, Lecture Notes in Engineering 54, Springer, 1989 (v6.97 used)",
              "doi": "10.1007/978-3-642-84010-4_1"},
}

# ------------------------------------------------------------------------------------------------ cases
DA_COMMON_NOTE = (
    "3D-printed UIUC research propeller (Deters 2014). Blades are inserted in interchangeable hubs: three extra hubs "
    "change the blade angle by -10, -5 and +5 deg (Deters et al. 2014 Sec. III, Table 2), which is how the 'pitched' "
    "data sets were produced; their offset is the experiment's pitch_deg (nominal hub value). ")
CASE_DEFS = [
    dict(id="uiuc_da4002_9in", name="UIUC DA4002 9-in propeller (2 blades, c/R = 0.18 rectangular, SDA1075)",
         family="da", vol=2, prop="da4002", D_in=9.05, B=2, scan="da4002_9x6.75_geom.txt", design="da4002_geom.txt",
         rect_tip=True, foils="sda1075",
         sets=[("9x6.75", 0.0, "da4002_9x6.75_geom.txt"), ("9x2.85", -10.0, "da4002_9x2.85_geom.txt"),
               ("9x4.76", -5.0, "da4002_9x4.76_geom.txt"), ("9x8.95", 5.0, "da4002_9x8.95_geom.txt")]),
    dict(id="uiuc_da4002_5in", name="UIUC DA4002 5-in propeller (2 blades, c/R = 0.18 rectangular, SDA1075)",
         family="da", vol=2, prop="da4002", D_in=5.0, B=2, scan="da4002_5x3.75_geom.txt", design="da4002_geom.txt",
         rect_tip=True, foils="sda1075",
         sets=[("5x3.75", 0.0, "da4002_5x3.75_geom.txt"), ("5x1.58", -10.0, "da4002_5x1.58_geom.txt"),
               ("5x2.65", -5.0, "da4002_5x2.65_geom.txt"), ("5x4.92", 5.0, "da4002_5x4.92_geom.txt")]),
    dict(id="uiuc_da4022_9in", name="UIUC DA4022 9-in propeller (2 blades, c/R = 0.23 rectangular, SDA1075)",
         family="da", vol=2, prop="da4022", D_in=9.1, B=2, scan="da4022_9x6.75_geom.txt", design="da4022_geom.txt",
         rect_tip=True, foils="sda1075", sets=[("9x6.75", 0.0, None)]),
    dict(id="uiuc_da4022_9in_3b", name="UIUC DA4022 9-in propeller, 3-blade hub (c/R = 0.23 rectangular, SDA1075)",
         family="da", vol=2, prop="da4022", D_in=9.1, B=3, scan="da4022_9x6.75_geom.txt", design="da4022_geom.txt",
         rect_tip=True, foils="sda1075", sets=[("9x6.75_3b", 0.0, None)]),
    dict(id="uiuc_da4022_9in_4b", name="UIUC DA4022 9-in propeller, 4-blade hub (c/R = 0.23 rectangular, SDA1075)",
         family="da", vol=2, prop="da4022", D_in=9.1, B=4, scan="da4022_9x6.75_geom.txt", design="da4022_geom.txt",
         rect_tip=True, foils="sda1075", sets=[("9x6.75_4b", 0.0, None)]),
    dict(id="uiuc_da4052_9in", name="UIUC DA4052 9-in propeller (2 blades, PROPID minimum-induced-loss design, "
                                    "SDA1045/SDA1075/SDA1100)",
         family="da", vol=2, prop="da4052", D_in=9.05, B=2, scan="da4052_9x6.75_geom.txt", design="da4052_geom.txt",
         rect_tip=False, foils="sda3",
         sets=[("9x6.75", 0.0, None), ("9x2.85", -10.0, None), ("9x4.76", -5.0, None), ("9x8.95", 5.0, None)]),
    dict(id="uiuc_da4052_9in_3b", name="UIUC DA4052 9-in propeller, 3-blade hub (SDA1045/SDA1075/SDA1100)",
         family="da", vol=2, prop="da4052", D_in=9.05, B=3, scan="da4052_9x6.75_geom.txt", design="da4052_geom.txt",
         rect_tip=False, foils="sda3", sets=[("9x6.75_3b", 0.0, None)]),
    dict(id="uiuc_apce_10x7", name="APC Thin Electric 10x7 (10x7E), UIUC Vol. 1", family="apc", vol=1,
         prop="apce", size="10x7", D_in=10.0, B=2, pe0="10x7E-PERF.PE0", scan="apce_10x7_geom.txt",
         blend=("e63", "naca4412"), anchors_rR=[0.22, 0.40, 0.56, 0.732]),
    dict(id="uiuc_apcsp_11x7", name="APC Sport 11x7, UIUC Vol. 1", family="apc", vol=1,
         prop="apcsp", size="11x7", D_in=11.0, B=2, pe0="11x7-PERF.PE0", scan="apcsp_11x7_geom.txt",
         blend=("clarky", "clarky"), anchors_rR=[0.15, 0.255, 0.45, 0.65, 0.85, 1.0]),
    dict(id="uiuc_apcsf_10x4p7", name="APC Slow Flyer 10x4.7 (10x4.7SF), UIUC Vol. 1", family="apc", vol=1,
         prop="apcsf", size="10x4.7", D_in=10.0, B=2, pe0="10x47SF-PERF.PE0", scan="apcsf_10x4.7_geom.txt",
         blend=("e63", "e63"), anchors_rR=[0.17, 0.26, 0.40]),
]


# ------------------------------------------------------------------------------------------------ helpers
def read_table(path):
    """UIUC data files: header line + whitespace columns; CR, LF or CRLF line ends."""
    txt = Path(path).read_bytes().decode("latin-1").replace("\r\n", "\n").replace("\r", "\n")
    lines = [l for l in txt.split("\n") if l.strip()]
    # a few files carry an extra first line (e.g. 'RPM average = 4042'): the header is the line before the numbers
    k = next(i for i, l in enumerate(lines) if re.match(r"^\s*-?\d", l))
    head = lines[k - 1].split()
    rows = [[float(v) for v in l.split()] for l in lines[k:]]
    return pd.DataFrame(rows, columns=head)


def write_csv(path, df, header_lines, float_fmt="%.6g"):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        for h in header_lines:
            f.write("# " + h + "\n")
        df.to_csv(f, index=False, float_format=float_fmt)


def md5(path):
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


# ------------------------------------------------------------------------------------------------ airfoils
def read_coords(path):
    rows = []
    lines = Path(path).read_bytes().decode("latin-1").replace("\r", "\n").split("\n")
    for line in lines[1:]:
        p = line.replace(",", " ").split()
        if len(p) >= 2:
            try:
                rows.append((float(p[0]), float(p[1])))
            except ValueError:
                pass
    xy = np.array(rows)
    if xy[:, 0].max() > 1.5:  # Lednicer
        nu, nl = int(xy[0, 0]), int(xy[0, 1])
        up, lo = xy[1:1 + nu], xy[1 + nu:1 + nu + nl]
        xy = np.vstack([up[::-1], lo[1:]])
    return xy


def normalise(xy):
    """LE = point of minimum x, TE = midpoint of the two end points; rotate chord to the x axis, unit chord."""
    xy = np.asarray(xy, float)
    ile = int(np.argmin(xy[:, 0]))
    le = xy[ile]
    te = 0.5 * (xy[0] + xy[-1])
    ang = math.atan2(te[1] - le[1], te[0] - le[0])
    c = math.hypot(te[0] - le[0], te[1] - le[1])
    rot = np.array([[math.cos(-ang), -math.sin(-ang)], [math.sin(-ang), math.cos(-ang)]])
    out = (xy - le) @ rot.T / c
    return out, math.degrees(ang)


def naca4(m, p, t, n=121):
    beta = np.linspace(0, math.pi, n)
    x = 0.5 * (1 - np.cos(beta))
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 - 0.1015 * x ** 4)
    yc = np.where(x < p, m / p ** 2 * (2 * p * x - x ** 2), m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x - x ** 2))
    dyc = np.where(x < p, 2 * m / p ** 2 * (p - x), 2 * m / (1 - p) ** 2 * (p - x))
    th = np.arctan(dyc)
    xu, yu = x - yt * np.sin(th), yc + yt * np.cos(th)
    xl, yl = x + yt * np.sin(th), yc - yt * np.cos(th)
    return np.vstack([np.c_[xu[::-1], yu[::-1]], np.c_[xl[1:], yl[1:]]])


def surfaces(xy, xg):
    """Upper and lower y at the x grid, from a dense parametric spline through the Selig coordinates."""
    tck, _ = splprep([xy[:, 0], xy[:, 1]], s=0, k=3)
    u = np.linspace(0, 1, 20001)
    x, y = splev(u, tck)
    i = int(np.argmin(x))
    xu, yu = x[:i + 1][::-1], y[:i + 1][::-1]
    xl, yl = x[i:], y[i:]
    return np.interp(xg, xu, yu), np.interp(xg, xl, yl)


def blend_section(xy1, xy2, f, tc, n=121):
    """APC-style section: camber and thickness shape blended (fraction f of airfoil 2), thickness scaled to tc."""
    xg = 0.5 * (1 - np.cos(np.linspace(0, math.pi, n)))
    u1, l1 = surfaces(xy1, xg)
    u2, l2 = surfaces(xy2, xg)
    c1, t1 = 0.5 * (u1 + l1), u1 - l1
    c2, t2 = 0.5 * (u2 + l2), u2 - l2
    cam = (1 - f) * c1 + f * c2
    shape = (1 - f) * t1 / t1.max() + f * t2 / t2.max()
    thick = shape / shape.max() * tc
    yu, yl = cam + thick / 2, cam - thick / 2
    return np.vstack([np.c_[xg[::-1], yu[::-1]], np.c_[xg[1:], yl[1:]]])


def write_dat(path, name, xy):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write(name + "\n")
        for x, y in xy:
            f.write("%10.6f %10.6f\n" % (x, y))


def parse_sda(paper_txt):
    t = Path(paper_txt).read_text(encoding="utf-8", errors="replace").replace("\u2212", "-")
    out = {}
    for name in ["SDA1045", "SDA1075", "SDA1100"]:
        i = t.find("%s coordinates" % name)
        up, lo = [], []
        for line in t[i:i + 6000].splitlines()[1:]:
            nums = re.findall(r"-?\d+\.\d+", line)
            if len(nums) == 4:
                up.append(nums[:2]); lo.append(nums[2:])
            elif len(nums) == 2 and up:
                up.append(nums)
                break
        out[name.lower()] = np.array(up + lo, float)
    return out


# ------------------------------------------------------------------------------------------------ xfoil
def xfoil_sweep(datfile, re_, ncrit, a0, a1, da, workdir):
    pol = Path(workdir) / ("p_%d_%d_%g_%g.pol" % (re_, ncrit, a0, a1))
    if pol.exists():
        pol.unlink()
    cmds = ["PLOP", "G F", "", "LOAD %s" % datfile, "PANE", "OPER", "VPAR", "N %g" % ncrit, "",
            "VISC %g" % re_, "ITER 250", "PACC", pol.name, "", "ASEQ %g %g %g" % (a0, a1, da), "", "QUIT"]
    try:
        r = subprocess.run([str(XFOIL)], input="\n".join(cmds) + "\n", text=True, capture_output=True,
                           cwd=workdir, timeout=180)
        log = r.stdout[-3000:]
    except subprocess.TimeoutExpired:
        log = "TIMEOUT"
    rows = []
    if pol.exists():
        for line in pol.read_text().splitlines():
            p = line.split()
            if len(p) >= 7:
                try:
                    rows.append([float(v) for v in p[:7]])
                except ValueError:
                    pass
    return rows, log


def xfoil_polar(datfile, re_, ncrit, cache_dir):
    """Merged XFOIL polar for one Re/Ncrit: sweeps 0..20 and 0..-10 from a fresh start, then fresh restarts at
    8 and 14 deg to fill points where the continuation sweep failed to converge. Cached by coordinate md5."""
    key = "%s_%s_re%d_n%g" % (Path(datfile).stem, md5(datfile)[:8], re_, ncrit)
    cf_ = Path(cache_dir) / (key + ".csv")
    if cf_.exists():
        return pd.read_csv(cf_)
    with tempfile.TemporaryDirectory() as td:
        shutil.copy(datfile, Path(td) / "foil.dat")
        allrows = {}

        def run(a0, a1, da):
            rows, _ = xfoil_sweep("foil.dat", re_, ncrit, a0, a1, da, td)
            for r in rows:
                a = round(r[0], 2)
                if a not in allrows:
                    allrows[a] = r

        run(0, 20, 0.5)
        run(-0.5, -10, -0.5)
        # fresh restarts only where the continuation sweep left gaps
        for a0, a1, da in [(8, 20, 0.5), (14, 20, 0.5), (-5, -10, -0.5)]:
            want = np.arange(a0, a1 + da / 2, da)
            if sum(round(a, 2) not in allrows for a in want) > 2:
                run(a0, a1, da)
    df = pd.DataFrame(sorted(allrows.values()), columns=["alpha_deg", "cl", "cd", "cdp", "cm", "xtr_top", "xtr_bot"])
    df = df[(df.cd > 0) & (df.cd < 0.5)]
    df.to_csv(cf_, index=False)
    return df


def build_polars(datfile, label, cache_dir, ncrit):
    jobs = {}
    with cf.ThreadPoolExecutor(max_workers=10) as ex:
        for re_ in RE_LIST:
            jobs[re_] = ex.submit(xfoil_polar, datfile, re_, ncrit, cache_dir)
    parts, kept = [], []
    for re_, j in jobs.items():
        d = j.result()
        grid = np.arange(math.ceil(d.alpha_deg.min()), math.floor(d.alpha_deg.max()) + 1) if len(d) else []
        if len(d) < 12 or len(grid) < 8 or not (d.alpha_deg.min() < 0 < d.alpha_deg.max()):
            print("   drop %s Re=%g Ncrit=%g: %d points" % (label, re_, ncrit, len(d)))
            continue
        d = d.copy()
        d.insert(0, "re", re_)
        d["ncrit"] = ncrit
        d["source"] = "XFOIL 6.97 Ncrit=%g" % ncrit
        parts.append(d[["re", "alpha_deg", "cl", "cd", "cm", "ncrit", "source"]])
        kept.append(re_)
    return pd.concat(parts, ignore_index=True), kept


def clarky_lsat(raw):
    """Measured Clark-Y (B) polar, UIUC LSAT Vol 3, drag-polar runs (alpha, cl, cd)."""
    txt = (raw / "lsat/v03/volume03/DRAG03.TXT").read_bytes().decode("latin-1").replace("\r", "")
    i = txt.find("CLARKYB.DRG\n")
    j = txt.find("CLARKYBT.DRG")
    seg = txt[i:j]
    rows = []
    blocks = seg.split("Average Reynolds #:")[1:]
    for b in blocks:
        lines = b.strip().split("\n")
        re_ = float(lines[0])
        for l in lines:
            p = l.split()
            if len(p) >= 7:
                try:
                    a, cl, cd = float(p[0]), float(p[1]), float(p[2])
                except ValueError:
                    continue
                rows.append([re_, a, cl, cd])
    df = pd.DataFrame(rows, columns=["re_run", "alpha_deg", "cl", "cd"])
    # runs at ~60k and ~100k are split in two files each (drag polar + 1-2 high-alpha points): group by nominal Re
    df["re"] = (df.re_run / 1e4).round().astype(int) * 1e4
    df.loc[(df.re_run > 55e3) & (df.re_run < 65e3), "re"] = 60e3
    df["source"] = "measured, UIUC LSAT Vol. 3 (Lyon et al. 1998), Clark-Y (B), DRAG03.TXT, run Re=" + \
        df.re_run.astype(int).astype(str)
    return df[["re", "alpha_deg", "cl", "cd", "source"]].sort_values(["re", "alpha_deg"])


# ------------------------------------------------------------------------------------------------ APC PE0
def read_pe0(path):
    txt = Path(path).read_text(errors="replace")
    rows = []
    for line in txt.splitlines():
        p = line.split()
        if len(p) == 14:
            try:
                rows.append([float(v) for v in p])
            except ValueError:
                pass
    cols = ["station_in", "chord_in", "pitch_quoted_in", "pitch_lete_in", "pitch_prather_in", "sweep_in", "rake_in",
            "thickness_ratio", "twist_deg", "max_thick_in", "area_in2", "zhigh_in", "cgy_in", "cgz_in"]
    df = pd.DataFrame(rows, columns=cols)
    meta = {}
    for k in ["RADIUS", "HUBRAD", "HUBTRA", "BLADES"]:
        m = re.search(r"%s:\s+([\d.]+)" % k, txt)
        meta[k] = float(m.group(1))
    af = re.findall(r"AIRFOIL(\d):\s+([\d.]+),\s+(\S+)", txt)
    meta["airfoils"] = [(float(r), n) for _, r, n in af]
    meta["header"] = txt.splitlines()[0].strip() + " / " + txt.splitlines()[1].strip()
    return df, meta


# ------------------------------------------------------------------------------------------------ experiments
def find_runs(datadir, base):
    """Performance runs and static file of one configuration, e.g. base='da4002_9x6.75' or 'apce_10x7'."""
    runs, static = [], None
    pr = re.compile(r"^%s_([0-9a-z]{6})_(\d{4,5})\.txt$" % re.escape(base))
    ps = re.compile(r"^%s_static_([0-9a-z]{6})\.txt$" % re.escape(base))
    for f in sorted(os.listdir(datadir)):
        m = pr.match(f)
        if m:
            runs.append((f, m.group(1), int(m.group(2))))
        elif ps.match(f):
            static = f
    return runs, static


# ------------------------------------------------------------------------------------------------ build
def tip_fix_da(scan, design, rect):
    """PropellerScanner tip artefacts (chord and twist collapse at r/R = 0.95-1.00 although the 3D-printed blade has
    a square / designed tip, see photos). Returns corrected table + list of notes."""
    s = scan.copy()
    notes = []
    d = design.set_index("r/R")
    rows = [0.95, 1.00] if rect else [1.00]
    for rr in rows:
        i = s.index[np.isclose(s["r/R"], rr)][0]
        if rect:
            ref = 0.90
            c_new = float(s.loc[np.isclose(s["r/R"], ref), "c/R"].iloc[0])
        else:
            ref = 0.95
            c_new = float(d.loc[rr, "c/R"])
        b_ref = float(s.loc[np.isclose(s["r/R"], ref), "beta"].iloc[0])
        b_new = b_ref + float(d.loc[rr, "beta"] - d.loc[ref, "beta"])
        notes.append("r/R=%.2f: scan c/R=%.4f beta=%.3f replaced by c/R=%.4f (%s) and beta=%.3f (scan beta at "
                     "r/R=%.2f + design twist change %.2f->%.2f)" % (
                         rr, s.loc[i, "c/R"], s.loc[i, "beta"], c_new,
                         "scan c/R at r/R=0.90, rectangular blade" if rect else "design c/R",
                         b_new, ref, ref, rr))
        s.loc[i, "c/R"] = c_new
        s.loc[i, "beta"] = b_new
    return s, notes


def da_airfoil_column(rR, foils):
    if foils == "sda1075":
        return ["sda1075"] * len(rR)
    out = []
    for r in rR:
        if r <= 0.45 + 1e-9:
            out.append("sda1045")
        elif abs(r - 0.75) < 1e-9:
            out.append("sda1075")
        elif r >= 1.0 - 1e-9:
            out.append("sda1100")
        else:
            out.append("transition")
    return out


def build_case(cd, raw, foil_lib, args):
    case = CASES / cd["id"]
    for sub in ["airfoils", "experiment", "extraction", "sources"]:
        (case / sub).mkdir(parents=True, exist_ok=True)
    vol = cd["vol"]
    datadir = raw / "UIUC-propDB" / ("volume-%d" % vol) / "data"
    D = cd["D_in"] * IN2M
    R = D / 2
    meta_notes = []
    experiments = []
    airfoils_used = {}
    src_keys = ["uiucpropdb"]

    # ---------------- geometry
    if cd["family"] == "da":
        scan = read_table(datadir / cd["scan"])
        design = read_table(datadir / cd["design"])
        geo, tipnotes = tip_fix_da(scan, design, cd["rect_tip"])
        foils = da_airfoil_column(geo["r/R"].to_numpy(), cd["foils"])
        g = pd.DataFrame({"r_m": geo["r/R"] * R, "chord_m": geo["c/R"] * R, "twist_deg": geo["beta"],
                          "airfoil": foils, "thickness_ratio": 0.12})
        geo_hdr = [
            "%s geometry: UIUC Propeller Database Vol. 2, file %s ('Geometry (tested)': PropellerScanner "
            "measurement of the 3D-printed propeller, average of both blades; tabulated r/R, c/R, beta)." % (
                cd["name"], cd["scan"]),
            "Converted: r_m = r/R * R, chord_m = c/R * R with R = %.4f m (true diameter %.2f in from the database "
            "geometry plot); twist_deg = beta (blade angle from the plane of rotation, chord line) unchanged." % (
                R, cd["D_in"]),
            "Tip correction (PropellerScanner artefact, square/designed tip visible in the photos): " + "; ".join(tipnotes),
            "Uncorrected scan: extraction/geometry_scan_raw.csv; design geometry (Deters 2014 Fig. %s): "
            "extraction/geometry_design.csv." % ("26" if cd["prop"] == "da4002" else "22"),
            "Airfoil: %s (Deters, Ananda & Selig 2014, Sec. III); t/c = 0.12 for all SDA airfoils." % (
                "SDA1075 along the whole blade" if cd["foils"] == "sda1075" else
                "SDA1045 at r/R <= 0.45, SDA1075 at 0.75, SDA1100 at the tip, linear transitions in between"),
        ]
        write_csv(case / "geometry.csv", g, geo_hdr)
        raw_g = pd.DataFrame({"r_m": scan["r/R"] * R, "chord_m": scan["c/R"] * R, "twist_deg": scan["beta"],
                              "airfoil": da_airfoil_column(scan["r/R"].to_numpy(), cd["foils"]), "thickness_ratio": 0.12})
        write_csv(case / "extraction" / "geometry_scan_raw.csv", raw_g,
                  ["Raw PropellerScanner geometry, UIUC propDB Vol. 2 file %s, R = %.4f m, no tip correction" % (cd["scan"], R)])
        des_g = pd.DataFrame({"r_m": design["r/R"] * R, "chord_m": design["c/R"] * R, "twist_deg": design["beta"],
                              "airfoil": da_airfoil_column(design["r/R"].to_numpy(), cd["foils"]),
                              "thickness_ratio": 0.12})
        write_csv(case / "extraction" / "geometry_design.csv", des_g,
                  ["Design geometry, UIUC propDB Vol. 2 file %s ('Geometry (design)'), scaled to R = %.4f m" % (
                      cd["design"], R)])
        Rhub = float(g.r_m.iloc[0])
        names = ["sda1075"] if cd["foils"] == "sda1075" else ["sda1045", "sda1075", "sda1100"]
        for n in names:
            airfoils_used[n] = foil_lib[n]
        src_keys += ["deters2014", "deters2014phd"]
    else:
        pe0, pm = read_pe0(raw / "apc" / "PE0-FILES_WEB" / cd["pe0"])
        Rpe0 = pm["RADIUS"] * IN2M
        assert abs(Rpe0 - R) < 1e-6, (Rpe0, R)
        sel = pe0[pe0.station_in >= 0.15 * pm["RADIUS"] - 1e-9].reset_index(drop=True)
        rR = (sel.station_in / pm["RADIUS"]).to_numpy()
        # airfoil column: anchor stations (nearest PE0 station to each target r/R) named, transition in between
        anchor_idx = [int(np.argmin(np.abs(rR - a))) for a in cd["anchors_rR"]]
        col = []
        for i in range(len(sel)):
            if i in anchor_idx:
                col.append(foil_lib["apc_anchor"][cd["id"]][anchor_idx.index(i)])
            elif i < anchor_idx[0]:
                col.append(foil_lib["apc_anchor"][cd["id"]][0])
            elif i > anchor_idx[-1]:
                col.append(foil_lib["apc_anchor"][cd["id"]][-1])
            else:
                col.append("transition")
        g = pd.DataFrame({"r_m": sel.station_in * IN2M, "chord_m": sel.chord_in * IN2M, "twist_deg": sel.twist_deg,
                          "airfoil": col, "thickness_ratio": sel.thickness_ratio})
        a1, a2 = pm["airfoils"][0], pm["airfoils"][1]
        geo_hdr = [
            "%s geometry: APC manufacturer geometry file %s (%s), from PE0-FILES_WEB-202602.zipx; tabulated." % (
                cd["name"], cd["pe0"], pm["header"]),
            "Columns used: STATION (in) -> r_m, CHORD (in) -> chord_m (x 0.0254), TWIST (deg) -> twist_deg, THICKNESS "
            "RATIO -> thickness_ratio. APC defines TWIST and CHORD with the mold LE and TE parting lines, i.e. the LE-TE "
            "chord line; TWIST = atan(quoted pitch / 2 pi r) = blade angle from the plane of rotation (solver convention).",
            "Stations inboard of r = 0.15 R dropped (hub; APC HUBRAD = %.2f in, HUBTRA (hub transition) = %.2f in)." % (
                pm["HUBRAD"], pm["HUBTRA"]),
            "APC airfoil sections: AIRFOIL1 %s at %.2f in (transition start), AIRFOIL2 %s at %.2f in (transition end), "
            "'scaled based on thickness ratio'; APC12 = NACA 4412 (APC note in the file)." % (a1[1], a1[0], a2[1], a2[0]),
            "Airfoil column: named XFOIL sections at anchor stations r/R = %s (blend fraction and t/c of the station, "
            "see case.json airfoils), 'transition' rows in between (solver interpolates cl/cd)." % (
                ", ".join("%.3f" % rR[i] for i in anchor_idx)),
            "UIUC PropellerScanner geometry of the tested specimen (Vol. 1 file %s): extraction/geometry_uiuc_scan.csv" % cd["scan"],
        ]
        write_csv(case / "geometry.csv", g, geo_hdr)
        scan = read_table(datadir / cd["scan"])
        sg = pd.DataFrame({"r_m": scan["r/R"] * R, "chord_m": scan["c/R"] * R, "twist_deg": scan["beta"],
                           "airfoil": [col[int(np.argmin(np.abs(rR - x)))] if x >= rR[0] else col[0]
                                       for x in scan["r/R"]],
                           "thickness_ratio": np.interp(scan["r/R"], rR, sel.thickness_ratio)})
        # airfoil column of the scan: anchors only at the PE0 anchor radii -> name nearest, rest transition
        sg_air = []
        for x in scan["r/R"]:
            j = int(np.argmin(np.abs(np.array(cd["anchors_rR"]) - x)))
            if abs(cd["anchors_rR"][j] - x) < 0.026:
                sg_air.append(foil_lib["apc_anchor"][cd["id"]][j])
            elif x < cd["anchors_rR"][0]:
                sg_air.append(foil_lib["apc_anchor"][cd["id"]][0])
            elif x > cd["anchors_rR"][-1]:
                sg_air.append(foil_lib["apc_anchor"][cd["id"]][-1])
            else:
                sg_air.append("transition")
        sg["airfoil"] = sg_air
        write_csv(case / "extraction" / "geometry_uiuc_scan.csv", sg,
                  ["UIUC propDB Vol. 1 PropellerScanner geometry, file %s (tabulated r/R, c/R, beta), R = %.4f m;"
                   " thickness_ratio interpolated from the APC PE0 file; airfoil names as in geometry.csv" % (cd["scan"], R)])
        pe0.to_csv(case / "extraction" / ("apc_%s.csv" % Path(cd["pe0"]).stem), index=False)
        Rhub = float(g.r_m.iloc[0])
        for n in foil_lib["apc_anchor"][cd["id"]]:
            airfoils_used[n] = foil_lib[n]
        if cd["id"] == "uiuc_apcsp_11x7":
            airfoils_used["clarky_uiuc_lsat"] = foil_lib["clarky_uiuc_lsat"]
        src_keys += ["apcpe0", "apceng", "brandt2011", "brandt2005", "selig2015corr", "uiucads"]

    # ---------------- experiments
    if cd["family"] == "da":
        sets = cd["sets"]
    else:
        sets = [(cd["size"], 0.0, None)]
    scan0 = read_table(datadir / cd["scan"])
    for (cfg, pitch, pscanfile) in sets:
        base = "%s_%s" % (cd["prop"], cfg)
        runs, static = find_runs(datadir, base)
        if not runs:
            raise RuntimeError("no runs for " + base)
        tag = "pitch%+03d" % pitch if pitch != 0 else "pitch0"
        tag = tag.replace("+", "p").replace("-", "m")
        pscan = None
        if pscanfile:
            ps = read_table(datadir / pscanfile)
            m = (ps["r/R"] >= 0.40 - 1e-9) & (ps["r/R"] <= 0.90 + 1e-9)
            pscan = float((ps.loc[m, "beta"].to_numpy() - scan0.loc[m, "beta"].to_numpy()).mean())
        parts = []
        for f, rid, rpm in runs:
            d = read_table(datadir / f)
            d.columns = [c.strip() for c in d.columns]
            out = pd.DataFrame({"J": d["J"], "ct": d["CT"], "cp": d["CP"], "eta": d["eta"]})
            out["rpm"] = rpm
            out["v_ms"] = out.J * rpm / 60.0 * D
            out["pitch_deg"] = pitch
            if pscan is not None:
                out["pitch_deg_scan"] = round(pscan, 2)
            out["run"] = "%s" % f
            parts.append(out)
        perf = pd.concat(parts, ignore_index=True)
        bl = "B=%d" % cd["B"]
        hdr = [
            "%s, %s, blade-angle offset %+g deg (%s), performance vs J at constant rpm (wind tunnel speed varied)." % (
                cd["name"], bl, pitch, base),
            "Source: UIUC Propeller Database Vol. %d (%s), data files %s; tabulated (J, CT, CP, eta), "
            "copied without change, no digitising." % (
                vol, "Deters 2014 PhD / Deters et al. 2014" if vol == 2 else "Brandt 2005 MS thesis + later tests, "
                "ver. 3 2015-02-03, with the 2015 motor-fixture drag correction", ", ".join(r[0] for r in runs)),
            "rpm = nominal value in the file name (the run's set speed); v_ms = J * rpm/60 * D with D = %.4f m "
            "(true diameter %.2f in, coefficients in the database are based on it)." % (D, cd["D_in"]),
            "Coefficients: CT = T/(rho n^2 D^4), CP = P/(rho n^3 D^5), eta = J CT/CP; tunnel corrections already applied "
            "by UIUC (%s)." % ("fairing blockage via source panels + Glauert closed-section correction, Deters et al. "
                               "2014 Sec. II.E" if vol == 2 else "Glauert/blockage/pitot-inflow velocity corrections, "
                               "Brandt & Selig 2011 Sec. III, and motor-fixture drag correction Selig & Ananda 2015"),
            "pitch_deg = blade-angle offset relative to geometry.csv (the geometry holds the actual blade angles of the "
            "0-deg configuration)." + (" pitch_deg_scan = mean (r/R 0.40-0.90) difference between the PropellerScanner "
                                       "geometry of this configuration (%s) and of the 0-deg one: an alternative estimate "
                                       "of the real offset." % pscanfile if pscanfile else ""),
            "Uncertainty: not quantified by the source; repeatability tests (Deters et al. 2014 Figs. 3-4) and the "
            "UIUC/OSU comparison suggest roughly +-3-5 % in CT and CP; low-J points of each run carry the largest "
            "velocity uncertainty.",
        ]
        fn = "experiment/performance_%s.csv" % tag if len(sets) > 1 or pitch != 0 else "experiment/performance.csv"
        write_csv(case / fn, perf, hdr)
        experiments.append({
            "file": fn, "kind": "performance", "x": "J", "y": ["ct", "cp", "eta"],
            "defaults": {"pitch_deg": pitch, "rho": FLUID["rho"], "kin_viscosity": FLUID["kin_viscosity"]},
            "conditions": "UIUC 2.8x4.0 ft low-turbulence tunnel (Tu < 0.1 %%), %s, constant rpm per run (%s rpm), "
                          "tunnel speed varied; blade-angle offset %+g deg; tunnel corrections applied by UIUC." % (
                              bl, ", ".join(str(r[2]) for r in runs), pitch),
            "source": "UIUC Propeller Database Vol. %d, files %s (tabulated)" % (vol, ", ".join(r[0] for r in runs)),
            "uncertainty": "not given by the source; ~+-3-5 % in CT/CP estimated from the published repeatability"})
        if static:
            s = read_table(datadir / static)
            st = pd.DataFrame({"rpm": s["RPM"], "J": 0.0, "v_ms": 0.0, "ct": s["CT"], "cp": s["CP"]})
            st["pitch_deg"] = pitch
            sfn = "experiment/static_%s.csv" % tag if len(sets) > 1 or pitch != 0 else "experiment/static.csv"
            write_csv(case / sfn, st, [
                "%s, %s, blade-angle offset %+g deg: static performance (J = 0, tunnel off), CT and CP vs rpm." % (
                    cd["name"], bl, pitch),
                "Source: UIUC Propeller Database Vol. %d, file %s; tabulated (RPM, CT, CP), copied without change." % (
                    vol, static),
                "Re in the database static plots is based on c(0.75R) and 0.75*Omega*R. No tunnel corrections needed.",
                "Uncertainty: not quantified by the source; scatter between neighbouring rpm points ~+-0.002 in CT.",
            ], float_fmt="%.6g")
            experiments.append({
                "file": sfn, "kind": "performance", "x": "rpm", "y": ["ct", "cp"],
                "defaults": {"pitch_deg": pitch, "v_ms": 0.0, "J": 0.0, "rho": FLUID["rho"],
                             "kin_viscosity": FLUID["kin_viscosity"]},
                "conditions": "static (J = 0, tunnel off), %s, rpm swept %d-%d; blade-angle offset %+g deg" % (
                    bl, st.rpm.min(), st.rpm.max(), pitch),
                "source": "UIUC Propeller Database Vol. %d, file %s (tabulated)" % (vol, static),
                "uncertainty": "not given by the source"})

    # OSU data for the APC 10x7E (McCrink & Gregory), from the MATLAB .fig on the UIUC site
    if cd["id"] == "uiuc_apce_10x7":
        osu = osu_10x7(raw)
        osu["rpm"] = 6500
        osu["v_ms"] = osu.J * 6500 / 60.0 * D
        osu["pitch_deg"] = 0.0
        write_csv(case / "experiment/performance_osu.csv", osu, [
            "APC Thin Electric 10x7E measured at Ohio State University (McCrink & Gregory), as plotted in the UIUC/OSU "
            "comparison on the UIUC Propeller Data Site (June 2015).",
            "Source: MATLAB figure files CT_OSU_UIUC.fig and CP_OSU_UIUC.fig (series '10x7_{OSU}') from "
            "https://m-selig.ae.illinois.edu/props/comparisons/OSU_UIUC_PROP_DATA ; the plotted XData/YData were read "
            "exactly from the .fig (no raster digitising); eta = J*ct/cp computed here.",
            "rpm = 6500 is an ASSUMPTION: the figure compares with the UIUC 6531/6519-rpm runs and McCrink & Gregory "
            "(2017) state the comparison runs were taken at constant rpm at conditions similar to Brandt & Selig; the "
            "exact OSU rpm is not stated. OSU data: multi-axis load cell on the motor, fairing blockage corrected.",
            "Uncertainty: McCrink & Gregory 2017 show 95 % confidence bars (Monte Carlo); not given numerically.",
        ])
        experiments.append({
            "file": "experiment/performance_osu.csv", "kind": "performance", "x": "J", "y": ["ct", "cp", "eta"],
            "defaults": {"pitch_deg": 0.0, "rho": FLUID["rho"], "kin_viscosity": FLUID["kin_viscosity"]},
            "conditions": "Ohio State 3x5 ft tunnel, constant rpm (assumed 6500, not stated), speed varied",
            "source": "McCrink & Gregory (OSU) data in the UIUC/OSU comparison .fig files (exact plotted values)",
            "uncertainty": "rpm assumed; data 95 % CI shown graphically in McCrink & Gregory 2017"})
        src_keys.append("mccrink2017")

    # ---------------- airfoil files
    af_meta = {}
    for name, info in airfoils_used.items():
        dst_dat = case / "airfoils" / (name + ".dat")
        if info.get("dat"):
            shutil.copy(info["dat"], dst_dat)
        dst_pol = case / "airfoils" / (name + "_polars.csv")
        write_csv(dst_pol, info["polars"], info["polar_header"])
        entry = {"coordinates": "airfoils/%s.dat" % name if info.get("dat") else None,
                 "polars": "airfoils/%s_polars.csv" % name, "polar_source": info["polar_source"]}
        if info.get("thickness_ratio"):
            entry["thickness_ratio"] = info["thickness_ratio"]
        if info.get("alt") is not None:
            alt = case / "airfoils" / (name + "_ncrit5_polars.csv")
            write_csv(alt, info["alt"], info["alt_header"])
            entry["alternative_polars"] = {"ncrit5": "airfoils/%s_ncrit5_polars.csv" % name}
        if info.get("note"):
            entry["note"] = info["note"]
        if info.get("used") is False:
            entry["used_in_geometry"] = False
        if entry["coordinates"] is None:
            del entry["coordinates"]
        af_meta[name] = entry
        if info.get("xfoil_dir"):
            pass
    src_keys.append("xfoil")
    if cd["id"] == "uiuc_apcsp_11x7":
        src_keys.append("lsat3")
    if cd["family"] == "apc":
        src_keys.append("mccrink2017")

    # ---------------- case.json
    if cd["family"] == "da":
        notes = (DA_COMMON_NOTE +
                 "Geometry = UIUC 'tested' PropellerScanner geometry (both blades averaged) with the tip artefact "
                 "corrected (see geometry.csv header); PropellerScanner may underestimate twist (Uhlig, cited by Deters "
                 "et al. 2014), and the scanned twist is 1-2 deg above design over the outer blade, so geometry_design.csv "
                 "(in extraction/) is a useful sensitivity check. Root: r/R 0.15-0.25 is the blade neck between hub block "
                 "and blade (chord drops to c/R ~0.11-0.14); the solver blade starts at the first station (Rhub = "
                 "0.15 R). Blade material: 3D-printed plastic, flexible at higher rpm (Deters et al. 2014 say the 12 %% "
                 "airfoils were chosen to stiffen the blades). Polars: XFOIL Ncrit = 9 (tunnel Tu < 0.1 %%); Ncrit = 5 "
                 "alternative files are listed in each airfoil entry. Chord Re at 0.75R ranges ~%s. Fluid: UIUC lab "
                 "air assumed rho = 1.17 kg/m3, nu = 1.55e-5 m2/s (Urbana, ~225 m altitude, ~22 C; the database "
                 "gives coefficients only)." % re_range(cd, g, R))
    else:
        notes = (
            "Commercial APC propeller tested in UIUC Vol. 1 (Brandt 2005 / Brandt & Selig 2011, data ver. 3 with the "
            "2015 motor-fixture drag correction). Geometry = APC's own design file (current release v2025-1001); APC "
            "revises designs, so the specimen tested in 2005-2008 may differ slightly: the UIUC PropellerScanner geometry "
            "of the tested specimen is in extraction/geometry_uiuc_scan.csv and is compared in "
            "extraction/check_%s.png. Airfoils are KNOWN from the APC file (not assumed): %s. ASSUMPTIONS in the "
            "airfoil model: (1) APC 'scales airfoils based on thickness ratio' - taken as thickness-only scaling "
            "(camber line kept); (2) the AIRFOIL1->AIRFOIL2 transition is a linear blend in r of camber line and "
            "thickness shape; (3) sections are modelled by a few XFOIL anchor sections with linear cl/cd transitions "
            "between them. Published BEM validations of APC thin-electric props used NACA 4412 outboard (McCrink & "
            "Gregory 2017: Clark-Y 0-5 %% span, NACA 4412 elsewhere, from cutting the blades), consistent with APC's "
            "APC12 = NACA 4412. Polars: XFOIL Ncrit = 9 (UIUC tunnel Tu < 0.1 %%), Ncrit = 5 alternatives listed. "
            "Rhub = first geometry station (0.15 R). Fluid: rho = 1.17 kg/m3, nu = 1.55e-5 m2/s assumed (UIUC lab, "
            "the database gives coefficients only). Chord Re at 0.75R ranges ~%s." % (
                cd["id"], foil_lib["apc_desc"][cd["id"]], re_range(cd, g, R)))
    if cd["B"] > 2:
        notes += " This folder is the %d-blade configuration of the same blades (separate case because B is a " \
                 "case-level value in the runner)." % cd["B"]
    meta = {
        "id": cd["id"], "name": cd["name"], "turbine_type": 1, "B": cd["B"], "R": round(R, 6), "Rhub": round(Rhub, 6),
        "fluid": FLUID, "geometry_file": "geometry.csv",
        "alternative_geometry": ({"scan_raw": "extraction/geometry_scan_raw.csv",
                                  "design": "extraction/geometry_design.csv"} if cd["family"] == "da" else
                                 {"uiuc_scan": "extraction/geometry_uiuc_scan.csv"}),
        "airfoils": af_meta, "experiments": experiments,
        "sources": [SOURCES[k] for k in dict.fromkeys(src_keys)],
        "notes": notes,
    }
    with open(case / "case.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)
        f.write("\n")
    return case


def re_range(cd, g, R):
    c75 = float(np.interp(0.75 * R, g.r_m, g.chord_m))
    rpms = []
    vol = cd["vol"]
    return "%.0fk-%.0fk (c75 = %.1f mm, 0.75*Omega*R at the lowest/highest tested rpm)" % tuple(
        [x / 1e3 for x in cd["_re75"]] + [c75 * 1e3]) if "_re75" in cd else "see experiment files"


def osu_10x7(raw):
    import scipy.io as sio
    base = raw / "UIUC-propDB/comparisons/OSU_UIUC_PROP_DATA/x/00_all_OSU_UIUC_PROP_DATA"
    out = {}
    for q, f in (("ct", "CT_OSU_UIUC.fig"), ("cp", "CP_OSU_UIUC.fig")):
        d = sio.loadmat(base / f, squeeze_me=True, struct_as_record=False)

        def walk(h, acc):
            if getattr(h, "type", None) == "graph2d.lineseries":
                acc.append((h.properties.DisplayName, np.atleast_1d(h.properties.XData),
                            np.atleast_1d(h.properties.YData)))
            ch = getattr(h, "children", None)
            if ch is not None:
                for c in np.atleast_1d(ch):
                    walk(c, acc)
            return acc
        for n, x, y in walk(d["hgS_070000"], []):
            if n == "10x7_{OSU}":
                m = ~np.isnan(y)
                out[q] = pd.Series(y[m], index=np.round(x[m], 4))
    df = pd.concat(out, axis=1).dropna().sort_index()
    df.index.name = "J"
    df = df.reset_index()
    df["eta"] = df.J * df.ct / df.cp
    return df


# ------------------------------------------------------------------------------------------------ airfoil library
def build_foil_lib(raw, work, skip_xfoil=False):
    lib = {"apc_anchor": {}, "apc_desc": {}}
    cache = work / "xfoil_cache"
    cache.mkdir(parents=True, exist_ok=True)
    datdir = work / "dat"
    datdir.mkdir(parents=True, exist_ok=True)

    def add(name, xy, title, srcnote, tc=None, note=None):
        dat = datdir / (name + ".dat")
        write_dat(dat, title, xy)
        lib[name] = {"dat": dat, "title": title, "srcnote": srcnote, "thickness_ratio": tc, "note": note}

    sda = parse_sda(raw / "UIUC-propDB/papers/DetersAnandaSelig-2014-AIAA-2014-2151.txt")
    for n, xy in sda.items():
        xyn, ang = normalise(xy)
        add(n, xyn, "%s  (Deters, Ananda & Selig 2014, AIAA 2014-2151, Appendix A Table %d; normalised: LE at "
                    "origin, chord rotated %.3f deg to x axis)" % (n.upper(), {"sda1045": 3, "sda1075": 4, "sda1100": 5}[n], -ang),
            "coordinates tabulated in Deters et al. 2014 Appendix A", tc=0.12)

    e63, ang_e = normalise(read_coords(raw / "foils/e63.dat"))
    clk, ang_c = normalise(read_coords(raw / "foils/clarky.dat"))
    n4412 = naca4(0.04, 0.4, 0.12)
    base = {"e63": e63, "clarky": clk, "naca4412": n4412}
    base_t = {k: thickness(v) for k, v in base.items()}

    for cd in CASE_DEFS:
        if cd["family"] != "apc":
            continue
        pe0, pm = read_pe0(raw / "apc" / "PE0-FILES_WEB" / cd["pe0"])
        rR_all = (pe0.station_in / pm["RADIUS"]).to_numpy()
        sel = rR_all >= 0.15 - 1e-9
        rR = rR_all[sel]
        tcs = pe0.thickness_ratio.to_numpy()[sel]
        (r1, n1), (r2, n2) = pm["airfoils"][0], pm["airfoils"][1]
        names, desc = [], []
        for a in cd["anchors_rR"]:
            i = int(np.argmin(np.abs(rR - a)))
            r_in = rR[i] * pm["RADIUS"]
            f = 0.0 if cd["blend"][0] == cd["blend"][1] else float(np.clip((r_in - r1) / (r2 - r1), 0, 1))
            tc = float(tcs[i])
            b1, b2 = cd["blend"]
            if f == 0.0 and abs(tc - base_t[b1]) < 0.0015:
                name = b1
                xy = base[b1]
                title = {"e63": "E63 (4.25%%) UIUC coordinate database e63.dat, normalised (rotated %.3f deg)" % -ang_e,
                         "clarky": "CLARK Y UIUC coordinate database clarky.dat, normalised"}[b1]
            else:
                if f == 0.0:
                    name = "%s_t%s" % (b1, ("%.1f" % (tc * 100)).replace(".", "p"))
                    title = "%s thickness-scaled to t/c = %.4f (camber kept), APC %s r/R = %.3f" % (
                        b1.upper(), tc, cd["pe0"], rR[i])
                else:
                    name = "%s%s_f%s_t%s" % (b1, b2, ("%.2f" % f).replace(".", "p"), ("%.1f" % (tc * 100)).replace(".", "p"))
                    title = "blend %s->%s fraction %.3f, t/c = %.4f, APC %s r/R = %.3f" % (
                        b1.upper(), "APC12(NACA4412)" if b2 == "naca4412" else b2.upper(), f, tc, cd["pe0"], rR[i])
                xy = blend_section(base[b1], base[b2], f, tc)
            if name not in lib:
                add(name, xy, title, "generated: " + title, tc=round(tc, 4),
                    note="APC section at r/R = %.3f of %s: %s" % (rR[i], cd["pe0"], title))
            names.append(name)
            desc.append("%s at r/R=%.3f (f=%.2f, t/c=%.3f)" % (name, rR[i], f, tc))
        lib["apc_anchor"][cd["id"]] = names
        lib["apc_desc"][cd["id"]] = "APC file: AIRFOIL1 %s at %.2f in, AIRFOIL2 %s at %.2f in; anchors: %s" % (
            n1, r1, n2, r2, "; ".join(desc))
    # native Clark-Y 11.7 % for the XFOIL-vs-measurement check
    if "clarky" not in lib:
        add("clarky", clk, "CLARK Y UIUC coordinate database clarky.dat, normalised", "UIUC coordinate database",
            tc=round(thickness(clk), 4))

    # polars
    names = [k for k in lib if k not in ("apc_anchor", "apc_desc")]
    if not skip_xfoil:
        # warm the cache with one pool over every foil / Re / Ncrit (Ncrit 9 first)
        with cf.ThreadPoolExecutor(max_workers=XFOIL_WORKERS) as ex:
            futs = [ex.submit(xfoil_polar, lib[n]["dat"], re_, nc, cache)
                    for nc in NCRITS for n in names for re_ in RE_LIST]
            for i, f in enumerate(cf.as_completed(futs)):
                f.result()
                if i % 20 == 0:
                    print("  xfoil %d/%d" % (i, len(futs)), flush=True)
    for n in names:
        info = lib[n]
        print(" polars", n, flush=True)
        if skip_xfoil:
            continue
        p9, kept9 = build_polars(info["dat"], n, cache, 9.0)
        p5, kept5 = build_polars(info["dat"], n, cache, 5.0)
        hdr = ["%s polars: XFOIL 6.97, viscous, Mach 0, free transition, Ncrit = 9, PANE default repanelling (160 "
               "panels), sweeps 0..+20 and -0.5..-10 deg in 0.5 deg steps from a fresh boundary layer, plus fresh "
               "restarts at 8, 14 and -5 deg where the continuation sweep left gaps; only converged points kept." % n] + [
            "Coordinates: airfoils/%s.dat (%s)." % (n, info["title"]),
            "Ncrit = 9 chosen for the UIUC low-turbulence tunnel (Tu < 0.1 %, Deters et al. 2014; Brandt & Selig 2011).",
            "Re kept: %s. XFOIL is least reliable at Re < 40k (large laminar separation bubbles, hysteresis)." % (
                ", ".join("%g" % r for r in kept9))]
        info["polars"] = p9
        info["polar_header"] = hdr
        info["polar_source"] = "XFOIL 6.97, Ncrit=9, Re %s (computed; no measured polar exists for this section)" % (
            "/".join("%gk" % (r / 1e3) for r in kept9))
        info["alt"] = p5
        info["alt_header"] = [h.replace("Ncrit = 9", "Ncrit = 5") for h in hdr[:1]] + hdr[1:2] + [
            "Alternative (higher-disturbance) transition setting for sensitivity studies.",
            "Re kept: %s." % ", ".join("%g" % r for r in kept5)]
    lsat = clarky_lsat(raw)
    lib["clarky_uiuc_lsat"] = {
        "dat": lib["clarky"]["dat"], "polars": lsat, "used": False,
        "polar_header": ["Clark-Y measured polar: UIUC Low-Speed Airfoil Tests, Summary of Low-Speed Airfoil Data Vol. 3 "
                         "(Lyon et al. 1998), data file DRAG03.TXT, airfoil 'Clark-Y (B)' (clean, model by J. Robertson); "
                         "tabulated alpha, cl, cd (drag-polar runs, cd from wake rake, spanwise average).",
                         "Runs at Re 59.5k/60.7k are grouped as re=60000 and 99k/100k as re=100000; measured range only "
                         "(no post-stall). The coordinates file is the nominal 11.7 % Clark-Y (UIUC clarky.dat).",
                         "Not referenced by geometry.csv: provided to check/replace the XFOIL Clark-Y-family polars."],
        "polar_source": "measured, UIUC LSAT Vol. 3 (Lyon et al. 1998), Clark-Y (B), t/c 11.7 %, Re 60k-400k",
        "thickness_ratio": 0.117,
        "note": "measured reference polar (nominal Clark-Y 11.7 %); geometry.csv uses thickness-scaled XFOIL sections"}
    return lib


def thickness(xy):
    xg = 0.5 * (1 - np.cos(np.linspace(0, math.pi, 201)))
    u, l = surfaces(xy, xg)
    return float((u - l).max())


# ------------------------------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default=os.environ.get("UIUC_RAW", "/tmp/uiuc"))
    ap.add_argument("--work", default=os.environ.get("UIUC_WORK", "/tmp/uiuc/build"))
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--skip-xfoil", action="store_true")
    a = ap.parse_args()
    raw, work = Path(a.raw), Path(a.work)
    lib = build_foil_lib(raw, work, a.skip_xfoil)
    for cd in CASE_DEFS:
        if a.only and cd["id"] not in a.only:
            continue
        cd["_re75"] = re75(cd, raw)
        print("case", cd["id"], flush=True)
        build_case(cd, raw, lib, a)


def re75(cd, raw):
    """Chord Reynolds number at 0.75R based on 0.75*Omega*R, lowest and highest tested rpm (excl. static)."""
    vol = cd["vol"]
    datadir = raw / "UIUC-propDB" / ("volume-%d" % vol) / "data"
    D = cd["D_in"] * IN2M
    R = D / 2
    sc = read_table(datadir / cd["scan"])
    c75 = float(np.interp(0.75, sc["r/R"], sc["c/R"])) * R
    rpms = []
    sets = cd.get("sets") or [(cd["size"], 0, None)]
    for cfg, _, _ in sets:
        runs, _ = find_runs(datadir, "%s_%s" % (cd["prop"], cfg))
        rpms += [r[2] for r in runs]
    lo, hi = min(rpms), max(rpms)
    f = lambda n: 0.75 * n * 2 * math.pi / 60 * R * c75 / FLUID["kin_viscosity"]
    return (f(lo), f(hi))


if __name__ == "__main__":
    main()

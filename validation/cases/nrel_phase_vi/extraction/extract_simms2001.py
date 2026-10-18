"""Extract the MEASURED curves from the vector figures of Simms et al. 2001 (NREL/TP-500-29494),
Figures 8-21 (pages 27-40): upwind, 0 deg yaw, sequence S, wind speeds 7, 10, 13, 15, 20, 25 m/s.

Each figure has four panels (performance / aeroelastic / wake / CFD codes); the measured data are drawn in
every panel as the thick (1.32 pt) polyline with filled black diamonds. The figures are vector graphics
(Excel charts), so the polyline vertices are read directly with PyMuPDF (page.get_drawings()) - no raster
digitising. Each panel is calibrated from its own tick marks (short 0.12-pt strokes on the axes) and the
numeric tick labels next to them. The four copies of the measured curve are extracted independently and
compared; their spread is the extraction uncertainty.

Run:  python3 extract_simms2001.py   (from this folder; reads ../sources/29494.pdf)
Writes simms2001_measured_raw.csv (every panel) and prints the per-figure consensus + spread.
"""
import csv
import os
import re

import fitz  # PyMuPDF
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "..", "sources", "29494.pdf")

FIGURES = {  # page (1-based) -> (figure, quantity)
    27: ("Fig. 8", "lss_torque_nm"),
    28: ("Fig. 9", "root_flap_moment_nm"),
    29: ("Fig. 10", "cn_095"),
    30: ("Fig. 11", "ct_095"),
    31: ("Fig. 12", "cm_095"),
    32: ("Fig. 13", "q_local_pa_095"),
    33: ("Fig. 14", "cn_030"),
    34: ("Fig. 15", "ct_030"),
    35: ("Fig. 16", "cm_030"),
    36: ("Fig. 17", "q_local_pa_030"),
    37: ("Fig. 18", "cn_063"),
    38: ("Fig. 19", "ct_063"),
    39: ("Fig. 20", "cm_063"),
    40: ("Fig. 21", "q_local_pa_063"),
}

NUM = re.compile(r"^-?\d+(\.\d+)?$")


def segs(drawing):
    out = []
    for it in drawing["items"]:
        if it[0] == "l":
            out.append((it[1].x, it[1].y, it[2].x, it[2].y))
    return out


def fit(pos, val):
    A = np.vstack([np.ones(len(pos)), pos]).T
    coef, *_ = np.linalg.lstsq(A, np.asarray(val, float), rcond=None)
    resid = np.asarray(val) - A @ coef
    return coef, float(np.max(np.abs(resid))) if len(resid) else 0.0


def calibrate(grid, words):
    """grid: the 0.12-pt drawing of one panel. Returns maps (wind(y), value(x)) in unrotated page coords."""
    s = segs(grid)
    r = grid["rect"]
    # short ticks
    vticks = [(x1, y1, x2, y2) for x1, y1, x2, y2 in s if abs(x1 - x2) < 0.01 and 0.5 < abs(y1 - y2) < 4]
    hticks = [(x1, y1, x2, y2) for x1, y1, x2, y2 in s if abs(y1 - y2) < 0.01 and 0.5 < abs(x1 - x2) < 4]
    # value axis: ticks are vertical strokes below the plot at x positions; labels below them
    vx = sorted({round(t[0], 3) for t in vticks})
    vy = max(max(t[1], t[3]) for t in vticks)
    val_pos, val_val = [], []
    for x in vx:
        cands = [w for w in words if NUM.match(w[4]) and w[1] > vy - 0.5 and w[1] < vy + 12
                 and abs(0.5 * (w[0] + w[2]) - x) < 4]
        if len(cands) == 1:
            val_pos.append(x)
            val_val.append(float(cands[0][4]))
    # wind axis: horizontal strokes at y positions on the right of the plot; labels right of them
    wy = sorted({round(t[1], 3) for t in hticks})
    wx = max(max(t[0], t[2]) for t in hticks)
    w_pos, w_val = [], []
    for y in wy:
        cands = [w for w in words if NUM.match(w[4]) and w[0] > wx - 0.5 and w[0] < wx + 12
                 and abs(0.5 * (w[1] + w[3]) - y) < 4]
        if len(cands) == 1:
            w_pos.append(y)
            w_val.append(float(cands[0][4]))
    cv, ev = fit(val_pos, val_val)
    cw, ew = fit(w_pos, w_val)
    return cw, cv, (ew, ev, len(w_pos), len(val_pos))


def main():
    doc = fitz.open(PDF)
    rows = []
    summary = {}
    for page_no, (fig, qty) in FIGURES.items():
        page = doc[page_no - 1]
        dr = page.get_drawings()
        words = page.get_text("words")
        grids = [x for x in dr if abs((x.get("width") or 0) - 0.12) < 0.01 and len(x["items"]) > 10]
        thick = [x for x in dr if (x.get("width") or 0) > 1.0]
        diamonds = [x for x in dr if x["type"] == "f" and x.get("fill") == (0.0, 0.0, 0.0) and len(x["items"]) == 4]
        panels = []
        for k, t in enumerate(sorted(thick, key=lambda x: (round(x["rect"].x0 / 200), x["rect"].y0))):
            g = [gg for gg in grids if gg["rect"].contains(t["rect"]) or gg["rect"].intersects(t["rect"])]
            g = min(g, key=lambda gg: gg["rect"].get_area())
            cw, cv, err = calibrate(g, words)
            pts = []
            for it in t["items"]:
                for pt in (it[1], it[2]):
                    if not any(abs(pt.x - a) < 1e-3 and abs(pt.y - b) < 1e-3 for a, b in pts):
                        pts.append((pt.x, pt.y))
            # check against filled diamonds inside this panel
            dcent = [((x["rect"].x0 + x["rect"].x1) / 2, (x["rect"].y0 + x["rect"].y1) / 2)
                     for x in diamonds if g["rect"].contains(x["rect"])]
            res = []
            for (x, y) in pts:
                v = cw[0] + cw[1] * y
                q = cv[0] + cv[1] * x
                dmin = min((np.hypot(x - a, y - b) for a, b in dcent), default=np.nan)
                res.append((v, q, dmin))
                rows.append(dict(page=page_no, figure=fig, quantity=qty, panel=k, v_ms=round(v, 3),
                                 value=round(q, 5), marker_offset_pt=round(dmin, 2),
                                 calib_resid_wind=round(err[0], 4), calib_resid_value=round(err[1], 5),
                                 n_wind_ticks=err[2], n_value_ticks=err[3],
                                 value_per_pt=round(abs(cv[1]), 5)))
            panels.append(sorted(res))
        arr = np.array([[r[1] for r in p] for p in panels])
        varr = np.array([[r[0] for r in p] for p in panels])
        summary[qty] = (fig, varr.mean(0), arr.mean(0), arr.max(0) - arr.min(0), abs(cv[1]))
    with open(os.path.join(HERE, "simms2001_measured_raw.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    for qty, (fig, v, m, spread, scale) in summary.items():
        print(f"{fig:8s} {qty:22s} v=" + " ".join(f"{a:6.2f}" for a in v))
        print(f"{'':31s} mean=" + " ".join(f"{a:8.4f}" for a in m))
        print(f"{'':31s} spread(4 panels)=" + " ".join(f"{a:7.4f}" for a in spread) + f"   value/pt={scale:.4g}")
    return summary


if __name__ == "__main__":
    main()

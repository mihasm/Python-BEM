"""Extract the MEASURED data (mean +- one standard deviation) from the vector figures of
Sørensen, Michelsen & Schreck 2002, Wind Energy 5:151-169, doi:10.1002/we.64:

  Fig. 3 (p. 6)  low-speed shaft torque vs wind speed         measured: '+' marker + solid error bar
  Fig. 4 (p. 7)  root flap moment vs wind speed               measured: '+' marker + solid error bar
  Fig. 5 (p. 7)  root edge moment vs wind speed               measured: '+' marker + solid error bar
  Fig. 6 (p. 8)  Cn vs r/R, cases S0700000 ... S2500000       measured: circle + dotted error bar
  Fig. 7 (p. 9)  Ct vs r/R, same cases                        measured: circle + dotted error bar

The figures are vector graphics; marker centres and error-bar end points are read with PyMuPDF and mapped
to data units with each panel's own tick marks and tick labels (linear least-squares fit, residual printed).
Computed CFD curves (triangles / squares) are ignored. Legend samples are excluded.

Run: python3 extract_sorensen2002.py   (reads ../sources/Sorensen2002_NS_predictions_PhaseVI.pdf)
Writes sorensen2002_measured_raw.csv.
"""
import csv
import os
import re

import fitz
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "..", "sources", "Sorensen2002_NS_predictions_PhaseVI.pdf")
NUM = re.compile(r"^-?\d+(\.\d+)?$")


def num(s):
    s = s.replace("\u2212", "-").replace("\u2013", "-")
    return float(s) if NUM.match(s) else None


def lines(dr):
    out = []
    for d in dr:
        for it in d["items"]:
            if it[0] == "l":
                out.append(dict(x1=it[1].x, y1=it[1].y, x2=it[2].x, y2=it[2].y, w=d.get("width") or 0,
                                dashed=bool(d.get("dashes") and d["dashes"] not in ("[] 0",)), d=d))
    return out


def fit(pos, val):
    A = np.vstack([np.ones(len(pos)), pos]).T
    c, *_ = np.linalg.lstsq(A, np.asarray(val, float), rcond=None)
    return c, float(np.max(np.abs(np.asarray(val) - A @ c)))


def calibrate(frame, L, words):
    x0, y0, x1, y1 = frame.x0, frame.y0, frame.x1, frame.y1
    yt, xt = set(), set()
    for s in L:
        if abs(s["y1"] - s["y2"]) < 1e-3 and 0.5 < abs(s["x1"] - s["x2"]) < 3.5:  # horizontal tick
            if min(abs(s["x1"] - x0), abs(s["x2"] - x0)) < 0.05 and y0 - 0.1 <= s["y1"] <= y1 + 0.1:
                yt.add(round(s["y1"], 3))
        if abs(s["x1"] - s["x2"]) < 1e-3 and 0.5 < abs(s["y1"] - s["y2"]) < 3.5:  # vertical tick
            if min(abs(s["y1"] - y1), abs(s["y2"] - y1)) < 0.05 and x0 - 0.1 <= s["x1"] <= x1 + 0.1:
                xt.add(round(s["x1"], 3))
    yp, yv, xp, xv = [], [], [], []
    for y in sorted(yt):
        c = [w for w in words if num(w[4]) is not None and w[2] < x0 and w[2] > x0 - 25
             and abs(0.5 * (w[1] + w[3]) - y) < 3.5]
        if len(c) == 1:
            yp.append(y); yv.append(num(c[0][4]))
    for x in sorted(xt):
        c = [w for w in words if num(w[4]) is not None and w[1] > y1 and w[1] < y1 + 12
             and abs(0.5 * (w[0] + w[2]) - x) < 4]
        if len(c) == 1:
            xp.append(x); xv.append(num(c[0][4]))
    cx, ex = fit(xp, xv)
    cy, ey = fit(yp, yv)
    return cx, cy, ex, ey, len(xp), len(yp)


def crosses(L):
    """'+' markers: a vertical and a horizontal stroke of equal length (~4.7 pt) with a common midpoint."""
    v = [s for s in L if abs(s["x1"] - s["x2"]) < 1e-3 and 4.0 < abs(s["y1"] - s["y2"]) < 5.5]
    h = [s for s in L if abs(s["y1"] - s["y2"]) < 1e-3 and 4.0 < abs(s["x1"] - s["x2"]) < 5.5]
    out = []
    for a in v:
        mx, my = a["x1"], 0.5 * (a["y1"] + a["y2"])
        for b in h:
            if abs(0.5 * (b["x1"] + b["x2"]) - mx) < 0.05 and abs(b["y1"] - my) < 0.05:
                out.append((mx, my))
    return out


def circles(dr):
    out = []
    for d in dr:
        its = d["items"]
        if len(its) == 4 and all(i[0] == "c" for i in its):
            r = d["rect"]
            if 3.5 < r.width < 5.5:
                out.append(((r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2))
    return out


def errorbar(L, xm, ym, dashed):
    """vertical stroke through the marker centre (std bar); returns (top_y, bottom_y) or None."""
    best = None
    for s in L:
        if abs(s["x1"] - s["x2"]) > 1e-3 or abs(s["x1"] - xm) > 0.06:
            continue
        ln = abs(s["y1"] - s["y2"])
        if ln < 0.3 or (4.0 < ln < 5.5 and abs(0.5 * (s["y1"] + s["y2"]) - ym) < 0.05 and not s["dashed"]):
            continue  # the '+' marker itself
        if dashed is not None and s["dashed"] != dashed:
            continue
        lo, hi = min(s["y1"], s["y2"]), max(s["y1"], s["y2"])
        if not (lo - 0.05 <= ym <= hi + 0.05):
            continue
        asym = abs(0.5 * (lo + hi) - ym)
        if best is None or asym < best[2]:
            best = (lo, hi, asym)
    return best


def main():
    doc = fitz.open(PDF)
    rows = []
    jobs = [(6, "Fig. 3", ["lss_torque_nm"], "cross"),
            (7, "Fig. 4/5", ["root_flap_moment_nm", "root_edge_moment_nm"], "cross"),
            (8, "Fig. 6", ["cn"] * 6, "circle"),
            (9, "Fig. 7", ["ct"] * 6, "circle")]
    cases = ["S0700000", "S1000000", "S1300000", "S1500000", "S2000000", "S2500000"]
    for page_no, fig, qtys, kind in jobs:
        page = doc[page_no - 1]
        dr = page.get_drawings()
        words = page.get_text("words")
        L = lines(dr)
        frames = [d["rect"] for d in dr if len(d["items"]) == 1 and d["items"][0][0] == "qu"
                  and d["rect"].width > 100]
        # Fig. 6 panel (f) has its frame drawn as plain line segments inside one 33-item path
        frames += [d["rect"] for d in dr if len(d["items"]) > 20 and d["items"][0][0] == "l"
                   and d["rect"].width > 180 and 70 < d["rect"].height < 76]
        frames = sorted(frames, key=lambda r: r.y0)
        # legend sample lines (long horizontal strokes with no data meaning) -> exclusion zones
        legend = [s for s in L if abs(s["y1"] - s["y2"]) < 1e-3 and 30 < abs(s["x1"] - s["x2"]) < 45]
        markers = crosses(L) if kind == "cross" else circles(dr)
        for k, fr in enumerate(frames):
            cx, cy, ex, ey, nx, ny = calibrate(fr, L, words)
            qty = qtys[k]
            for (xm, ym) in sorted(markers):
                if not (fr.x0 - 1 <= xm <= fr.x1 + 1 and fr.y0 - 1 <= ym <= fr.y1 + 1):
                    continue
                if any(abs(ym - s["y1"]) < 1.0 and min(s["x1"], s["x2"]) - 3 <= xm <= max(s["x1"], s["x2"]) + 3
                       for s in legend):
                    continue
                eb = errorbar(L, xm, ym, dashed=True if kind == "circle" else False)
                xval = cx[0] + cx[1] * xm
                yval = cy[0] + cy[1] * ym
                if eb:
                    ytop = cy[0] + cy[1] * eb[0]
                    ybot = cy[0] + cy[1] * eb[1]
                    std = 0.5 * abs(ytop - ybot)
                    asym = cy[1] * eb[2]
                else:
                    std, asym = np.nan, np.nan
                rows.append(dict(page=page_no, figure=fig, panel=k,
                                 case=cases[k] if kind == "circle" else "",
                                 quantity=qty, x=round(xval, 4), value=round(yval, 5),
                                 std=round(std, 5) if std == std else "",
                                 bar_asym=round(abs(asym), 5) if asym == asym else "",
                                 calib_resid_x=round(ex, 5), calib_resid_y=round(ey, 5), n_xticks=nx, n_yticks=ny,
                                 value_per_pt=round(abs(cy[1]), 6)))
    with open(os.path.join(HERE, "sorensen2002_measured_raw.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(r["figure"], r["panel"], r["case"], r["quantity"], r["x"], r["value"], r["std"], "asym", r["bar_asym"],
              "res", r["calib_resid_x"], r["calib_resid_y"], r["n_xticks"], r["n_yticks"])


if __name__ == "__main__":
    main()

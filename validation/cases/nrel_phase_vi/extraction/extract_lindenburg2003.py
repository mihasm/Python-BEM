"""Extract the MEASURED diamonds from Lindenburg 2003, ECN-C--03-025, page 63:
  Figure 6.2  Stationary shaft torque [kNm] vs wind speed  (15 measured points, 6-25 m/s)
  Figure 6.3  Blade root flapping moment [kNm] vs wind speed (6 measured points)
The figures are vector graphics: the measured points are 'qu' (quad = diamond) items of one path; their
centres are mapped to data units with the major tick marks of each chart and the tick labels.

Run: python3 extract_lindenburg2003.py   (reads ../sources/Lindenburg2003_ECN-C-03-025.pdf)
Writes lindenburg2003_measured_raw.csv
"""
import csv
import os

import fitz
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "..", "sources", "Lindenburg2003_ECN-C-03-025.pdf")


def num(s):
    try:
        return float(s.replace("\u2212", "-"))
    except ValueError:
        return None


def main():
    page = fitz.open(PDF)[62]
    dr = page.get_drawings()
    words = page.get_text("words")
    charts = [("Fig. 6.2", "lss_torque_knm", 0, 280), ("Fig. 6.3", "root_flap_moment_knm", 300, 560)]
    rows = []
    for fig, qty, ylo, yhi in charts:
        segs = []
        for d in dr:
            if not (ylo <= d["rect"].y0 and d["rect"].y1 <= yhi):
                continue
            for it in d["items"]:
                if it[0] == "l":
                    segs.append((it[1].x, it[1].y, it[2].x, it[2].y))
        # major ticks: ~3.9 pt long, perpendicular to the frame
        xt = sorted({round(x1, 2) for x1, y1, x2, y2 in segs if abs(x1 - x2) < 1e-3 and 3.5 < abs(y1 - y2) < 4.5})
        yt = sorted({round(y1, 2) for x1, y1, x2, y2 in segs if abs(y1 - y2) < 1e-3 and 3.5 < abs(x1 - x2) < 4.5})
        ybot = max(max(y1, y2) for x1, y1, x2, y2 in segs if abs(x1 - x2) < 1e-3 and 3.5 < abs(y1 - y2) < 4.5)
        xleft = min(min(x1, x2) for x1, y1, x2, y2 in segs if abs(y1 - y2) < 1e-3 and 3.5 < abs(x1 - x2) < 4.5)
        xp, xv, yp, yv = [], [], [], []
        for x in xt:
            c = [w for w in words if num(w[4]) is not None and ybot < w[1] < ybot + 6 and abs((w[0] + w[2]) / 2 - x) < 4]
            if len(c) == 1:
                xp.append(x); xv.append(num(c[0][4]))
        for y in yt:
            c = [w for w in words if num(w[4]) is not None and xleft - 25 < w[0] < xleft and abs((w[1] + w[3]) / 2 - y) < 4]
            if len(c) == 1:
                yp.append(y); yv.append(num(c[0][4]))
        cx = np.polyfit(xp, xv, 1)
        cy = np.polyfit(yp, yv, 1)
        rx = np.max(np.abs(np.polyval(cx, xp) - xv))
        ry = np.max(np.abs(np.polyval(cy, yp) - yv))
        for d in dr:
            if not (ylo <= d["rect"].y0 and d["rect"].y1 <= yhi):
                continue
            for it in d["items"]:
                if it[0] != "qu":
                    continue
                q = it[1]
                xc = (q.ul.x + q.ur.x + q.ll.x + q.lr.x) / 4
                yc = (q.ul.y + q.ur.y + q.ll.y + q.lr.y) / 4
                v = float(np.polyval(cx, xc))
                val = float(np.polyval(cy, yc))
                # the legend sample sits next to the word 'Measured'
                leg = [w for w in words if w[4] == "Measured" and ylo <= w[1] <= yhi]
                is_legend = any(abs(yc - (w[1] + w[3]) / 2) < 3 and 0 < w[0] - xc < 15 for w in leg)
                rows.append(dict(figure=fig, quantity=qty, v_ms=round(v, 3), value=round(val, 5),
                                 legend=is_legend, calib_resid_x=round(rx, 4), calib_resid_y=round(ry, 5),
                                 n_xticks=len(xp), n_yticks=len(yp), value_per_pt=round(abs(cy[0]), 6)))
    with open(os.path.join(HERE, "lindenburg2003_measured_raw.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()

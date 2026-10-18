"""Exact extraction of Krogstad & Adaramola (2012), Wind Energy 15:743-756, doi:10.1002/we.502, Fig. 5
(performance at yaw angles 0..50 deg, pitch 0, U ~ 10 m/s) from the VECTOR drawing in the PDF (PyMuPDF).

PDF: sources/krogstad_adaramola_2012_we502.pdf, page index 5 (journal page 6 of the preprint layout).
NOTE: the figure caption swaps the panel names; the panel labels/axis titles are authoritative:
 (a) C_P vs lambda, (b) C_q vs lambda, (c) C_T vs lambda, (d) K_P vs J. We use (a) and (c).

Axis calibration = the plot frame lines, which coincide with the end ticks (checked against the tick-label
centres, within 0.4 pt):
 (a) x: 107.8 pt -> 0, 285.1 pt -> 12 ;  y: 205.0 pt -> 0.0, 72.1 pt -> 0.5
 (c) x: 331.0 pt -> 0, 509.4 pt -> 12 ;  y: 205.0 pt -> 0.0, 72.3 pt -> 1.4
Marker centre = centre of the marker path's bounding box. Series by marker style (legend):
 0 deg  black filled circle | 10 deg open star | 20 deg dark-red filled circle | 30 deg blue open square
 40 deg green open circle   | 50 deg red filled star
Precision: marker centres are exact to ~0.1 pt = 0.007 in lambda, 0.0004 in CP, 0.001 in CT.
"""
import os
import fitz
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "..", "sources", "krogstad_adaramola_2012_we502.pdf")
doc = fitz.open(PDF)
page = doc[5]

PANELS = {
    "cp": dict(x0=107.8, x1=285.1, xv=(0, 12), y0=205.0, y1=72.1, yv=(0.0, 0.5), legend=(110, 70, 150, 130)),
    "ct": dict(x0=331.0, x1=509.4, xv=(0, 12), y0=205.0, y1=72.3, yv=(0.0, 1.4), legend=(460, 140, 510, 200)),
}


def close(c, ref, tol=0.05):
    return c is not None and all(abs(a - b) < tol for a, b in zip(c, ref))


def style(d):
    col, fill, n = d.get("color"), d.get("fill"), len(d["items"])
    blk = (0.14, 0.12, 0.13)
    if d["type"] == "fs" and n == 4 and close(fill, (0, 0, 0)):
        return 0
    if d["type"] == "fs" and n == 10 and close(fill, (1, 1, 1)) and (close(col, (0, 0, 0)) or close(col, blk)):
        return 10
    if d["type"] == "fs" and n == 4 and close(fill, (0.5, 0, 0)):
        return 20
    if d["type"] == "s" and n == 4 and close(col, (0, 0, 1)) and d["items"][0][0] == "l" or \
            (d["type"] == "s" and close(col, (0, 0, 1)) and d["items"][0][0] == "re"):
        return 30
    if d["type"] == "s" and n == 4 and close(col, (0, 1, 0)) and d["items"][0][0] == "c":
        return 40
    if d["type"] == "fs" and n == 10 and close(fill, (1, 0, 0)):
        return 50
    return None


rows = []
drawings = page.get_drawings()
for name, P in PANELS.items():
    pts = {}
    for d in drawings:
        r = d["rect"]
        if not (2.0 < r.width < 5.0 and 2.0 < r.height < 5.0):
            continue
        cx, cy = (r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2
        if not (P["x0"] - 3 <= cx <= P["x1"] + 3 and P["y1"] - 3 <= cy <= P["y0"] + 3):
            continue
        lx0, ly0, lx1, ly1 = P["legend"]
        if lx0 <= cx <= lx1 and ly0 <= cy <= ly1:
            continue
        s = style(d)
        if s is None:
            continue
        x = P["xv"][0] + (cx - P["x0"]) / (P["x1"] - P["x0"]) * (P["xv"][1] - P["xv"][0])
        y = P["yv"][0] + (cy - P["y0"]) / (P["y1"] - P["y0"]) * (P["yv"][1] - P["yv"][0])
        pts.setdefault(s, set()).add((round(x, 3), round(y, 4)))
    for s, v in pts.items():
        for x, y in sorted(v):
            rows.append(dict(panel=name, yaw_deg=s, tsr=x, value=y))

df = pd.DataFrame(rows)
print(df.groupby(["panel", "yaw_deg"]).size())
df.to_csv(os.path.join(HERE, "krogstad_adaramola2012_fig5_raw.csv"), index=False)

# overlay
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
pix = page.get_pixmap(dpi=200, clip=fitz.Rect(70, 60, 520, 220))
img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)
s = 200 / 72
fig = plt.figure(figsize=(pix.width / 100, pix.height / 100), dpi=100)
ax = fig.add_axes([0, 0, 1, 1]); ax.imshow(img); ax.axis("off")
for name, P in PANELS.items():
    d = df[df.panel == name]
    px = P["x0"] + (d.tsr - P["xv"][0]) / (P["xv"][1] - P["xv"][0]) * (P["x1"] - P["x0"])
    py = P["y0"] + (d.value - P["yv"][0]) / (P["yv"][1] - P["yv"][0]) * (P["y1"] - P["y0"])
    ax.plot((px - 70) * s, (py - 60) * s, "+", color="magenta", ms=7, mew=1.0)
fig.savefig(os.path.join(HERE, "overlay_krogstad_adaramola2012_fig5.png"))

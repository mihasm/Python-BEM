"""Exact extraction of the MEASURED points of Krogstad & Eriksen (2013), Renewable Energy 50:325-333,
doi:10.1016/j.renene.2012.06.044, Fig. 3 (BT1 turbine T2, D = 0.894 m, U = 10 m/s, pitch 0) from the
vector drawing in the PDF (PyMuPDF), page index 3.

'Measured' = large black filled circles (4-segment Bezier paths, filled black; a second 11-segment path with
the same centre is the outline and is ignored). (a) CP markers are 3.7 pt, (b) CT markers 4.4 pt wide.
The legend marker ("Measured" entry, at px 410.7,165.7 in (a) and 460.1,360.1 in (b)) is removed by position.
Axis calibration = plot frame (coincides with the end ticks 0 and 12, -0.1 and 0.6, 0 and 1.4):
 (a) x 333.92 -> 0, 544.45 -> 12 ; y 224.65 -> -0.1, 75.22 -> 0.6
 (b) x 334.20 -> 0, 544.72 -> 12 ; y 422.89 -> 0.0, 273.45 -> 1.4
Precision ~0.05 pt: < 0.005 in TSR, < 0.001 in CP/CT.
"""
import os
import fitz
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
page = fitz.open(os.path.join(HERE, "..", "sources", "krogstad_eriksen_2013_BT1.pdf"))[3]
PAN = {
    "cp": dict(x0=333.92, x1=544.45, y0=224.65, y1=75.22, v0=-0.1, v1=0.6, size=(3.4, 4.0), legend=(405, 160, 416, 171)),
    "ct": dict(x0=334.20, x1=544.72, y0=422.89, y1=273.45, v0=0.0, v1=1.4, size=(4.1, 4.7), legend=(455, 355, 465, 365)),
}
rows = []
for name, P in PAN.items():
    for d in page.get_drawings():
        r = d["rect"]
        if len(d["items"]) != 4 or not d.get("fill") or max(d["fill"]) > 0.2:
            continue
        if not (P["size"][0] < r.width < P["size"][1] and P["size"][0] < r.height < P["size"][1]):
            continue
        cx, cy = (r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2
        if not (P["x0"] - 2 < cx < P["x1"] + 2 and P["y1"] - 2 < cy < P["y0"] + 2):
            continue
        lx0, ly0, lx1, ly1 = P["legend"]
        if lx0 < cx < lx1 and ly0 < cy < ly1:
            continue
        tsr = (cx - P["x0"]) / (P["x1"] - P["x0"]) * 12
        v = P["v0"] + (cy - P["y0"]) / (P["y1"] - P["y0"]) * (P["v1"] - P["v0"])
        rows.append(dict(panel=name, tsr=round(tsr, 3), value=round(v, 4), px=cx, py=cy))
df = pd.DataFrame(rows).sort_values(["panel", "tsr"])
print(df.to_string())
df.drop(columns=["px", "py"]).to_csv(os.path.join(HERE, "krogstad_eriksen2013_fig3_raw.csv"), index=False)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
clip = fitz.Rect(300, 60, 560, 450)
pix = page.get_pixmap(dpi=220, clip=clip)
img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)
k = 220 / 72
fig = plt.figure(figsize=(pix.width / 100, pix.height / 100), dpi=100)
ax = fig.add_axes([0, 0, 1, 1]); ax.imshow(img); ax.axis("off")
ax.plot((df.px - clip.x0) * k, (df.py - clip.y0) * k, "x", color="red", ms=8, mew=1.5)
fig.savefig(os.path.join(HERE, "overlay_krogstad_eriksen2013_fig3.png"))

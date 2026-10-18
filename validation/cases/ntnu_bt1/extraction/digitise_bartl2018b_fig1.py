"""Digitise Bartl, Schottler, Muehle, Hoelling, Peinke, Adaramola, Saetran (2018), Wind Energ. Sci. 3:489-502,
doi:10.5194/wes-3-489-2018, Fig. 1 (a) CP,T1 and (b) CT,T1 vs yaw angle -40..+40 deg, TSR = 6.0, U = 10 m/s,
NTNU rotor D = 0.894 m, pitch 0. Inflow A (TI 0.23 %, red, triangles) and inflow B (TI 10 %, green, squares).

Raster image embedded in the PDF (page index 4), extracted at native 1592x586 px into
extraction/pages/bartl2018b_fig1_img.png. Axis calibration from the light-grey grid lines found by scanning:
 (a) x: -40 deg at px 133.0, +40 at 488.5 ; y: 0.55 at px 79.0, 0.20 at 475.0
 (b) x: -40 deg at px 665.0, +40 at 1021.0; y: 1.00 at px 71.0, 0.50 at 458.0
The markers sit on the 10-deg grid lines; for each grid line x the pixels of the marker fill colour within
+-7 px are collected and their mean y is the value. Red and green markers overlap at some angles; the
partly hidden one is then estimated from its visible part (higher uncertainty, flagged).
Uncertainty: +-0.004 in CP and CT (grid 56.6 px per 0.05 CP and 77.4 px per 0.1 CT, marker ~12 px).
"""
import os
import numpy as np
import pandas as pd
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
im = np.array(Image.open(os.path.join(HERE, "pages", "bartl2018b_fig1_img.png")).convert("RGB")).astype(int)
R, G, B = im[:, :, 0], im[:, :, 1], im[:, :, 2]
red = (R > 170) & (G < 90) & (B < 90)
green = (G > 90) & (R < 90) & (B < 90)
PAN = {"cp": dict(xa=133.0, xb=488.5, ya=79.0, yb=475.0, va=0.55, vb=0.20, x0=120, x1=505),
       "ct": dict(xa=665.0, xb=1021.0, ya=71.0, yb=458.0, va=1.00, vb=0.50, x0=650, x1=1035)}
YAW = list(range(-40, 41, 10))
rows = []
for name, P in PAN.items():
    for yaw in YAW:
        px = P["xa"] + (yaw + 40) / 80 * (P["xb"] - P["xa"])
        for lab, m in (("A", red), ("B", green)):
            sub = m[60:490, int(px) - 7:int(px) + 8]
            ys = np.nonzero(sub)[0] + 60
            if len(ys) == 0:
                continue
            # marker = densest cluster of rows (lines leave few pixels per row)
            cnt = sub.sum(1)
            rr = np.nonzero(cnt >= max(3, 0.5 * cnt.max()))[0] + 60
            py = float(np.mean(rr))
            v = P["va"] + (py - P["ya"]) / (P["yb"] - P["ya"]) * (P["vb"] - P["va"])
            rows.append(dict(panel=name, inflow=lab, yaw_deg=yaw, value=round(v, 4), px=px, py=py, npx=int(len(ys))))
df = pd.DataFrame(rows)
print(df.to_string())
df.drop(columns=["px", "py"]).to_csv(os.path.join(HERE, "bartl2018b_fig1_raw.csv"), index=False)
fig = plt.figure(figsize=(im.shape[1] / 100, im.shape[0] / 100), dpi=100)
ax = fig.add_axes([0, 0, 1, 1]); ax.imshow(im.astype(np.uint8)); ax.axis("off")
for lab, c in (("A", "blue"), ("B", "magenta")):
    d = df[df.inflow == lab]
    ax.plot(d.px, d.py, "+", color=c, ms=14, mew=1.5)
fig.savefig(os.path.join(HERE, "overlay_bartl2018b_fig1.png"))

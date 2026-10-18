"""Digitise Fig. 2 of Nanos et al. (2022) (G06 chord c/R and twist vs r/R) from the embedded raster
(nanos2022_fig2_original.png, 2067x1400 px, image 0 of PDF page 6).
Calibration from gridlines found by pixel scan: x: r/R 0.2/0.4/0.6/0.8 at px 570.5/891.0/1211.5/1532.0;
y: c/R 0.05/0.10/0.15/0.20 (right twist axis 0..30 deg spans the same height as c/R 0..0.25, i.e. twist = 120*(c/R value)) at px 952.5/719.5/487.5/254.5,
axis line y=0 at px 1185. Curve y per column = mean of pixels of the curve colour (lines ~8 px thick).
"""
import os
import numpy as np
import cv2
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
img = cv2.imread(os.path.join(HERE, "nanos2022_fig2_original.png"))[:, :, ::-1].astype(int)
X = (570.5, 0.2, 1532.0, 0.8)
Y = (952.5, 0.05, 254.5, 0.20)
lin = lambda p, c: c[1] + (p - c[0]) * (c[3] - c[1]) / (c[2] - c[0])
inv = lambda v, c: c[0] + (v - c[1]) * (c[2] - c[0]) / (c[3] - c[1])
COL = {"chord": (230, 153, 0), "twist": (89, 179, 230)}
masks = {k: np.sqrt(((img - np.array(v)) ** 2).sum(-1)) < 40 for k, v in COL.items()}
for m in masks.values():
    m[:230, 1300:] = False     # legend
    m[1180:, :] = False        # axis
    m[:, :245] = False; m[:, 1853:] = False   # tick labels of both axes
masks['twist'][:1090, 1836:] = False   # blue right-axis line/label pixels above the tip part of the twist curve
curves = {}
for k, m in masks.items():
    xs, ys = [], []
    for x in range(250, 1853):
        yy = np.nonzero(m[:, x])[0]
        if len(yy):
            xs.append(x); ys.append(yy.mean())
    curves[k] = (np.array(xs), np.array(ys))
rR = np.round(np.concatenate([[0.06, 0.07, 0.08, 0.09], np.arange(0.10, 0.951, 0.025), [0.96, 0.97, 0.98, 0.99, 0.995]]), 4)
out = pd.DataFrame({"r_R": rR})
for k, (xs, ys) in curves.items():
    px = inv(rR, X)
    out[k + "_px_y"] = np.interp(px, xs, ys)
out["c_R"] = np.round(lin(out.chord_px_y, Y), 5)
out["twist_deg"] = np.round(lin(out.twist_px_y, Y) * 120, 3)    # right axis 0..30 deg spans left 0..0.25 -> x120
out.to_csv(os.path.join(HERE, "fig2_digitised.csv"), index=False)
# root part (for the notes): chord/twist near r/R 0..0.06
for k, (xs, ys) in curves.items():
    sel = xs < inv(0.07, X)
    print(k, "root samples r/R,value:", [(round(lin(x, X), 3), round(lin(y, Y) * (120 if k == "twist" else 1), 3)) for x, y in zip(xs[sel][::8], ys[sel][::8])])
ov = img[:, :, ::-1].astype(np.uint8).copy()
for _, r in out.iterrows():
    x = int(round(inv(r.r_R, X)))
    cv2.drawMarker(ov, (x, int(round(r.chord_px_y))), (0, 0, 0), cv2.MARKER_CROSS, 18, 2)
    cv2.drawMarker(ov, (x, int(round(r.twist_px_y))), (0, 0, 255), cv2.MARKER_CROSS, 18, 2)
cv2.imwrite(os.path.join(HERE, "fig2_overlay.png"), ov)
print(out[["r_R", "c_R", "twist_deg"]].to_string())

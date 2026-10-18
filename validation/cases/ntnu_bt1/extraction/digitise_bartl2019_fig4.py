"""Digitise Bartl, Sagmo, Bracchi & Saetran (2019), Eur. J. Mech. B/Fluids 75:180-192, Fig. 4
(NTNU measured S826 lift and drag, Re = 0.5e5 ... 6.0e5).

The figure is a raster image embedded in the PDF (sources/bartl_2019_s826_polars.pdf, page 18,
image extracted at its native 2248x1294 px into extraction/pages/bartl2019_p17_img.png).
Each Reynolds number is drawn in its own exact RGB colour (read from the legend), so the curves are
separated by colour. For every measured angle of attack (the marker positions, found from the filled
markers: -15,-10,-8,-6,-4,-2,0,2,4,6,8,9,...,18,20,25,30 deg) the pixels of the curve's colour in a
narrow column around the marker x position are collected and their median y is taken.
Where a curve is completely hidden under curves drawn later (high Re curves lie on top of each other),
the value is taken as the mean of the visible later-drawn curves at that angle and flagged 'occluded'.

Drag: the full lines are wake-rake (momentum deficit) values, the dashed lines force-balance values.
Following the paper (p.16-17: "At the onset of stall at alpha = 13 deg and beyond, drag results from
direct force measurements as indicated by dashed lines, are deemed to give a more accurate drag
estimation") cd is the wake-rake value for alpha < 13 deg and the force-balance value for alpha >= 13.
In a column with both, the force-balance (dashed) value is the upper (larger cd) cluster.

Axis calibration from the grid lines found by scanning (see find_grid below):
  (a) x: alpha=0 at px 481.0, alpha=10 at 662.25 ... -> 18.125 px/deg;  y: CL=1.6 at 165.5, CL=-0.4 at 1082.5
  (b) x: alpha=0 at px 1603.0, alpha=30 at 2147.5 -> 18.15 px/deg;     y: CD=0 at 1151, CD=0.7 at 126
"""
import numpy as np
import pandas as pd
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

HERE = os.path.dirname(os.path.abspath(__file__))
IMG = os.path.join(HERE, "pages", "bartl2019_p17_img.png")
im = np.array(Image.open(IMG).convert("RGB")).astype(int)

RE = [0.5e5, 0.7e5, 1.0e5, 2.0e5, 3.0e5, 4.0e5, 5.0e5, 6.0e5]
COL = [(0, 166, 255), (0, 133, 255), (0, 99, 255), (0, 66, 255), (0, 33, 255), (0, 0, 255), (0, 0, 177), (0, 0, 0)]
ALPHA = [-15, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 25, 30]

# calibration
ax_a = dict(x0=481.0, sx=18.125, y0=899.0, sy=-458.5)        # px = x0 + sx*alpha ; py = y0 + sy*CL
ax_b = dict(x0=1603.0, sx=18.15, y0=1151.0, sy=-(1151 - 126) / 0.7)


def mask(c, tol=30):
    return (np.abs(im - np.array(c)).sum(2) <= tol)


def column_values(m, px, half=4, ylim=(100, 1150), xlim=None):
    x0, x1 = int(round(px - half)), int(round(px + half)) + 1
    sub = m[ylim[0]:ylim[1], x0:x1]
    ys = np.nonzero(sub.any(axis=1))[0] + ylim[0]
    return ys


def marker_y(m, px, ys, half=9):
    """y of the marker centre inside a column cluster: rows where the colour run around px is widest
    (a marker is wider than the connecting lines, which matters on steep line segments)."""
    x0, x1 = int(round(px - half)), int(round(px + half)) + 1
    w = np.array([m[y, x0:x1].sum() for y in ys])
    sel = ys[w >= 0.7 * w.max()]
    return float(np.mean(sel))


def clusters(ys, gap=12):
    if len(ys) == 0:
        return []
    out, cur = [], [ys[0]]
    for y in ys[1:]:
        if y - cur[-1] > gap:
            out.append(cur)
            cur = [y]
        else:
            cur.append(y)
    out.append(cur)
    return [np.array(c) for c in out]


# legend boxes are excluded
legend_a = (650, 1130, 725, 1020)  # y0,y1,x0,x1
rows = []
masks = []
for c in COL:
    m = mask(c)
    if c == (0, 0, 0):  # black: exclude axes, ticks, text -> only keep inside plot interiors
        pass
    y0, y1, x0, x1 = legend_a
    m[y0:y1, x0:x1] = False
    m[:, :195] = False
    m[:, 1040:1318] = False
    m[:, 2162:] = False
    m[:100, :] = False
    m[1148:, :] = False
    masks.append(m)

res = {}
FORCE_MIN = {20: 0.33, 25: 0.45, 30: 0.55}   # solid (wake) lines end below these values, see figure
for i in reversed(range(len(RE))):
    m = masks[i]
    for a in ALPHA:
        # lift
        px = ax_a["x0"] + ax_a["sx"] * a
        cl = clusters(column_values(m, px))
        if cl:
            # choose the biggest cluster (marker), median y
            best = max(cl, key=len)
            cl_v, flag_l = (marker_y(m, px, best) - ax_a["y0"]) / ax_a["sy"], "visible"
        else:
            vis = [res[(j, a)]["cl"] for j in range(i + 1, len(RE)) if res[(j, a)]["flag_cl"] == "visible"]
            cl_v, flag_l = (np.mean(vis) if vis else np.nan), "occluded"
        # drag
        px = ax_b["x0"] + ax_b["sx"] * a
        cd_cl = clusters(column_values(m, px))
        cd_cl = [c for c in cd_cl if len(c) >= 3]
        if cd_cl:
            vals = [(marker_y(m, px, c) - ax_b["y0"]) / ax_b["sy"] for c in cd_cl]
            # beyond 20 deg the dashed (force) lines lie far above the solid (wake) lines: keep only those
            if a in FORCE_MIN:
                vals = [v for v in vals if v > FORCE_MIN[a]]
                cd_cl = [c for c, v in zip(cd_cl, [(marker_y(m, px, c) - ax_b["y0"]) / ax_b["sy"] for c in cd_cl]) if v > FORCE_MIN[a]]
        if cd_cl:
            if a >= 13 and len(vals) >= 2:
                cd_v, meth = max(vals), "force balance (dashed)"
            elif a >= 13:
                cd_v, meth = vals[0], "single cluster (wake=force or one hidden)"
            else:
                # below 13 deg only one line; pick biggest cluster
                cd_v, meth = vals[int(np.argmax([len(c) for c in cd_cl]))], "wake rake"
            flag_d = "visible"
        else:
            vis = [res[(j, a)]["cd"] for j in range(i + 1, len(RE)) if res[(j, a)]["flag_cd"] == "visible"]
            cd_v, meth, flag_d = (np.mean(vis) if vis else np.nan), ("force balance" if a >= 13 else "wake rake"), "occluded"
        res[(i, a)] = dict(cl=cl_v, cd=cd_v, flag_cl=flag_l, flag_cd=flag_d, meth=meth)

for i, re in enumerate(RE):
    for a in ALPHA:
        r = res[(i, a)]
        rows.append(dict(re=re, alpha_deg=a, cl=round(float(r["cl"]), 3), cd=round(float(r["cd"]), 4),
                         cl_digitised=r["flag_cl"], cd_digitised=r["flag_cd"], cd_method=r["meth"]))
df = pd.DataFrame(rows)
df.to_csv(os.path.join(HERE, "bartl2019_fig4_raw.csv"), index=False)
print(df.to_string())

# overlay
fig = plt.figure(figsize=(22.48, 12.94), dpi=100)
axo = fig.add_axes([0, 0, 1, 1])
axo.imshow(im.astype(np.uint8))
for i, re in enumerate(RE):
    d = df[df.re == re]
    xa = ax_a["x0"] + ax_a["sx"] * d.alpha_deg
    ya = ax_a["y0"] + ax_a["sy"] * d.cl
    xb = ax_b["x0"] + ax_b["sx"] * d.alpha_deg
    yb = ax_b["y0"] + ax_b["sy"] * d.cd
    oc = d.cl_digitised == "occluded"
    axo.plot(xa[~oc], ya[~oc], "o", mfc="none", mec="red", ms=9, mew=1.2)
    axo.plot(xa[oc], ya[oc], "s", mfc="none", mec="orange", ms=9, mew=1.2)
    oc = d.cd_digitised == "occluded"
    axo.plot(xb[~oc], yb[~oc], "o", mfc="none", mec="red", ms=9, mew=1.2)
    axo.plot(xb[oc], yb[oc], "s", mfc="none", mec="orange", ms=9, mew=1.2)
axo.axis("off")
fig.savefig(os.path.join(HERE, "overlay_bartl2019_fig4.png"))

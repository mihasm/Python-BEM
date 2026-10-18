"""Digitise Fig. 8 of Nanos et al. (2022), WES 7:1263 (G06 performance) from the embedded raster image
(nanos2022_fig8_original.png = image 0 of PDF page 12, 2067 x 1279 px, extracted losslessly with PyMuPDF).

Axis calibration from gridline pixel positions found by scanning for light-gray lines (see calib below),
markers found per curve colour. find_markers() was a first attempt (box-filter density, not used for the data);
the final points are run_method4() = method 2 (markers along each curve) + method 3 (hollow-outline search)
+ method 4 (hidden markers read from the curve line). Output: fig8_markers_final.csv, fig8_overlay_final.png
(black = method 2, magenta = 3, blue = 4).
"""
import os
import numpy as np
import cv2
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
img = cv2.imread(os.path.join(HERE, "nanos2022_fig8_original.png"))[:, :, ::-1].astype(int)  # RGB

COL = {  # RGB of the legend lines (sampled in the legend of panel a)
    "c0": (229, 147, 0),     # beta 0 / 6 m/s   (orange, circles)
    "c1": (90, 179, 230),    # beta 2 / 7 m/s   (light blue, diamonds)
    "c2": (0, 153, 128),     # beta 3 / 8 m/s   (green, stars)
    "c3": (20, 125, 220),    # beta 4 / 9 m/s   (blue, squares)
    "c4": (234, 95, 0),      # beta 6 / 10 m/s  (red-orange, triangles)
}
# pixel <-> data calibration: (px1, v1, px2, v2) from gridlines (light-gray line scan, see README of this folder)
CAL = {
    "a": dict(x=(136.0, 4.0, 599.5, 8.0), y=(416.5, 0.0, 91.5, 0.4), box=(136, 12, 657, 497), qty="cp"),
    "b": dict(x=(949.0, 5.0, 1296.5, 8.0), y=(476.5, 0.2, 88.0, 0.7), box=(833, 12, 1354, 497), qty="ct"),
    "c": dict(x=(1645.0, 5.0, 1991.5, 8.0), y=(495.0, 0.0, 91.5, 0.05), box=(1529, 12, 2050, 497), qty="cq"),
    "d": dict(x=(136.0, 4.0, 483.0, 8.0), y=(1176.0, 0.10, 707.0, 0.40), box=(136, 688, 657, 1176), qty="cp"),
    "e": dict(x=(1005.5, 6.0, 1180.0, 8.0), y=(1176.0, 0.4, 808.5, 0.7), box=(832, 688, 1354, 1176), qty="ct"),
    "f": dict(x=(1702.0, 6.0, 1875.5, 8.0), y=(1176.0, 0.02, 809.0, 0.05), box=(1529, 688, 2050, 1176), qty="cq"),
}
LEGEND_BOX = {"a": (243, 252, 400, 497), "d": (262, 930, 445, 1176)}  # exclude legend areas


def lin(p, cal):
    p1, v1, p2, v2 = cal
    return v1 + (p - p1) * (v2 - v1) / (p2 - p1)


def inv(v, cal):
    p1, v1, p2, v2 = cal
    return p1 + (v - v1) * (p2 - p1) / (v2 - v1)


def find_markers(panel, ckey, tol=34, min_density=0.22):
    cal = CAL[panel]
    x0, y0, x1, y1 = cal["box"]
    oy, ox = max(y0 - 15, 0), max(x0 - 15, 0)
    sub = img[oy:y1 + 15, ox:x1 + 15]
    d = np.sqrt(((sub - np.array(COL[ckey])) ** 2).sum(-1))
    mask = (d < tol).astype(np.uint8)
    if panel in LEGEND_BOX:
        lx0, ly0, lx1, ly1 = LEGEND_BOX[panel]
        mask[ly0 - oy:ly1 - oy, lx0 - ox:lx1 - ox] = 0
    dens = cv2.boxFilter(mask.astype(float), -1, (25, 25), normalize=True)
    pts = []
    work = dens.copy()
    while True:
        iy, ix = np.unravel_index(np.argmax(work), work.shape)
        if work[iy, ix] < min_density:
            break
        # refine: marker outline extents in a 31x31 window, using rows/cols through the peak
        wy0, wy1 = max(iy - 15, 0), iy + 16
        wx0, wx1 = max(ix - 15, 0), ix + 16
        win = mask[wy0:wy1, wx0:wx1]
        ys, xs = np.nonzero(win)
        cy = wy0 + 0.5 * (ys.min() + ys.max())
        cx = wx0 + xs.mean()
        pts.append((cx + ox, cy + oy, work[iy, ix]))
        work[max(iy - 18, 0):iy + 19, max(ix - 18, 0):ix + 19] = 0
    pts = sorted(pts)
    return [(lin(px, cal["x"]), lin(py, cal["y"]), px, py, dn) for px, py, dn in pts]


# ---------------------------------------------------------------------------------------------
# Method 2 (used for the final data): walk along each curve. For every pixel column the curve's y is the
# median of its colour pixels; the number of colour pixels in a 23x23 window centred on (x, y_curve) gives
# a 1-D profile along x whose peaks are the hollow markers (outline adds pixels to the line).
# Marker centre: x = peak (refined by symmetric extent), y = midpoint of the outline extent in the window.
def walk_curve(panel, ckey, tol=30, half=11, min_sep=30, min_gain=1.6):
    cal = CAL[panel]
    x0, y0, x1, y1 = cal["box"]
    oy, ox = max(y0 - 15, 0), max(x0 - 15, 0)
    sub = img[oy:y1 + 15, ox:x1 + 15]
    mask = (np.sqrt(((sub - np.array(COL[ckey])) ** 2).sum(-1)) < tol)
    if panel in LEGEND_BOX:
        lx0, ly0, lx1, ly1 = LEGEND_BOX[panel]
        mask[max(ly0 - oy, 0):ly1 - oy, max(lx0 - ox, 0):lx1 - ox] = False
    H, W = mask.shape
    ycur = np.full(W, np.nan)
    for x in range(W):
        ys = np.nonzero(mask[:, x])[0]
        if len(ys):
            ycur[x] = np.median(ys)
    prof = np.zeros(W)
    for x in range(W):
        if np.isnan(ycur[x]):
            continue
        yc = int(ycur[x])
        prof[x] = mask[max(yc - half, 0):yc + half + 1, max(x - half, 0):x + half + 1].sum()
    base = np.median(prof[prof > 0]) if (prof > 0).any() else 0
    pts = []
    work = prof.copy()
    while work.max() > min_gain * base:
        x = int(np.argmax(work))
        yc = int(ycur[x])
        win = mask[max(yc - half - 2, 0):yc + half + 3, max(x - half - 2, 0):x + half + 3]
        ys, xs = np.nonzero(win)
        cy = max(yc - half - 2, 0) + 0.5 * (ys.min() + ys.max())
        cx = max(x - half - 2, 0) + 0.5 * (xs.min() + xs.max())
        pts.append((cx + ox, cy + oy))
        work[max(x - min_sep, 0):x + min_sep + 1] = 0
    pts.sort()
    return [(lin(px, cal["x"]), lin(py, cal["y"]), px, py) for px, py in pts]


def run_method2(panels=("a", "b", "d", "e")):
    rows = []
    for panel in panels:
        for ck in COL:
            for lam, val, px, py in walk_curve(panel, ck):
                rows.append(dict(panel=panel, curve=ck, tsr=round(lam, 4), value=round(val, 5),
                                 px=round(px, 1), py=round(py, 1)))
    df = pd.DataFrame(rows)
    ov = img[:, :, ::-1].astype(np.uint8).copy()
    for _, r in df.iterrows():
        cv2.drawMarker(ov, (int(round(r.px)), int(round(r.py))), (0, 0, 0), cv2.MARKER_CROSS, 16, 2)
    cv2.imwrite(os.path.join(HERE, "fig8_overlay_method2.png"), ov)
    df.to_csv(os.path.join(HERE, "fig8_markers_method2.csv"), index=False)
    return df


# ---------------------------------------------------------------------------------------------
# Method 3: complete method 2 where curves overlap. Slide along the curve every 2 px (curve y = median of
# its colour pixels in the column). A position is 'marker-like' if, through it, BOTH the column and the
# row of this colour's mask span >= 14 px (window +-13 px) AND each is either >= 2 separate runs (hollow
# outline, with or without the line inside) or one run >= 13 px (filled marker / line along that axis).
# A bare line fails one of the two directions (it is only ~4-7 px thick across). Consecutive marker-like
# positions are grouped; the group's marker centre = extent midpoint of the colour pixels in the window.
# Union with method 2 (de-duplicated within 18 px).
def _runs(v):
    idx = np.nonzero(v)[0]
    if len(idx) == 0:
        return 0, 0, 0
    nr = 1 + int((np.diff(idx) > 1).sum())
    longest = max(len(g) for g in np.split(idx, np.where(np.diff(idx) > 1)[0] + 1))
    return nr, idx.max() - idx.min() + 1, longest


def ring_search(panel, ckey, tol=30, half=13):
    cal = CAL[panel]
    x0, y0, x1, y1 = cal["box"]
    m = (np.sqrt(((img - np.array(COL[ckey])) ** 2).sum(-1)) < tol)
    if panel in LEGEND_BOX:
        lx0, ly0, lx1, ly1 = LEGEND_BOX[panel]
        m[ly0:ly1, lx0:lx1] = False
    ok = []
    for xp in range(x0 - 5, x1 + 5, 2):
        ys = np.nonzero(m[y0 - 10:y1 + 10, xp])[0]
        if len(ys) == 0:
            continue
        for g in np.split(ys, np.where(np.diff(ys) > 25)[0] + 1):
            yc = int(np.median(g)) + y0 - 10
            col = m[yc - half:yc + half + 1, xp]
            row = m[yc, xp - half:xp + half + 1]
            good = True
            for v in (col, row):
                nr, span, longest = _runs(v)
                if not (span >= 14 and (nr >= 2 or longest >= 13)):
                    good = False
            if good:
                ok.append((xp, yc))
    out = []
    ok.sort()
    groups = []
    for xp, yc in ok:
        if groups and xp - groups[-1][-1][0] <= 4 and abs(yc - groups[-1][-1][1]) < 12:
            groups[-1].append((xp, yc))
        else:
            groups.append([(xp, yc)])
    for g in groups:
        if len(g) < 3:          # a marker is ~20 px wide -> several consecutive hits
            continue
        xp = int(round(np.mean([a for a, b in g]))); yc = int(round(np.mean([b for a, b in g])))
        win = m[yc - half:yc + half + 1, xp - half:xp + half + 1]
        wy, wx = np.nonzero(win)
        out.append((xp - half + 0.5 * (wx.min() + wx.max()), yc - half + 0.5 * (wy.min() + wy.max())))
    return out


def run_method3(panels=("a", "b", "d", "e")):
    m2 = run_method2(panels)
    rows = [dict(r) for _, r in m2.iterrows()]
    for r in rows:
        r["method"] = 2
    for panel in panels:
        cal = CAL[panel]
        for ck in COL:
            have = [(r["px"], r["py"]) for r in rows if r["panel"] == panel and r["curve"] == ck]
            for cx, cy in ring_search(panel, ck):
                if all(np.hypot(cx - hx, cy - hy) > 18 for hx, hy in have):
                    have.append((cx, cy))
                    rows.append(dict(panel=panel, curve=ck, tsr=round(lin(cx, cal["x"]), 4),
                                     value=round(lin(cy, cal["y"]), 5), px=round(cx, 1), py=round(cy, 1), method=3))
    df = pd.DataFrame(rows).sort_values(["panel", "curve", "tsr"])
    ov = img[:, :, ::-1].astype(np.uint8).copy()
    for _, r in df.iterrows():
        colr = (0, 0, 0) if r.method == 2 else (255, 0, 255)
        cv2.drawMarker(ov, (int(round(r.px)), int(round(r.py))), colr, cv2.MARKER_CROSS, 16, 2)
    cv2.imwrite(os.path.join(HERE, "fig8_overlay_final.png"), ov)
    df.to_csv(os.path.join(HERE, "fig8_markers_final.csv"), index=False)
    return df


# ---------------------------------------------------------------------------------------------
# Method 4 (fill-in for markers hidden under other curves): the tests were run at tsr set-points on a
# 0.5 grid (all automatically detected markers sit within +-0.06 of 4.0, 4.5, ... except beta=4 at 7.65-7.69,
# 8.07). For a grid set-point inside a curve's x-range with no detected marker within 0.15, the value is
# read from the curve's own polyline at that x (median y of its colour pixels in columns x-1..x+1, the
# cluster closest to the linear interpolation of the neighbouring detected markers). Flagged method = 4.
def run_method4(panels=("a", "b", "d", "e"), tol=30):
    df = run_method3(panels)
    rows = [dict(r) for _, r in df.iterrows()]
    for panel in panels:
        cal = CAL[panel]
        x0, y0, x1, y1 = cal["box"]
        for ck in COL:
            m = (np.sqrt(((img - np.array(COL[ck])) ** 2).sum(-1)) < tol)
            if panel in LEGEND_BOX:
                lx0, ly0, lx1, ly1 = LEGEND_BOX[panel]
                m[ly0:ly1, lx0:lx1] = False
            m[:y0 - 12] = False; m[y1 + 12:] = False
            m[:, :x0 - 12] = False; m[:, x1 + 12:] = False
            cols = np.nonzero(m.any(0))[0]
            det = sorted([(r["tsr"], r["value"]) for r in rows if r["panel"] == panel and r["curve"] == ck])
            if len(det) < 2:
                continue
            lo, hi = det[0][0], det[-1][0]          # only fill between detected end markers
            for lam in np.arange(4.0, 10.01, 0.5):
                if lam < lo - 0.01 or lam > hi + 0.01:
                    continue
                if min(abs(lam - t) for t, v in det) < 0.15:
                    continue
                xp = int(round(inv(lam, cal["x"])))
                ys = np.nonzero(m[:, xp - 1:xp + 2].any(1))[0]
                if len(ys) == 0:
                    continue
                guess = inv(np.interp(lam, [t for t, v in det], [v for t, v in det]), cal["y"])
                groups = np.split(ys, np.where(np.diff(ys) > 6)[0] + 1)
                g = min(groups, key=lambda g: abs(np.median(g) - guess))
                if abs(np.median(g) - guess) > 40:
                    continue
                py = float(np.median(g))
                rows.append(dict(panel=panel, curve=ck, tsr=round(lam, 4), value=round(lin(py, cal["y"]), 5),
                                 px=float(xp), py=py, method=4))
    df = pd.DataFrame(rows).sort_values(["panel", "curve", "tsr"])
    ov = img[:, :, ::-1].astype(np.uint8).copy()
    for _, r in df.iterrows():
        colr = {2: (0, 0, 0), 3: (255, 0, 255), 4: (255, 0, 0)}[int(r.method)]
        cv2.drawMarker(ov, (int(round(r.px)), int(round(r.py))), colr, cv2.MARKER_CROSS, 16, 2)
    cv2.imwrite(os.path.join(HERE, "fig8_overlay_final.png"), ov)
    df.to_csv(os.path.join(HERE, "fig8_markers_final.csv"), index=False)
    return df


if __name__ == "__main__":
    final = run_method4()
    print(final.groupby(["panel", "curve"]).size().unstack())

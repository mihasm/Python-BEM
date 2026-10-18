"""Digitise marker plots that exist only as raster images (Excel charts) in the New MEXICO reports.

Images were taken out of the PDFs unchanged (PyMuPDF Pixmap of the embedded image, no resampling):
  e14048_fig22_0.png : ECN-E--14-048 Fig. 22(a)  Cdax vs TSR, balance, MEXICO (black squares) / New MEXICO config 0 (red dots)
  e14048_fig22_1.png : ECN-E--14-048 Fig. 22(b)  Cdax vs TSR, pressures, same series
  e14048_fig30_0.png : ECN-E--14-048 Fig. 30(a)  Cdax vs TSR, balance, New MEXICO rough (config 0) / partially clean (config 3)
  e14048_fig30_1.png : ECN-E--14-048 Fig. 30(b)  same, pressures
  m16063_fig4_0.png  : Boorsma & Schepers 2016 (ECN-M--16-063, J. Phys. Conf. Ser. 753 022004) Fig. 4(b), Cp vs TSR,
                       rough/clean, from pressures (filled) and from generator torque (open)
Axis calibration: the light-grey gridlines are found as image columns/rows with many grey pixels, the axis
lines as long dark lines; their values are read from the tick labels (see CAL below, checked by eye on the
image). Markers: colour mask (black / red), then
  - filled markers: distance transform of the mask >= 2.6 px (connecting lines are thinner) -> blob centroids;
  - open markers: holes of the mask (regions not connected to the image border) of marker size -> centroids.
Legend areas are excluded by a box per image. The result is written to <image>_points.csv and an overlay
<image>_overlay.png; the overlay was inspected for every image.
"""
import os, sys, json
import numpy as np, cv2, pandas as pd
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))

# pixel -> value: two (pixel, value) pairs per axis taken from gridlines/axis lines
CAL = {
    'e14048_fig22_0.png': dict(x=[(59, 2.0), (691, 16.0)], y=[(595, 0.0), (10, 1.4)], legend=(380, 200, 716, 330),
                               series={'black_filled': 'MEXICO 2006 balance', 'red_filled': 'New MEXICO balance config 0'}),
    'e14048_fig22_1.png': dict(x=[(51, 2.0), (695, 16.0)], y=[(595, 0.0), (10, 1.4)], legend=(380, 200, 716, 330),
                               series={'black_filled': 'MEXICO 2006 pressure', 'red_filled': 'New MEXICO pressure config 0'}),
    'e14048_fig30_0.png': dict(x=[(53, 4.0), (689, 16.0)], y=[(594, 0.0), (9, 1.4)], legend=(380, 230, 716, 380),
                               series={'black_filled': 'New MEXICO balance config 0 (rough)',
                                       'red_filled': 'New MEXICO balance config 3 (clean outboard)'}),
    'e14048_fig30_1.png': dict(x=[(54, 4.0), (690, 16.0)], y=[(595, 0.0), (10, 1.4)], legend=(380, 230, 716, 380),
                               series={'black_filled': 'New MEXICO pressure config 0 (rough)',
                                       'red_filled': 'New MEXICO pressure config 3 (clean outboard)'}),
    'm16063_fig4_0.png': dict(x=[(76, 4.0), (509, 10.0)], y=[(595, 0.1), (13, 0.5)], legend=(195, 350, 470, 580),
                              series={'black_filled': 'rough (pressure)', 'red_filled': 'clean (pressure)',
                                      'black_open': 'rough (generator)', 'red_open': 'clean (generator)'}),
}


def lin(pairs):
    (p0, v0), (p1, v1) = pairs
    k = (v1 - v0) / (p1 - p0)
    return lambda p: v0 + k * (p - p0)


def masks(img):
    b, g, r = [img[:, :, i].astype(int) for i in range(3)]
    red = (r > 170) & (g < 110) & (b < 110)
    black = (r < 90) & (g < 90) & (b < 90)
    return {'red': red, 'black': black}


def filled_centres(m, dmin=2.6, amin=6):
    d = ndimage.distance_transform_edt(m)
    core = d >= dmin
    lab, n = ndimage.label(core)
    out = []
    for i in range(1, n + 1):
        ys, xs = np.nonzero(lab == i)
        if len(xs) >= amin:
            out.append((xs.mean(), ys.mean(), len(xs)))
    return out


def open_centres(m, amin=6, amax=120):
    filled = ndimage.binary_fill_holes(m)
    holes = filled & ~m
    lab, n = ndimage.label(holes)
    out = []
    for i in range(1, n + 1):
        ys, xs = np.nonzero(lab == i)
        if amin <= len(xs) <= amax and (xs.max() - xs.min()) < 14 and (ys.max() - ys.min()) < 14:
            out.append((xs.mean(), ys.mean(), len(xs)))
    return out


def run(name):
    c = CAL[name]
    img = cv2.imread(os.path.join(HERE, name))
    H, W = img.shape[:2]
    fx, fy = lin(c['x']), lin(c['y'])
    M = masks(img)
    rows = []
    for key, label in c['series'].items():
        col, kind = key.split('_')
        pts = filled_centres(M[col]) if kind == 'filled' else open_centres(M[col])
        for x, y, a in pts:
            # plot area only: right of the y axis, above the x axis, outside the legend
            if x < c['x'][0][0] + 3 or y > c['y'][0][0] - 3:
                continue
            L = c['legend']
            if L and L[0] <= x <= L[2] and L[1] <= y <= L[3]:
                continue
            # an open marker's ring also produces a 'filled' core if it is thick: filled squares of the
            # same colour sit on solid lines; we keep both kinds separate by construction
            rows.append(dict(series=label, x_value=fx(x), y_value=fy(y), px=x, py=y, area=a))
    df = pd.DataFrame(rows).sort_values(['series', 'x_value'])
    df.to_csv(os.path.join(HERE, name.replace('.png', '_points.csv')), index=False)
    ov = img.copy()
    colours = {}
    for i, (s, g) in enumerate(df.groupby('series')):
        cc = [(255, 0, 0), (0, 160, 0), (255, 0, 255), (0, 200, 255)][i % 4]
        for _, r in g.iterrows():
            cv2.drawMarker(ov, (int(round(r.px)), int(round(r.py))), cc, cv2.MARKER_CROSS, 14, 1)
        cv2.putText(ov, s, (60, 20 + 16 * i), cv2.FONT_HERSHEY_SIMPLEX, 0.4, cc, 1)
    # calibration check: draw the calibrated gridlines
    for p, v in c['x']:
        cv2.line(ov, (p, 0), (p, H - 1), (0, 200, 0), 1)
    for p, v in c['y']:
        cv2.line(ov, (0, p), (W - 1, p), (0, 200, 0), 1)
    big = cv2.resize(ov, None, fx=1.6, fy=1.6, interpolation=cv2.INTER_NEAREST)
    cv2.imwrite(os.path.join(HERE, name.replace('.png', '_overlay.png')), big)
    print(name)
    print(df[['series', 'x_value', 'y_value', 'area']].round(4).to_string(index=False))


if __name__ == '__main__':
    for n in (sys.argv[1:] or CAL):
        run(n)

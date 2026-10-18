"""Vector extraction of the data markers and curve fits from Bahaj et al. 2007 (Renewable Energy 32:407-426),
Figs. 7-11 (pages 11-14 of the publisher PDF, sources/bahaj2007_power_thrust_mct.pdf).

The figures in the PDF are vector graphics, so the marker positions are read exactly from the drawing
commands with PyMuPDF (no raster digitising). Each marker is reduced to its centre:
  - '+'  two short crossing lines (horizontal + vertical) sharing a midpoint
  - 'x'  two short crossing diagonal lines
  - '*'  asterisk: horizontal + vertical + two diagonals (or any >= 3 crossing lines)
  - 'diamond' a 'qu' (quad) item
  - 'square'  a small 're' item
  - 'circle'  four Bezier 'c' items forming a closed circle
Long connected polylines (solid or dashed) are the authors' curve fits; they are kept to verify the
"scale adjusted" offsets of the upper panels of Figs. 8 and 9.

Axes are calibrated from the gridlines of each panel (the thin-line drawing with 13-14 lines: the
vertical gridlines are the labelled TSR ticks, the horizontal ones the labelled C ticks).

Run: python3 vector_extract.py   -> writes extraction/raw_markers.json
"""
import json
import math
import os
from collections import defaultdict

import fitz

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, '..', 'sources', 'bahaj2007_power_thrust_mct.pdf')

# panel definitions: page index (0-based), approximate frame rect (x0,y0,x1,y1) in PDF points,
# gridline values along x and y (first to last), legend rectangles to exclude.
PANELS = {
    # Fig. 7, page 417 (pdf page 11)
    'fig7a_cp': dict(page=10, rect=(114.0, 70.9, 383.4, 313.7), xv=[2, 4, 6, 8, 10, 12],
                     yv=[0, .1, .2, .3, .4, .5, .6], legends=[(322.2, 75.3, 379.1, 123.5)]),
    'fig7b_ct': dict(page=10, rect=(113.9, 342.0, 383.6, 584.9), xv=[2, 4, 6, 8, 10, 12],
                     yv=[0, .2, .4, .6, .8, 1.0, 1.2], legends=[(322.3, 532.2, 379.2, 580.5)]),
    # Fig. 8, page 418 (pdf page 12): upper panels are "scale adjusted" (raw frame scale given here)
    'fig8_cp_upper': dict(page=11, rect=(105.1, 71.0, 224.2, 222.5), xv=[4, 5, 6, 7, 8, 9, 10],
                          yv=[.2, .3, .4, .5, .6, .7, .8], legends=[]),
    'fig8_cp_lower': dict(page=11, rect=(105.1, 226.3, 224.2, 377.8), xv=[4, 5, 6, 7, 8, 9, 10],
                          yv=[0, .1, .2, .3, .4, .5, .6], legends=[]),
    'fig8_ct_upper': dict(page=11, rect=(258.2, 71.0, 377.4, 222.5), xv=[4, 5, 6, 7, 8, 9, 10],
                          yv=[.4, .6, .8, 1.0, 1.2, 1.4, 1.6], legends=[]),
    'fig8_ct_lower': dict(page=11, rect=(258.2, 226.3, 377.4, 377.8), xv=[4, 5, 6, 7, 8, 9, 10],
                          yv=[0, .2, .4, .6, .8, 1.0, 1.2], legends=[]),
    # Fig. 9, page 419 (pdf page 13)
    'fig9_cp_upper': dict(page=12, rect=(114.1, 70.9, 231.4, 220.1), xv=[3, 4, 5, 6, 7, 8, 9],
                          yv=[.1, .2, .3, .4, .5, .6, .7], legends=[]),
    'fig9_cp_lower': dict(page=12, rect=(114.1, 223.9, 231.4, 373.1), xv=[3, 4, 5, 6, 7, 8, 9],
                          yv=[0, .1, .2, .3, .4, .5, .6], legends=[]),
    'fig9_ct_upper': dict(page=12, rect=(265.0, 70.9, 382.3, 220.1), xv=[3, 4, 5, 6, 7, 8, 9],
                          yv=[0, .2, .4, .6, .8, 1.0, 1.2], legends=[]),
    'fig9_ct_lower': dict(page=12, rect=(265.0, 223.9, 382.3, 373.1), xv=[3, 4, 5, 6, 7, 8, 9],
                          yv=[0, .2, .4, .6, .8, 1.0, 1.2], legends=[]),
    # Fig. 10 (20 deg hub pitch, yaw) and Fig. 11 (25 deg), page 420 (pdf page 14)
    'fig10a_cp': dict(page=13, rect=(68.7, 71.0, 223.7, 233.8), xv=[2, 4, 6, 8, 10, 12],
                      yv=[0, .1, .2, .3, .4, .5, .6], legends=[(167.9, 75.4, 219.3, 114.5)]),
    'fig10b_ct': dict(page=13, rect=(259.8, 71.0, 414.8, 233.8), xv=[2, 4, 6, 8, 10, 12],
                      yv=[0, .2, .4, .6, .8, 1.0, 1.2], legends=[(359.0, 190.3, 410.4, 229.4)]),
    'fig11a_cp': dict(page=13, rect=(68.6, 311.3, 223.7, 474.3), xv=[2, 4, 6, 8, 10, 12],
                      yv=[0, .1, .2, .3, .4, .5, .6], legends=[(173.2, 315.7, 219.3, 345.6)]),
    'fig11b_ct': dict(page=13, rect=(259.9, 311.3, 415.0, 474.3), xv=[2, 4, 6, 8, 10, 12],
                      yv=[0, .2, .4, .6, .8, 1.0, 1.2], legends=[(364.4, 315.7, 410.6, 345.6)]),
}


def inside(pt, rect, pad=0.0):
    x, y = pt
    return rect[0] - pad <= x <= rect[2] + pad and rect[1] - pad <= y <= rect[3] + pad


def gridlines(drawings, rect):
    """Return sorted x of vertical and y of horizontal full-length lines spanning the panel."""
    w, h = rect[2] - rect[0], rect[3] - rect[1]
    xs, ys = set(), set()
    for dr in drawings:
        for it in dr['items']:
            if it[0] != 'l':
                continue
            p, q = it[1], it[2]
            if not (inside((p.x, p.y), rect, 1.5) and inside((q.x, q.y), rect, 1.5)):
                continue
            if abs(p.x - q.x) < 0.05 and abs(p.y - q.y) > 0.95 * h:
                xs.add(round(p.x, 2))
            if abs(p.y - q.y) < 0.05 and abs(p.x - q.x) > 0.95 * w:
                ys.add(round(p.y, 2))
    return sorted(xs), sorted(ys)


def calib(pix, vals):
    """Linear least-squares map pixel -> value from matched gridline lists."""
    import numpy as np
    pix = np.asarray(pix, float)
    vals = np.asarray(vals, float)
    a, b = np.polyfit(pix, vals, 1)
    resid = np.max(np.abs(a * pix + b - vals))
    return a, b, resid


def markers_and_curves(drawings, rect, legends):
    segs = []      # short straight segments (candidate marker strokes)
    markers = []   # (kind, x, y)
    curves = []    # list of polylines (list of points) per drawing
    for di, dr in enumerate(drawings):
        dashed = bool(dr.get('dashes')) and dr.get('dashes') not in ('[] 0', '[] 0.0')
        items = dr['items']
        circ_pts = []
        poly = []
        for it in items:
            k = it[0]
            if k == 'qu':
                q = it[1]
                pts = [q.ul, q.ur, q.ll, q.lr]
                cx = sum(p.x for p in pts) / 4
                cy = sum(p.y for p in pts) / 4
                wq = max(p.x for p in pts) - min(p.x for p in pts)
                if wq < 8 and inside((cx, cy), rect, 3):
                    markers.append(('diamond', cx, cy, di))
            elif k == 're':
                r = it[1]
                if r.width < 8 and r.height < 8 and inside(((r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2), rect, 3):
                    markers.append(('square', (r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2, di))
            elif k == 'c':
                circ_pts.extend([it[1], it[4]])
                if len(circ_pts) == 8:
                    xs_ = [p.x for p in circ_pts]
                    ys_ = [p.y for p in circ_pts]
                    cx, cy = (max(xs_) + min(xs_)) / 2, (max(ys_) + min(ys_)) / 2
                    if max(xs_) - min(xs_) < 8 and inside((cx, cy), rect, 3):
                        markers.append(('circle', cx, cy, di))
                    circ_pts = []
            elif k == 'l':
                p, q = it[1], it[2]
                L = math.hypot(q.x - p.x, q.y - p.y)
                # curve polylines: consecutive connected segments
                if poly and abs(poly[-1][0] - p.x) < 1e-3 and abs(poly[-1][1] - p.y) < 1e-3:
                    poly.append((q.x, q.y))
                else:
                    if len(poly) > 15:
                        curves.append(dict(pts=poly, dashed=dashed, drawing=di))
                    poly = [(p.x, p.y), (q.x, q.y)]
                if L < 7.5:
                    segs.append(((p.x + q.x) / 2, (p.y + q.y) / 2, math.degrees(math.atan2(q.y - p.y, q.x - p.x)) % 180, L, di))
        if len(poly) > 15:
            curves.append(dict(pts=poly, dashed=dashed, drawing=di))
        # a marker drawn as 4 separate short lines forming a closed rhombus/square (seen once: Fig. 8a diamond
        # at TSR 7.09). Its strokes are removed from the crossing-stroke candidates below.
        ls = [it for it in items if it[0] == 'l']
        for a in range(len(ls) - 3):
            quad = ls[a:a + 4]
            P = [(q[1].x, q[1].y, q[2].x, q[2].y) for q in quad]
            xs_ = [v for p in P for v in (p[0], p[2])]
            ys_ = [v for p in P for v in (p[1], p[3])]
            if max(xs_) - min(xs_) > 4 or max(ys_) - min(ys_) > 4 or max(xs_) - min(xs_) < 1:
                continue
            ends = {}
            for p in P:
                for e in ((round(p[0], 2), round(p[1], 2)), (round(p[2], 2), round(p[3], 2))):
                    ends[e] = ends.get(e, 0) + 1
            if len(ends) == 4 and all(v == 2 for v in ends.values()):
                diag = all(abs(p[0] - p[2]) > 0.3 and abs(p[1] - p[3]) > 0.3 for p in P)
                cx, cy = (max(xs_) + min(xs_)) / 2, (max(ys_) + min(ys_)) / 2
                if inside((cx, cy), rect, 3):
                    markers.append(('diamond' if diag else 'square', cx, cy, di))
                    for p in P:
                        segs[:] = [s for s in segs if not (abs(s[0] - (p[0] + p[2]) / 2) < 1e-3
                                                             and abs(s[1] - (p[1] + p[3]) / 2) < 1e-3)]
    # group short segments by common midpoint -> crossing markers
    used = [False] * len(segs)
    for i, s in enumerate(segs):
        if used[i]:
            continue
        grp = [i]
        for j in range(i + 1, len(segs)):
            if not used[j] and abs(segs[j][0] - s[0]) < 0.35 and abs(segs[j][1] - s[1]) < 0.35 and segs[j][4] == s[4]:
                grp.append(j)
        if len(grp) < 2:
            continue
        # orientation classes: 0 = horizontal, 90 = vertical, 45/135 = diagonals
        angset = sorted({int(round(segs[g][2] / 45.0)) * 45 % 180 for g in grp})
        if len(angset) < 2:
            continue  # e.g. a tick mark drawn twice: not a marker (markers are crossing strokes)
        cx = sum(segs[g][0] for g in grp) / len(grp)
        cy = sum(segs[g][1] for g in grp) / len(grp)
        if not inside((cx, cy), rect, 3):
            continue
        for g in grp:
            used[g] = True
        if angset == [0, 45, 90, 135]:
            kind = '*'
        elif angset == [0, 90]:
            kind = '+'
        elif angset == [45, 135]:
            kind = 'x'
        else:
            kind = '?' + str(angset)
        # two coincident markers of the same kind (e.g. a repeated point) give twice the strokes
        nrep = max(1, round(len(grp) / len(angset)) if kind != '*' else round(len(grp) / 4))
        for _ in range(nrep):
            markers.append((kind, cx, cy, s[4]))
    # drop legend samples
    markers = [m for m in markers if not any(inside((m[1], m[2]), L, 0.5) for L in legends)]
    curves = [c for c in curves if sum(inside(p, rect, 1) for p in c['pts']) > 0.8 * len(c['pts'])]
    return markers, curves


def main():
    doc = fitz.open(PDF)
    out = {}
    cache = {}
    for name, P in PANELS.items():
        if P['page'] not in cache:
            cache[P['page']] = doc[P['page']].get_drawings()
        drs = cache[P['page']]
        gx, gy = gridlines(drs, P['rect'])
        assert len(gx) == len(P['xv']), (name, gx)
        assert len(gy) == len(P['yv']), (name, gy)
        ax, bx, rx = calib(gx, P['xv'])
        ay, by, ry = calib(gy, sorted(P['yv'], reverse=True))
        mk, cv = markers_and_curves(drs, P['rect'], P['legends'])
        out[name] = dict(
            page=P['page'] + 1, gridx=gx, gridy=gy, calib=dict(ax=ax, bx=bx, ay=ay, by=by, resid_x=rx, resid_y=ry),
            markers=[dict(kind=k, px=x, py=y, drawing=d, x=ax * x + bx, y=ay * y + by) for k, x, y, d in mk],
            curves=[dict(dashed=c['dashed'], drawing=c['drawing'],
                         x=[ax * p[0] + bx for p in c['pts']], y=[ay * p[1] + by for p in c['pts']],
                         px=[p[0] for p in c['pts']], py=[p[1] for p in c['pts']]) for c in cv])
        kinds = defaultdict(int)
        for m in mk:
            kinds[(m[0], m[3])] += 1
        print(f"{name}: grid {len(gx)}x{len(gy)} resid {rx:.2e},{ry:.2e}; markers {dict(kinds)}; "
              f"curves {[(c['drawing'], c['dashed'], len(c['pts'])) for c in cv]}")
    with open(os.path.join(HERE, 'raw_markers.json'), 'w') as f:
        json.dump(out, f, indent=1)


if __name__ == '__main__':
    main()

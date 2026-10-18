"""Generic per-figure digitiser for the NACA R-640 performance charts (CT, CP, eta vs V/nD)."""
import json, numpy as np
from digilib import *
from calib import grids, grids2
from scipy import ndimage as ndi

COLS = [(255, 0, 0), (0, 150, 0), (0, 0, 255), (200, 0, 200), (255, 128, 0), (0, 170, 170), (120, 60, 0)]


def trace_curve(d, seed, heading, box, xmin, xmax, ymax, max_turn=25, turn_pen=0.012, **kw):
    def stop(x, y):
        return x < xmin or x > xmax or y > ymax
    f = trace(d, seed, heading, box, stop=stop, max_turn=max_turn, turn_pen=turn_pen, **kw)
    b = trace(d, seed, heading + 180, box, stop=stop, max_turn=max_turn, turn_pen=turn_pen, **kw)
    p = np.vstack([b[::-1], f[1:]])
    p[:, 0] = np.clip(p[:, 0], xmin, xmax); p[:, 1] = np.minimum(p[:, 1], ymax)
    return p


def dash_mask(d, p, thr=0.12):
    """True where the traced line is dashed: gaps (low darkness) within +-12 px along the path."""
    seg = np.r_[0, np.cumsum(np.hypot(*np.diff(p, axis=0).T))]
    s = np.arange(0, seg[-1], 1.0)
    xs = np.interp(s, seg, p[:, 0]); ys = np.interp(s, seg, p[:, 1])
    # max darkness across +-2 px normal to the path (approximate with a 5x5 max filter)
    dm = ndi.maximum_filter(d, size=5)
    v = ndi.map_coordinates(dm, [ys, xs], order=1)
    gap = v < thr
    gapn = np.convolve(gap.astype(float), np.ones(25), 'same') > 0
    return np.interp(seg, s, gapn.astype(float)) > 0.5


def digitize(name, page, gbox, xcal, ycal, seeds, tbox, ylim_px, xlim_px, extra_overlay_box=None,
             vthr=200, ycal_from='h', guides=None, tol=7.0, norm=False, minlen=120, tkw=None, dkw=None, extend=True):
    """xcal=(value_of_first_vline, step) ycal=(value_of_last_hline(bottom), step) ; seeds {label: (x,y,heading,max_turn)}"""
    g, vl, hl = grids2(page, gbox)
    # calibration: vertical lines -> x, horizontal lines -> y (values from index, assumes no line missed;
    # missed lines are detected by spacing and re-indexed)
    def index(lines):
        lines = np.asarray(lines); dd = np.abs(np.diff(lines))
        sp = np.median(dd[dd < 1.4 * np.percentile(dd, 10)])
        idx = np.round(np.abs(lines - lines[0]) / sp).astype(int)
        return idx, sp
    ix, spx = index(vl); iy, spy = index(hl[::-1])
    # optional explicit anchors: xcal/ycal = (value, step, anchor_px) -> index relative to the anchor line
    if len(xcal) == 3:
        ix = np.round((np.asarray(vl) - xcal[2]) / spx).astype(int)
    if len(ycal) == 3:
        iy = np.round((ycal[2] - np.asarray(hl[::-1])) / spy).astype(int)
    def robust(pairs, step):
        A = Axis(pairs)
        keep = [pr for pr, r in zip(pairs, A.resid) if abs(r) < 0.06 * abs(step)]
        return Axis(keep), len(pairs) - len(keep)
    X, dropx = robust([(p, xcal[0] + xcal[1] * i) for p, i in zip(vl, ix)], xcal[1])
    Y, dropy = robust([(p, ycal[0] + ycal[1] * i) for p, i in zip(hl[::-1], iy)], ycal[1])
    if dropx or dropy:
        print('  calibration: dropped %d v / %d h gridlines as outliers' % (dropx, dropy))
    d = darkness_map(g, vl, hl, **(dkw or {}))
    res = {}
    for lab, (sx, sy, hd, mt) in seeds.items():
        p = trace_curve(d, (sx, sy), hd, tbox, xlim_px[0], xlim_px[1], ylim_px, max_turn=mt, **(tkw or {}))
        if guides and lab in guides:
            # hand waypoints (page px, picked on zoomed coordinate crops) replace the tracker where the
            # tracker is ambiguous (curve crossings / tangencies); the line itself is still located
            # automatically within +-tol px of the waypoint polyline
            wp = np.asarray(guides[lab], float)
            gp, ok, _ = guided_trace(d, wp, tol=tol)
            lo, hi = wp[:, 0].min(), wp[:, 0].max()
            keep = (p[:, 0] < lo) | (p[:, 0] > hi)
            p = np.vstack([p[keep & (p[:, 0] < lo)], gp, p[keep & (p[:, 0] > hi)]])
            p = p[np.argsort(p[:, 0])]
        if extend:
            p = extend_left(d, p, xlim_px[0], ylim_px, tbox)
        dm = dash_mask(d, p)
        res[lab] = (p, dm)
    ob = extra_overlay_box or tbox
    curves = []
    for k, (lab, (p, dm)) in enumerate(res.items()):
        curves.append((p, COLS[k % len(COLS)]))
    overlay(g, ob, curves, 'ov_%s.png' % name, marks=[(p[dm][i, 0], p[dm][i, 1], (255, 200, 0)) for (p, dm) in res.values() for i in range(0, dm.sum(), 6)])
    out = {}
    for lab, (p, dm) in res.items():
        out[lab] = {'x': X(p[:, 0]).tolist(), 'y': Y(p[:, 1]).tolist(), 'dashed': dm.tolist(), 'px': p.tolist()}
    cal = {'vlines': vl, 'hlines': hl, 'x_resid_max': float(np.abs(X.resid).max()), 'y_resid_max': float(np.abs(Y.resid).max()),
           'x_px_per_unit': 1 / X.k, 'y_px_per_unit': 1 / Y.k, 'nv': len(vl), 'nh': len(hl)}
    json.dump({'name': name, 'page': page, 'calibration': cal, 'curves': out}, open('dig_%s.json' % name, 'w'))
    print(name, 'cal resid x %.4f y %.5f, nv %d nh %d' % (cal['x_resid_max'], cal['y_resid_max'], len(vl), len(hl)))
    for lab, c in out.items():
        x = np.array(c['x']); y = np.array(c['y']); dm = np.array(c['dashed'])
        pass
        print('  %-4s n=%3d  x %.3f..%.3f  y %.4f..%.4f  dashed up to x=%s' % (lab, len(x), x.min(), x.max(), y.min(), y.max(),
              ('%.2f' % x[dm].max()) if dm.any() else '-'))
    return out


def auto_seeds(page, gbox, y, x0, x1, thr=150, vthr=200, axis='h', norm=False):
    """Curve crossings of row y (axis 'h') or column (axis 'v', y = column x, x0..x1 = rows)."""
    g, vl, hl = grids(page, gbox, thr=vthr, norm=norm)
    d = darkness_map(g, vl, hl, atten=0.0)
    if axis == 'h':
        prof = d[int(y) - 1:int(y) + 2, x0:x1].mean(0)
    else:
        prof = d[x0:x1, int(y) - 1:int(y) + 2].mean(1)
    on = prof > 0.3
    runs, i = [], 0
    while i < len(on):
        if on[i]:
            j = i
            while j < len(on) and on[j]: j += 1
            runs.append((x0 + (i + j - 1) / 2.0, j - i)); i = j
        else:
            i += 1
    return runs


def row_peaks(page, gbox, y, x0, x1, vthr=238, norm=True, thr=0.12, atten=0.6):
    g, vl, hl = grids2(page, gbox)
    d = darkness_map(g, vl, hl, atten=atten)
    prof = d[int(y) - 1:int(y) + 2, x0:x1].mean(0)
    return [i + x0 for i in range(1, len(prof) - 1) if prof[i] >= prof[i - 1] and prof[i] >= prof[i + 1] and prof[i] > thr], d


def probe(page, gbox, y, x0, x1, heading, tbox, xlim, ylim, **kw):
    pk, d = row_peaks(page, gbox, y, x0, x1)
    for x in pk:
        p = trace_curve(d, (x, y), heading, tbox, xlim[0], xlim[1], ylim, max_turn=12, min_score=0.04, max_gap=16)
        print(x, len(p), np.round(p[0]).tolist(), np.round(p[-1]).tolist())


def extend_left(d, p, xmin, ymax, box, gap=40, turn=8):
    """Continue a curve leftwards from its min-x end (dashed low-J part the tracker lost): heading from the
    last 25 px, small turn limit, long gap tolerance."""
    p = p[np.argsort(p[:, 0])]
    if p[0, 0] <= xmin + 3:
        return p
    a = p[0]; b = p[min(len(p) - 1, 8)]
    hd = np.degrees(np.arctan2(a[1] - b[1], a[0] - b[0]))
    stop = lambda x, y: x < xmin or y > ymax
    e = trace(d, tuple(a), hd, box, stop=stop, max_turn=turn, turn_pen=0.012, min_score=0.03, max_gap=gap)
    e[:, 0] = np.clip(e[:, 0], xmin, None)
    return np.vstack([e[::-1][:-1], p])

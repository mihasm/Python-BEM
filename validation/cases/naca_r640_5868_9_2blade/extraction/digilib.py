"""Small curve-tracing digitiser for scanned NACA charts (raster).

Workflow per figure: (1) find gridline pixel positions (long dark rows/cols),
(2) calibrate pixel->data from labelled gridlines, (3) trace each curve from a
seed point by following the dark ridge (gridlines suppressed), (4) overlay.
"""
import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage as ndi


def load_gray(path):
    return np.array(Image.open(path).convert('L')).astype(float)


def gridlines(gray, box, axis, thr=180, frac=0.35, merge=4):
    """Return centres of gridlines inside box=(x0,y0,x1,y1). axis='v' or 'h'."""
    x0, y0, x1, y1 = box
    reg = gray[y0:y1, x0:x1] < thr
    if axis == 'v':
        prof = reg.sum(0) / reg.shape[0]; off = x0
    else:
        prof = reg.sum(1) / reg.shape[1]; off = y0
    idx = np.where(prof > frac)[0]
    groups = []
    for i in idx:
        if groups and i - groups[-1][-1] <= merge:
            groups[-1].append(i)
        else:
            groups.append([i])
    out = []
    for g in groups:
        w = prof[g]
        out.append(off + float(np.sum(np.array(g) * w) / np.sum(w)))
    return out


def darkness_map(gray, vlines=(), hlines=(), half=2, blur=1.0, atten=0.35):
    """Darkness 0..1 with gridlines blanked (so the tracker is not pulled along them)."""
    d = np.clip((235.0 - gray) / 170.0, 0, 1)
    d = ndi.gaussian_filter(d, blur)
    for x in vlines:
        xi = int(round(x)); d[:, max(0, xi - half):xi + half + 1] *= atten
    for y in hlines:
        yi = int(round(y)); d[max(0, yi - half):yi + half + 1, :] *= atten
    return d


def _sample(d, xs, ys):
    return ndi.map_coordinates(d, [ys, xs], order=1, mode='constant', cval=0.0)


def trace(d, seed, angle_deg, box, step=3.0, look=20.0, max_turn=30, turn_pen=0.003,
          min_score=0.07, max_gap=8, max_steps=4000, stop=None, recentre=3):
    """Follow a dark line from seed (x,y) heading angle_deg (image coords, 0=+x, 90=+y down).
    stop(x,y)->bool optional. Returns array of (x,y) points."""
    x, y = map(float, seed)
    th = np.deg2rad(angle_deg)
    pts = [(x, y)]
    gap = 0
    x0, y0, x1, y1 = box
    t = np.arange(2.0, look + 0.01, 1.0)
    for _ in range(max_steps):
        best = None
        for dd in np.arange(-max_turn, max_turn + 0.1, 2.0):
            a = th + np.deg2rad(dd)
            sc = _sample(d, x + t * np.cos(a), y + t * np.sin(a)).mean() - turn_pen * abs(dd)
            if best is None or sc > best[0]:
                best = (sc, a)
        sc, a = best
        if sc < min_score:
            gap += 1
            if gap > max_gap:
                break
            a = th  # coast straight through a gap (gridline crossing, dash gap)
        else:
            gap = 0
        nx, ny = x + step * np.cos(a), y + step * np.sin(a)
        if recentre and gap == 0:
            # move to the darkness centroid across the line
            nrm = a + np.pi / 2
            s = np.arange(-recentre, recentre + 0.01, 0.5)
            w = _sample(d, nx + s * np.cos(nrm), ny + s * np.sin(nrm))
            if w.sum() > 0:
                off = float((w * s).sum() / w.sum())
                nx += off * np.cos(nrm); ny += off * np.sin(nrm)
        # update heading smoothly
        h = np.arctan2(ny - y, nx - x)
        th = np.arctan2(0.6 * np.sin(th) + 0.4 * np.sin(h), 0.6 * np.cos(th) + 0.4 * np.cos(h)) if gap == 0 else th
        x, y = nx, ny
        if not (x0 <= x <= x1 and y0 <= y <= y1):
            break
        if stop is not None and stop(x, y):
            pts.append((x, y)); break
        pts.append((x, y))
    # drop trailing coasted points
    if gap:
        pts = pts[:-gap] if gap < len(pts) else pts
    return np.array(pts)


class Axis:
    """Linear pixel->value map from two (pixel, value) pairs, optionally least squares over many."""
    def __init__(self, pairs):
        p = np.array(pairs, float)
        A = np.vstack([p[:, 0], np.ones(len(p))]).T
        self.k, self.c = np.linalg.lstsq(A, p[:, 1], rcond=None)[0]
        self.resid = p[:, 1] - (self.k * p[:, 0] + self.c)

    def __call__(self, px):
        return self.k * np.asarray(px) + self.c

    def inv(self, v):
        return (np.asarray(v) - self.c) / self.k


def overlay(gray, box, curves, path, marks=None, scale=1.0):
    """curves: list of (pts_px Nx2, colour). marks: list of (x,y,colour) points drawn as circles."""
    x0, y0, x1, y1 = box
    im = Image.fromarray(gray[y0:y1, x0:x1].astype(np.uint8)).convert('RGB')
    dr = ImageDraw.Draw(im)
    for pts, col in curves:
        q = [(px - x0, py - y0) for px, py in pts]
        if len(q) > 1:
            dr.line(q, fill=col, width=1)
    if marks:
        for mx, my, col in marks:
            dr.ellipse([mx - x0 - 4, my - y0 - 4, mx - x0 + 4, my - y0 + 4], outline=col, width=2)
    if scale != 1.0:
        im = im.resize((int(im.width * scale), int(im.height * scale)))
    im.save(path)


def coord_crop(gray, box, path, tick=50, scale=1.0):
    """Crop with a light pixel-coordinate ruler, for picking seeds by eye."""
    x0, y0, x1, y1 = box
    im = Image.fromarray(gray[y0:y1, x0:x1].astype(np.uint8)).convert('RGB')
    dr = ImageDraw.Draw(im)
    for X in range((x0 // tick + 1) * tick, x1, tick):
        dr.line([(X - x0, 0), (X - x0, 8)], fill=(255, 0, 0)); dr.text((X - x0 + 2, 8), str(X), fill=(255, 0, 0))
    for Y in range((y0 // tick + 1) * tick, y1, tick):
        dr.line([(0, Y - y0), (8, Y - y0)], fill=(0, 0, 255)); dr.text((10, Y - y0 - 5), str(Y), fill=(0, 0, 255))
    if scale != 1.0:
        im = im.resize((int(im.width * scale), int(im.height * scale)))
    im.save(path)


def trace_both(d, seed, angle_deg, box, xmin=None, xmax=None, **kw):
    """Trace forwards (angle) and backwards (angle+180) from an interior seed; x-limits stop."""
    def stop(x, y):
        return (xmin is not None and x < xmin) or (xmax is not None and x > xmax)
    f = trace(d, seed, angle_deg, box, stop=stop, **kw)
    b = trace(d, seed, angle_deg + 180, box, stop=stop, **kw)
    return np.vstack([b[::-1], f[1:]])


def line_maps(gray, thr=170, minlen=250):
    """Binary masks of long horizontal / vertical lines (morphological opening)."""
    import cv2
    bw = (gray < thr).astype(np.uint8)
    hk = cv2.getStructuringElement(cv2.MORPH_RECT, (minlen, 1))
    vk = cv2.getStructuringElement(cv2.MORPH_RECT, (1, minlen))
    return cv2.morphologyEx(bw, cv2.MORPH_OPEN, hk), cv2.morphologyEx(bw, cv2.MORPH_OPEN, vk)


def grid_positions(mask, axis, lo, hi, span_lo, span_hi, minfrac=0.5, merge=4):
    """Centres of long lines: axis 'h' -> rows within [lo,hi] whose mask covers > minfrac of columns span."""
    if axis == 'h':
        prof = mask[lo:hi, span_lo:span_hi].mean(1)
    else:
        prof = mask[span_lo:span_hi, lo:hi].mean(0)
    idx = np.where(prof > minfrac)[0]
    groups = []
    for i in idx:
        if groups and i - groups[-1][-1] <= merge:
            groups[-1].append(i)
        else:
            groups.append([i])
    return [lo + float(np.average(g, weights=prof[g])) for g in groups]


def refine_cross(gray, xg, yg, win=14, half=45):
    """Refine an eye-picked gridline crossing: darkest column profile in a band above/below the
    crossing (avoiding the crossing itself), and likewise for the row."""
    xg, yg = int(round(xg)), int(round(yg))
    d = 255.0 - gray
    band_v = np.vstack([d[yg - half:yg - 6, xg - win:xg + win + 1], d[yg + 6:yg + half, xg - win:xg + win + 1]])
    pv = band_v.mean(0)
    band_h = np.hstack([d[yg - win:yg + win + 1, xg - half:xg - 6], d[yg - win:yg + win + 1, xg + 6:xg + half]])
    ph = band_h.mean(1)
    def centre(p, off):
        i = int(np.argmax(p)); lo, hi = max(0, i - 2), min(len(p), i + 3)
        w = p[lo:hi] - np.median(p)
        w = np.clip(w, 0, None)
        return off + (np.sum(np.arange(lo, hi) * w) / w.sum() if w.sum() > 0 else i)
    return centre(pv, xg - win), centre(ph, yg - win)


class Affine:
    """data = A @ [px, py, 1]; fitted on refined gridline crossings."""
    def __init__(self, px, data):
        px = np.asarray(px, float); data = np.asarray(data, float)
        M = np.column_stack([px, np.ones(len(px))])
        self.A, *_ = np.linalg.lstsq(M, data, rcond=None)
        self.resid = data - M @ self.A
        Mi = np.column_stack([data, np.ones(len(data))])
        self.Ai, *_ = np.linalg.lstsq(Mi, px, rcond=None)

    def __call__(self, pts):
        pts = np.atleast_2d(pts)
        return np.column_stack([pts, np.ones(len(pts))]) @ self.A

    def inv(self, dat):
        dat = np.atleast_2d(dat)
        return np.column_stack([dat, np.ones(len(dat))]) @ self.Ai


def guided_trace(d, guide, step=2.0, tol=9.0, thr=0.25, smooth=0.5):
    """Follow a dark line along a hand-given guide polyline (Nx2 page px, any orientation).
    At each arc-length station the darkness profile along the guide normal (+-tol px) is searched for
    peaks; the peak closest to the previous offset (continuity) is taken. Stations without a peak are
    left out (gaps of dashed lines / label clutter). Returns (pts Mx2, ok-mask over stations, offsets)."""
    guide = np.asarray(guide, float)
    seg = np.r_[0, np.cumsum(np.hypot(*np.diff(guide, axis=0).T))]
    s = np.arange(0, seg[-1] + 1e-9, step)
    gx = np.interp(s, seg, guide[:, 0]); gy = np.interp(s, seg, guide[:, 1])
    tx = np.gradient(gx); ty = np.gradient(gy); n = np.hypot(tx, ty); tx /= n; ty /= n
    nx, ny = -ty, tx
    o = np.arange(-tol, tol + 0.01, 0.5)
    pts, offs, ok = [], [], []
    prev = 0.0
    for i in range(len(s)):
        xs = gx[i] + o * nx[i]; ys = gy[i] + o * ny[i]
        v = ndi.map_coordinates(d, [ys, xs], order=1)
        # local maxima above thr
        pk = [j for j in range(1, len(v) - 1) if v[j] >= v[j - 1] and v[j] >= v[j + 1] and v[j] > thr]
        if not pk:
            ok.append(False); offs.append(np.nan); continue
        # sub-sample centroid of each peak
        cands = []
        for j in pk:
            lo, hi = max(0, j - 3), min(len(v), j + 4)
            w = np.clip(v[lo:hi] - thr * 0.5, 0, None)
            cands.append(float((o[lo:hi] * w).sum() / w.sum()))
        c = min(cands, key=lambda q: abs(q - prev))
        if abs(c - prev) > 6.0 and len(pts) > 3:
            ok.append(False); offs.append(np.nan); continue
        prev = smooth * prev + (1 - smooth) * c if np.isfinite(prev) else c
        prev = c
        pts.append((gx[i] + c * nx[i], gy[i] + c * ny[i])); ok.append(True); offs.append(c)
    return np.array(pts), np.array(ok), np.array(offs)

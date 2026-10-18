"""Extract the New MEXICO measured points from Mexnext-III final report (ECN-E--18-003),
Figure A.4 (p. 149: sectional Fn, Ft for cases 1.1/1.2/1.3) and Figure A.5 (p. 150: Cdax, Cp).

The figures are vector graphics: the measured points are black filled diamonds and the
repeatability band (stdev of datapoint averages) is a grey (0.74) filled polygon. We take the
exact path coordinates with PyMuPDF and calibrate the axes from the tick marks (short black
lines) matched to the numeric tick labels. Output: extraction/mexnext3_figA4_raw.csv,
extraction/mexnext3_figA5_raw.csv, and overlay PNGs for checking.
"""
import fitz, numpy as np, csv, os

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, '..', 'sources', 'ecn_e18003_mexnext3.pdf')
doc = fitz.open(PDF)


def words_num(page):
    out = []
    for w in page.get_text('words'):
        try:
            out.append((float(w[4]), (w[0] + w[2]) / 2, (w[1] + w[3]) / 2, w))
        except ValueError:
            pass
    return out


def ticks(page):
    xt, yt = [], []
    seen = set()
    for x in page.get_drawings():
        if x['type'] == 's' and x.get('color') == (0.0, 0.0, 0.0) and len(x['items']) == 1:
            it = x['items'][0]
            if it[0] != 'l':
                continue
            a, b = it[1], it[2]
            key = (round(a.x, 1), round(a.y, 1), round(b.x, 1), round(b.y, 1))
            if key in seen:
                continue
            seen.add(key)
            if abs(a.x - b.x) < 0.01 and 1 < abs(a.y - b.y) < 5:
                xt.append((a.x, max(a.y, b.y)))       # vertical tick -> x axis
            elif abs(a.y - b.y) < 0.01 and 1 < abs(a.x - b.x) < 5:
                yt.append((min(a.x, b.x), a.y))       # horizontal tick -> y axis
    return xt, yt


def calib(tick_list, labels, axis):
    """match each tick to the nearest numeric label, fit a line value = k*pix + n"""
    pairs = []
    for t in tick_list:
        best = None
        for v, cx, cy, w in labels:
            if axis == 'x':
                dist = abs(cx - t[0]) + 0.5 * abs(cy - t[1] - 5)
                if cy < t[1] or cy - t[1] > 12:
                    continue
            else:
                dist = abs(cy - t[1]) + 0.5 * abs(t[0] - cx - 6)
                if cx > t[0] or t[0] - cx > 16:
                    continue
            if best is None or dist < best[0]:
                best = (dist, v)
        if best and best[0] < 4:
            pairs.append(((t[0] if axis == 'x' else t[1]), best[1]))
    pairs = sorted(set(pairs))
    P = np.array(pairs)
    k, n = np.polyfit(P[:, 0], P[:, 1], 1)
    resid = np.max(np.abs(P[:, 0] * k + n - P[:, 1]))
    return k, n, pairs, resid


def diamonds(page):
    pts = []
    for x in page.get_drawings():
        if x['type'] == 'f' and x.get('fill') == (0.0, 0.0, 0.0) and len(x['items']) == 3:
            r = x['rect']
            pts.append(((r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2))
    return sorted(set((round(a, 3), round(b, 3)) for a, b in pts))


def bands(page):
    out = []
    for x in page.get_drawings():
        if x['type'] == 'f' and x.get('fill') and abs(x['fill'][0] - 0.74) < 0.01 and len(x['items']) >= 3:
            poly = []
            for it in x['items']:
                if it[0] == 'l':
                    poly += [(it[1].x, it[1].y), (it[2].x, it[2].y)]
            out.append((x['rect'], poly))
    return out


# ---------------- Figure A.4, p.149 ----------------
page = doc[148]
W = words_num(page)
xt, yt = ticks(page)
# panel layout: left column x<260 (Fn), right column (Ft); rows by y
rows = [(100, 290), (300, 505), (520, 720)]
cases = ['1.1', '1.2', '1.3']
cols = {'fn': (40, 260), 'ft': (262, 470)}
D = diamonds(page)
B = bands(page)
res = []
cal_used = {}
for (y0, y1), case in zip(rows, cases):
    for q, (x0, x1) in cols.items():
        lab = [w for w in W if x0 <= w[1] <= x1 and y0 <= w[2] <= y1]
        xtp = [t for t in xt if x0 <= t[0] <= x1 and y0 <= t[1] <= y1]
        ytp = [t for t in yt if x0 <= t[0] <= x1 + 20 and y0 <= t[1] <= y1]
        kx, nx, px, rx = calib(xtp, lab, 'x')
        ky, ny, py, ry = calib(ytp, lab, 'y')
        cal_used[(case, q)] = (kx, nx, ky, ny, px, py, rx, ry)
        # plot frame: points must be inside x-range of the ticks area (exclude legend diamonds at x>470)
        dd = [p for p in D if x0 + 15 <= p[0] <= x1 and y0 <= p[1] <= y1 - 10]
        bb = [b for b in B if x0 <= b[0].x0 <= x1 and y0 <= b[0].y0 <= y1]
        for (u, v) in dd:
            r_m = kx * u + nx
            val = ky * v + ny
            # band: vertical extent of the grey polygon at this x
            lo = hi = np.nan
            for rect, poly in bb:
                P = np.array(poly)
                near = P[np.abs(P[:, 0] - u) < 0.6]
                if len(near):
                    vals = ky * near[:, 1] + ny
                    lo, hi = vals.min(), vals.max()
            res.append(dict(case=case, quantity=q, r_m=r_m, value=val, band_lo=lo, band_hi=hi, px=u, py=v))
        print(case, q, 'xcal pts', len(px), 'resid', rx, 'ycal pts', len(py), 'resid', ry, 'n points', len(dd), 'bands', len(bb))

with open(os.path.join(HERE, 'mexnext3_figA4_raw.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(res[0].keys()))
    w.writeheader()
    for r in res:
        w.writerow(r)

# overlay
pix = page.get_pixmap(dpi=200)
import cv2
img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)[:, :, :3].copy()
img = cv2.cvtColor(img, cv2.COLOR_RGB2BGR)
s = 200 / 72
for r in res:
    case, q = r['case'], r['quantity']
    kx, nx, ky, ny = cal_used[(case, q)][:4]
    u = (r['r_m'] - nx) / kx
    v = (r['value'] - ny) / ky
    cv2.circle(img, (int(u * s), int(v * s)), 9, (0, 0, 255), 2)
    for b in (r['band_lo'], r['band_hi']):
        if not np.isnan(b):
            vb = (b - ny) / ky
            cv2.line(img, (int(u * s) - 8, int(vb * s)), (int(u * s) + 8, int(vb * s)), (255, 0, 0), 2)
cv2.imwrite(os.path.join(HERE, 'overlay_mexnext3_figA4.png'), img)

# ---------------- Figure A.5, p.150 ----------------
page = doc[149]
W = words_num(page)
xt, yt = ticks(page)
D = diamonds(page)
B = bands(page)
res5 = []
for q, (x0, x1) in {'cdax': (40, 260), 'cp': (262, 470)}.items():
    lab = [w for w in W if x0 <= w[1] <= x1 and 100 <= w[2] <= 400]
    ytp = [t for t in yt if x0 <= t[0] <= x1 + 20 and 100 <= t[1] <= 400]
    ky, ny, py, ry = calib(ytp, lab, 'y')
    dd = [p for p in D if x0 + 15 <= p[0] <= x1 and 100 <= p[1] <= 380]
    print('A5', q, 'ycal', py, 'resid', ry, 'pts', dd)
    # x positions: 24 m/s left, 15 middle, 10 right
    xs = sorted(set(round(p[0], 1) for p in dd))
    labels = {0: '1.3', 1: '1.2', 2: '1.1'}
    for (u, v) in dd:
        idx = xs.index(round(u, 1))
        res5.append(dict(case=labels[idx], quantity=q, value=ky * v + ny, px=u, py=v))
with open(os.path.join(HERE, 'mexnext3_figA5_raw.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(res5[0].keys()))
    w.writeheader()
    for r in res5:
        w.writerow(r)
for r in res5:
    print(r)

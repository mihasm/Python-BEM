"""Extract the Risø-A1-21 polar with leading-edge roughness (zigzag tape, run005), Re = 1.6e6, from
Fuglsang, Dahl & Antoniou 1999, Risø-R-1112(EN), Figure 6-8 (report page 38 = PDF page 39).

The figure is vector graphics. In it: filled circles = 'LER measurement', open circles = smooth
measurement, filled squares = EllipSys2D. The filled circles are filled Bezier paths ('c' items);
we take their centres. Axes are calibrated from the dotted grid lines (exact), matched to the tick
labels. The three panels (CL, CD, CM vs alpha) share alpha values point by point, so the points of
the three panels are matched by their alpha.

Output: extraction/risoe_A1_21_LER_raw.csv and overlay_risoe_fig6-8.png
"""
import os, csv
import fitz, numpy as np, cv2

HERE = os.path.dirname(os.path.abspath(__file__))
doc = fitz.open(os.path.join(HERE, '..', 'sources', 'risoe_r1112_A1_airfoils.pdf'))
page = doc[38]

nums = []
for w in page.get_text('words'):
    try:
        nums.append((float(w[4]), (w[0] + w[2]) / 2, (w[1] + w[3]) / 2))
    except ValueError:
        pass

# grid lines (dotted) and frame
vlines, hlines = set(), set()
for x in page.get_drawings():
    if x['type'] != 's':
        continue
    for it in x['items']:
        if it[0] == 'l':
            a, b = it[1], it[2]
            if abs(a.x - b.x) < 0.05 and abs(a.y - b.y) > 50:
                vlines.add(round(a.x, 2))
            if abs(a.y - b.y) < 0.05 and abs(a.x - b.x) > 50:
                hlines.add(round(a.y, 2))
        elif it[0] == 're':
            r = it[1]
            vlines.update([round(r.x0, 2), round(r.x1, 2)])
            hlines.update([round(r.y0, 2), round(r.y1, 2)])

panels = {'cl': (40, 245), 'cd': (250, 455), 'cm': (460, 665)}


def fit(pairs):
    P = np.array(pairs)
    k, n = np.polyfit(P[:, 0], P[:, 1], 1)
    return k, n, float(np.max(np.abs(P[:, 0] * k + n - P[:, 1])))


cal = {}
for q, (y0, y1) in panels.items():
    xl = [(v, cx, cy) for v, cx, cy in nums if y0 < cy < y1 and cy > y1 - 25 and cx > 220]
    yl = [(v, cx, cy) for v, cx, cy in nums if y0 < cy < y1 and cx < 230 and cy < y1 - 15]
    xp = []
    for v, cx, cy in xl:
        g = min(vlines, key=lambda u: abs(u - cx))
        if abs(g - cx) < 2:
            xp.append((g, v))
    yp = []
    for v, cx, cy in yl:
        g = min(hlines, key=lambda u: abs(u - cy))
        if abs(g - cy) < 2:
            yp.append((g, v))
    kx, nx, rx = fit(sorted(set(xp)))
    ky, ny, ry = fit(sorted(set(yp)))
    cal[q] = (kx, nx, ky, ny)
    print(q, 'x pairs', sorted(set(xp)), 'resid', rx, '| y pairs', sorted(set(yp)), 'resid', ry)

pts = {q: [] for q in panels}
for x in page.get_drawings():
    if x['type'] == 'f' and all(i[0] == 'c' for i in x['items']):
        r = x['rect']
        cx, cy = (r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2
        for q, (y0, y1) in panels.items():
            if y0 < cy < y1:
                kx, nx, ky, ny = cal[q]
                pts[q].append((kx * cx + nx, ky * cy + ny, cx, cy))

# drop the legend symbols: they lie right of the text "LER measurement" - detect by alpha > 35 or
# by being outside the plot frame; legend symbols sit at fixed screen positions in each panel
for q in pts:
    ys = sorted(set(round(p[3], 1) for p in pts[q]))
    pts[q] = [p for p in pts[q] if -10 <= p[0] <= 35]
    print(q, len(pts[q]), 'points')

# legend symbol in CL panel sits at alpha ~ 33..35 (inside the frame); remove points whose
# alpha does not occur in the other panels
def key(a):
    return round(a, 1)

alist = {q: [p[0] for p in pts[q]] for q in pts}
rows = []
for a, cl, cx, cy in sorted(pts['cl']):
    match = {}
    for q in ('cd', 'cm'):
        cand = [p for p in pts[q] if abs(p[0] - a) < 0.15]
        if len(cand) == 1:
            match[q] = cand[0][1]
    if len(match) == 2:
        rows.append(dict(alpha_deg=round(a, 3), cl=round(cl, 4), cd=round(match['cd'], 5), cm=round(match['cm'], 4)))
    else:
        print('unmatched CL point', a, cl, match)
print(len(rows), 'matched rows')
with open(os.path.join(HERE, 'risoe_A1_21_LER_raw.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

# overlay
dpi = 200
s = dpi / 72
pix = page.get_pixmap(dpi=dpi)
img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)[:, :, :3].copy()
img = cv2.cvtColor(img, cv2.COLOR_RGB2BGR)
for r in rows:
    for q in ('cl', 'cd', 'cm'):
        kx, nx, ky, ny = cal[q]
        u = (r['alpha_deg'] - nx) / kx
        v = (r[q] - ny) / ky
        cv2.circle(img, (int(u * s), int(v * s)), 7, (0, 0, 255), 1)
cv2.imwrite(os.path.join(HERE, 'overlay_risoe_fig6-8.png'), img)
for r in rows:
    print(r)

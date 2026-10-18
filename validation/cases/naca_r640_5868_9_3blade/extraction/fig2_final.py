import sys, json, numpy as np
sys.path.insert(0, '.')
from digilib import *
g = load_gray('../p6-000.png')
box = (378, 1596, 1176, 2736)
vl = gridlines(g, box, 'v'); hl = gridlines(g, box, 'h')
d = darkness_map(g, vl, hl)
B = (380, 1600, 1172, 2730)
paths = {
 'P/D35': trace_both(d, (503, 2140), -40, B, xmin=455, xmax=1162, max_turn=10),
 'b/D':   trace_both(d, (503, 2203), -57, B, xmin=402, xmax=1162, max_turn=30),
 'P/D25': trace_both(d, (503, 2318), -40, B, xmin=402, xmax=1162, max_turn=30),
 'P/D15': trace_both(d, (503, 2443), -25, B, xmin=402, xmax=1162, max_turn=30),
 'h/b':   trace_both(d, (800, 2480.5), 2, B, xmin=402, xmax=1162, max_turn=14),
}
# root part of the 35 deg curve (it crosses the dashed b/D curve of 37-3647 at r/R~0.245): separate seed
_root = trace_both(d, (425, 2235), -38, B, xmin=402, xmax=458, max_turn=12)
paths['P/D35'] = np.vstack([_root[_root[:, 0] < 455], paths['P/D35']])
# x calibration: vertical gridlines 402.8 ... 1153 = r/R 0.20 ... 0.95 (every 0.05)
vx = [v for v in vl if 400 < v < 1160]
X = Axis([(p, 0.20 + 0.05 * i) for i, p in enumerate(vx)])
# y calibration: horizontal gridlines bottom 2734.4 = 0, each line up = 0.02 h/b (0.005 b/D, 0.1 P/D)
hy = sorted([h for h in hl if h > 1600], reverse=True)
Yhb = Axis([(p, 0.02 * i) for i, p in enumerate(hy)])
print('x resid', np.abs(X.resid).max(), 'y resid', np.abs(Yhb.resid).max(), 'n', len(vx), len(hy))
scale = {'h/b': 1.0, 'b/D': 0.25, 'P/D35': 5.0, 'P/D25': 5.0, 'P/D15': 5.0}
data = {}
for k, p in paths.items():
    x = X(p[:, 0]); y = Yhb(p[:, 1]) * scale[k]
    o = np.argsort(x); x, y = x[o], y[o]
    data[k] = (x, y)
    print(k, 'r/R %.3f..%.3f' % (x.min(), x.max()))
np.save('fig2_data.npy', {k: np.vstack(v) for k, v in data.items()}, allow_pickle=True)
cols = {'P/D35': (255, 0, 0), 'b/D': (0, 160, 0), 'P/D25': (0, 0, 255), 'P/D15': (200, 0, 200), 'h/b': (255, 128, 0)}
overlay(g, (300, 1550, 1250, 2760), [(v, cols[k]) for k, v in paths.items()], 'fig2_overlay.png')
# stations
rr = np.round(np.arange(0.20, 0.951, 0.05), 3)
out = []
for r in rr:
    row = {'r_R': r}
    for k, (x, y) in data.items():
        row[k] = float(np.interp(r, x, y)) if x.min() - 1e-3 <= r <= x.max() + 1e-3 else np.nan
    for s in (15, 25, 35):
        pd = row['P/D%d' % s]
        row['beta%d' % s] = np.degrees(np.arctan(pd / (np.pi * r)))
        row['dbeta%d' % s] = row['beta%d' % s] - s  # blade angle relative to setting at 0.75R
    out.append(row)
import pandas as pd
df = pd.DataFrame(out)
pd.set_option('display.width', 200)
print(df.round(4).to_string())
df.to_csv('fig2_stations.csv', index=False)

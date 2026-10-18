import sys, json, numpy as np
sys.path.insert(0, '.')
from digilib import *
g = load_gray('../p6-000.png')
box = (378, 1596, 1176, 2736)
vl = gridlines(g, box, 'v'); hl = gridlines(g, box, 'h')
d = darkness_map(g, vl, hl)
# interior seeds picked by eye on coordinate crops (page pixel coords), heading in image degrees
seeds = {'P/D35': ((503, 2140), -40, 10), 'b/D': ((503, 2203), -57, 30), 'P/D25': ((503, 2318), -40, 30),
         'P/D15': ((503, 2443), -25, 30), 'h/b':   ((652, 2441), 15, 30)}
cols = {'P/D35': (255, 0, 0), 'b/D': (0, 160, 0), 'P/D25': (0, 0, 255), 'P/D15': (200, 0, 200), 'h/b': (255, 128, 0)}
res = {}
for k, (s, a, mt) in seeds.items():
    res[k] = trace_both(d, s, a, (380, 1600, 1172, 2730), xmin=402, xmax=1162, max_turn=mt)
    print(k, len(res[k]), np.round(res[k][0]), np.round(res[k][-1]))
overlay(g, (300, 1550, 1250, 2760), [(v, cols[k]) for k, v in res.items()], 'fig2_overlay_try.png')
np.save('fig2_paths.npy', res, allow_pickle=True)

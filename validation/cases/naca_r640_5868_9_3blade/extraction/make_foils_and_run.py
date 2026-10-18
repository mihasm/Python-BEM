"""Clark Y family by simple ordinate (y) scaling of the UIUC Clark Y (11.7 %), then XFOIL polars."""
import numpy as np, subprocess, os, sys
from concurrent.futures import ThreadPoolExecutor
# UIUC clarky.dat is in Lednicer format: title, '61. 61.', 61 upper points LE->TE, blank, 61 lower LE->TE
L = open('clarky_uiuc.dat').read().splitlines()
pts = [list(map(float, l.split())) for l in L[2:] if len(l.split()) == 2]
upper = np.array(pts[:61]); lower = np.array(pts[61:122])
xy = np.vstack([upper[::-1], lower[1:]])   # Selig order: TE -> upper -> LE -> lower -> TE
t0 = 0.0
# thickness of the original (max over x of upper - lower)
iLE = np.argmin(xy[:, 0]); up = xy[:iLE + 1][::-1]; lo = xy[iLE:]
xs = np.linspace(0, 1, 2001)
t0 = (np.interp(xs, up[:, 0], up[:, 1]) - np.interp(xs, lo[:, 0], lo[:, 1])).max()
print('original t/c', t0)
TCS = [0.08, 0.09, 0.10, 0.12, 0.14, 0.18, 0.24]
XF = '/Users/miha/Documents/dev/Python-BEM/xfoil_executables/xfoil'
RES = [5e5, 1e6, 2e6, 3e6]
def write(tc):
    n = 'clarky_t%02d' % round(tc * 100)
    with open(n + '.dat', 'w') as f:
        f.write('%s Clark Y scaled to t/c=%.2f (all y-ordinates x %.4f of UIUC clarky.dat, t/c %.4f)\n' % (n, tc, tc / t0, t0))
        for x, y in xy: f.write('%.6f %.6f\n' % (x, y * tc / t0))
    return n
def run(n, re):
    out = '%s_re%.0e.pol' % (n, re)
    if os.path.exists(out): os.remove(out)
    cmds = ['PLOP', 'G F', '', 'LOAD %s.dat' % n, 'PPAR', 'N 200', '', '', 'OPER', 'VISC %g' % re, 'VPAR', 'N 9', '',
            'ITER 200', 'PACC', out, '', 'ASEQ 0 -10 -0.5', 'INIT', 'ASEQ 0.5 20 0.5', 'PACC', '', 'QUIT']
    try:
        subprocess.run([XF], input='\n'.join(cmds) + '\n', text=True, capture_output=True, timeout=600)
    except subprocess.TimeoutExpired:
        print('timeout', n, re)
    return out
names = [write(tc) for tc in TCS]
jobs = [(n, re) for n in names for re in RES]
with ThreadPoolExecutor(8) as ex:
    for o in ex.map(lambda a: run(*a), jobs):
        print(o, sum(1 for l in open(o) if l[:9].strip().replace('-', '').replace('.', '').isdigit()) if os.path.exists(o) else 'MISSING', flush=True)

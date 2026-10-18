"""Sanity check of a marine propeller case: loads every CSV with pandas, checks units/ranges, plots
geometry, sections, polars and the experiment into extraction/check_<case>.png.  Run: python3 check_case.py"""
import json, os, glob
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
CASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
meta = json.load(open(os.path.join(CASE, 'case.json')))
g = pd.read_csv(os.path.join(CASE, meta['geometry_file']), comment='#')
R = meta['R']
assert abs(g.r_m.iloc[-1] - R) < 1e-6, 'radius must end at R'
assert g.chord_m.max() < R, 'chord in m'
assert (np.diff(g.twist_deg) < 0).all(), 'propeller blade angle must decrease outboard'
for n, a in meta['airfoils'].items():
    p = pd.read_csv(os.path.join(CASE, a['polars']), comment='#')
    assert {'re', 'alpha_deg', 'cl', 'cd'} <= set(p.columns), n
    assert (p.cd > 0).all(), n
fig, ax = plt.subplots(2, 3, figsize=(16, 9))
ax[0, 0].plot(g.r_m / R, g.chord_m / (2 * R), 'o-', label='c/D'); ax[0, 0].plot(g.r_m / R, g.thickness_ratio, 's-', label='t/c')
ax[0, 0].set_xlabel('r/R'); ax[0, 0].legend(); ax[0, 0].set_title(meta['id'] + ' geometry')
ax[0, 1].plot(g.r_m / R, g.twist_deg, 'o-'); ax[0, 1].set_xlabel('r/R'); ax[0, 1].set_ylabel('twist (beta) [deg]')
for n, a in meta['airfoils'].items():
    xy = np.loadtxt(os.path.join(CASE, a['coordinates']), skiprows=1)
    ax[0, 2].plot(xy[:, 0], xy[:, 1], lw=0.8, label=n[-5:])
ax[0, 2].set_aspect('equal'); ax[0, 2].legend(fontsize=6, ncol=3); ax[0, 2].set_title('sections')
names = list(meta['airfoils'])
for n in [names[0], names[len(names) // 2], names[-1]]:
    p = pd.read_csv(os.path.join(CASE, meta['airfoils'][n]['polars']), comment='#')
    for re, gg in p.groupby('re'):
        ax[1, 0].plot(gg.alpha_deg, gg.cl, lw=0.8); ax[1, 1].plot(gg.cd, gg.cl, lw=0.8, label=f'{n[-5:]} {re:.0e}')
ax[1, 0].set_xlabel('alpha'); ax[1, 0].set_ylabel('cl'); ax[1, 1].set_xlabel('cd'); ax[1, 1].set_xlim(0, 0.05)
ax[1, 1].legend(fontsize=5, ncol=2); ax[1, 1].set_title('XFOIL polars (root/mid/tip sections)')
for e in meta['experiments']:
    d = pd.read_csv(os.path.join(CASE, e['file']), comment='#')
    assert (d.eta[d.ct > 0] < 1).all() and (d.cp > 0).all(), e['file']
    lab = os.path.basename(e['file'])[12:-4]
    ax[1, 2].plot(d.J, d.ct, 'o-', ms=3, label='ct ' + lab); ax[1, 2].plot(d.J, d.cp, 's--', ms=3, label='cp ' + lab)
    ax[1, 2].plot(d.J, d.eta, '^:', ms=3, label='eta ' + lab)
ax[1, 2].set_ylim(-0.05, 1); ax[1, 2].set_xlabel('J'); ax[1, 2].legend(fontsize=6); ax[1, 2].grid(alpha=.3)
fig.tight_layout(); out = os.path.join(CASE, 'extraction', 'check_%s.png' % meta['id']); fig.savefig(out, dpi=90)
print('OK', meta['id'], 'stations', len(g), 'airfoils', len(names), '->', out)

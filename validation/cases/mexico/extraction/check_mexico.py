"""Sanity check of the MEXICO case files: units, ranges, plots -> extraction/check_mexico.png"""
import os, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

C = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
meta = json.load(open(os.path.join(C, 'case.json')))
rd = lambda p: pd.read_csv(os.path.join(C, p), comment='#')
g = rd('geometry.csv')
R = meta['R']
assert abs(g.r_m.iloc[-1] - R) < 1e-9 and g.r_m.is_monotonic_increasing
assert g.chord_m.max() < 0.3 and (np.diff(g.twist_deg[5:]) <= 0).all(), 'twist must decrease outboard'
assert set(g.airfoil) <= set(meta['airfoils']) | {'transition'}
for e in meta['experiments']:
    d = rd(e['file'])
    for y in e['y']:
        assert y in d.columns, (e['file'], y)
    if 'cp' in d:
        assert d.cp.max() < 16 / 27, 'Cp above Betz'
    print('ok', e['file'], len(d), 'rows')
for n, a in meta['airfoils'].items():
    if 'polars' in a:
        p = rd(a['polars'])
        print(n, 'Re', sorted(p.re.unique()), 'alpha %.1f..%.1f' % (p.alpha_deg.min(), p.alpha_deg.max()),
              'cl max %.2f cd min %.4f' % (p.cl.max(), p.cd.min()))

fig, ax = plt.subplots(3, 3, figsize=(16, 13))
a = ax[0, 0]; a.plot(g.r_m, g.chord_m, 'o-'); a.set_ylabel('chord [m]'); a.set_xlabel('r [m]')
for i, r in g.iterrows():
    a.annotate(r.airfoil[:4], (r.r_m, r.chord_m), fontsize=6)
a = ax[0, 1]; a.plot(g.r_m, g.twist_deg, 'o-'); a.set_ylabel('twist [deg] (pitch -2.3 not incl.)'); a.set_xlabel('r [m]')
a = ax[0, 2]
for n in ('du91_w2_250', 'risoe_a1_21', 'naca64_418'):
    xy = np.loadtxt(os.path.join(C, 'airfoils', n + '.dat'), skiprows=1)
    a.plot(xy[:, 0], xy[:, 1], label=n)
a.set_aspect('equal'); a.legend(fontsize=7)
a = ax[1, 0]; b = ax[1, 1]
for n, st in (('du91_w2_250', '-'), ('risoe_a1_21', 'o-'), ('naca64_418', '--')):
    p = rd(meta['airfoils'][n]['polars'])
    for re, q in p.groupby('re'):
        if re in (5e5, 7e5, 1.6e6):
            a.plot(q.alpha_deg, q.cl, st, ms=3, label='%s Re %.1e' % (n, re))
            b.plot(q.alpha_deg, q.cd, st, ms=3)
a.set_ylabel('cl'); b.set_ylabel('cd'); a.legend(fontsize=6); a.set_xlabel('alpha'); b.set_xlabel('alpha')
s = rd('experiment/newmexico2014_spanwise_axial.csv')
s6 = rd('experiment/mexico2006_spanwise_axial.csv')
a = ax[1, 2]
for (v, q), col in zip(s.groupby('v_ms'), 'bgr'):
    a.errorbar(q.r_m, q.fn_n_per_m, q.fn_n_per_m_err, fmt='o-', color=col, label='2014 Fn %.1f m/s' % v)
    a.plot(q.r_m, q.ft_n_per_m, 's--', color=col)
for (v, q), col in zip(s6.groupby('v_ms'), 'bgr'):
    a.plot(q.r_m, q.fn_n_per_m, 'x:', color=col, label='2006 Fn %.2f m/s' % v)
a.set_xlabel('r [m]'); a.set_ylabel('Fn (o) / Ft (s) [N/m], rotor axes'); a.legend(fontsize=6)
a = ax[2, 0]; b = ax[2, 1]
for f, mk in (('newmexico2014_tsr_sweep_rough', 'o'), ('newmexico2014_tsr_sweep_clean', 's'), ('mexico2006_tsr_sweep', 'x')):
    d = rd('experiment/%s.csv' % f)
    a.plot(d.tsr, d.ct, mk, label=f + ' ct (balance)')
    a.plot(d.tsr, d.ct_pressure, mk, mfc='none', label=f + ' ct_pressure')
    if 'cp' in d:
        b.plot(d.tsr, d.cp, mk, label=f + ' cp (pressure)')
        b.plot(d.tsr, d.cp_generator, mk, mfc='none', label=f + ' cp_generator')
r14 = rd('experiment/newmexico2014_axial_rotor_loads.csv')
a.plot(r14.tsr, r14.ct, 'k*', ms=12, label='Mexnext-III cases (pressure)'); b.plot(r14.tsr, r14.cp, 'k*', ms=12)
a.set_xlabel('tsr'); a.set_ylabel('Ct (Cdax)'); a.legend(fontsize=6)
b.set_xlabel('tsr'); b.set_ylabel('Cp'); b.legend(fontsize=6); b.axhline(16 / 27, color='k', lw=.5)
a = ax[2, 2]
m6 = rd('experiment/mexico2006_balance_thrust.csv')
a.plot(r14.v_ms, r14.thrust_n, 'k*', label='2014 Fax (pressure, weighting factors)')
a.plot(m6.v_ms, m6.thrust_n, 'ro', label='2006 balance - tower drag (Bechmann)')
a.set_xlabel('v [m/s]'); a.set_ylabel('thrust [N]'); a.legend(fontsize=7)
for x in ax.flat:
    x.grid(alpha=.3)
fig.tight_layout()
fig.savefig(os.path.join(C, 'extraction', 'check_mexico.png'), dpi=110)
print('wrote check_mexico.png')

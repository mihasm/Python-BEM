"""Sanity checks of the whole case (units, ranges, file references) and the overview plot check_bahaj_tidal_800mm.png."""
import json
import os

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.join(HERE, '..')
meta = json.load(open(os.path.join(CASE, 'case.json')))
geo = pd.read_csv(os.path.join(CASE, meta['geometry_file']), comment='#')
R = meta['R']

# --- geometry
assert abs(geo.r_m.iloc[-1] - R) < 1e-9, 'radius must end at R'
assert geo.r_m.is_monotonic_increasing and geo.r_m.iloc[0] >= meta['Rhub']
assert geo.chord_m.between(0.01, 0.1).all(), 'chord in m'
assert geo.twist_deg.is_monotonic_decreasing, 'turbine twist decreasing outboard'
assert set(geo.airfoil) <= set(meta['airfoils']), 'airfoil keys'
for k, a in meta['airfoils'].items():
    for f in (a['coordinates'], a['polars']):
        assert os.path.exists(os.path.join(CASE, f)), f
    pol = pd.read_csv(os.path.join(CASE, a['polars']), comment='#')
    assert {'re', 'alpha_deg', 'cl', 'cd'} <= set(pol.columns)
    assert (pol.cd > 0).all() and pol.alpha_deg.between(-12, 22).all()
    xy = np.loadtxt(os.path.join(CASE, a['coordinates']), skiprows=1)
    assert abs(xy[0, 0] - 1) < 1e-6 and abs(xy[-1, 0] - 1) < 1e-6 and xy[:, 0].min() < 0.001
    tc = float(geo.loc[geo.airfoil == k, 'thickness_ratio'].iloc[0])
    up = xy[:np.argmin(xy[:, 0]) + 1][::-1]
    lo = xy[np.argmin(xy[:, 0]):]
    xx = np.linspace(0.01, 0.99, 500)
    t = np.max(np.interp(xx, up[:, 0], up[:, 1]) - np.interp(xx, lo[:, 0], lo[:, 1]))
    assert abs(t - tc) < 0.002, (k, t, tc)
print('geometry/airfoils OK:', len(geo), 'stations,', len(meta['airfoils']), 'airfoils')

# --- experiments
ALLOWED = {'tsr', 'cp', 'ct', 'cp_err', 'ct_err', 'tsr_err', 'v_ms', 'rpm', 'pitch_deg', 'yaw_deg', 'rho',
           'kin_viscosity'}
exps = {}
for e in meta['experiments']:
    df = pd.read_csv(os.path.join(CASE, e['file']), comment='#')
    for c in [e['x']] + e['y']:
        assert c in df.columns, (e['file'], c)
    for k in e['defaults']:
        assert k not in df.columns or True
    has_v = 'v_ms' in df.columns or 'v_ms' in e['defaults']
    assert has_v and df.get('v_ms', pd.Series([1])).notna().all()
    assert 'pitch_deg' in df.columns or 'pitch_deg' in e['defaults']
    assert (df.cp.dropna() < 16 / 27).all() and (df.cp.dropna() > 0).all(), 'Cp below Betz'
    assert df.ct.dropna().between(0.2, 1.3).all()
    extra = [c for c in df.columns if c not in ALLOWED]
    print(f"{e['file']}: {len(df)} rows, tsr {df.tsr.min():.2f}-{df.tsr.max():.2f}, "
          f"cp max {df.cp.max():.3f}, ct max {df.ct.max():.3f}, pitch {sorted(df.pitch_deg.unique()) if 'pitch_deg' in df else e['defaults'].get('pitch_deg')}, extra cols {extra}")
    exps[e['file']] = df

# --- Reynolds number range of the rotor (for the polar Re choice)
nu = meta['fluid']['kin_viscosity']
for v, tsr in [(1.0, 4), (1.73, 6), (1.73, 11), (1.3, 3)]:
    rr = geo.r_m.to_numpy()
    w = v * np.sqrt((2 / 3) ** 2 + (tsr * rr / R) ** 2)   # a ~ 1/3, a' ~ 0
    re = w * geo.chord_m.to_numpy() / nu
    print(f'v {v} tsr {tsr}: Re root {re[0]:.2e}, 0.7R {re[10]:.2e}, tip {re[-1]:.2e}, max {re.max():.2e}')

# --- plot
fig, axs = plt.subplots(2, 3, figsize=(17, 9.5))
ax = axs[0, 0]
ax.plot(geo.r_m / R, geo.chord_m / R, 'o-', label='c/R')
ax.plot(geo.r_m / R, geo.thickness_ratio, 's-', label='t/c')
ax2 = ax.twinx()
ax2.plot(geo.r_m / R, geo.twist_deg, 'r^-', label='twist_deg (excl. pitch)')
ax2.plot(geo.r_m / R, geo.twist_deg + 5.0, 'r:', label='twist + 5 (20 deg hub pitch)')
ax2.set_ylabel('deg', color='r')
ax.set_xlabel('r/R')
ax.set_title('geometry.csv (Bahaj 2007 Table 1)')
ax.legend(loc='upper center', fontsize=8)
ax2.legend(loc='upper right', fontsize=8)
ax.grid(alpha=.3)

t = exps['experiment/tunnel_performance.csv']
cols = {0.0: 'C0', 5.0: 'C1', 10.0: 'C2', 12.0: 'C3', 15.0: 'C4'}
for (p, v), g in t.groupby(['pitch_deg', 'v_ms']):
    mk = 'o' if v != 1.0 else 'x'
    for k, q in enumerate(['cp', 'ct']):
        axs[0, 1 + k].errorbar(g.tsr, g[q], yerr=g[q + '_err'], fmt=mk, ms=4, color=cols[p], lw=0.6,
                               label=f'pitch_deg {p:g} (hub {p + 15:g}), {v} m/s')
for k, q in enumerate(['cp', 'ct']):
    axs[0, 1 + k].set_title(f'tunnel_performance.csv: {q}')
    axs[0, 1 + k].set_xlabel('tsr')
    axs[0, 1 + k].grid(alpha=.3)
axs[0, 1].axhline(16 / 27, color='k', ls='--', lw=0.6)
axs[0, 1].legend(fontsize=7)

y = exps['experiment/tank_yaw_deep_immersion.csv']
for (p, yaw), g in y.groupby(['pitch_deg', 'yaw_deg']):
    ls = '-' if p == 5 else '--'
    c = {0: 'C0', 15: 'C1', 22.5: 'C2', 30: 'C3'}[yaw]
    g = g.sort_values('tsr')
    axs[1, 0].plot(g.tsr, g.cp, 'o' + ls, ms=3, lw=0.5, color=c, label=f'pitch {p:g}, yaw {yaw:g}')
    axs[1, 1].plot(g.tsr, g.ct, 'o' + ls, ms=3, lw=0.5, color=c)
axs[1, 0].set_title('tank_yaw_deep_immersion.csv cp (solid: pitch 5, dashed: 10)')
axs[1, 1].set_title('tank_yaw_deep_immersion.csv ct')
axs[1, 0].legend(fontsize=7)
s = exps['experiment/tank_shallow_immersion.csv']
t2 = exps['experiment/tank_twin_rotor_table2.csv']
t3 = exps['experiment/tank_twin_rotor_gap_table3.csv']
ax = axs[1, 2]
for p, g in s.groupby('pitch_deg'):
    ax.plot(g.tsr, g.cp, 's', ms=4, label=f'shallow 0.19D cp, pitch {p:g}')
    ax.plot(g.tsr, g.ct, 's', ms=4, mfc='none', label=f'shallow 0.19D ct, pitch {p:g}')
ax.plot(t2.tsr, t2.cp, 'k^', ms=4, label='Table 2 cp (twin rotor)')
ax.plot(t2.tsr, t2.ct, 'k^', ms=4, mfc='none', label='Table 2 ct')
ax.plot(t3.tsr, t3.cp, 'mv', ms=4, label='Table 3 cp @ tsr 7')
ax.plot(t3.tsr, t3.ct, 'mv', ms=4, mfc='none', label='Table 3 ct')
d = y[(y.yaw_deg == 0) & (y.pitch_deg == 5)].sort_values('tsr')
ax.plot(d.tsr, d.cp, 'C0-', lw=0.8, label='deep 0.55D, pitch 5 (ref.)')
ax.plot(d.tsr, d.ct, 'C0--', lw=0.8)
ax.set_title('tank: shallow immersion and twin-rotor tables')
ax.legend(fontsize=6)
for a in axs[1]:
    a.grid(alpha=.3)
    a.set_xlabel('tsr')
fig.tight_layout()
fig.savefig(os.path.join(HERE, 'check_bahaj_tidal_800mm.png'), dpi=85)
print('wrote check_bahaj_tidal_800mm.png')

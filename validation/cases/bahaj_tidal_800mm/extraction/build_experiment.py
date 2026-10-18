"""Build the experiment CSVs of the case from the vector-extracted markers (raw_markers.json, made by
vector_extract.py) of Bahaj et al. 2007, Renewable Energy 32:407-426, Figs. 7-11, plus Tables 2-3 (typed in).

Pitch convention (see case.json notes): the paper's "hub pitch angle" is the nose-tail angle of the section at
r = 80 mm (0.2R) measured from the rotor plane; Table 1 "pitch distribution" is 15.0 deg at 0.2R and 0 at the tip,
so it is the twist relative to the tip. geometry.csv uses that distribution as twist_deg, and the global
pitch_deg written here is the tip angle = hub_pitch - 15.0 deg.

"Scale adjusted" upper panels of Figs. 8 and 9: the authors shifted each data set vertically for clarity.
The shift of each set is recovered by matching its curve fit in the upper panel against the same curve fit
(same dash pattern) in the lower, unshifted panel; the shift is constant along the curve (std < 0.01) and is
a round multiple of 0.1, which is then removed from the markers.

Run: python3 build_experiment.py  -> ../experiment/*.csv, overlay_*.png, crosscheck printout
"""
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.join(HERE, '..', 'experiment')
R = json.load(open(os.path.join(HERE, 'raw_markers.json')))
TIP_FROM_HUB = 15.0   # Table 1: pitch distribution at 0.2R is 15.0 deg (tip = 0)


def kind_of(m):
    k = m['kind']
    if k.startswith('?') and '45' in k and '135' in k:
        return 'x'   # an 'x' marker that touched a frame/tick stroke
    return k


def pts(panel, kind):
    a = [(m['x'], m['y']) for m in R[panel]['markers'] if kind_of(m) == kind]
    return sorted(a)


# --- offsets of the scale-adjusted upper panels (curve drawing ids: upper, lower) -------------------
CURVE_PAIRS = {
    ('fig8', 'cp'): {'*': (10, 35), 'x': (12, 36), 'diamond': (14, 37), 'square': (16, 38)},
    ('fig8', 'ct'): {'*': (50, 75), 'x': (52, 76), 'diamond': (54, 77), 'square': (56, 78)},
    ('fig9', 'cp'): {'x': (9, 31), 'diamond': (11, 32), 'square': (13, 33)},
    ('fig9', 'ct'): {'x': (44, 63), 'diamond': (45, 64)},
}
OFFSETS = {}
for (fig, q), d in CURVE_PAIRS.items():
    for kind, (du, dl) in d.items():
        cu = [c for c in R[f'{fig}_{q}_upper']['curves'] if c['drawing'] == du][0]
        cl = [c for c in R[f'{fig}_{q}_lower']['curves'] if c['drawing'] == dl][0]
        xu, yu = np.array(cu['x']), np.array(cu['y'])
        xl, yl = np.array(cl['x']), np.array(cl['y'])
        o = np.argsort(xl)
        m = (xu >= xl.min()) & (xu <= xl.max())
        dy = yu[m] - np.interp(xu[m], xl[o], yl[o])
        off = round(float(np.median(dy)) / 0.1) * 0.1
        assert abs(np.median(dy) - off) < 0.01 and dy.std() < 0.01, (fig, q, kind, np.median(dy), dy.std())
        OFFSETS[(fig, q, kind)] = off
        print(f'offset {fig} {q} {kind:8s}: median {np.median(dy):+.4f} std {dy.std():.4f} -> {off:+.1f}')


def pair_cp_ct(cp, ct, tol=0.0045):
    """Markers of the same test point have the same TSR in the Cp and the Ct panel (checked: identical to
    <0.003). Pair them greedily by smallest TSR difference; a point present in only one panel keeps NaN in the
    other column."""
    cand = sorted((abs(xc - xp), i, j) for i, (xp, _) in enumerate(cp) for j, (xc, _) in enumerate(ct)
                  if abs(xc - xp) < tol)
    ui, uj, rows = set(), set(), []
    for d, i, j in cand:
        if i in ui or j in uj:
            continue
        ui.add(i)
        uj.add(j)
        rows.append(((cp[i][0] + ct[j][0]) / 2, cp[i][1], ct[j][1]))
    rows += [(x, y, np.nan) for i, (x, y) in enumerate(cp) if i not in ui]
    rows += [(x, np.nan, y) for j, (x, y) in enumerate(ct) if j not in uj]
    return sorted(rows)


def frame(rows, **cols):
    df = pd.DataFrame(rows, columns=['tsr', 'cp', 'ct'])
    for k, v in cols.items():
        df[k] = v
    return df


def write(df, name, header):
    path = os.path.join(EXP, name)
    with open(path, 'w') as f:
        for line in header.strip('\n').split('\n'):
            f.write('# ' + line + '\n')
        df.to_csv(f, index=False, float_format='%.4f')
    print('wrote', name, len(df), 'rows')


RHO, NU = 999.1, 1.14e-6   # fresh water at ~15 deg C ASSUMED (temperature not reported)
AREA, RADIUS = np.pi * 0.4 ** 2, 0.4
EXTR_ERR = 0.003            # extraction uncertainty of Cp/Ct (vector extraction; marker centre vs value)
DQ, DT = 0.25, 5.0          # dynamometer accuracy (paper Sec. 4.3): torque +-0.25 N m, thrust +-5 N


def add_errors(df):
    """cp_err / ct_err: root-sum-square of the dynamometer accuracy converted at the row's speed and TSR, and the
    extraction uncertainty. Blockage-correction uncertainty is NOT included (not quantified by the authors)."""
    q = 0.5 * RHO * df['v_ms'] ** 2 * AREA
    omega = df['tsr'] * df['v_ms'] / RADIUS
    cp_meas = DQ * omega / (q * df['v_ms'])
    ct_meas = DT / q
    df['cp_err'] = np.sqrt(cp_meas ** 2 + EXTR_ERR ** 2).round(4)
    df['ct_err'] = np.sqrt(ct_meas ** 2 + EXTR_ERR ** 2).round(4)
    df.loc[df['cp'].isna(), 'cp_err'] = np.nan
    df.loc[df['ct'].isna(), 'ct_err'] = np.nan
    return df


COMMON = """Bahaj AS, Molland AF, Chaplin JR, Batten WMJ (2007) Power and thrust measurements of marine current turbines
under various hydrodynamic flow conditions in a cavitation tunnel and a towing tank. Renewable Energy 32(3):407-426,
doi:10.1016/j.renene.2006.01.012.
Extraction: VECTOR extraction of the marker centres from the publisher PDF drawing commands (PyMuPDF), axes calibrated
on the gridlines (residual < 0.003 TSR, < 0.0003 in C). Scripts: extraction/vector_extract.py + build_experiment.py;
overlays extraction/overlay_fig*.png. Every marker of the figure is one row (repeated TSRs = separate measurements).
Extraction uncertainty: +-0.005 in TSR, +-0.003 in Cp/Ct.
cp_err/ct_err = RSS of the dynamometer accuracy (paper Sec. 4.3: torque +-0.25 N m, thrust +-5 N, converted at the
row's v_ms and tsr) and the extraction uncertainty; blockage-correction uncertainty not included.
Coefficients: Cp = Q*Omega/(0.5 rho U^3 A), Ct = T/(0.5 rho U^2 A), A = pi*0.4^2, TSR = Omega R / U.
BLOCKAGE-CORRECTED values (paper Sec. 6.2 and App. A: Barnsley & Wellicome actuator-disc wake-blockage
correction to the equivalent free-stream speed U_F; tunnel 2.4 x 1.2 m -> A/C = 0.175; up to -18 % Cp, -11 % Ct
in the tunnel and -8 % / -5 % in the tank at Ct = 0.8). Compare with the open-water (unbounded) BEM directly.
pitch_deg = solver global pitch = TIP chord angle = hub_pitch_deg - 15.0, where hub_pitch_deg (extra column, the
paper's 'hub pitch angle') is the nose-tail angle at r = 80 mm = 0.2R from the rotor plane and geometry.csv
twist_deg is the Table 1 'pitch distribution' (15.0 deg at 0.2R ... 0 at the tip).
v_ms = nominal tunnel / carriage speed of the figure legend or caption (the blockage-corrected U_F is a few % lower
and not given; v_ms only sets the Reynolds number here). Torque/thrust from the hub dynamometer, ahead of bearings
and seals (no mechanical losses in Cp). Thrust = rotor thrust only (dynamometer between rotor and shaft).
Extra columns ignored by the runner: hub_pitch_deg, facility, tip_immersion_D, figure."""

# ===================================================================== Fig. 7 (+ Fig. 8 1.0 m/s): tunnel
FIG7 = [('+', 15, 1.40), ('diamond', 20, 1.73), ('circle', 25, 1.54), ('square', 27, 1.30), ('*', 30, 1.54)]
dfs = []
for kind, hub, v in FIG7:
    rows = pair_cp_ct(pts('fig7a_cp', kind), pts('fig7b_ct', kind))
    dfs.append(frame(rows, v_ms=v, pitch_deg=hub - TIP_FROM_HUB, hub_pitch_deg=float(hub),
                     facility='cavitation_tunnel', figure='7'))
fig7 = pd.concat(dfs, ignore_index=True)


def fig_upper(fig, kind):
    cp = [(x, y - OFFSETS[(fig, 'cp', kind)]) for x, y in pts(f'{fig}_cp_upper', kind)]
    ct = [(x, y - OFFSETS[(fig, 'ct', kind)]) for x, y in pts(f'{fig}_ct_upper', kind)] \
        if (fig, 'ct', kind) in OFFSETS else []
    return pair_cp_ct(cp, ct)


def f89(fig, kind, v, hub, fac, h):
    return frame(fig_upper(fig, kind), v_ms=v, pitch_deg=hub - TIP_FROM_HUB, hub_pitch_deg=float(hub),
                 facility=fac, tip_immersion_D=h, figure=fig[3:])


# Fig. 8 / 9 series (the 1.7 and 1.5 m/s tunnel sets and the deep-tank sets duplicate Fig. 7 and Figs. 10/11
# point by point - checked in crosscheck.txt - so only the sets that exist nowhere else are written)
tun10 = f89('fig8', '*', 1.0, 20, 'cavitation_tunnel', np.nan)
tun17 = f89('fig8', 'x', 1.7, 20, 'cavitation_tunnel', np.nan)
tun15 = f89('fig9', 'x', 1.5, 25, 'cavitation_tunnel', np.nan)
deep8 = f89('fig8', 'diamond', 1.4, 20, 'towing_tank', 0.55)
deep9 = f89('fig9', 'diamond', 1.4, 25, 'towing_tank', 0.55)
shal8 = f89('fig8', 'square', 1.5, 20, 'towing_tank', 0.19)
shal9 = f89('fig9', 'square', 1.2, 25, 'towing_tank', 0.19)

tunnel = add_errors(pd.concat([fig7, tun10], ignore_index=True).drop(columns=['tip_immersion_D'], errors='ignore'))
tunnel = tunnel[['tsr', 'cp', 'ct', 'cp_err', 'ct_err', 'v_ms', 'pitch_deg', 'hub_pitch_deg', 'facility', 'figure']]
write(tunnel, 'tunnel_performance.csv', COMMON + """
SOURCE: Fig. 7 (p. 417) "Comparisons of hub pitch angles at the cavitation tunnel - zero yaw", (a) Cp, (b) Ct
  (figure = 7): hub pitch 15 deg at 1.40 m/s ('+'), 20 deg at 1.73 m/s (diamond), 25 deg at 1.54 m/s (circle),
  27 deg at 1.30 m/s (square), 30 deg at 1.54 m/s (asterisk);
plus Fig. 8 (p. 418) upper-left/right panels, cavitation tunnel 1.0 m/s, 20 deg hub pitch (asterisk; figure = 8).
  Fig. 8's upper panels are 'scale adjusted' (each data set shifted vertically by the authors); the shift (+0.3 in
  Cp, +0.6 in Ct for this set) was recovered by matching the set's curve fit with the unshifted fit in the lower
  panel (constant to < 0.002) and removed.
QinetiQ Haslar cavitation tunnel, working section 2.4 m x 1.2 m x 5 m, rotor near the tunnel centre, zero yaw,
non-cavitating pressure. Cp and Ct markers paired by identical TSR; a point present in only one panel has the other
column empty (25 deg: Ct-only points at TSR 7.32 and 7.95). The Fig. 8 1.7 m/s and Fig. 9 1.5 m/s tunnel sets are the
same measurements as the 20 deg and 25 deg sets of Fig. 7 (identical TSRs, |dCp|,|dCt| < 0.003) and are not repeated.
The 1.0 vs 1.7 m/s sets at 20 deg show the (small) Reynolds-number effect.""")

yawdfs = []
for figno, hub, series in [('10', 20, [('diamond', 0.0), ('square', 15.0), ('circle', 22.5), ('+', 30.0)]),
                           ('11', 25, [('diamond', 0.0), ('square', 15.0), ('+', 30.0)])]:
    for kind, yaw_ in series:
        yawdfs.append(frame(pair_cp_ct(pts(f'fig{figno}a_cp', kind), pts(f'fig{figno}b_ct', kind)), v_ms=1.4,
                            pitch_deg=hub - TIP_FROM_HUB, yaw_deg=yaw_, hub_pitch_deg=float(hub),
                            facility='towing_tank', tip_immersion_D=0.55, figure=figno))
yaw = add_errors(pd.concat(yawdfs, ignore_index=True))
yaw = yaw[['tsr', 'cp', 'ct', 'cp_err', 'ct_err', 'v_ms', 'pitch_deg', 'yaw_deg', 'hub_pitch_deg', 'facility',
           'tip_immersion_D', 'figure']]
write(yaw, 'tank_yaw_deep_immersion.csv', COMMON + """
SOURCE: Fig. 10 (p. 420) "Effect of yaw with 20 deg hub pitch at towed speed of 1.4 m/s" (yaw 0, 15, 22.5, 30 deg)
and Fig. 11 "Effect of yaw with 25 deg hub pitch at towed speed of 1.4 m/s" (yaw 0, 15, 30 deg); (a) Cp, (b) Ct.
Markers: 0 deg diamond, 15 deg square, 22.5 deg circle, 30 deg '+'.
Towing tank (Southampton Institute) 60 m x 3.7 m x 1.8 m deep, A/C = 0.0755, carriage 1.4 m/s.
yaw_deg = angle between the rotor axis and the towing direction. Cp and Ct are referred to the full carriage speed
and the full disc area (Eqs. 1-3), not to U*cos(yaw); thrust is along the rotor axis (shaft dynamometer).
tip_immersion_D = 0.55 (deep) is INFERRED, not stated: the 0-deg-yaw points of Fig. 10 / Fig. 11 are point by point
the deep-immersion tank points of Fig. 8 / Fig. 9 (|dCp|, |dCt| < 0.001; crosscheck.txt). Fig. 8's legend gives
1.5 m/s for that deep-tank set, Fig. 10's caption 1.4 m/s: 1.4 is used (affects only the Reynolds number).
(The brief mentioned yaw 0/7.5/15 deg: the paper only has 0, 15, 22.5, 30 deg (20 deg) and 0, 15, 30 deg (25 deg).)
The yaw-0 rows are the deep-tank vs tunnel comparison of Figs. 8/9: tank results slightly lower at low TSR.""")

shallow = add_errors(pd.concat([shal8, shal9], ignore_index=True))
shallow = shallow[['tsr', 'cp', 'ct', 'cp_err', 'ct_err', 'v_ms', 'pitch_deg', 'hub_pitch_deg', 'facility',
                   'tip_immersion_D', 'figure']]
write(shallow, 'tank_shallow_immersion.csv', COMMON + """
SOURCE: Fig. 8 (p. 418, 20 deg hub pitch, 1.5 m/s) and Fig. 9 (p. 419, 25 deg hub pitch, 1.2 m/s), upper panels,
square markers: towing tank, SHALLOW tip immersion h = 0.19D (blade tip 150 mm below the free surface at its top
position). The upper panels are 'scale adjusted'; this set has zero shift (checked against the lower panel fits).
Fig. 9's Ct panel has no shallow-immersion data: the 25 deg rows have Cp only.
FREE-SURFACE EFFECT: Cp 10-15 % and Ct ~5 % below the deep-immersion results (paper Sec. 7.1); an unbounded BEM is
expected to over-predict these rows - they document the free-surface effect, not the rotor aerodynamics.""")

# ===================================================================== Tables 2 and 3 (tabulated)
t2 = pd.DataFrame([
    (0.25, 100, 6.6, 0.0, 0.42, 0.84), (0.25, 100, 6.6, 6.5, 0.43, 0.83), (0.25, 100, 6.6, 7.8, 0.44, 0.84),
    (0.25, 100, 9.4, 0.0, 0.32, 0.94), (0.25, 100, 9.4, 7.1, 0.30, 0.95), (0.25, 100, 9.4, 8.6, 0.30, 0.95),
    (0.25, 100, 9.4, 10.3, 0.30, 0.94),
    (0.5, 200, 6.8, 0.0, 0.43, 0.80), (0.5, 200, 6.8, 5.2, 0.41, 0.82), (0.5, 200, 6.8, 7.6, 0.42, 0.81),
    (0.5, 200, 6.8, 9.74, 0.41, 0.85), (0.5, 200, 8.8, 0.0, 0.32, 0.93), (0.5, 200, 8.8, 7.7, 0.27, 0.94),
    (0.5, 200, 8.8, 8.4, 0.26, 0.95), (0.5, 200, 8.8, 9.7, 0.26, 0.95),
    (1.0, 400, 7.4, 0.0, 0.41, 0.86), (1.0, 400, 7.4, 6.2, 0.42, 0.86), (1.0, 400, 7.4, 7.1, 0.42, 0.88),
    (1.0, 400, 7.4, 9.1, 0.42, 0.87),
], columns=['gap_over_R', 'gap_mm', 'tsr', 'tsr_dummy', 'cp', 'ct'])
t2['v_ms'] = 1.4
t2['pitch_deg'] = 20.0 - TIP_FROM_HUB
t2 = add_errors(t2)
t2['cp_err'] = np.sqrt(t2['cp_err'] ** 2 - EXTR_ERR ** 2 + 0.005 ** 2).round(4)   # rounding to 2 decimals
t2['ct_err'] = np.sqrt(t2['ct_err'] ** 2 - EXTR_ERR ** 2 + 0.005 ** 2).round(4)
t2['tsr_err'] = 0.1
t2['dummy_rotating'] = t2['tsr_dummy'] > 0
write(t2[['tsr', 'cp', 'ct', 'tsr_err', 'cp_err', 'ct_err', 'v_ms', 'pitch_deg', 'gap_over_R', 'gap_mm',
          'dummy_rotating', 'tsr_dummy']],
      'tank_twin_rotor_table2.csv', """Bahaj et al. 2007, Renewable Energy 32:407-426, Table 2 (p. 421) "Effect of dummy TSR on power and thrust
coefficient". TABULATED (typed from the paper). tsr given as e.g. '6.6 +- 0.1' -> tsr_err = 0.1; Cp, Ct to 2 decimals
(cp_err/ct_err = RSS of dynamometer accuracy and the +-0.005 rounding).
Towing tank, deep tip immersion h = 0.55D; the instrumented rotor runs beside an identical non-instrumented dummy rotor
with tip-to-tip gap g (gap_mm, gap_over_R); tsr_dummy = 0: dummy static. Extra blockage corrections were applied when
the dummy rotated. Blockage-corrected like all data of the paper.
ASSUMED (not stated for the twin-rotor tests): hub pitch 20 deg (pitch_deg = 5.0) - the static-dummy values match the
20 deg deep-tank curve of Fig. 8 within ~0.03 (crosscheck.txt); carriage speed 1.4 m/s (as the 20 deg deep-tank
runs; the tank tests used 0.8-1.5 m/s) - it only sets the Reynolds number.
The paper concludes the dummy rotor has no measurable influence, so these are effectively single-rotor points.""")

t3 = pd.DataFrame([
    (0.025, 10, 0.42, 0.40, 0.88, 0.86), (0.0625, 25, 0.42, 0.42, 0.88, 0.88), (0.125, 50, 0.41, 0.40, 0.88, 0.86),
    (0.25, 100, 0.42, 0.41, 0.89, 0.88), (0.375, 150, 0.42, 0.41, 0.88, 0.87), (0.5, 200, 0.42, 0.41, 0.86, 0.88),
    (1.0, 400, 0.41, 0.42, 0.88, 0.89)],
    columns=['gap_over_R', 'gap_mm', 'cp_st', 'cp_rot', 'ct_st', 'ct_rot'])
t3l = pd.concat([
    pd.DataFrame({'tsr': 7.0, 'cp': t3.cp_st, 'ct': t3.ct_st, 'gap_over_R': t3.gap_over_R, 'gap_mm': t3.gap_mm,
                  'dummy_rotating': False}),
    pd.DataFrame({'tsr': 7.0, 'cp': t3.cp_rot, 'ct': t3.ct_rot, 'gap_over_R': t3.gap_over_R, 'gap_mm': t3.gap_mm,
                  'dummy_rotating': True})], ignore_index=True)
t3l['v_ms'] = 1.4
t3l['pitch_deg'] = 20.0 - TIP_FROM_HUB
t3l = add_errors(t3l)
t3l['cp_err'] = np.sqrt(t3l['cp_err'] ** 2 - EXTR_ERR ** 2 + 0.005 ** 2).round(4)
t3l['ct_err'] = np.sqrt(t3l['ct_err'] ** 2 - EXTR_ERR ** 2 + 0.005 ** 2).round(4)
write(t3l[['tsr', 'cp', 'ct', 'cp_err', 'ct_err', 'v_ms', 'pitch_deg', 'gap_over_R', 'gap_mm', 'dummy_rotating']],
      'tank_twin_rotor_gap_table3.csv', """Bahaj et al. 2007, Renewable Energy 32:407-426, Table 3 (p. 421) "Effect of gap width on power and thrust
coefficient interpolated at TSR = 7". TABULATED (typed); the table's columns CP/CT (static) and CP/CT (rotating) are
written as rows with dummy_rotating = False / True. Values interpolated by the authors to TSR = 7 (tsr = 7.0 exact).
Towing tank, deep immersion. ASSUMED as for Table 2: 20 deg hub pitch (pitch_deg = 5.0), 1.4 m/s.
cp_err/ct_err = RSS of dynamometer accuracy and the +-0.005 rounding (the interpolation error is not known).""")

# keep the old names used by the cross-checks below
fig8 = pd.concat([tun17, deep8], ignore_index=True)
fig9 = pd.concat([tun15, deep9], ignore_index=True)
fig8['v_ms'] = fig8['v_ms'].where(fig8.facility != 'cavitation_tunnel', 1.7)
# ===================================================================== cross-checks
lines = []


def near(df, x, tol=0.01):
    d = (df['tsr'] - x).abs()
    return df.loc[d.idxmin()] if d.min() < tol else None


def compare(a, b, label):
    dcp, dct = [], []
    for _, r in a.iterrows():
        s = near(b, r['tsr'])
        if s is None:
            continue
        if np.isfinite(r['cp']) and np.isfinite(s['cp']):
            dcp.append(r['cp'] - s['cp'])
        if np.isfinite(r['ct']) and np.isfinite(s['ct']):
            dct.append(r['ct'] - s['ct'])
    lines.append(f'{label}: {len(dcp)} Cp pairs, max |dCp| = {np.max(np.abs(dcp)):.4f}; '
                 f'{len(dct)} Ct pairs, max |dCt| = {np.max(np.abs(dct)):.4f}')


compare(fig7[fig7.hub_pitch_deg == 20], fig8[(fig8.v_ms == 1.7)], 'Fig7 20deg vs Fig8 tunnel 1.7 m/s (upper, de-shifted)')
compare(fig7[fig7.hub_pitch_deg == 25], fig9[(fig9.facility == 'cavitation_tunnel')], 'Fig7 25deg vs Fig9 tunnel 1.5 m/s')
compare(yaw[(yaw.figure == '10') & (yaw.yaw_deg == 0)], fig8[fig8.facility == 'towing_tank'], 'Fig10 0deg yaw vs Fig8 deep tank')
compare(yaw[(yaw.figure == '11') & (yaw.yaw_deg == 0)], fig9[fig9.facility == 'towing_tank'], 'Fig11 0deg yaw vs Fig9 deep tank')
deep = fig8[fig8.facility == 'towing_tank'].dropna(subset=['cp', 'ct']).sort_values('tsr')
for _, r in t2[~t2.dummy_rotating.astype(bool)].iterrows():
    lines.append(f"Table2 static-dummy tsr {r.tsr}: Cp {r.cp} Ct {r.ct}  vs Fig.8 deep 20deg (interp) "
                 f"Cp {np.interp(r.tsr, deep.tsr, deep.cp):.3f} Ct {np.interp(r.tsr, deep.tsr, deep.ct):.3f}")
for (fig, q, kind), off in OFFSETS.items():
    lines.append(f'upper-panel shift {fig} {q} {kind}: {off:+.1f}')
open(os.path.join(HERE, 'crosscheck.txt'), 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))

# audit file: all de-shifted Fig. 8/9 sets (including the ones that duplicate Figs. 7/10/11), for overlay.py
pd.concat([tun10, tun17, deep8, shal8, tun15, deep9, shal9], ignore_index=True).to_csv(
    os.path.join(HERE, 'fig8_fig9_deshifted_all_sets.csv'), index=False, float_format='%.4f')

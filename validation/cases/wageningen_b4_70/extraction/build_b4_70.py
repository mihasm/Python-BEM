"""Builds validation/cases/wageningen_b4_70 (Wageningen B-screw series B4-70, P/D = 1.0) from the series definition.

Sources (copies in ../sources/):
  Oosterveld M.W.C., van Oossanen P. (1975). Further computer-analyzed data of the Wageningen B-screw series.
    International Shipbuilding Progress 22(251), 251-262. doi:10.3233/ISP-1975-2225102.
    Table 1 (blade contour c Z/(D AE/A0), a/c, b/c, t/D = Ar - Br Z), Tables 2/3 (section functions V1, V2),
    Table 5 (KT, KQ polynomials at Rn = 2e6), eq. 7 (Rn at 0.75R).
  Kuiper G. (1992). The Wageningen Propeller Series. MARIN Publ. 92-001: sect. 4.2 + Fig. 4.1 (four-bladed pitch
    reduction 82.2/88.7/95/99.2 % at 0.2/0.3/0.4/0.5R), sect. 4.5 (pitch defined on the face), 4.7 (rake 15 deg),
    4.8 (hub d/D = 1/6), Tables 4.7-4.12 (same section data, cross-check).
Run: python3 build_b4_70.py [--xfoil]
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from marine_common import write_selig, run_xfoil_set, write_csv, propeller_beta_deg  # noqa

CASE = os.path.dirname(HERE)
REPO = os.path.abspath(os.path.join(CASE, '..', '..', '..'))
XFOIL = os.path.join(REPO, 'xfoil_executables', 'xfoil')
Z, EAR, PD = 4, 0.70, 1.0
D = 0.240            # ASSUMPTION: typical NSMB B-series model diameter; results are made Re-equivalent anyway
R = D / 2
RN = 2e6
RHO, NU = 1000.0, 1.139e-6
RES = [2.5e5, 5e5, 1e6, 2e6, 3e6]
NCRIT = 9

# Oosterveld & van Oossanen 1975 Table 1 (four- to seven-bladed): r/R, c Z/(D AE/A0), b/c, Ar, Br
T1 = [(0.2, 1.662, 0.350, 0.0526, 0.0040), (0.3, 1.882, 0.350, 0.0464, 0.0035), (0.4, 2.050, 0.351, 0.0402, 0.0030),
      (0.5, 2.152, 0.355, 0.0340, 0.0025), (0.6, 2.187, 0.389, 0.0278, 0.0020), (0.7, 2.144, 0.443, 0.0216, 0.0015),
      (0.8, 1.970, 0.479, 0.0154, 0.0010), (0.9, 1.582, 0.500, 0.0092, 0.0005), (1.0, 0.000, 0.500, 0.0030, 0.0)]
PITCH_RED = {0.2: 0.822, 0.3: 0.887, 0.4: 0.950, 0.5: 0.992}   # Kuiper Fig. 4.1, B4 only
# BB-series contour (Oosterveld Table 4), used only for the SHAPE of the tip between 0.9R and 1.0R
BB = {0.9: 1.798, 0.95: 1.434, 0.975: 1.122}

# V1, V2 (Oosterveld Tables 2, 3; Kuiper Tables 4.7-4.10). P = +1 LE ... 0 max thickness ... -1 TE
P_LE = [1.0, .95, .9, .85, .8, .7, .6, .5, .4, .2, 0.0]
P_TE = [-1.0, -.95, -.9, -.8, -.7, -.6, -.5, -.4, -.2, 0.0]
R_ROWS = [0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 1.0]
V1_LE = {0.15: [.3860, .3150, .2642, .2230, .1870, .1320, .0920, .0615, .0384, .0096, 0],
         0.2: [.3560, .2821, .2353, .2000, .1685, .1180, .0804, .0520, .0304, .0049, 0],
         0.25: [.3256, .2513, .2068, .1747, .1465, .1008, .0669, .0417, .0224, .0031, 0],
         0.3: [.2923, .2186, .1760, .1445, .1191, .0790, .0503, .0300, .0148, .0027, 0],
         0.4: [.2181, .1467, .1088, .0833, .0637, .0357, .0189, .0090, .0033, 0, 0],
         0.5: [.1278, .0778, .0500, .0328, .0211, .0085, .0034, .0008, 0, 0, 0],
         0.6: [.0382, .0169, .0067, .0022, .0006, 0, 0, 0, 0, 0, 0]}
V1_TE = {0.15: [.3000, .2824, .2650, .2300, .1950, .1610, .1280, .0955, .0365, 0],
         0.2: [.2826, .2630, .2400, .1967, .1570, .1207, .0880, .0592, .0172, 0],
         0.25: [.2598, .2372, .2115, .1651, .1246, .0899, .0579, .0350, .0084, 0],
         0.3: [.2306, .2040, .1790, .1333, .0943, .0623, .0376, .0202, .0033, 0],
         0.4: [.1467, .1200, .0972, .0630, .0395, .0214, .0116, .0044, 0, 0],
         0.5: [.0522, .0420, .0330, .0190, .0100, .0040, .0012, 0, 0, 0]}
V2_LE = {0.15: [0, .1300, .2600, .3665, .4520, .5995, .7105, .8055, .8825, .9760, 1],
         0.2: [0, .1560, .2840, .3905, .4777, .6190, .7277, .8170, .8875, .9750, 1],
         0.25: [0, .1758, .3042, .4108, .4982, .6359, .7415, .8259, .8899, .9751, 1],
         0.3: [0, .1890, .3197, .4265, .5130, .6505, .7520, .8315, .8920, .9750, 1],
         0.4: [0, .1935, .3235, .4335, .5220, .6590, .7593, .8345, .8933, .9725, 1],
         0.5: [0, .1750, .3056, .4135, .5039, .6430, .7478, .8275, .8880, .9710, 1],
         0.6: [0, .1485, .2720, .3775, .4620, .6060, .7200, .8090, .8790, .9690, 1],
         0.7: [0, .1240, .2337, .3300, .4140, .5615, .6840, .7850, .8660, .9675, 1],
         0.8: [0, .1050, .2028, .2925, .3765, .5265, .6545, .7635, .8520, .9635, 1],
         0.85: [0, .1000, .1950, .2830, .3660, .5160, .6455, .7550, .8450, .9615, 1],
         0.9: [0, .0975, .1900, .2775, .3600, .5100, .6400, .7500, .8400, .9600, 1]}
V2_TE = {0.15: [0, .0540, .1325, .2870, .4280, .5585, .6770, .7805, .9360, 1],
         0.2: [0, .0640, .1455, .3060, .4535, .5842, .6995, .7984, .9446, 1],
         0.25: [0, .0725, .1567, .3228, .4740, .6050, .7184, .8139, .9519, 1],
         0.3: [0, .0800, .1670, .3360, .4885, .6195, .7335, .8265, .9583, 1],
         0.4: [0, .0905, .1810, .3500, .5040, .6353, .7525, .8415, .9645, 1],
         0.5: [0, .0950, .1865, .3569, .5140, .6439, .7580, .8456, .9639, 1],
         0.6: [0, .0965, .1885, .3585, .5110, .6415, .7530, .8426, .9613, 1],
         0.7: [0, .0975, .19, .36, .51, .64, .75, .84, .96, 1]}

# Oosterveld & van Oossanen 1975 Table 5: (C, s, t, u, v), KT/KQ = sum C J^s (P/D)^t (AE/A0)^u Z^v, Rn = 2e6
KT_C = [(0.00880496, 0, 0, 0, 0), (-0.204554, 1, 0, 0, 0), (0.166351, 0, 1, 0, 0), (0.158114, 0, 2, 0, 0),
        (-0.147581, 2, 0, 1, 0), (-0.481497, 1, 1, 1, 0), (0.415437, 0, 2, 1, 0), (0.0144043, 0, 0, 0, 1),
        (-0.0530054, 2, 0, 0, 1), (0.0143481, 0, 1, 0, 1), (0.0606826, 1, 1, 0, 1), (-0.0125894, 0, 0, 1, 1),
        (0.0109689, 1, 0, 1, 1), (-0.133698, 0, 3, 0, 0), (0.00638407, 0, 6, 0, 0), (-0.00132718, 2, 6, 0, 0),
        (0.168496, 3, 0, 1, 0), (-0.0507214, 0, 0, 2, 0), (0.0854559, 2, 0, 2, 0), (-0.0504475, 3, 0, 2, 0),
        (0.010465, 1, 6, 2, 0), (-0.00648272, 2, 6, 2, 0), (-0.00841728, 0, 3, 0, 1), (0.0168424, 1, 3, 0, 1),
        (-0.00102296, 3, 3, 0, 1), (-0.0317791, 0, 3, 1, 1), (0.018604, 1, 0, 2, 1), (-0.00410798, 0, 2, 2, 1),
        (-0.000606848, 0, 0, 0, 2), (-0.0049819, 1, 0, 0, 2), (0.0025983, 2, 0, 0, 2), (-0.000560528, 3, 0, 0, 2),
        (-0.00163652, 1, 2, 0, 2), (-0.000328787, 1, 6, 0, 2), (0.000116502, 2, 6, 0, 2), (0.000690904, 0, 0, 1, 2),
        (0.00421749, 0, 3, 1, 2), (0.0000565229, 3, 6, 1, 2), (-0.00146564, 0, 3, 2, 2)]
KQ_C = [(0.00379368, 0, 0, 0, 0), (0.00886523, 2, 0, 0, 0), (-0.032241, 1, 1, 0, 0), (0.00344778, 0, 2, 0, 0),
        (-0.0408811, 0, 1, 1, 0), (-0.108009, 1, 1, 1, 0), (-0.0885381, 2, 1, 1, 0), (0.188561, 0, 2, 1, 0),
        (-0.00370871, 1, 0, 0, 1), (0.00513696, 0, 1, 0, 1), (0.0209449, 1, 1, 0, 1), (0.00474319, 2, 1, 0, 1),
        (-0.00723408, 2, 0, 1, 1), (0.00438388, 1, 1, 1, 1), (-0.0269403, 0, 2, 1, 1), (0.0558082, 3, 0, 1, 0),
        (0.0161886, 0, 3, 1, 0), (0.00318086, 1, 3, 1, 0), (0.015896, 0, 0, 2, 0), (0.0471729, 1, 0, 2, 0),
        (0.0196283, 3, 0, 2, 0), (-0.0502782, 0, 1, 2, 0), (-0.030055, 3, 1, 2, 0), (0.0417122, 2, 2, 2, 0),
        (-0.0397722, 0, 3, 2, 0), (-0.00350024, 0, 6, 2, 0), (-0.0106854, 3, 0, 0, 1), (0.00110903, 3, 3, 0, 1),
        (-0.000313912, 0, 6, 0, 1), (0.0035985, 3, 0, 1, 1), (-0.00142121, 0, 6, 1, 1), (-0.00383637, 1, 0, 2, 1),
        (0.0126803, 0, 2, 2, 1), (-0.00318278, 2, 3, 2, 1), (0.00334268, 0, 6, 2, 1), (-0.00183491, 1, 1, 0, 2),
        (0.000112451, 3, 2, 0, 2), (-0.0000297228, 3, 6, 0, 2), (0.000269551, 1, 0, 1, 2), (0.00083265, 2, 0, 1, 2),
        (0.00155334, 0, 2, 1, 2), (0.000302683, 0, 6, 1, 2), (-0.0001843, 0, 0, 2, 2), (-0.000425399, 0, 3, 2, 2),
        (0.0000869243, 3, 3, 2, 2), (-0.0004659, 0, 6, 2, 2), (0.0000554194, 1, 6, 2, 2)]


def poly(C, J, pd, ear, z):
    return sum(c * J ** s * pd ** t * ear ** u * z ** v for c, s, t, u, v in C)


def _row(table, r, default_zero_above):
    """Linear interpolation in r/R between tabulated rows (rows missing above = zero / the last row)."""
    keys = sorted(table)
    if r > keys[-1] + 1e-9:
        return np.zeros(len(table[keys[0]])) if default_zero_above else np.array(table[keys[-1]])
    i = np.searchsorted(keys, r) - 1
    i = max(i, 0)
    k0, k1 = keys[i], keys[i + 1]
    w = (r - k0) / (k1 - k0)
    return (1 - w) * np.array(table[k0]) + w * np.array(table[k1])


def bseries_section(r, b_c, n=81):
    """Unit-chord B-series section at r/R with max thickness 1 (scaled later), edges t_le = t_te = 0.
    Returns x (LE->TE), y_face, t (thickness) in units of t_max."""
    x = 0.5 * (1 - np.cos(np.linspace(0, np.pi, n)))
    P = np.where(x <= b_c, (b_c - x) / b_c, -(x - b_c) / (1 - b_c))
    v1le, v2le = _row(V1_LE, r, True), _row(V2_LE, r, False)
    v1te, v2te = _row(V1_TE, r, True), _row(V2_TE, r, False)
    f1le, f2le = PchipInterpolator(P_LE[::-1], v1le[::-1]), PchipInterpolator(P_LE[::-1], v2le[::-1])
    f1te, f2te = PchipInterpolator(P_TE, v1te), PchipInterpolator(P_TE, v2te)
    yf = np.where(P >= 0, f1le(np.clip(P, 0, 1)), f1te(np.clip(P, -1, 0)))
    t = np.where(P >= 0, f2le(np.clip(P, 0, 1)), f2te(np.clip(P, -1, 0)))
    # rounded nose: the tabulated thickness ends in a wedge (t_le = 0); within 5 % of the LE->max-thickness
    # distance it is replaced by t(x0) sqrt(x/x0) (elliptic-type nose), which XFOIL needs for a sane LE.
    x0 = 0.05 * b_c
    t0 = float(np.interp(x0, x, t))
    t = np.where(x < x0, t0 * np.sqrt(np.clip(x / x0, 0, 1)), t)
    return x, yf, t


def make_section(r, t_c, b_c):
    x, yf, t = bseries_section(r, b_c)
    yf, t = yf * t_c, t * t_c
    yb = yf + t
    # nose-tail line: LE (x=0, y=yf[0]) to TE (x=1, y=yf[-1]); rotate so that it is the x axis
    delta = np.arctan2(yf[0] - yf[-1], 1.0)   # > 0: LE above the face reference line -> chord pitched up
    def rot(xx, yy):
        # rotate about the LE by +delta (brings the TE up to the LE height)
        xr, yr = xx, yy - yf[0]
        X = xr * np.cos(delta) - yr * np.sin(delta)
        Y = xr * np.sin(delta) + yr * np.cos(delta)
        return X, Y
    Xb, Yb = rot(x, yb)
    Xf, Yf = rot(x, yf)
    L = Xf[-1]
    Xb, Yb, Xf, Yf = Xb / L, Yb / L, Xf / L, Yf / L
    X = np.concatenate([Xb[::-1], Xf[1:]])
    Y = np.concatenate([Yb[::-1], Yf[1:]])
    camber = np.max(0.5 * (Yb + Yf))
    return X, Y, np.degrees(delta), camber


def fname(x):
    return f'b4_70_r{int(round(x * 1000)):04d}'


def main(do_xfoil):
    os.makedirs(os.path.join(CASE, 'airfoils'), exist_ok=True)
    os.makedirs(os.path.join(CASE, 'experiment'), exist_ok=True)
    t1 = pd.DataFrame(T1, columns=['r_R', 'K', 'b_c', 'Ar', 'Br'])
    stations = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975, 1.0]
    rows, foils = [], {}
    for x in stations:
        if x in BB:
            K = 1.582 * BB[x] / BB[0.9] if x > 0.9 else float(t1.K[np.isclose(t1.r_R, x)].iloc[0])
        else:
            K = float(np.interp(x, t1.r_R, t1.K))
        c_D = K * EAR / Z
        t_D = float(np.interp(x, t1.r_R, t1.Ar)) - float(np.interp(x, t1.r_R, t1.Br)) * Z
        b_c = float(np.interp(x, t1.r_R, t1.b_c))
        pd_face = PD * PITCH_RED.get(x, 1.0)
        beta_face = float(propeller_beta_deg(pd_face, x))
        if x <= 0.95:
            t_c = t_D / c_D
            X, Y, delta, camber = make_section(x, t_c, b_c)
            name = fname(x)
            p = os.path.join(CASE, 'airfoils', name + '.dat')
            write_selig(p, f'Wageningen B-series section r/R={x:.3f} t/c={t_c:.5f} (face-referenced, chord-rotated)', X, Y)
            foils[name] = (p, t_c, delta, camber)
            foil, dl, tc_out = name, delta, t_c
        else:
            foil = fname(0.95)
            dl = foils[fname(0.95)][2]
            tc_out = foils[fname(0.95)][1]
        rows.append(dict(r_m=x * R, chord_m=c_D * D, twist_deg=beta_face + dl, airfoil=foil, thickness_ratio=tc_out,
                         r_R=x, c_D=c_D, t_D=t_D, P_D_face=pd_face, beta_face_deg=beta_face, delta_chord_deg=dl))
    geo = pd.DataFrame(rows)
    write_csv(os.path.join(CASE, 'geometry.csv'), [
        'Wageningen B-screw series B4-70, P/D = 1.0 (Z = 4, AE/A0 = 0.70), GENERATED from the series definition:',
        '  Oosterveld & van Oossanen (1975) ISP 22(251) Table 1: c/D = K(r) (AE/A0)/Z, t/D = Ar - Br Z, b/c; Tables 2-3: V1, V2.',
        '  Kuiper (1992) MARIN 92-001 Fig. 4.1/sect. 4.2: B4 pitch reduction to 82.2/88.7/95.0/99.2 % at 0.2/0.3/0.4/0.5R',
        '  (the figure label at 0.5R reads "92.2%", an evident misprint for 99.2% - the curve is monotonic).',
        f'D = {D} m is an ASSUMED model size (dimensionless case; the experiment rows set rpm for Rn(0.75R) = 2e6).',
        'Hub d/D = 1/6 (Kuiper 4.8) -> Rhub; blade data start at 0.2R. Rake 15 deg and the small B-series skew are ignored.',
        'Tip: B-series contour is tabulated only to 0.9R (K = 1.582) and 0 at 1.0R; 0.95R/0.975R chords use the SHAPE of',
        '  the BB-series contour (Oosterveld Table 4: K 1.798/1.434/1.122 at 0.9/0.95/0.975R) scaled to the B value at 0.9R.',
        'twist_deg = beta_face + delta_chord: beta_face = atan(P_face/D / (pi r/R)) with the pitch defined on the FACE',
        '  (pressure side reference line, Kuiper 4.5); delta_chord = angle of the nose-tail line relative to the face line',
        '  of the generated section (non-zero only inside 0.6R, as Kuiper eq. 4.5 says; <= ~1.5 deg at 0.2R).',
        'airfoil: B-series sections generated at each station (edge thickness 0, rounded nose, see extraction/build_b4_70.py);',
        '  0.975R and 1.0R use the 0.95R section. Extra columns are documentation only.'], geo, float_format='%.6f')

    if do_xfoil:
        pol = run_xfoil_set(XFOIL, os.path.join(HERE, 'xfoil_raw'), {k: v[0] for k, v in foils.items()}, RES, ncrit=NCRIT)
        for name, df in pol.items():
            df['source'] = f'XFOIL 6.99 Ncrit={NCRIT}, 200 panels, free transition'
            p, tc, dl, cam = foils[name]
            write_csv(os.path.join(CASE, 'airfoils', name + '_polars.csv'), [
                f'{name}: Wageningen B-series section, t/c = {tc:.5f}, max camber of the generated section {cam:.5f} c,',
                f'  nose-tail line rotated by {dl:.3f} deg from the face reference line (included in twist_deg).',
                'COMPUTED polar: XFOIL 6.99 (repository binary), viscous, free transition, Ncrit = 9, 200 panels, ITER 150-300,',
                '  sweeps 0 -> +18 and 0 -> -10 deg (0.5 deg); only converged points. No measured polars of these sections exist.',
                '  The B-series polynomials are faired to Rn = 2e6 assuming turbulent friction (ITTC-57 line; Oosterveld eq. 8-10);',
                '  free-transition XFOIL gives lower drag than that. Script: extraction/build_b4_70.py --xfoil.'], df,
                float_format='%.5f')
            print(name, df.groupby('re').alpha_deg.agg(['count', 'min', 'max']).to_dict('index'))

    # performance from the Rn = 2e6 polynomials, rpm per row so that Rn(0.75R) = 2e6 (Oosterveld eq. 7)
    c075 = float(np.interp(0.75, geo.r_R, geo.chord_m))
    J = np.round(np.arange(0.0, 1.30001, 0.05), 3)
    kt = np.array([poly(KT_C, j, PD, EAR, Z) for j in J])
    kq = np.array([poly(KQ_C, j, PD, EAR, Z) for j in J])
    keep = kt > 0
    J, kt, kq = J[keep], kt[keep], kq[keep]
    nD = RN * NU / (c075 * np.sqrt(J ** 2 + (0.75 * np.pi) ** 2))
    rpm = nD / D * 60
    perf = pd.DataFrame({'J': J, 'ct': kt, 'cp': 2 * np.pi * kq, 'eta': J * kt / (2 * np.pi * kq), 'kt': kt, 'kq': kq,
                         'kq10': 10 * kq, 'rpm': rpm, 'v_ms': J * nD})
    write_csv(os.path.join(CASE, 'experiment', 'performance_polynomial_rn2e6.csv'), [
        'Wageningen B4-70, P/D = 1.0: open-water characteristics from the Oosterveld & van Oossanen (1975) Table 5',
        '  polynomials (regression of the NSMB open-water tests of 120 B-series models, corrected by Lerbs equivalent-profile',
        '  method to Rn = 2e6). NOT raw measurements but the faired series data (published regression; residual vs the',
        '  original tests of the order of +-2-3 % in KT/KQ). Coefficients copied from the paper (checked against the scan of',
        '  Table 5 and against Bernitsas et al. 1981 / Molland et al. A4.1 as reproduced in public code).',
        'ct = KT, cp = 2*pi*KQ, eta = J KT/(2 pi KQ); kt, kq, kq10 = polynomial values. Rows with KT <= 0 dropped.',
        f'rpm per row: n D = Rn nu / (c0.75 sqrt(J^2 + (0.75 pi)^2)) with Rn = 2e6, nu = {NU}, c0.75 = {c075:.5f} m (eq. 7),',
        '  so the runner reproduces the polynomial Reynolds number; v_ms = J n D.'], perf, float_format='%.6f')
    print(perf.iloc[::4])
    # cross-check of the implementation: MIT 2.611 notes (2006): B5-75, P/D 1.0 at J = 0.64 -> KT = 0.544 J^2 = 0.223
    print('check B5-75 P/D 1.0 J 0.64: KT =', round(poly(KT_C, 0.64, 1.0, 0.75, 5), 4), '(MIT notes: 0.223)')


if __name__ == '__main__':
    main('--xfoil' in sys.argv)

"""Builds validation/cases/dtmb_4119 (DTMB/DTRC propeller 4119, 3 blades, D = 1 ft) from tabulated sources.

Sources (copies in ../sources/):
  Geometry : Koyama K. (1993), Comparative calculations of propellers by surface panel method - workshop organized
             by the 20th ITTC Propulsor Committee. Papers of Ship Research Institute, Suppl. No. 15, Appendix
             'Input data sent', Table 1.1.1(a) 'Blade geometry of DTRC4119' (p. 65), tabulated (DTRC 'xyzprop' print).
             https://www.nmri.go.jp/service/repository_data/PNM19150001-03.pdf
  Open water: same appendix, Table 1.3.1 'Open water test results for DTRC4119' (p. 72), Jessup's experiment,
             DTMB towing basin; Fig. 1.3.1 (p. 75) says Rn = 1 x 10^6.
  Sections : NACA 66 (TMB modified nose and tail) + a = 0.8 mean line, Brockett (1966) DTMB Report 1780 Tables 1/3.
Run: python3 build_dtmb4119.py [--xfoil]   (--xfoil runs XFOIL for missing combinations; without it the polars
     are rebuilt from the .pol files in extraction/xfoil_raw)
"""
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from marine_common import (naca66mod_a08, write_selig, run_xfoil_set, write_csv, propeller_beta_deg,  # noqa
                           complete_polars, polar_header_lines)

CASE = os.path.dirname(HERE)
REPO = os.path.abspath(os.path.join(CASE, '..', '..', '..'))
XFOIL = os.path.join(REPO, 'xfoil_executables', 'xfoil')
D = 12.0 * 0.0254          # 0.3048 m
R = D / 2
B = 3

# Koyama 1993 Table 1.1.1(a): r/R, c/D, P/D, t/c (T/CHORD), f/c (FM/CHORD); rows reordered root -> tip
GEOM = [
    (0.200, 0.32000, 1.10500, 0.20550, 0.01429),
    (0.250, 0.34200, 1.10370, 0.17870, 0.01985),
    (0.300, 0.36350, 1.10220, 0.15530, 0.02318),
    (0.400, 0.40480, 1.09830, 0.11800, 0.02303),
    (0.500, 0.43920, 1.09320, 0.09016, 0.02182),
    (0.600, 0.46100, 1.08790, 0.06960, 0.02072),
    (0.700, 0.46220, 1.08390, 0.05418, 0.02003),
    (0.800, 0.43470, 1.08110, 0.04206, 0.01967),
    (0.900, 0.36130, 1.07850, 0.03321, 0.01817),
    (0.925, 0.32621, 1.07779, 0.03252, 0.01744),
    (0.950, 0.27750, 1.07700, 0.03228, 0.01631),
    (0.975, 0.20449, 1.07609, 0.03211, 0.01450),
    (0.990, 0.13276, 1.07546, 0.03187, 0.01298),
    (0.995, 0.09479, 1.07523, 0.03175, 0.01239),
    (1.000, 0.00000, 1.07500, 0.03160, 0.01175),
]
# airfoils are built at these stations; the others are 'transition' (or the 0.95 section towards the tip)
FOIL_STATIONS = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
RES = [2e5, 4e5, 7e5, 1e6, 1.5e6]
NCRIT = 9
RHO, NU = 998.2, 1.0e-6     # fresh water ~20 C (not stated; Rn = 1e6 at 0.7R with n = 10 rps implies nu ~ 1.0e-6)


def fname(x):
    return f'dtmb4119_r{int(round(x * 1000)):04d}'


def main(do_xfoil):
    os.makedirs(os.path.join(CASE, 'airfoils'), exist_ok=True)
    os.makedirs(os.path.join(CASE, 'experiment'), exist_ok=True)
    g = pd.DataFrame(GEOM, columns=['r_R', 'c_D', 'P_D', 't_c', 'f_c'])
    rows = []
    for _, s in g.iterrows():
        x = s.r_R
        if any(abs(x - fs) < 1e-6 for fs in FOIL_STATIONS):
            foil = fname(x)
        elif x > 0.95:
            foil = fname(0.95)
        else:
            foil = 'transition'
        rows.append(dict(r_m=x * R, chord_m=s.c_D * D, twist_deg=float(propeller_beta_deg(s.P_D, x)), airfoil=foil,
                         thickness_ratio=s.t_c, r_R=x, P_D=s.P_D, camber_ratio=s.f_c))
    geo = pd.DataFrame(rows)
    write_csv(os.path.join(CASE, 'geometry.csv'), [
        'DTMB (DTRC) propeller 4119 - 3 blades, D = 12.000 in = 0.3048 m, no skew, no rake, hub d/D = 0.2',
        'Source: Koyama (1993) Papers of Ship Research Institute Suppl. 15, Table 1.1.1(a) "Blade geometry of DTRC4119"',
        '  (p. 65, DTRC xyzprop print; TABULATED, copied unchanged; rows re-sorted root->tip). Same table: 22nd ITTC',
        '  Propulsion Committee RANS/panel workshop (Jessup 1998), Energy Equip. Sys. 6 (2018) Table 2 (agrees).',
        'twist_deg = blade angle beta from the plane of rotation = atan( (P/D) / (pi r/R) ) of the nose-tail-line pitch;',
        '  matches the table PITCHAN column (e.g. 0.7R: 0.45794 rad = 26.24 deg). No global pitch (fixed-pitch propeller).',
        'chord_m = (c/D) * 0.3048; tip row 1.0 has zero chord (table). thickness_ratio = T/CHORD, camber_ratio = FM/CHORD.',
        "airfoil: NACA 66 (TMB mod.) + a=0.8 sections built at 0.2..0.95R; 'transition' rows interpolate between the",
        '  neighbouring named sections; rows 0.975..1.0 use the 0.95R section (t/c 0.032 is nearly constant there).',
        'Extra columns r_R, P_D, camber_ratio are documentation only.'], geo, float_format='%.6f')

    foils, shapes = {}, {}
    for x in FOIL_STATIONS:
        s = g[np.isclose(g.r_R, x)].iloc[0]
        X, Y = naca66mod_a08(s.t_c, s.f_c)
        p = os.path.join(CASE, 'airfoils', fname(x) + '.dat')
        write_selig(p, f'DTMB4119 r/R={x:.3f} NACA66mod t/c={s.t_c:.5f} a=0.8 f/c={s.f_c:.5f}', X, Y)
        foils[fname(x)] = p
        shapes[fname(x)] = (X, Y, s.t_c, s.f_c)

    pol = run_xfoil_set(XFOIL, os.path.join(HERE, 'xfoil_raw'), foils, RES, ncrit=NCRIT, run=do_xfoil)
    for name, xf in pol.items():
        X, Y, t_c, f_c = shapes[name]
        # transitional friction: Jessup reports laminar flow over much of the untripped model blades at Rn ~1e6
        df, kinds = complete_polars(xf, X, Y, t_c, RES, NCRIT, friction='transitional')
        write_csv(os.path.join(CASE, 'airfoils', name + '_polars.csv'), [
            f'{name}: NACA 66 (TMB modified) thickness t/c = {t_c:.5f}, NACA a=0.8 mean line f/c = {f_c:.5f}',
            '  No measured 2D polars of these sections at Re 2e5-1.5e6 exist. Coordinates: airfoils/' + name + '.dat',
            '  (Brockett 1966 ordinates, see extraction/marine_common.py). Script: extraction/build_dtmb4119.py.',
            'XFOIL rows: XFOIL 6.99 (repository binary xfoil_executables/xfoil), viscous, free transition,',
            f'  Ncrit = {NCRIT} (Jessup reports laminar flow over much of the untripped model blades at Rn ~1e6), 200 panels,',
            '  ITER 150-300, sweeps 0 -> +18 / -10 deg (0.5 deg) plus 0 -> +12 / -6 deg (0.25 deg); only converged points.',
            *polar_header_lines(kinds, 'transitional')], df, float_format='%.5f')
        print(name, kinds)

    # Open water, Koyama 1993 Table 1.3.1
    ow = pd.DataFrame({'J': [0.5, 0.7, 0.833, 0.9, 1.1], 'kt': [0.285, 0.200, 0.146, 0.120, 0.034],
                       'kq10': [0.477, 0.360, 0.280, 0.239, 0.106], 'eta_table': [0.489, 0.632, 0.692, 0.725, 0.575]})
    ow['kq'] = ow.kq10 / 10
    ow['ct'] = ow.kt
    ow['cp'] = 2 * np.pi * ow.kq
    ow['eta'] = ow.J * ow.ct / ow.cp
    ow = ow[['J', 'ct', 'cp', 'eta', 'kt', 'kq', 'kq10', 'eta_table']]
    write_csv(os.path.join(CASE, 'experiment', 'performance_open_water.csv'), [
        'DTMB propeller 4119, open-water test (Jessup), DTMB towing basin; uniform inflow, uncavitating.',
        'Source: Koyama (1993) Papers of Ship Research Institute Suppl. 15, Table 1.3.1 "Open water test results for',
        '  DTRC4119" (p. 72), TABULATED (KT, 10KQ, eta0 at 5 J). Fig. 1.3.1 (p. 75) shows the faired curves, Rn = 1e6.',
        'Conversion to the program convention: ct = KT, cp = 2*pi*KQ, eta = J*ct/cp. kt, kq, kq10, eta_table = source values.',
        '  NOTE: the source eta0 column is NOT consistent with its own KT/KQ at J = 0.5, 0.9, 1.1 (0.489/0.725/0.575 vs',
        '  0.475/0.719/0.562 from J KT/(2 pi KQ)); it agrees at J 0.7/0.833 - probably read off the faired curve. Use eta.',
        'Extraction uncertainty: none (table, 3 decimals: +-0.0005 KT, +-0.0005 10KQ). Measurement uncertainty not stated.',
        'Design point (lifting-surface design, not measured): J = 0.833, KT = 0.154, 10KQ = 0.290, eta = 0.706.'],
        ow, float_format='%.5f')
    print(ow)


if __name__ == '__main__':
    main('--xfoil' in sys.argv)

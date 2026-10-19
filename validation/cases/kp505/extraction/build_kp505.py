"""Builds validation/cases/kp505 (KRISO KP505 propeller of the KCS container ship, 5 blades, model D = 250 mm).

Sources (copies in ../sources/):
  Geometry  : KRISO/MOERI 'propgeom.txt' (r/R, P/D, rake, skew, C/D, fo/C, to/D), in the Tokyo 2015 CFD workshop
              archive KCS_prop_MOERI.zip (https://t2015.nmri.go.jp/file/Geometry_IGES_files/kcs/KCS_prop_MOERI.zip,
              'Same with SIMMAN2014 / CFDWS2010'); identical table in SIMMAN 2008 'KCS propgeom NMRI.zip'.
  Open water: (a) MOERI model D = 250 mm: 'openwater.txt' in the same archive (J, Kt, 10Kq, Eta0, 21 points).
              (b) NMRI model D = 104.6 mm: SIMMAN 2008 'KCS_POT_NMRI.zip' -> KCS_POT_NMRI.txt (32 rps, Rn = nD^2/nu =
                  3.68e5, 22.2 C, NMRI No.3 towing tank, 2005-09-12).
Run: python3 build_kp505.py [--xfoil]   (--xfoil runs XFOIL for missing combinations; without it the polars are
     rebuilt from the .pol files in extraction/xfoil_raw)
"""
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
D = 0.250
R = D / 2
# propgeom.txt: r/R, P/D, rake, skew(deg), C/D, fo/C, to/D
GEOM = [
    (0.180, 0.834700, 0.0, -4.720, 0.231300, 0.028448, 0.045850),
    (0.250, 0.891200, 0.0, -6.980, 0.261800, 0.029641, 0.040710),
    (0.300, 0.926900, 0.0, -7.820, 0.280900, 0.029477, 0.037120),
    (0.400, 0.978300, 0.0, -7.740, 0.313800, 0.026769, 0.030470),
    (0.500, 1.007900, 0.0, -5.560, 0.340300, 0.022010, 0.024590),
    (0.600, 1.013000, 0.0, -1.500, 0.357300, 0.017324, 0.019470),
    (0.700, 0.996700, 0.0, 4.110, 0.359000, 0.014039, 0.014920),
    (0.800, 0.956600, 0.0, 10.480, 0.337600, 0.011996, 0.010730),
    (0.900, 0.900600, 0.0, 17.170, 0.279700, 0.010440, 0.006930),
    (0.950, 0.868300, 0.0, 20.630, 0.222500, 0.010067, 0.005280),
    (1.000, 0.833100, 0.0, 24.180, 0.000100, 8.700000, 0.003690),
]
FOIL_STATIONS = [0.18, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
RES = [1e5, 2e5, 4e5, 7e5, 1e6, 1.5e6]
NCRIT = 9
MOERI_RPS = 18.0     # ASSUMPTION: rate of revolution of the MOERI open-water test is not given in the archive
RHO, NU = 998.0, 1.14e-6   # ASSUMPTION fresh water ~15 C for the MOERI test (not stated)
NMRI_NU = 32 * 0.1046 ** 2 / 3.68e5          # 9.51e-7 m^2/s, from the stated Rn and rps (22.2 C)
NMRI_RPM_EQ = 32 * 60 * (0.1046 / D) ** 2    # rpm of the 250 mm geometry giving the NMRI chord Reynolds numbers


def fname(x):
    return f'kp505_r{int(round(x * 1000)):04d}'


def openwater(path_txt, header, rows):
    d = pd.DataFrame(rows, columns=['J', 'kt', 'kq10'])
    d['kq'] = d.kq10 / 10
    d['ct'] = d.kt
    d['cp'] = 2 * np.pi * d.kq
    d['eta'] = np.where(d.J > 0, d.J * d.ct / d.cp, 0.0)
    return d[['J', 'ct', 'cp', 'eta', 'kt', 'kq', 'kq10']]


def main(do_xfoil):
    os.makedirs(os.path.join(CASE, 'airfoils'), exist_ok=True)
    os.makedirs(os.path.join(CASE, 'experiment'), exist_ok=True)
    g = pd.DataFrame(GEOM, columns=['r_R', 'P_D', 'rake', 'skew_deg', 'c_D', 'f_c', 't_D'])
    g['t_c'] = g.t_D / g.c_D
    rows = []
    for _, s in g.iterrows():
        x = s.r_R
        if any(abs(x - fs) < 1e-6 for fs in FOIL_STATIONS):
            foil = fname(x)
        elif x > 0.95:
            foil = fname(0.95)
        else:
            foil = 'transition'
        tc = s.t_c if x < 0.999 else g.t_c.iloc[-2]
        rows.append(dict(r_m=x * R, chord_m=s.c_D * D, twist_deg=float(propeller_beta_deg(s.P_D, x)), airfoil=foil,
                         thickness_ratio=tc, r_R=x, P_D=s.P_D, skew_deg=s.skew_deg, camber_ratio=min(s.f_c, 0.05)))
    geo = pd.DataFrame(rows)
    write_csv(os.path.join(CASE, 'geometry.csv'), [
        'KRISO KP505 (propeller of the KRISO Container Ship KCS) - 5 blades, model D = 0.250 m, hub d/D = 0.180,',
        '  skew 32 deg (balanced, -7.8 deg at 0.3R .. +24.2 deg at the tip), no rake, AE/A0 = 0.800.',
        'Source: KRISO/MOERI propgeom.txt (Tokyo 2015 / SIMMAN 2008 & 2014 / Gothenburg 2010 workshop archives),',
        '  TABULATED, copied unchanged (r/R, P/D, skew, C/D, fo/C, to/D); identical in SIMMAN 2008 KCS propgeom NMRI.zip.',
        'twist_deg = blade angle beta from the plane of rotation = atan( (P/D) / (pi r/R) ); no global pitch (FPP).',
        'chord_m = (C/D) * 0.250 (expanded chord). thickness_ratio = (to/D)/(C/D). Tip row: C/D = 0.0001 and fo/C = 8.7',
        '  in the source (meaningless at zero chord): airfoil and t/c of 0.95R used there, camber_ratio clipped.',
        'Skew is NOT representable in BEM (sections are taken on the cylinder as if unskewed); it is kept as a column.',
        "airfoil: NACA 66 thickness + a=0.8 mean line sections built at the named stations; 'transition' (0.25R) between.",
        'Extra columns r_R, P_D, skew_deg, camber_ratio are documentation only.'], geo, float_format='%.6f')

    foils, shapes = {}, {}
    for x in FOIL_STATIONS:
        s = g[np.isclose(g.r_R, x)].iloc[0]
        X, Y = naca66mod_a08(s.t_c, s.f_c)
        p = os.path.join(CASE, 'airfoils', fname(x) + '.dat')
        write_selig(p, f'KP505 r/R={x:.3f} NACA66mod t/c={s.t_c:.5f} a=0.8 f/c={s.f_c:.5f}', X, Y)
        foils[fname(x)] = p
        shapes[fname(x)] = (X, Y, s.t_c, s.f_c)
    pol = run_xfoil_set(XFOIL, os.path.join(HERE, 'xfoil_raw'), foils, RES, ncrit=NCRIT, run=do_xfoil)
    for name, xf in pol.items():
        X, Y, t_c, f_c = shapes[name]
        df, kinds = complete_polars(xf, X, Y, t_c, RES, NCRIT, friction='transitional')
        write_csv(os.path.join(CASE, 'airfoils', name + '_polars.csv'), [
            f'{name}: NACA 66 thickness t/c = {t_c:.5f} (KRISO: "NACA66 thickness form"; built with the DTMB-modified',
            f'  ordinates of Brockett 1966), NACA a=0.8 mean line f/c = {f_c:.5f}. No measured polars exist.',
            'XFOIL rows: XFOIL 6.99 (repository binary), viscous, free transition, Ncrit = 9, 200 panels, ITER 150-300,',
            '  sweeps 0 -> +18 and 0 -> -10 deg (0.5 deg) plus 0 -> +12 / -6 deg (0.25 deg); only converged points',
            '  (raw output in extraction/xfoil_raw). Script: extraction/build_kp505.py (helpers extraction/marine_common.py).',
            *polar_header_lines(kinds, 'transitional')], df, float_format='%.5f')
        print(name, kinds)

    moeri = [(0.000, 0.518, 0.725), (0.100, 0.482, 0.677), (0.150, 0.458, 0.646), (0.200, 0.435, 0.622),
             (0.250, 0.412, 0.589), (0.300, 0.387, 0.557), (0.350, 0.361, 0.531), (0.400, 0.336, 0.497),
             (0.450, 0.310, 0.466), (0.500, 0.285, 0.437), (0.550, 0.259, 0.405), (0.600, 0.235, 0.376),
             (0.650, 0.209, 0.343), (0.700, 0.185, 0.311), (0.750, 0.161, 0.278), (0.800, 0.137, 0.247),
             (0.850, 0.111, 0.215), (0.900, 0.083, 0.181), (0.950, 0.053, 0.140), (1.000, 0.022, 0.096),
             (1.050, -0.012, 0.046)]
    d = openwater(None, None, moeri)
    d['eta_table'] = [0.000, 0.113, 0.169, 0.223, 0.278, 0.332, 0.379, 0.431, 0.477, 0.519, 0.561, 0.597, 0.631,
                      0.665, 0.691, 0.705, 0.697, 0.654, 0.575, 0.373, -0.437]
    write_csv(os.path.join(CASE, 'experiment', 'performance_moeri_250mm.csv'), [
        'KP505 open water, MOERI (KRISO) model D = 250 mm (scale 31.6), towing tank. TABULATED: openwater.txt in',
        '  KCS_prop_MOERI.zip (Tokyo 2015 workshop, https://t2015.nmri.go.jp/kcs_gc.html; also CFDWS2010/SIMMAN2014).',
        'ct = KT, cp = 2*pi*KQ, eta = J*ct/cp; kt, kq, kq10, eta_table = source values (eta_table agrees to rounding).',
        'Rate of revolution / Reynolds number NOT stated in the archive: defaults rpm = 1080 (18 rps) is an ASSUMPTION',
        '  (Re(0.7R) ~ 8e5 at J 0.7). Extraction uncertainty: none (table, 3 decimals).'], d, float_format='%.5f')

    nmri = [(0.000, 0.5327, 0.7517), (0.100, 0.4937, 0.7058), (0.150, 0.4719, 0.6813), (0.200, 0.4469, 0.6538),
            (0.250, 0.4208, 0.6232), (0.300, 0.3922, 0.5895), (0.350, 0.3657, 0.5589), (0.400, 0.3425, 0.5314),
            (0.450, 0.3143, 0.5008), (0.500, 0.2895, 0.4702), (0.550, 0.2647, 0.4396), (0.600, 0.2407, 0.4090),
            (0.650, 0.2162, 0.3784), (0.700, 0.1931, 0.3478), (0.700, 0.1943, 0.3478), (0.750, 0.1688, 0.3172),
            (0.800, 0.1414, 0.2805), (0.850, 0.1148, 0.2468), (0.900, 0.0870, 0.2132), (0.950, 0.0581, 0.1704),
            (1.000, 0.0293, 0.1275), (1.050, -0.0033, 0.0786)]
    d2 = openwater(None, None, nmri)
    write_csv(os.path.join(CASE, 'experiment', 'performance_nmri_105mm.csv'), [
        'KP505 open water, NMRI model D = 104.6 mm (scale 75.5), NMRI No.3 towing tank, sting dynamometer, 2005-09-12,',
        '  n = 32 rps, Tw = 22.2 C, Rn = n D^2 / nu = 3.68e5, immersion I/D = 1.5. TABULATED: KCS_POT_NMRI.txt',
        '  (SIMMAN 2008, http://www.simman2008.dk/PDF/KCS_POT_NMRI.zip). Source note: "KQ data are only as a guide due',
        '  to the small size propeller". J = 0.7 appears twice in the source (repeat point) - both kept.',
        'ct = KT, cp = 2*pi*KQ, eta = J*ct/cp. kt, kq, kq10 = source values.',
        f'REYNOLDS-EQUIVALENT CONDITIONS for the 250 mm case geometry: defaults rpm = {NMRI_RPM_EQ:.1f} (= 32 rps x',
        f'  (0.1046/0.250)^2) and kin_viscosity = {NMRI_NU:.3e} (from the stated Rn) reproduce the chord Reynolds numbers',
        '  of the 104.6 mm model (Re(0.7R) ~ 3e5); the coefficients are size-independent.'], d2, float_format='%.5f')
    print(d.tail(3), d2.tail(3))


if __name__ == '__main__':
    main('--xfoil' in sys.argv)

"""geometry.csv of the Bahaj et al. (2007) 800 mm tidal rotor from Table 1 (tabulated).

Table 1 of Bahaj et al. 2007 (p. 410) = Table 1 of Batten et al. 2008 (p. 1087; 9 of the 17 stations, identical
values) = Table 1 of Batten et al. 2007 Ocean Eng. Columns: r/R, radius (mm), c/R, pitch distribution (deg), t/c (%).
"""
import os

import pandas as pd

from make_airfoils import TC_TABLE1, foil_name

HERE = os.path.dirname(os.path.abspath(__file__))
R = 0.4
T1 = [  # r/R, r_mm, c/R, pitch_dist_deg, t/c %
    (0.20, 80, 0.1250, 15.0, 24.0), (0.25, 100, 0.1203, 12.1, 22.5), (0.30, 120, 0.1156, 9.5, 20.7),
    (0.35, 140, 0.1109, 7.6, 19.5), (0.40, 160, 0.1063, 6.1, 18.7), (0.45, 180, 0.1016, 4.9, 18.1),
    (0.50, 200, 0.0969, 3.9, 17.6), (0.55, 220, 0.0922, 3.1, 17.1), (0.60, 240, 0.0875, 2.4, 16.6),
    (0.65, 260, 0.0828, 1.9, 16.1), (0.70, 280, 0.0781, 1.5, 15.6), (0.75, 300, 0.0734, 1.2, 15.1),
    (0.80, 320, 0.0688, 0.9, 14.6), (0.85, 340, 0.0641, 0.6, 14.1), (0.90, 360, 0.0594, 0.4, 13.6),
    (0.95, 380, 0.0547, 0.2, 13.1), (1.00, 400, 0.0500, 0.0, 12.6)]
# Batten et al. 2008 Table 1 (cross-check): r/R, c/R, pitch, t/c
B08 = [(0.2, 0.1250, 15.0, 24.0), (0.3, 0.1156, 9.5, 20.7), (0.4, 0.1063, 6.1, 18.7), (0.5, 0.0969, 3.9, 17.6),
       (0.6, 0.0875, 2.4, 16.6), (0.7, 0.0781, 1.5, 15.6), (0.8, 0.0688, 0.9, 14.6), (0.9, 0.0594, 0.4, 13.6),
       (1.0, 0.0500, 0.0, 12.6)]

t1 = {round(a, 2): (c, p, t) for a, _, c, p, t in T1}
for a, c, p, t in B08:
    assert t1[round(a, 2)] == (c, p, t), a
for a, r_mm, *_ in T1:
    assert abs(a * 400 - r_mm) < 1e-9
assert [t for *_, t in T1] == TC_TABLE1

df = pd.DataFrame({
    'r_m': [a * R for a, *_ in T1],
    'chord_m': [round(c * R, 5) for _, _, c, _, _ in T1],
    'twist_deg': [p for _, _, _, p, _ in T1],
    'airfoil': [foil_name(t) for *_, t in T1],
    'thickness_ratio': [t / 100 for *_, t in T1],
})
hdr = """Bahaj et al. (2007) 800 mm marine current turbine, 3 blades, R = 0.400 m - blade geometry
Source: Bahaj, Molland, Chaplin & Batten 2007, Renewable Energy 32:407-426, Table 1 (p. 410), 17 stations,
  TABULATED, copied unchanged (r = r/R * 0.4 m, chord = c/R * 0.4 m, t/c). Cross-checked against Batten et al. 2008,
  Renewable Energy 33:1085-1096, Table 1 (9 stations, identical) and Batten et al. 2007 Ocean Eng. 34:1013, Table 1.
twist_deg = Table 1 'pitch distribution' = chord angle relative to the TIP chord (15.0 deg at 0.2R, 0 at the tip),
  i.e. the solver twist WITHOUT global pitch. The paper's 'hub pitch angle' (20, 25 deg, ...) is the nose-tail angle
  at r = 80 mm (0.2R) from the rotor plane, so the local angle to the rotor plane = twist_deg + (hub_pitch - 15);
  the experiment files carry pitch_deg = hub_pitch - 15 (5 deg for the 20 deg design setting).
  Checked: this solver with pitch_deg = 5 gives alpha = 18/14/10/7.4 deg (TSR 4) and 10/4.3/2.8/2.0 deg (TSR 6) at
  r/R 0.2/0.4/0.6/0.8, as the authors' BEM for 20 deg hub pitch (Batten 2008 Fig. 4a: ~16/14/10.5/7.5, ~7.8/4.5/3/2).
airfoil: NACA 63-8xx at the exact Table 1 t/c of the station (the blades were lofted from 63-812/815/818/821/824
  interpolated to these thicknesses, Bahaj 2007 Sec. 3.1); coordinates generated with naca456 (a = 1.0 mean line,
  cl_i = 0.8), see extraction/make_airfoils.py.
Blade starts at r = 0.08 m (0.2R): inboard of it the section blends to a circular root clamped in the hub (hub/boss
  diameter 100 mm -> Rhub = 0.05 m); that root is not tabulated and is not modelled (as in the authors' BEM, x =
  0.2-1.0). Tip: square tip, chord 0.020 m at r = R (Table 1)."""
with open(os.path.join(HERE, '..', 'geometry.csv'), 'w') as f:
    for L in hdr.split('\n'):
        f.write('# ' + L + '\n')
    df.to_csv(f, index=False, float_format='%.5g')
print(df)

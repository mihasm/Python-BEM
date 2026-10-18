"""Builds the Selig-format airfoil coordinate files of the MEXICO case.

du91_w2_250.dat   : DU 91-W2-250 coordinates (blunt TE, t/c=0.25) from the SoftwareX code repository
                    ElsevierSoftwareX/SOFTX-D-16-00046, AirfoilCoordinates/DU91W2250.dat (already Selig order).
                    Cross-checked against the Risø Wind Turbine Airfoil Catalogue profile (Bertagnolio et al.
                    2001, Risø-R-1280, data dir DU91-W2-250/prof.dat, stored reverse order) below.
naca64_418.dat    : NACA 64(3)-418 from the UIUC database via airfoiltools.com (naca643418-il).
risoe_a1_21.dat   : RISØ-A1-21. The coordinates are proprietary (Risø catalogue: "Airfoil shape: Not
                    available"). RECONSTRUCTED from the 23 pressure-sensor positions (x/c, y/c) of the MEXICO
                    blade at the 60%R section (the RISØ-A1-21 section), Mexnext-III final report ECN-E--18-003,
                    Table D.2 (p. 261), with a periodic smoothing spline per surface; trailing edge closed at
                    x/c=1 with the TE thickness taken as 0.4% (estimate). Shape accuracy ~0.3% c; adequate for
                    thickness and the runner's Montgomerie extrapolation, NOT for XFOIL.
The same Table D.2 sensor coordinates at 25%R (DU) and 82/92%R (NACA) are compared with the public DU/NACA
files as a check of the method; see check_airfoil_shapes.png.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'airfoils')
RAW = os.path.join(HERE, 'raw_airfoil_data')
os.makedirs(OUT, exist_ok=True)


def load_xy(path, skip=1):
    rows = []
    with open(path) as f:
        for line in f.read().splitlines()[skip:]:
            p = line.split()
            if len(p) >= 2:
                try:
                    rows.append((float(p[0]), float(p[1])))
                except ValueError:
                    pass
    return np.array(rows)


def write_selig(path, name, xy, comment):
    with open(path, 'w') as f:
        f.write(name + '\n')
        for x, y in xy:
            f.write('%10.6f %10.6f\n' % (x, y))
    print('wrote', path, len(xy), 'points; t/c = %.4f' % thickness(xy), comment)


def thickness(xy):
    i = np.argmin(xy[:, 0])
    up, lo = xy[:i + 1][::-1], xy[i:]
    xs = np.linspace(0.01, 0.99, 400)
    yu = np.interp(xs, up[:, 0], up[:, 1])
    yl = np.interp(xs, lo[:, 0], lo[:, 1])
    return float(np.max(yu - yl))


# ---- DU 91-W2-250
du = load_xy(os.path.join(RAW, 'DU91W2250_softx.dat'))
write_selig(os.path.join(OUT, 'du91_w2_250.dat'), 'DU 91-W2-250 (SoftwareX SOFTX-D-16-00046 repository)', du, '')
cat = load_xy(os.path.join(RAW, 'du91w2250_prof_catalogue.dat'))   # TE -> lower -> LE -> upper
cat = cat[::-1]                                                     # -> Selig order

# ---- NACA 64(3)-418
na = load_xy(os.path.join(RAW, 'naca643418_airfoiltools.dat'))
write_selig(os.path.join(OUT, 'naca64_418.dat'), 'NACA 64(3)-418 (UIUC / airfoiltools naca643418-il)', na, '')

# ---- sensor coordinates, Mexnext-III Table D.2 (x/c, y/c in %), per section
sens = {
    25: [(62.04, 9.49), (50.66, 11.53), (38.58, 12.79), (29.67, 12.68), (22.49, 11.73), (15.64, 10.10),
         (10.57, 8.36), (5.30, 5.80), (2.54, 3.90), (0.20, 1.16), (0.00, 0.00), (0.19, -1.22), (0.95, -2.71),
         (3.33, -5.13), (8.05, -7.91), (15.73, -10.58), (24.35, -12.01), (40.74, -11.42), (50.59, -9.31),
         (60.90, -6.05), (69.73, -2.96), (80.54, 0.02), (89.65, 0.88)],
    60: [(83.61, 3.24), (73.55, 5.09), (63.29, 7.10), (52.80, 9.45), (42.59, 11.64), (32.54, 12.50),
         (22.74, 11.90), (14.57, 10.21), (6.34, 6.63), (2.04, 3.29), (0.44, 1.31), (0.00, 0.00), (0.32, -0.85),
         (1.67, -1.76), (4.79, -3.04), (13.79, -6.18), (24.55, -8.23), (34.44, -8.41), (42.81, -7.55),
         (51.56, -5.93), (60.62, -3.72), (69.73, -1.54), (78.96, 0.01), (89.79, 0.60)],
    82: [(73.98, 5.83), (51.90, 9.89), (39.90, 10.85), (27.97, 10.28), (17.47, 8.61), (7.66, 5.76), (3.47, 3.80),
         (1.29, 2.28), (0.32, 1.16), (0.00, 0.00), (0.20, -0.89), (1.10, -1.98), (3.24, -3.13), (18.99, -6.35),
         (39.24, -7.08), (49.44, -6.17), (60.01, -4.65), (69.23, -3.12), (80.11, -1.34), (89.53, -0.14)],
}
for k in sens:
    sens[k] = np.array(sens[k]) / 100.0

# ---- RISØ-A1-21 reconstruction from the 60%R sensors
s60 = sens[60]
iL = np.argmin(s60[:, 0] + np.abs(s60[:, 1]))
up = s60[:iL + 1][::-1]          # LE -> TE, upper (y>0)
lo = s60[iL:]                    # LE -> TE, lower
te_half = 0.002
# parametrise each surface by sqrt(x) (clusters near LE) and interpolate monotonically (PCHIP)
def surf(pts, y_te):
    x = np.concatenate([pts[:, 0], [1.0]])
    y = np.concatenate([pts[:, 1], [y_te]])
    s = np.sqrt(x)
    f = PchipInterpolator(s, y)
    sg = np.linspace(0, 1, 81) ** 1.0
    return sg ** 2, f(sg)
xu, yu = surf(up, te_half + 0.0)
xl, yl = surf(lo, -te_half)
ris = np.vstack([np.column_stack([xu, yu])[::-1], np.column_stack([xl, yl])[1:]])
write_selig(os.path.join(OUT, 'risoe_a1_21.dat'),
            'RISO-A1-21 reconstructed from MEXICO 60%R sensor positions (Mexnext-III Table D.2)', ris,
            '(design t/c = 0.21)')

# ---- check plot
fig, axs = plt.subplots(3, 1, figsize=(10, 9))
axs[0].plot(du[:, 0], du[:, 1], 'k-', lw=1, label='DU91-W2-250 SoftX (used)')
axs[0].plot(cat[:, 0], cat[:, 1], 'g--', lw=1, label='Risø catalogue prof.dat')
axs[0].plot(sens[25][:, 0], sens[25][:, 1], 'ro', ms=4, label='MEXICO 25%R sensors (Table D.2)')
axs[1].plot(ris[:, 0], ris[:, 1], 'k-', lw=1, label='RISØ-A1-21 reconstruction (used)')
axs[1].plot(sens[60][:, 0], sens[60][:, 1], 'ro', ms=4, label='MEXICO 60%R sensors (Table D.2)')
axs[2].plot(na[:, 0], na[:, 1], 'k-', lw=1, label='NACA 64(3)-418 (used)')
axs[2].plot(sens[82][:, 0], sens[82][:, 1], 'ro', ms=4, label='MEXICO 82%R sensors (Table D.2)')
for a in axs:
    a.set_aspect('equal'); a.grid(alpha=.3); a.legend(fontsize=8, loc='lower right')
fig.tight_layout()
fig.savefig(os.path.join(HERE, 'check_airfoil_shapes.png'), dpi=130)

"""
Generates the reference data of the case `structural_analytical`:

  experiment/section_properties.csv   area, centroid, second moments of area, extreme-fibre distances of
                                      solid and hollow (constant wall thickness) sections
  experiment/beam_closed_form.csv     cantilever bending moment / shear / bending stress / deflection under
                                      closed-form spanwise loads
  experiment/centrifugal_closed_form.csv  centrifugal force and stress of rotating bars
  airfoils/*.dat                      the exact shapes (Selig format) used, so a test feeds the program the same outline

Independent of the solver code: polygon integrals (Green's theorem) on very dense outlines, the exact inward offset of
the outline (shapely/GEOS buffer) for the hollow sections, closed-form integrals where they exist, and XFOIL's BEND
command as a second, independent implementation.

Run (needs numpy, scipy, pandas, matplotlib, shapely; XFOIL from the repository for the BEND cross-check):
    python3 -m venv /tmp/venv && PIP_USER=0 /tmp/venv/bin/pip install numpy scipy pandas matplotlib shapely
    /tmp/venv/bin/python validation/cases/structural_analytical/extraction/make_reference.py
"""
import math
import os
import re
import subprocess

import numpy as np
import pandas as pd
from scipy import integrate
from shapely.geometry import Polygon
from shapely.geometry.polygon import orient

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.dirname(HERE)
REPO = os.path.abspath(os.path.join(CASE, "..", "..", ".."))
XFOIL = os.path.join(REPO, "xfoil_executables", "xfoil")
pi = math.pi


# ----------------------------------------------------------------------------------------------------------------------
# geometry helpers
# ----------------------------------------------------------------------------------------------------------------------
def ring_integrals(x, y):
    """Signed area, first and second moments of a closed polygon about the ORIGIN (Green's theorem, exact for the
    polygon). Returns A, Sx=∫x dA, Sy=∫y dA, Iyy_o=∫y² dA, Ixx_o=∫x² dA, Ixy_o=∫xy dA."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if x[0] == x[-1] and y[0] == y[-1]:
        x, y = x[:-1], y[:-1]
    xj, yj = np.roll(x, -1), np.roll(y, -1)
    cr = x * yj - xj * y
    A = cr.sum() / 2
    Sx = ((x + xj) * cr).sum() / 6
    Sy = ((y + yj) * cr).sum() / 6
    Iy2 = (cr * (y * y + y * yj + yj * yj)).sum() / 12  # ∫y² dA
    Ix2 = (cr * (x * x + x * xj + xj * xj)).sum() / 12  # ∫x² dA
    Ixy = (cr * (x * yj + 2 * x * y + 2 * xj * yj + xj * y)).sum() / 24  # ∫xy dA
    return np.array([A, Sx, Sy, Iy2, Ix2, Ixy])


def section_props(outer, inners=()):
    """Centroidal properties of a region bounded by `outer` minus `inners` (lists of (x, y) arrays).
    Ix = ∫(y-yc)² dA (about the centroidal axis PARALLEL TO x = the chord: flapwise),
    Iy = ∫(x-xc)² dA (edgewise), Ixy = ∫(x-xc)(y-yc) dA."""
    r = ring_integrals(*outer)
    tot = r * np.sign(r[0])
    for xi, yi in inners:
        r = ring_integrals(xi, yi)
        tot = tot - r * np.sign(r[0])
    A, Sx, Sy, Iy2, Ix2, Ixy = tot
    xc, yc = Sx / A, Sy / A
    allx = np.asarray(outer[0])
    ally = np.asarray(outer[1])
    return dict(A_m2=A, xc_m=xc, yc_m=yc,
                Ix_m4=Iy2 - A * yc ** 2, Iy_m4=Ix2 - A * xc ** 2, Ixy_m4=Ixy - A * xc * yc,
                y_up_m=ally.max() - yc, y_low_m=yc - ally.min(),
                x_le_m=xc - allx.min(), x_te_m=allx.max() - xc)


def naca4(code, n=2001, c=1.0, closed_te=False):
    """NACA 4-digit outline (Abbott & von Doenhoff 1959, §6.4; thickness applied perpendicular to the camber line).
    Returns Selig-ordered x, y (TE -> upper -> LE -> lower -> TE), n points per side, cosine spacing."""
    m = int(code[0]) / 100.0
    p = int(code[1]) / 10.0
    t = int(code[2:]) / 100.0
    beta = np.linspace(0, pi, n)
    x = 0.5 * (1 - np.cos(beta))
    a4 = -0.1036 if closed_te else -0.1015
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 + a4 * x ** 4)
    if m == 0:
        yc = np.zeros_like(x)
        dyc = np.zeros_like(x)
    else:
        yc = np.where(x < p, m / p ** 2 * (2 * p * x - x ** 2), m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x - x ** 2))
        dyc = np.where(x < p, 2 * m / p ** 2 * (p - x), 2 * m / (1 - p) ** 2 * (p - x))
    th = np.arctan(dyc)
    xu, yu = x - yt * np.sin(th), yc + yt * np.cos(th)
    xl, yl = x + yt * np.sin(th), yc - yt * np.cos(th)
    xs = np.concatenate([xu[::-1], xl[1:]]) * c
    ys = np.concatenate([yu[::-1], yl[1:]]) * c
    return xs, ys


def ccw(x, y):
    """Counter-clockwise copy (positive area) of a closed outline."""
    if ring_integrals(x, y)[0] < 0:
        return np.asarray(x)[::-1], np.asarray(y)[::-1]
    return np.asarray(x), np.asarray(y)


def hollow_exact(x, y, wall):
    """Exact inner boundary of a constant-thickness wall: inward offset of the outline by `wall`
    (Minkowski/buffer, GEOS). Returns the inner outline (x, y) (or None when nothing is left)."""
    P = Polygon(list(zip(x, y)))
    H = P.buffer(-wall, quad_segs=256)
    if H.is_empty:
        return None
    if H.geom_type != "Polygon":  # e.g. TE region cut off into a separate piece: keep the biggest
        H = max(H.geoms, key=lambda g: g.area)
    H = orient(H, 1.0)
    xi, yi = H.exterior.xy
    return np.asarray(xi), np.asarray(yi)


def midline_bredt(x, y, wall):
    """Bredt-Batho torsion constant of a thin closed single-cell section: J = 4 A_m² t / S_m,
    A_m and S_m of the wall mid-line (inward offset by wall/2)."""
    m = hollow_exact(x, y, wall / 2)
    if m is None:
        return float("nan")
    P = Polygon(list(zip(*m)))
    return 4 * P.area ** 2 * wall / P.length


def rotate(x, y, theta_deg):
    """Counter-clockwise rotation of the outline by theta about the origin (x right, y up)."""
    t = math.radians(theta_deg)
    x, y = np.asarray(x), np.asarray(y)
    return x * math.cos(t) - y * math.sin(t), x * math.sin(t) + y * math.cos(t)


def write_selig(path, name, x, y, fmt="%.8f"):
    with open(path, "w") as f:
        f.write(name + "\n")
        for a, b in zip(x, y):
            f.write((fmt + " " + fmt + "\n") % (a, b))


def xfoil_bend(code):
    """XFOIL 6.x BEND on the airfoils/naca<code>.dat outline, re-splined and paneled by XFOIL with 300 nodes.
    (Not XFOIL's own NACA generator: that one evidently adds the thickness VERTICALLY to the camber line - its A, Xc and
    Iyy of NACA 4412 equal those of NACA 0012 to 6 digits - while the NACA definition adds it normal to the camber
    line, which gives +0.35 % area and +0.7 % edgewise inertia for NACA 4412.) Returns dict or None."""
    if not os.path.exists(XFOIL):
        return None
    dat = "naca%s.dat" % code  # relative: XFOIL truncates long file names; run inside airfoils/
    cmd = "PLOP\nG F\n\nLOAD %s\nPPAR\nN 300\n\n\nBEND\n\nQUIT\n" % dat
    try:
        out = subprocess.run([XFOIL], input=cmd, capture_output=True, text=True, timeout=60,
                             cwd=os.path.join(CASE, "airfoils")).stdout
    except Exception as e:  # noqa
        print("xfoil failed", e)
        return None
    with open(os.path.join(HERE, "xfoil_bend", "naca%s_bend.txt" % code), "w") as f:
        f.write(out)

    def grab(label, start=0):
        m = re.search(re.escape(label) + r"\s*=\s*([-+0-9.Ee]+)", out[start:])
        return float(m.group(1)) if m else float("nan")

    skin = out.find("(skin)")
    return dict(A=grab("Area"), S=grab("Slen"), Iy=grab("Iyy"), Ix=grab("Ixx"), J=grab("J   "),
                xc=grab("Xc"), Iy_t=grab("Iyy/t", skin), Ix_t=grab("Ixx/t", skin), J_t=grab("J/t", skin))


# ----------------------------------------------------------------------------------------------------------------------
# 1. section properties
# ----------------------------------------------------------------------------------------------------------------------
rows = []
os.makedirs(os.path.join(HERE, "xfoil_bend"), exist_ok=True)
SRC_NACA = ("NACA 4-digit outline (Abbott & von Doenhoff 1959, Theory of Wing Sections, §6.4, a4=-0.1015, open TE "
            "closed by a straight segment). Properties: polygon integrals on 2x40001 cosine-spaced points "
            "(converged, see conv_rel); generated by extraction/make_reference.py")


def add(row):
    rows.append(row)


naca_list = ["0006", "0009", "0012", "0015", "0018", "0021", "0024", "2412", "4412", "4415", "4418", "6412"]
for code in naca_list:
    t_c = int(code[2:]) / 100
    camber = int(code[0]) / 100
    x, y = ccw(*naca4(code, n=40001))
    p = section_props((x, y))
    x2, y2 = ccw(*naca4(code, n=20001))
    p2 = section_props((x2, y2))
    conv = max(abs(p[k] / p2[k] - 1) for k in ("A_m2", "Ix_m4", "Iy_m4"))
    # the file the program can load (Selig, 161 points per side)
    xf, yf = naca4(code, n=161)
    write_selig(os.path.join(CASE, "airfoils", "naca%s.dat" % code), "NACA %s" % code, xf, yf)
    pf = section_props(ccw(xf, yf))
    xb = xfoil_bend(code)
    row = dict(case_id="naca%s_solid" % code, shape="naca%s" % code, airfoil_file="airfoils/naca%s.dat" % code,
               construction="solid", c_m=1.0, t_c=t_c, camber_c=camber, wall_m=0.0, theta_deg=0.0, **p,
               K_A=p["A_m2"] / t_c, K_I=p["Ix_m4"] / (t_c * (t_c ** 2 + camber ** 2)),
               J_m4=float("nan"),
               A_datfile_rel_err=pf["A_m2"] / p["A_m2"] - 1, Ix_datfile_rel_err=pf["Ix_m4"] / p["Ix_m4"] - 1,
               Iy_datfile_rel_err=pf["Iy_m4"] / p["Iy_m4"] - 1,
               conv_rel=conv,
               xfoil_A_m2=xb["A"] if xb else np.nan, xfoil_Ix_m4=xb["Ix"] if xb else np.nan,
               xfoil_Iy_m4=xb["Iy"] if xb else np.nan,
               method="exact (polygon integral, dense outline)", source=SRC_NACA, tolerance_rel=1e-3)
    if camber == 0:
        # closed-form check of the symmetric section: A = 2∫yt dx, Ix = 2/3∫yt³ dx (substitution x = u²)
        f = lambda u: 5 * t_c * (0.2969 * u - 0.1260 * u ** 2 - 0.3516 * u ** 4 + 0.2843 * u ** 6 - 0.1015 * u ** 8)
        A_q = 2 * integrate.quad(lambda u: f(u) * 2 * u, 0, 1, epsabs=1e-14, epsrel=1e-13)[0]
        Ix_q = 2 / 3 * integrate.quad(lambda u: f(u) ** 3 * 2 * u, 0, 1, epsabs=1e-16, epsrel=1e-13)[0]
        row["closed_form_A_m2"] = A_q
        row["closed_form_Ix_m4"] = Ix_q
    add(row)
    # Drela rule of thumb (MIT 16.01 notes): A ≈ 0.60 c t, I ≈ 0.036 c t (t² + h²)
    add(dict(case_id="naca%s_drela_rule" % code, shape="naca%s" % code, airfoil_file="airfoils/naca%s.dat" % code,
             construction="solid", c_m=1.0, t_c=t_c, camber_c=camber, wall_m=0.0, theta_deg=0.0,
             A_m2=0.60 * t_c, Ix_m4=0.036 * t_c * (t_c ** 2 + camber ** 2), K_A=0.60, K_I=0.036,
             method="approximation A=K_A c t, I=K_I c t (t^2+h^2)",
             source="Drela, M. 'Area and Bending Inertia of Airfoil Sections', MIT 16.01 Unified Engineering notes "
                    "(OCW Fall 2005), eqs. (6),(7),(10),(11): K_A~0.60, K_I~0.036 'for most common airfoils'",
             tolerance_rel=0.15))
    # hollow, constant wall thickness (exact inward offset) for several walls
    if code in ("0012", "0018", "2412", "4412", "4415"):
        for wall in (0.002, 0.005, 0.010, 0.020):
            inner = hollow_exact(x, y, wall)
            ph = section_props((x, y), [inner])
            J = midline_bredt(x, y, wall)
            r = dict(case_id="naca%s_hollow_w%g" % (code, wall), shape="naca%s" % code,
                     airfoil_file="airfoils/naca%s.dat" % code, construction="hollow", c_m=1.0, t_c=t_c,
                     camber_c=camber, wall_m=wall, theta_deg=0.0, **ph, J_m4=J,
                     method="exact: outline minus its inward offset by wall (GEOS buffer, 256 segments per quarter "
                            "circle); J by Bredt-Batho on the wall mid-line",
                     source="geometry as naca*_solid; Bredt-Batho J = 4 A_m^2 t / S_m (e.g. Megson, Aircraft "
                            "Structures for Engineering Students, ch. 18)",
                     tolerance_rel=2e-3)
            if xb:  # thin-wall estimate from XFOIL's skin integrals (outer contour, valid for wall << t)
                r["thinwall_A_m2"] = xb["S"] * wall
                r["thinwall_Ix_m4"] = xb["Ix_t"] * wall
                r["thinwall_Iy_m4"] = xb["Iy_t"] * wall
            add(r)
    # rotated (pitched/twisted) section: properties in the rotated frame, about the centroid
    if code in ("4412",):
        for theta in (10.0, 30.0, 90.0):
            xr, yr = rotate(x, y, theta)
            add(dict(case_id="naca%s_solid_rot%g" % (code, theta), shape="naca%s" % code,
                     airfoil_file="airfoils/naca%s.dat" % code, construction="solid", c_m=1.0, t_c=t_c,
                     camber_c=camber, wall_m=0.0, theta_deg=theta, **section_props((xr, yr)),
                     method="exact (polygon integral), outline rotated counter-clockwise by theta about the origin "
                            "(x from LE to TE, y to the upper side); Ix/Iy about the rotated frame's axes",
                     source=SRC_NACA, tolerance_rel=1e-3))
            inner = hollow_exact(x, y, 0.005)
            xr2, yr2 = rotate(*inner, theta)
            add(dict(case_id="naca%s_hollow_w0.005_rot%g" % (code, theta), shape="naca%s" % code,
                     airfoil_file="airfoils/naca%s.dat" % code, construction="hollow", c_m=1.0, t_c=t_c,
                     camber_c=camber, wall_m=0.005, theta_deg=theta, **section_props((xr, yr), [(xr2, yr2)]),
                     method="as naca4412_hollow_w0.005, rotated counter-clockwise by theta",
                     source=SRC_NACA, tolerance_rel=2e-3))

# small-blade scale: the program's default wall (0.001 m) on a 0.1 m chord
for code in ("0012", "4412"):
    c = 0.1
    x, y = ccw(*naca4(code, n=40001, c=c))
    for wall in (0.0, 0.001, 0.002):
        t_c = int(code[2:]) / 100
        if wall == 0:
            p = section_props((x, y))
            cons = "solid"
        else:
            p = section_props((x, y), [hollow_exact(x, y, wall)])
            cons = "hollow"
        add(dict(case_id="naca%s_c0.1_%s%s" % (code, cons, "" if wall == 0 else "_w%g" % wall), shape="naca%s" % code,
                 airfoil_file="airfoils/naca%s.dat" % code, construction=cons, c_m=c, t_c=t_c,
                 camber_c=int(code[0]) / 100, wall_m=wall, theta_deg=0.0, **p,
                 J_m4=midline_bredt(x, y, wall) if wall else float("nan"),
                 method="exact (polygon integral / exact inward offset)", source=SRC_NACA,
                 tolerance_rel=2e-3))

# circle (the program's cylindrical root, airfoil 'circle') and tube
Rc = 0.05
th = np.linspace(0, 2 * pi, 4001)
xcir, ycir = Rc * np.cos(th) + Rc, Rc * np.sin(th)  # LE at x=0 like an airfoil of chord 2R
write_selig(os.path.join(CASE, "airfoils", "circle.dat"), "circle (unit chord)",
            *(lambda t: (0.5 + 0.5 * np.cos(t), 0.5 * np.sin(t)))(np.linspace(0, 2 * pi, 201)))
add(dict(case_id="circle_solid", shape="circle", airfoil_file="airfoils/circle.dat", construction="solid",
         c_m=2 * Rc, t_c=1.0, camber_c=0, wall_m=0.0, theta_deg=0.0, A_m2=pi * Rc ** 2, xc_m=Rc, yc_m=0.0,
         Ix_m4=pi * Rc ** 4 / 4, Iy_m4=pi * Rc ** 4 / 4, Ixy_m4=0.0, y_up_m=Rc, y_low_m=Rc, x_le_m=Rc, x_te_m=Rc,
         J_m4=pi * Rc ** 4 / 2, method="closed form A=pi R^2, I=pi R^4/4", source="textbook (e.g. Roark's Formulas, "
         "7th ed., Table A.1 case 15)", tolerance_rel=1e-3))
for wall in (0.002, 0.005, 0.010):
    Ri = Rc - wall
    add(dict(case_id="circle_tube_w%g" % wall, shape="circle", airfoil_file="airfoils/circle.dat",
             construction="hollow", c_m=2 * Rc, t_c=1.0, camber_c=0, wall_m=wall, theta_deg=0.0,
             A_m2=pi * (Rc ** 2 - Ri ** 2), xc_m=Rc, yc_m=0.0, Ix_m4=pi * (Rc ** 4 - Ri ** 4) / 4,
             Iy_m4=pi * (Rc ** 4 - Ri ** 4) / 4, Ixy_m4=0.0, y_up_m=Rc, y_low_m=Rc, x_le_m=Rc, x_te_m=Rc,
             J_m4=pi * (Rc ** 4 - Ri ** 4) / 2, method="closed form annulus",
             source="textbook (Roark's Formulas, 7th ed., Table A.1 case 16)", tolerance_rel=1e-3))
# Chen, Yu & Capellaro (2010) Table I: aluminium tube R=0.3 m (outer), E=73 GPa, rho=2800 kg/m3, elasticity solution
E_al, rho_al = 73e9, 2800.0
chen_tab1 = {"1/15": (5.137e9, 2.024e8, 1.970e2), "1/7.5": (9.540e9, 3.301e8, 3.659e2),
             "1/5": (1.321e10, 4.042e8, 5.067e2), "2/7.5": (1.615e10, 4.424e8, 6.193e2),
             "1/3": (1.835e10, 4.587e8, 7.037e2)}
for key, (EA, EI, mu) in chen_tab1.items():
    ratio = eval(key)
    wall = ratio * 0.6
    Ri = 0.3 - wall
    add(dict(case_id="chen2010_tube_t2R_%s" % key.replace("/", "_"), shape="circle",
             airfoil_file="airfoils/circle.dat", construction="hollow", c_m=0.6, t_c=1.0, camber_c=0,
             wall_m=wall, theta_deg=0.0, A_m2=EA / E_al, xc_m=0.3, yc_m=0.0, Ix_m4=EI / E_al, Iy_m4=EI / E_al,
             Ixy_m4=0.0, y_up_m=0.3, y_low_m=0.3, x_le_m=0.3, x_te_m=0.3,
             closed_form_A_m2=pi * (0.09 - Ri ** 2), closed_form_Ix_m4=pi * (0.3 ** 4 - Ri ** 4) / 4,
             mass_per_length_kg_m=mu, rho_kg_m3=rho_al, E_pa=E_al,
             method="tabulated: A = EA/E, I = EI/E, Elasticity-theory row (mu = rho A)",
             source="Chen, Yu & Capellaro 2010, Wind Energy 13:497-516, doi:10.1002/we.372, Table I (p. 505); "
                    "R = 0.3 m outer radius, t/2R as in case_id", tolerance_rel=2e-3))

# ellipse (a = half chord, b = half thickness), solid: closed form
a_e, b_e = 0.5, 0.06
te = np.linspace(0, 2 * pi, 201)
write_selig(os.path.join(CASE, "airfoils", "ellipse12.dat"), "ellipse t/c=0.12", a_e + a_e * np.cos(te),
            b_e * np.sin(te))
add(dict(case_id="ellipse12_solid", shape="ellipse12", airfoil_file="airfoils/ellipse12.dat", construction="solid",
         c_m=1.0, t_c=0.12, camber_c=0, wall_m=0, theta_deg=0, A_m2=pi * a_e * b_e, xc_m=a_e, yc_m=0,
         Ix_m4=pi * a_e * b_e ** 3 / 4, Iy_m4=pi * a_e ** 3 * b_e / 4, Ixy_m4=0, y_up_m=b_e, y_low_m=b_e,
         x_le_m=a_e, x_te_m=a_e, J_m4=pi * a_e ** 3 * b_e ** 3 / (a_e ** 2 + b_e ** 2),
         method="closed form ellipse (J: Saint-Venant solution)", source="Roark's Formulas, 7th ed., Table A.1 "
         "case 17 and Table 10.7 case 2", tolerance_rel=1e-3))
xe, ye = ccw(a_e + a_e * np.cos(np.linspace(0, 2 * pi, 40001)), b_e * np.sin(np.linspace(0, 2 * pi, 40001)))
for wall in (0.002, 0.005):
    add(dict(case_id="ellipse12_hollow_w%g" % wall, shape="ellipse12", airfoil_file="airfoils/ellipse12.dat",
             construction="hollow", c_m=1.0, t_c=0.12, camber_c=0, wall_m=wall, theta_deg=0,
             **section_props((xe, ye), [hollow_exact(xe, ye, wall)]), J_m4=midline_bredt(xe, ye, wall),
             method="exact inward offset (GEOS buffer); the inner boundary of a constant-thickness ellipse wall is "
                    "NOT an ellipse", source="geometry closed form; this script", tolerance_rel=2e-3))

# Chen et al. (2010) isotropic blade-like section (Fig. 6, Table IV): leading-edge arc R=12.7 mm centred at the
# origin, straight edges tangent to the arc meeting at x2 = 38.1 mm, wall 1.524 mm, E=206.843 GPa, rho=1068.69 kg/m3
Rb, Lb, wb = 12.7e-3, 38.1e-3, 1.524e-3
alpha = math.asin(Rb / Lb)
ths = np.linspace(pi / 2 - alpha, 3 * pi / 2 + alpha, 20001)
xb_ = np.append(Rb * np.cos(ths), Lb)
yb_ = np.append(Rb * np.sin(ths), 0.0)
xb_, yb_ = ccw(xb_, yb_)
inner = hollow_exact(xb_, yb_, wb)
pb = section_props((xb_, yb_), [inner])
# Selig file of the outline: LE at x = 0, chord = R + L. The straight edges are sampled densely too (80 points each):
# an algorithm that offsets the outline POINTS (like the program's generate_hollow_foil) needs points along them.
ths_f = np.linspace(pi / 2 - alpha, 3 * pi / 2 + alpha, 121)
xa, ya = Rb * np.cos(ths_f), Rb * np.sin(ths_f)
s_edge = np.linspace(0, 1, 81)[:-1]
xf = np.concatenate([Lb + (xa[0] - Lb) * s_edge, xa, xa[-1] + (Lb - xa[-1]) * s_edge[1:], [Lb]]) + Rb
yf = np.concatenate([(ya[0]) * s_edge, ya, ya[-1] + (0 - ya[-1]) * s_edge[1:], [0.0]])
cb = Rb + Lb
write_selig(os.path.join(CASE, "airfoils", "chen2010_bladelike.dat"),
            "Chen et al. 2010 blade-like section (unit chord; real chord 50.8 mm)", xf / cb, yf / cb)
Eb, rhob = 206.843e9, 1068.69
tab4 = dict(VABS=(2.101e3, 1.050e4, 3.566e7, -3.046e5, 1.843e-1, 1.085e-5, 4.080e-5, 9.516e-3),
            Analytical=(2.101e3, 1.110e4, 3.567e7, -3.379e5, 1.843e-1, 1.085e-5, 4.081e-5, 9.513e-3),
            ANSYS=(2.101e3, 1.051e4, 3.567e7, -3.051e5, 1.843e-1, 1.085e-5, 4.081e-5, 9.513e-3),
            PreComp=(2.178e3, 9.100e3, 3.794e7, -3.238e-2, 1.960e-1, 1.125e-5, 4.702e-5, 1.000e-2))
pb_shift = dict(pb)
pb_shift["xc_m"] = pb["xc_m"] + Rb  # x measured from the leading edge in the row
add(dict(case_id="chen2010_bladelike_hollow_exact", shape="chen2010_bladelike",
         airfoil_file="airfoils/chen2010_bladelike.dat", construction="hollow", c_m=cb, t_c=2 * Rb / cb,
         camber_c=0, wall_m=wb, theta_deg=0, **pb_shift, J_m4=midline_bredt(xb_, yb_, wb),
         method="exact geometry of Chen et al. Fig. 6 (outline minus inward offset), this script",
         source="Chen, Yu & Capellaro 2010, Fig. 6 (p. 509); agrees with their Table IV to <0.02% (rows below)",
         tolerance_rel=2e-3))
for tool, (EI22, EI33, EA, S14, mu, i22, i33, xm2) in tab4.items():
    add(dict(case_id="chen2010_bladelike_table4_%s" % tool, shape="chen2010_bladelike",
             airfoil_file="airfoils/chen2010_bladelike.dat", construction="hollow", c_m=cb, t_c=2 * Rb / cb,
             camber_c=0, wall_m=wb, theta_deg=0, A_m2=EA / Eb, xc_m=xm2 + Rb, yc_m=0.0, Ix_m4=EI22 / Eb,
             Iy_m4=(EI33 - S14 ** 2 / EA) / Eb, Ixy_m4=0.0, mass_per_length_kg_m=mu, rho_kg_m3=rhob, E_pa=Eb,
             Ix_from_mass_m4=i22 / rhob, Iy_from_mass_m4=i33 / rhob,
             method="tabulated (%s column): A=EA/E, Ix=EI22/E, Iy=(EI33 - S14^2/EA)/E moved from the shear centre "
                    "to the centroid; Ix/Iy_from_mass = i22/rho, i33/rho (mass moments, about the mass centre)" % tool,
             source="Chen, Yu & Capellaro 2010, Table IV (p. 510)",
             tolerance_rel=0.002 if tool != "PreComp" else float("nan")))

sec = pd.DataFrame(rows)
cols = ["case_id", "shape", "airfoil_file", "construction", "c_m", "t_c", "camber_c", "wall_m", "theta_deg",
        "A_m2", "xc_m", "yc_m", "Ix_m4", "Iy_m4", "Ixy_m4", "y_up_m", "y_low_m", "x_le_m", "x_te_m", "J_m4",
        "K_A", "K_I", "tolerance_rel", "method", "source",
        "closed_form_A_m2", "closed_form_Ix_m4", "xfoil_A_m2", "xfoil_Ix_m4", "xfoil_Iy_m4",
        "thinwall_A_m2", "thinwall_Ix_m4", "thinwall_Iy_m4",
        "A_datfile_rel_err", "Ix_datfile_rel_err", "Iy_datfile_rel_err", "conv_rel",
        "mass_per_length_kg_m", "rho_kg_m3", "E_pa", "Ix_from_mass_m4", "Iy_from_mass_m4"]
sec = sec[cols]
HEADER_SEC = """# Reference section properties (kind: structural) - case structural_analytical
# Generated by extraction/make_reference.py (independent of the solver code). Tabulated literature values are rows whose
# `method` starts with "tabulated"; everything else is exact geometry/closed form (no digitising anywhere).
# Frame: x along the chord from the leading edge to the trailing edge, y normal to the chord towards the upper (suction)
#   side; the outline is NOT rotated unless theta_deg != 0 (then rotated counter-clockwise by theta about the origin).
#   Any other frame follows from the theta=0 row: for a frame turned by +theta (x'= x cos t + y sin t, y'= -x sin t + y cos t
#   i.e. the section rotated CLOCKWISE) Ix' = Ix cos^2 t + Iy sin^2 t - 2 Ixy sin t cos t,
#   Iy' = Iy cos^2 t + Ix sin^2 t + 2 Ixy sin t cos t, Ixy' = Ixy cos 2t + (Ix - Iy) sin t cos t (Ix=int y^2, Iy=int x^2,
#   Ixy=int xy); the rows rotated counter-clockwise here are the same formulas with t = -theta. A reflection (x -> -x or
#   y -> -y) only flips the sign of Ixy.
# Columns (SI units):
#   c_m chord; t_c thickness ratio; camber_c max camber/chord; wall_m constant wall thickness of a hollow section
#     (inward normal offset of the outline; 0 = solid); construction solid|hollow
#   A_m2 area; xc_m, yc_m centroid (x from the leading edge); Ix_m4 = int (y-yc)^2 dA (centroidal axis parallel to the
#     chord: FLAPWISE, the program's Ix at theta 0); Iy_m4 = int (x-xc)^2 dA (EDGEWISE, the program's Iy);
#     Ixy_m4 = int (x-xc)(y-yc) dA
#   y_up_m / y_low_m: distance from the centroid to the outermost upper / lower fibre; x_le_m / x_te_m: to the leading /
#     trailing edge (the program's norm_dist / tang_dist for the extreme-fibre bending stress)
#   J_m4 torsion constant (closed form, or Bredt-Batho 4 A_m^2 t / S_m for hollow sections) - for a future torsion check
#   K_A = A/(c^2 t/c), K_I = Ix/(c^4 tau (tau^2+eps^2)) - Drela's coefficients evaluated exactly
#   tolerance_rel: suggested pass band for a program result against A/Ix/Iy of the row
#   closed_form_*: independent closed-form integral (symmetric NACA: A=2 int yt dx, Ix=2/3 int yt^3 dx; tube: annulus)
#   xfoil_*: XFOIL 6.x BEND on XFOIL's own NACA outline (300 nodes) - an independent implementation (extraction/xfoil_bend/)
#   thinwall_*: thin-wall approximation from XFOIL's skin integrals (outer contour length x wall); valid for wall << t
#   *_datfile_rel_err: error of the same property computed from the coarse airfoils/<shape>.dat file (161 pts/side)
#     instead of the dense outline - the floor any program result on that file can reach
#   conv_rel: change between 2x20001 and 2x40001 outline points (numerical convergence of the reference)
#   mass_per_length_kg_m, rho_kg_m3, E_pa, Ix/Iy_from_mass_m4: the published stiffness/mass values these rows are from
"""
with open(os.path.join(CASE, "experiment", "section_properties.csv"), "w") as f:
    f.write(HEADER_SEC)
    sec.to_csv(f, index=False, float_format="%.8g")

# ----------------------------------------------------------------------------------------------------------------------
# 2. cantilever closed forms (bending moment, shear, stress, deflection, slope)
# ----------------------------------------------------------------------------------------------------------------------
r0, R = 0.1, 1.0
L = R - r0
E, rho = 70e9, 2700.0  # aluminium (assumed, typical 6061/6082: E 69-70 GPa, rho 2700 kg/m3)
c = 0.1
x0012, y0012 = ccw(*naca4("0012", n=40001, c=c))
p0012 = section_props((x0012, y0012))
I = p0012["Ix_m4"]
ymax = max(p0012["y_up_m"], p0012["y_low_m"])
EI = E * I
rs = np.round(np.linspace(r0, R, 19), 6)
beam = []
q, qt, P = 100.0, 200.0, 50.0
for r in rs:
    s = r - r0  # distance from the clamp
    u = R - r  # distance to the tip
    cases = {
        "B1_uniform_q": dict(desc="uniform load q=%g N/m on r0..R (normal to rotor plane)" % q,
                             V=q * u, M=q * u ** 2 / 2,
                             w=q * s ** 2 * (6 * L ** 2 - 4 * L * s + s ** 2) / (24 * EI),
                             slope=q * (L ** 3 - (L - s) ** 3) / (6 * EI)),
        "B2_linear_q": dict(desc="load rising linearly from 0 at r0 to qt=%g N/m at R" % qt,
                            V=qt / (2 * L) * (L ** 2 - s ** 2),
                            M=qt / (6 * L) * (2 * L ** 3 - 3 * L ** 2 * s + s ** 3),
                            w=qt * s ** 2 * (20 * L ** 3 - 10 * L ** 2 * s + s ** 3) / (120 * L * EI),
                            slope=qt * s * (8 * L ** 3 - 6 * L ** 2 * s + s ** 3) / (24 * L * EI)),
        "B3_tip_load": dict(desc="point load P=%g N at the tip R" % P, V=P, M=P * u,
                            w=P * s ** 2 * (3 * L - s) / (6 * EI), slope=P * s * (2 * L - s) / (2 * EI)),
    }
    for k, d in cases.items():
        beam.append(dict(case_id=k, load_case=d["desc"], r_m=r, r0_m=r0, R_m=R, shear_flap_n=d["V"], moment_flap_nm=d["M"],
                         stress_flap_pa=d["M"] * ymax / I, deflection_flap_m=d["w"], slope_flap_rad=d["slope"],
                         E_pa=E, Ix_m4=I, y_extreme_m=ymax))
beam = pd.DataFrame(beam)
# numeric self-check of the closed forms by direct integration of q(r)
for k, qf in (("B1_uniform_q", lambda s: q + 0 * s), ("B2_linear_q", lambda s: qt * s / L)):
    sub = beam[beam.case_id == k]
    ss = np.linspace(0, L, 200001)
    Mnum = np.array([np.trapezoid(qf(ss[ss >= s]) * (ss[ss >= s] - s), ss[ss >= s]) for s in sub.r_m - r0])
    assert np.allclose(Mnum, sub.moment_flap_nm, rtol=1e-6, atol=1e-6), k
# deflection/slope: integrate M/EI twice from the clamp (fine grid) and compare with the closed forms
ss = np.linspace(0, L, 200001)
for k, Mf in (("B1_uniform_q", lambda s: q * (L - s) ** 2 / 2),
              ("B2_linear_q", lambda s: qt / (6 * L) * (2 * L ** 3 - 3 * L ** 2 * s + s ** 3)),
              ("B3_tip_load", lambda s: P * (L - s))):
    sl = integrate.cumulative_trapezoid(Mf(ss) / EI, ss, initial=0)
    ww = integrate.cumulative_trapezoid(sl, ss, initial=0)
    sub = beam[beam.case_id == k]
    assert np.allclose(np.interp(sub.r_m - r0, ss, sl), sub.slope_flap_rad, rtol=1e-6, atol=1e-12), k
    assert np.allclose(np.interp(sub.r_m - r0, ss, ww), sub.deflection_flap_m, rtol=1e-6, atol=1e-12), k
HEADER_BEAM = """# Cantilever closed forms (kind: structural) - case structural_analytical
# Euler-Bernoulli cantilever clamped at r0_m, free at R_m, constant section: NACA 0012 solid, chord 0.1 m (Ix_m4 from
# section_properties.csv row naca0012_c0.1_solid), E = 70 GPa (aluminium, assumed). Loads act normal to the rotor plane.
# Formulas: Roark's Formulas for Stress and Strain (Young & Budynas, 7th ed. 2002), Table 8.1 cases 1a (tip load),
#   2a (uniform load) and 2a with linearly increasing load (w_a=0, w_l=qt); re-derived and checked numerically in
#   extraction/make_reference.py. s = r - r0, L = R - r0.
#   B1 uniform q:   V = q(R-r), M = q(R-r)^2/2, w = q s^2 (6L^2 - 4Ls + s^2)/(24EI), tip w = qL^4/(8EI)
#   B2 linear q(s)=qt s/L: M = qt/(6L)(2L^3 - 3L^2 s + s^3), w = qt s^2(20L^3 - 10L^2 s + s^3)/(120 L EI), tip 11qtL^4/(120EI)
#   B3 tip load P:  M = P(R-r), w = P s^2 (3L - s)/(6EI), tip PL^3/(3EI)
# Columns (structural column dictionary, see validation/structural_sources.md): shear_flap_n shear force, moment_flap_nm
#   bending moment out of the rotor plane about the section at r_m from all load outboard (the program's Ms_n),
#   stress_flap_pa = M*y_extreme/Ix (extreme-fibre flapwise bending stress, the program's stress_norm, here in Pa),
#   deflection_flap_m, slope_flap_rad (out-of-plane deflection and slope; for a future deflection solver).
# Note for comparing with the program: it sums discrete section forces dFn_j*(r_j - r_i); with a load given per element
#   the discrete sum differs from the integral by O(dr) - use fine sections or apply the load as the exact element force.
"""
with open(os.path.join(CASE, "experiment", "beam_closed_form.csv"), "w") as f:
    f.write(HEADER_BEAM)
    beam.to_csv(f, index=False, float_format="%.8g")

# ----------------------------------------------------------------------------------------------------------------------
# 3. centrifugal force and stress of rotating bars
# ----------------------------------------------------------------------------------------------------------------------
rpm = 600.0
om = rpm * 2 * pi / 60
K_A0012 = sec.loc[sec.case_id == "naca0012_solid", "A_m2"].iloc[0]  # A/c^2 of NACA 0012 (t/c 0.12)
cr, ct = 0.15, 0.05
cen = []
for r in rs:
    # C1: uniform cross-section
    A1 = K_A0012 * 0.1 ** 2
    F1 = rho * A1 * om ** 2 * (R ** 2 - r ** 2) / 2
    # C2: NACA 0012 solid, chord linear from cr at r0 to ct at R -> A(r) = K_A c(r)^2 (quadratic in r)
    cpoly = np.polynomial.Polynomial([cr - (ct - cr) / L * r0, (ct - cr) / L])
    Apoly = K_A0012 * cpoly ** 2
    integrand = (Apoly * np.polynomial.Polynomial([0, 1])).integ()
    F2 = rho * om ** 2 * (integrand(R) - integrand(r))
    mass2 = rho * (Apoly.integ()(R) - Apoly.integ()(r))
    cen.append(dict(case_id="C1_uniform_bar", description="uniform NACA0012 solid c=0.1 m, rho=%g, %g rpm" % (rho, rpm),
                    r_m=r, r0_m=r0, R_m=R, rpm=rpm, rho_kg_m3=rho, A_m2=A1, chord_m=0.1,
                    centrifugal_force_n=F1, stress_centrifugal_pa=F1 / A1, mass_outboard_kg=rho * A1 * (R - r)))
    cen.append(dict(case_id="C2_tapered_naca0012", description="NACA0012 solid, chord %g->%g m linear, rho=%g, %g rpm"
                    % (cr, ct, rho, rpm), r_m=r, r0_m=r0, R_m=R, rpm=rpm, rho_kg_m3=rho, A_m2=Apoly(r),
                    chord_m=cpoly(r), centrifugal_force_n=F2, stress_centrifugal_pa=F2 / Apoly(r),
                    mass_outboard_kg=mass2))
cen = pd.DataFrame(cen)
HEADER_CEN = """# Centrifugal force and stress of rotating bars (kind: structural) - case structural_analytical
# Exact integrals F(r) = rho omega^2 int_r^R A(xi) xi dxi, sigma(r) = F(r)/A(r) (axis of rotation at r = 0, blade along r,
#   no coning/pitch effect; the bar's own centrifugal load only). Uniform bar: sigma = rho omega^2 (R^2 - r^2)/2,
#   independent of A (standard result, e.g. Timoshenko & Goodier / any rotating-bar derivation; re-derived here).
# C2: solid NACA 0012 (A = K_A c^2 with K_A = 0.68504*0.12 from section_properties.csv), chord linear root->tip;
#   integrals of the polynomial A(r)*r done exactly (numpy.polynomial). rho = 2700 kg/m3 (aluminium, assumed).
# Columns: r_m station; A_m2 area there; centrifugal_force_n = axial force at r from the mass outboard (the program's dFc);
#   stress_centrifugal_pa (the program's stress_cent, in Pa); mass_outboard_kg = mass between r and R (r=r0 row: blade mass
#   from r0).
# Note for comparing with the program: it lumps each element's mass at its centre radius - a midpoint rule, error O(dr^2).
"""
with open(os.path.join(CASE, "experiment", "centrifugal_closed_form.csv"), "w") as f:
    f.write(HEADER_CEN)
    cen.to_csv(f, index=False, float_format="%.8g")
print(sec[["case_id", "A_m2", "Ix_m4", "Iy_m4", "xfoil_A_m2", "xfoil_Ix_m4", "xfoil_Iy_m4", "closed_form_Ix_m4"]]
      .to_string())

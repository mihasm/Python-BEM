"""Shared helpers for the marine-propeller validation cases (DTMB 4119, KP505, Wageningen B4-70).

A copy of this file lives in each case's extraction/ folder so every case is self-contained.

- naca66mod_a08(): NACA 66 (DTMB/TMB modified nose and tail) thickness + NACA a=0.8 mean line, built from the
  ordinates tabulated by Brockett (1966), DTMB Report 1780, Tables 1 and 3; leading-edge radius 0.448 T^2.
  Thickness laid off perpendicular to the mean line (as Brockett states for this section).
- run_xfoil_set(): XFOIL polars with the repository binary (xfoil_executables/xfoil, XFOIL 6.99, x86_64 via
  Rosetta), graphics off, sweeps from 0 deg (first pass 0->+18/0->-10 deg step 0.5, ITER 150; second pass for
  combinations with < 20 converged points 0->+12/0->-6 deg step 0.25, ITER 300), only converged points kept.
- complete_polars(): per Reynolds number the XFOIL polar where its converged points cover -4..+8 deg, else a
  parametric (XROTOR-type) polar from thin-airfoil theory and flat-plate friction (parametric_polar()). XFOIL does
  not converge on the thin (t/c 2-5 %) outer sections of these blades beyond a few degrees, free or tripped.
- write_selig(), write_polars(), write_csv(): file writers in the validation/README.md formats.
"""
import glob
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

# Brockett (1966) Table 3 (conventional stations) + Table 1 stations: x, Y/T (half thickness / total thickness)
_B_X = [0, .005, .0075, .0125, .025, .05, .075, .1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6, .65, .7, .75,
        .8, .85, .9, .95, .975, 1.0]
_B_YT = [0, .0665, .0812, .1044, .1466, .2066, .2525, .2907, .3521, .4000, .4363, .4637, .4832, .4952, .5000,
         .4962, .4846, .4653, .4383, .4035, .3612, .3110, .2532, .1877, .1143, .0748, .0333]
# a = 0.8 mean line, YC/F (Brockett Table 3)
_B_YC = [0, .0423, .0595, .0907, .1586, .2712, .3657, .4482, .5869, .6993, .7905, .8635, .9202, .9615, .9881,
         1.0, .9971, .9786, .9434, .8892, .8121, .7027, .5425, .3586, .1713, .0823, 0]


def cosine_x(n=90):
    b = np.linspace(0, np.pi, n)
    return 0.5 * (1 - np.cos(b))


def naca66mod_a08(t_c, f_c, n=90):
    """Returns Selig-ordered (x, y): upper surface TE->LE, lower surface LE->TE. t_c = T/c, f_c = camber/c."""
    x = cosine_x(n)
    s_nodes = np.sqrt(np.array(_B_X))
    yt_n = PchipInterpolator(s_nodes, _B_YT)(np.sqrt(x))
    # near the LE use the tabulated nose radius r_LE = 0.448 T^2  ->  Y/T = sqrt(0.896 x)  (matches .0665 at x=.005)
    nose = x < 0.005
    yt_n[nose] = np.sqrt(0.896 * x[nose])
    yt = t_c * yt_n
    yc_f = PchipInterpolator(_B_X, _B_YC)
    yc = f_c * yc_f(x)
    dyc = f_c * yc_f.derivative()(x)
    th = np.arctan(dyc)
    xu, yu = x - yt * np.sin(th), yc + yt * np.cos(th)
    xl, yl = x + yt * np.sin(th), yc - yt * np.cos(th)
    X = np.concatenate([xu[::-1], xl[1:]])
    Y = np.concatenate([yu[::-1], yl[1:]])
    return X, Y


def write_selig(path, name, X, Y):
    with open(path, 'w') as f:
        f.write(name + '\n')
        for a, b in zip(X, Y):
            f.write(f'{a:.6f} {b:.6f}\n')


def _read_pol(path):
    rows, started = [], False
    for L in open(path):
        if L.strip().startswith('------'):
            started = True
            continue
        if started and L.strip():
            p = L.split()
            try:
                rows.append([float(v) for v in p[:5]])
            except ValueError:
                pass
    return np.array(rows).reshape(-1, 5)


def _run_one(xfoil, raw, foil_dat, tag, re, ncrit, a0, a1, da):
    pol = os.path.join(raw, tag + '.pol')
    if os.path.exists(pol):
        os.remove(pol)
    dat = os.path.relpath(foil_dat, raw)
    cmds = '\n'.join(['PLOP', 'G F', '', f'LOAD {dat}', 'PPAR', 'N 200', '', '', 'OPER', 'VPAR', f'N {ncrit}', '',
                      f'VISC {re:.0f}', 'ITER 300', 'PACC', os.path.basename(pol), '', f'ASEQ {a0} {a1} {da}',
                      'PACC', '', 'QUIT', ''])
    try:
        r = subprocess.run([xfoil], input=cmds, text=True, capture_output=True, cwd=raw, timeout=900)
        log = r.stdout
    except subprocess.TimeoutExpired as e:
        log = f'TIMEOUT\n{e.stdout}'
    with open(os.path.join(raw, tag + '.log'), 'w') as f:
        f.write(log if isinstance(log, str) else str(log))
    return tag


def run_xfoil_set(xfoil, raw, foils, res, ncrit=9, sweeps=((0.0, 12.0, 0.25), (0.0, -6.0, -0.25)), first_tag=2,
                  min_points=20, workers=10, run=True):
    """foils: {name: dat_path}. Returns {name: DataFrame(re, alpha_deg, cl, cd, cm, ncrit)}.

    Runs the given sweeps (tagged s<first_tag>, s<first_tag+1>, ...) only for (foil, Re) combinations that have
    fewer than min_points converged points in the .pol files already in raw/ (an earlier pass used 0 -> +18 /
    0 -> -10 deg in 0.5 deg steps, ITER 150, tags s0/s1, which did not converge on the thin outer sections).
    All .pol files of a combination are merged (later passes win at equal alpha). run=False only reads raw/."""
    os.makedirs(raw, exist_ok=True)

    def merged(name, re):
        files = sorted(glob.glob(os.path.join(raw, f'{name}_re{int(re)}_n{ncrit}_s*.pol')), reverse=True)
        parts = [_read_pol(f) for f in files]
        d = np.vstack(parts) if parts else np.zeros((0, 5))
        if len(d):
            _, idx = np.unique(np.round(d[:, 0], 3), return_index=True)
            d = d[idx]
        return d

    jobs = []
    for name, dat in foils.items():
        for re in res:
            if len(merged(name, re)) >= min_points:
                continue
            for k, (a0, a1, da) in enumerate(sweeps):
                jobs.append((xfoil, raw, dat, f'{name}_re{int(re)}_n{ncrit}_s{first_tag + k}', re, ncrit, a0, a1, da))
    if not run:
        jobs = []
    print(f'XFOIL: {len(jobs)} runs', flush=True)
    with ThreadPoolExecutor(max_workers=workers) as ex:
        list(ex.map(lambda j: _run_one(*j), jobs))
    out = {}
    for name in foils:
        frames = []
        for re in res:
            d = merged(name, re)
            frames.append(pd.DataFrame({'re': float(re), 'alpha_deg': d[:, 0], 'cl': d[:, 1], 'cd': d[:, 2],
                                        'cm': d[:, 4], 'ncrit': ncrit}))
        out[name] = pd.concat(frames, ignore_index=True)
    return out


def _mean_line(X, Y, n=400):
    """Mean line z(x) of Selig coordinates (chord along x, LE at the smallest x), on theta = 0..pi."""
    X, Y = np.asarray(X, float), np.asarray(Y, float)
    i = int(np.argmin(X))
    xu, yu = X[:i + 1][::-1], Y[:i + 1][::-1]
    xl, yl = X[i:], Y[i:]
    th = np.linspace(0.0, np.pi, n)
    x = 0.5 * (1 - np.cos(th))
    z = 0.5 * (np.interp(x, xu, yu) + np.interp(x, xl, yl))
    return th, x, z


def thin_airfoil(X, Y):
    """Thin-airfoil theory of the section's mean line: zero-lift angle [deg] and cm about c/4."""
    th, x, z = _mean_line(X, Y)
    dz = np.gradient(z, x, edge_order=2)
    dz[[0, -1]] = dz[[1, -2]]
    a_l0 = -np.trapezoid(dz * (np.cos(th) - 1), th) / np.pi
    a1 = 2 / np.pi * np.trapezoid(dz * np.cos(th), th)
    a2 = 2 / np.pi * np.trapezoid(dz * np.cos(2 * th), th)
    return float(np.degrees(a_l0)), float(np.pi / 4 * (a2 - a1))


def skin_friction(re, friction):
    """Flat-plate Cf: 'turbulent' = ITTC-57 line; 'transitional' = Blasius below Re 5e5, Prandtl-Schlichting
    transitional law (0.074 Re^-0.2 - 1742/Re, transition at Re 5e5) above."""
    if friction == 'turbulent':
        return 0.075 / (np.log10(re) - 2) ** 2
    return 1.328 / np.sqrt(re) if re < 5e5 else 0.074 * re ** -0.2 - 1742.0 / re


POLAR_ALPHA = np.arange(-10.0, 14.01, 0.5)
CD_CL2 = 0.02          # drag due to lift, dCd/d(Cl - Cl_0)^2 (XROTOR-type quadratic bucket)
VISCOUS_SLOPE = 0.90   # viscous lift-slope factor on 2 pi (1 + 0.77 t/c); XFOIL gives 0.095-0.11 /deg at Re >= 7e5


def _smin(a, b, k=0.04):
    return -k * np.logaddexp(-a / k, -b / k)


def parametric_polar(X, Y, t_c, re, friction='transitional'):
    """XROTOR-type polar of a thin section, -10..+14 deg (the runner extends it to +-180 deg with Montgomerie).

    Lift: thin-airfoil zero-lift angle of the mean line, slope 0.90 * 2 pi (1 + 0.77 t/c) per rad, limited smoothly
    at cl_max = cl_0 + c(t/c) (Re/3e6)^0.1 and cl_min = cl_0 - c(t/c) (Re/3e6)^0.1, c = 0.60/0.75/0.85/1.0/1.1/1.2 at
    t/c 0/0.04/0.06/0.09/0.12/0.21 (estimate after the smooth NACA 6-series maxima of Abbott & von Doenhoff 1959),
    falling by 0.04 per deg past stall. Drag: 2 Cf (1 + 2 t/c + 60 (t/c)^4) (Hoerner form factor) + 0.02 (cl - cl_0)^2
    in the attached range, + 0.012 (deg past stall)^1.5 beyond it."""
    a_l0, cm = thin_airfoil(X, Y)
    slope = VISCOUS_SLOPE * 2 * np.pi * (1 + 0.77 * t_c) * np.pi / 180
    a = POLAR_ALPHA
    lin = slope * (a - a_l0)
    cl0 = slope * (0.0 - a_l0)
    dcl = float(np.interp(t_c, [0.0, 0.04, 0.06, 0.09, 0.12, 0.21], [0.60, 0.75, 0.85, 1.0, 1.1, 1.2]))
    dcl *= min(1.0, (re / 3e6) ** 0.1)
    cl_max, cl_min = cl0 + dcl, cl0 - dcl
    a_s, a_n = a_l0 + cl_max / slope, a_l0 + cl_min / slope
    past_s, past_n = np.clip(a - a_s, 0, None), np.clip(a_n - a, 0, None)
    cl = _smin(lin, cl_max - 0.04 * past_s)
    cl = -_smin(-cl, -(cl_min + 0.04 * past_n))
    cd0 = 2 * skin_friction(re, friction) * (1 + 2 * t_c + 60 * t_c ** 4)
    cd = cd0 + CD_CL2 * (np.clip(lin, cl_min, cl_max) - cl0) ** 2 + 0.012 * (past_s ** 1.5 + past_n ** 1.5)
    return pd.DataFrame({'re': float(re), 'alpha_deg': a, 'cl': cl, 'cd': cd, 'cm': cm})


def xfoil_covers(d, need=(-4.0, 8.0), max_gap=1.0):
    """True when the converged XFOIL points cover the `need` range without a gap wider than max_gap deg."""
    a = np.sort(d['alpha_deg'].to_numpy(float))
    if len(a) < 2 or a[0] > need[0] or a[-1] < need[1]:
        return False
    inside = a[(a >= need[0] - max_gap) & (a <= need[1] + max_gap)]
    return bool(np.all(np.diff(inside) <= max_gap + 1e-9))


def complete_polars(xfoil_df, X, Y, t_c, res, ncrit, friction='transitional'):
    """Per Re the XFOIL polar when it covers -4..+8 deg, else parametric_polar(). Returns (DataFrame, {re: kind})."""
    frames, kinds = [], {}
    for re in res:
        d = xfoil_df[np.isclose(xfoil_df['re'], re)]
        if xfoil_covers(d):
            d = d[['re', 'alpha_deg', 'cl', 'cd', 'cm']].copy()
            d['ncrit'] = ncrit
            d['source'] = f'XFOIL 6.99 Ncrit={ncrit} free transition'
            kinds[re] = 'xfoil'
        else:
            d = parametric_polar(X, Y, t_c, re, friction)
            d['ncrit'] = np.nan
            d['source'] = f'parametric (thin-airfoil lift, {friction} friction drag)'
            kinds[re] = 'parametric'
        frames.append(d)
    return pd.concat(frames, ignore_index=True), kinds


def polar_header_lines(kinds, friction):
    n_x = sum(1 for k in kinds.values() if k == 'xfoil')
    return [
        'Per Reynolds number (column source): XFOIL where its converged points cover -4..+8 deg without gaps',
        f'  ({n_x} of {len(kinds)} here), else a PARAMETRIC polar (XROTOR-type), because XFOIL does not converge on these',
        '  thin sections beyond a few degrees (free or tripped transition): lift from thin-airfoil theory of the mean line',
        '  (zero-lift angle; slope 0.90 * 2 pi (1 + 0.77 t/c) /rad), smooth stall at cl_0 +- c(t/c) (Re/3e6)^0.1 (estimate',
        '  after the smooth NACA 6-series maxima, Abbott & von Doenhoff 1959); drag 2 Cf (1 + 2 t/c + 60 (t/c)^4) +',
        f'  0.02 (cl - cl_0)^2 + post-stall rise, Cf {friction} flat plate (extraction/marine_common.py parametric_polar).',
        '  Parametric rows have no ncrit; cm there is the thin-airfoil value.']


def write_csv(path, header_lines, df, float_format='%.6g'):
    with open(path, 'w') as f:
        for h in header_lines:
            f.write('# ' + h + '\n')
        df.to_csv(f, index=False, float_format=float_format)


def propeller_beta_deg(p_d, r_R):
    """Blade angle from the plane of rotation for pitch ratio P/D at radius fraction r/R."""
    return np.degrees(np.arctan(np.asarray(p_d) / (np.pi * np.asarray(r_R))))

"""Shared helpers for the marine-propeller validation cases (DTMB 4119, KP505, Wageningen B4-70).

A copy of this file lives in each case's extraction/ folder so every case is self-contained.

- naca66mod_a08(): NACA 66 (DTMB/TMB modified nose and tail) thickness + NACA a=0.8 mean line, built from the
  ordinates tabulated by Brockett (1966), DTMB Report 1780, Tables 1 and 3; leading-edge radius 0.448 T^2.
  Thickness laid off perpendicular to the mean line (as Brockett states for this section).
- run_xfoil_set(): XFOIL polars with the repository binary (xfoil_executables/xfoil, XFOIL 6.99, x86_64 via
  Rosetta), graphics off, sweeps from 0 deg (first pass 0->+18/0->-10 deg step 0.5, ITER 150; second pass for
  combinations with < 20 converged points 0->+12/0->-6 deg step 0.25, ITER 300), only converged points kept.
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
                  min_points=20, workers=10):
    """foils: {name: dat_path}. Returns {name: DataFrame(re, alpha_deg, cl, cd, cm, ncrit)}.

    Runs the given sweeps (tagged s<first_tag>, s<first_tag+1>, ...) only for (foil, Re) combinations that have
    fewer than min_points converged points in the .pol files already in raw/ (an earlier pass used 0 -> +18 /
    0 -> -10 deg in 0.5 deg steps, ITER 150, tags s0/s1, which did not converge on the thin outer sections).
    All .pol files of a combination are merged (later passes win at equal alpha)."""
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


def write_csv(path, header_lines, df, float_format='%.6g'):
    with open(path, 'w') as f:
        for h in header_lines:
            f.write('# ' + h + '\n')
        df.to_csv(f, index=False, float_format=float_format)


def propeller_beta_deg(p_d, r_R):
    """Blade angle from the plane of rotation for pitch ratio P/D at radius fraction r/R."""
    return np.degrees(np.arctan(np.asarray(p_d) / (np.pi * np.asarray(r_R))))

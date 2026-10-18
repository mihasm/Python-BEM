"""XFOIL polars for the NACA 63-8xx sections of the Bahaj et al. (2007) 800 mm tidal rotor.

No measured polars exist at the rotor's chord Reynolds numbers (~0.5e5 - 4e5, see case.json notes): Abbott & von
Doenhoff give 63-series data at Re >= 3e6 and Molland et al. (2004) measured a 63-815 at Re = 0.8e6 only. The
authors themselves used XFOIL (Batten et al. 2008, Fig. 2), so XFOIL is used here too.

Binary: the repository's xfoil_executables/xfoil (XFOIL 6.97, x86_64, run via Rosetta), non-interactive:
  PLOP / G F (graphics off), LOAD <selig file>, PPAR N 240 (repanel, 240 panels), OPER, VPAR N <ncrit>,
  VISC <Re>, ITER 200, PACC <file>, ASEQ 0 -> +22 deg (step 0.5); a separate run ASEQ 0 -> -12 deg.
  Two sweeps from 0 so that each starts from a converged attached solution.
Ncrit 9 (XFOIL default = quiet tunnel; what Batten et al. used unless stated otherwise - they do not state it) is the
main polar set; Ncrit 5 (Tu ~ 0.3-0.5 %, Mack's correlation, a more turbulent water tunnel) is written as a
sensitivity set.

Output: ../airfoils/naca63-8xx_polars.csv (ncrit 9) and ../airfoils/naca63-8xx_polars_ncrit5.csv
Raw XFOIL polar files: extraction/xfoil_raw/ (console logs archived in sources/xfoil_run_logs.tar.xz)
"""
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))
XFOIL = os.path.join(REPO, 'xfoil_executables', 'xfoil')
AF = os.path.join(HERE, '..', 'airfoils')
RAW = os.path.join(HERE, 'xfoil_raw')
os.makedirs(RAW, exist_ok=True)

from make_airfoils import TC_TABLE1, foil_name  # noqa: E402
FOILS = [foil_name(tc) for tc in TC_TABLE1]   # one section per Table 1 station (exact t/c)
RES = [50000, 100000, 150000, 200000, 300000, 500000]
NCRITS = [9, 5]
SWEEPS = {'pos': (0.0, 22.0, 0.5), 'neg': (0.0, -12.0, -0.5)}
# retry sweeps, run only for (foil, Re, Ncrit) whose main sweeps left gaps on the 0.5-deg grid:
# a finer ascending sweep from a different start, a descending sweep from deep stall, and an ascending one from -12
RETRY = {'pos2': (0.25, 22.0, 0.25), 'down': (22.0, 0.0, -0.5), 'up': (-12.0, 0.0, 0.5)}
ALL_SWEEPS = {**SWEEPS, **RETRY}
GRID = np.arange(-12.0, 22.0001, 0.5)
PREFER = ['pos', 'neg', 'pos2', 'up', 'down']   # which sweep's value is kept when several converged


def run(foil, re, ncrit, sweep):
    a0, a1, da = ALL_SWEEPS[sweep]
    tag = f'{foil}_re{re}_n{ncrit}_{sweep}'
    pol = os.path.join(RAW, tag + '.pol')
    if os.path.exists(pol):
        os.remove(pol)
    # XFOIL wants short relative paths: work inside RAW
    dat = os.path.relpath(os.path.join(AF, foil + '.dat'), RAW)
    cmds = '\n'.join([
        'PLOP', 'G F', '',
        f'LOAD {dat}',
        'PPAR', 'N 240', '', '',
        'OPER', 'VPAR', f'N {ncrit}', '',
        f'VISC {re}', 'ITER 200',
        'PACC', os.path.basename(pol), '',
        f'ASEQ {a0} {a1} {da}',
        'PACC', '', 'QUIT', ''])
    try:
        r = subprocess.run([XFOIL], input=cmds, text=True, capture_output=True, cwd=RAW, timeout=600)
        log = r.stdout
    except subprocess.TimeoutExpired as e:
        log = f'TIMEOUT\n{e.stdout}'
    with open(os.path.join(RAW, tag + '.log'), 'w') as f:
        f.write(log)
    return tag


def read_pol(path):
    rows, started = [], False
    for L in open(path):
        if L.strip().startswith('------'):
            started = True
            continue
        if started and L.strip():
            p = L.split()
            rows.append([float(v) for v in p[:5]])  # alpha CL CD CDp CM
    return np.array(rows).reshape(-1, 5)


def main():
    jobs = [(f, re, n, s) for f in FOILS for re in RES for n in NCRITS for s in SWEEPS]
    # runs already done (polar file present) are skipped; per-run XFOIL console logs are written beside the polar
    # files but were archived to sources/xfoil_run_logs.tar.xz (git-ignored, 300 MB uncompressed)
    jobs = [j for j in jobs if not os.path.exists(os.path.join(RAW, f'{j[0]}_re{j[1]}_n{j[2]}_{j[3]}.pol'))]
    with ThreadPoolExecutor(max_workers=8) as ex:
        list(ex.map(lambda j: run(*j), jobs))
    # second pass: retry sweeps where the 0.5-deg grid has gaps
    retry = []
    for f in FOILS:
        for re in RES:
            for n in NCRITS:
                got = set()
                for s in SWEEPS:
                    p = os.path.join(RAW, f'{f}_re{re}_n{n}_{s}.pol')
                    if os.path.exists(p):
                        got |= set(np.round(read_pol(p)[:, 0], 2))
                if len(set(np.round(GRID, 2)) - got) > 0:
                    retry += [(f, re, n, s) for s in RETRY
                              if not os.path.exists(os.path.join(RAW, f'{f}_re{re}_n{n}_{s}.pol'))]
    print(len(retry), 'retry runs')
    with ThreadPoolExecutor(max_workers=8) as ex:
        list(ex.map(lambda j: run(*j), retry))
    for ncrit, suffix in [(9, ''), (5, '_ncrit5')]:
        for foil in FOILS:
            frames, dropped = [], []
            for re in RES:
                best = {}
                for s in reversed(PREFER):           # later (preferred) sweeps overwrite earlier ones
                    p = os.path.join(RAW, f'{foil}_re{re}_n{ncrit}_{s}.pol')
                    if os.path.exists(p):
                        for row in read_pol(p):
                            a = round(row[0], 2)
                            if abs(a * 2 - round(a * 2)) < 1e-6:    # keep the 0.5-deg grid only
                                best[a] = row
                d = np.array([best[a] for a in sorted(best)]).reshape(-1, 5)
                if len(d) == 0 or d[:, 0].max() < 12.0 or d[:, 0].min() > -5.0:
                    print(f'  DROPPED {foil} Re {re} Ncrit {ncrit}: converged only '
                          f'{d[:, 0].min() if len(d) else "-"}..{d[:, 0].max() if len(d) else "-"} deg')
                    dropped.append(f'Re {re}')
                    continue
                frames.append(pd.DataFrame({'re': re, 'alpha_deg': d[:, 0], 'cl': d[:, 1], 'cd': d[:, 2],
                                            'cm': d[:, 4], 'ncrit': ncrit}))
            df = pd.concat(frames, ignore_index=True)
            df['source'] = f'XFOIL 6.97 Ncrit={ncrit}, 240 panels, free transition'
            tc = int(foil[-3:]) / 1000
            path = os.path.join(AF, f'{foil}_polars{suffix}.csv')
            with open(path, 'w') as f:
                f.write(f'# {foil}: NACA 63-8xx at t/c = {tc:.3f} (63-series thickness, a=1.0 mean line, cl_i = 0.8) - COMPUTED polar\n')
                f.write(f'# XFOIL 6.97 (repository binary xfoil_executables/xfoil), viscous, free transition, Ncrit = {ncrit},\n')
                f.write('#   240 panels (PPAR N 240), ITER 200; sweeps from 0 deg to +22 and to -12 deg (step 0.5); where these left gaps,\n')
                f.write('#   retry sweeps (0.25->22 step 0.25, 22->0 descending, -12->0 ascending) fill them; 0.5-deg grid kept;\n')
                f.write('#   only converged points are kept (missing alphas = XFOIL did not converge).\n')
                f.write('# Coordinates: airfoils/' + foil + '.dat (naca456 generator, see extraction/make_airfoils.py).\n')
                f.write('# Re range chosen for the rotor: chord Re ~0.5e5 (root, low TSR) to ~4e5 (tip, TSR 11 at 1.73 m/s).\n')
                f.write('# No measured polars exist at these Re (Abbott & von Doenhoff: Re >= 3e6; Molland et al. 2004: 63-815 at 0.8e6).\n')
                f.write('# XFOIL under-predicts drag near and beyond stall (Molland et al. 2004, Fig. 6) - expect optimistic post-stall cl/cd.\n')
                if dropped:
                    f.write(f'# Omitted (XFOIL did not converge beyond a few degrees, set unusable): {", ".join(dropped)}.\n')
                f.write('# Re = 50000 sets are dominated by laminar separation bubbles and are noisy/unreliable in XFOIL.\n')
                f.write('# Script: extraction/run_xfoil.py; raw XFOIL output in extraction/xfoil_raw/.\n')
                df.to_csv(f, index=False, float_format='%.5f')
            cnt = df.groupby('re').alpha_deg.agg(['count', 'min', 'max'])
            print(foil, f'ncrit {ncrit}', ' | '.join(f'{int(r)}: {int(c["count"])} pts {c["min"]:.1f}..{c["max"]:.1f}'
                                                      for r, c in cnt.iterrows()))


if __name__ == '__main__':
    main()

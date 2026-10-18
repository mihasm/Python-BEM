"""XFOIL polars for the MEXICO airfoils (DU 91-W2-250 and NACA 64(3)-418).

XFOIL 6.97 binary of the repository (xfoil_executables/xfoil, x86_64 via Rosetta), graphics off,
PANE (160 panels), ITER 200, Ncrit = 9. Two transition settings:
  'trip'  : forced transition at x/c = 0.05 on both sides (XTR 0.05 0.05), mimicking the zig-zag tape at 5%
            chord used on the MEXICO blades in 2006 and on the inboard part (r/R<0.7) in New MEXICO 2014;
  'free'  : free transition, Ncrit 9 (clean outboard NACA part in New MEXICO configuration 3).
Each Re is run as two sweeps starting from 0 deg (up to +20, down to -10) so the boundary layer
solution continues smoothly; non-converged points are simply absent.

Output: extraction/xfoil/<airfoil>_<mode>_Re<re>.pol and merged CSVs (airfoils/ via build_case_files.py).
The .log files were truncated afterwards (head 400 + tail 200 lines): past stall XFOIL diverged to NaN and printed
~400 MB per run; NaN points are never written to the .pol files.
"""
import os, subprocess, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))
XFOIL = os.path.join(REPO, 'xfoil_executables', 'xfoil')
WORK = os.path.join(HERE, 'xfoil')
AIRF = os.path.join(HERE, '..', 'airfoils')
os.makedirs(WORK, exist_ok=True)

RES = [3e5, 5e5, 7e5, 1e6]
NCRIT = 9


def run(datfile, re, mode, a0, a1, da, tag):
    polname = '%s_%s_Re%d_%s.pol' % (os.path.basename(datfile)[:-4], mode, re, tag)
    pol = os.path.join(WORK, polname)
    if os.path.exists(pol):
        os.remove(pol)
    cmds = ['PLOP', 'G F', '', 'LOAD %s' % datfile, 'PPAR', 'N 160', '', '', 'PANE', 'OPER',
            'VISC %g' % re, 'ITER 200', 'VPAR', 'N %d' % NCRIT]
    if mode == 'trip':
        cmds += ['XTR 0.05 0.05']
    cmds += ['', 'PACC', polname, '', 'ASEQ %g %g %g' % (a0, a1, da), 'PACC', '', 'QUIT']
    p = subprocess.run([XFOIL], input='\n'.join(cmds) + '\n', capture_output=True, text=True, timeout=600,
                       cwd=WORK)
    with open(pol[:-4] + '.log', 'w') as f:
        f.write(p.stdout)
        f.write(p.stderr)
    return read_pol(pol)


def read_pol(pol):
    rows = []
    if not os.path.exists(pol):
        return pd.DataFrame(columns=['alpha_deg', 'cl', 'cd', 'cm'])
    started = False
    for line in open(pol):
        if line.strip().startswith('-----'):
            started = True
            continue
        if started:
            p = line.split()
            if len(p) >= 5:
                rows.append((float(p[0]), float(p[1]), float(p[2]), float(p[4])))
    return pd.DataFrame(rows, columns=['alpha_deg', 'cl', 'cd', 'cm'])


def polar(datfile, mode):
    out = []
    for re in RES:
        up = run(datfile, re, mode, 0, 20, 0.5, 'up')
        dn = run(datfile, re, mode, -0.5, -10, -0.5, 'dn')
        d = pd.concat([dn, up]).drop_duplicates('alpha_deg').sort_values('alpha_deg')
        d.insert(0, 're', re)
        d['ncrit'] = NCRIT
        d['source'] = 'XFOIL 6.97, Ncrit=9, %s' % ('forced transition x/c=0.05 both sides' if mode == 'trip'
                                                    else 'free transition')
        print(os.path.basename(datfile), mode, re, len(d), 'points, alpha %.1f..%.1f' %
              (d.alpha_deg.min(), d.alpha_deg.max()) if len(d) else 'NO POINTS', file=sys.stderr)
        out.append(d)
    return pd.concat(out, ignore_index=True)


if __name__ == '__main__':
    results = {}
    for name in ('du91_w2_250', 'naca64_418'):
        for mode in ('trip', 'free'):
            dat = os.path.join(WORK, name + '.dat')
            # XFOIL wants a short path without spaces: copy the coordinate file into the work dir
            with open(os.path.join(AIRF, name + '.dat')) as f, open(dat, 'w') as g:
                g.write(f.read())
            results[(name, mode)] = polar(name + '.dat', mode)
            results[(name, mode)].to_csv(os.path.join(WORK, '%s_%s_all.csv' % (name, mode)), index=False)

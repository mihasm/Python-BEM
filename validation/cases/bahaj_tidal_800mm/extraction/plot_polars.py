"""Plot the XFOIL polars (Ncrit 9 solid, Ncrit 5 dashed) for a visual check -> check_polars.png"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
AF = os.path.join(HERE, '..', 'airfoils')
from make_airfoils import TC_TABLE1, foil_name
foils = [foil_name(t) for t in (24.0, 19.5, 16.6, 14.6, 12.6)]   # 5 of the 17 sections
fig, axs = plt.subplots(3, 5, figsize=(22, 12))
cols = plt.cm.viridis
for j, f in enumerate(foils):
    for suf, ls in [('', '-'), ('_ncrit5', '--')]:
        d = pd.read_csv(os.path.join(AF, f'{f}_polars{suf}.csv'), comment='#')
        res = sorted(d.re.unique())
        for k, re in enumerate(res):
            g = d[d.re == re]
            c = cols(k / (len(res) - 1))
            lab = f'Re {re:.0f}' if suf == '' else None
            axs[0, j].plot(g.alpha_deg, g.cl, ls, color=c, lw=1, marker='.', ms=2, label=lab)
            axs[1, j].plot(g.alpha_deg, g.cd, ls, color=c, lw=1, marker='.', ms=2)
            axs[2, j].plot(g.alpha_deg, g.cl / g.cd, ls, color=c, lw=1)
    axs[0, j].set_title(f + ' (solid Ncrit 9, dashed Ncrit 5)')
    axs[1, j].set_ylim(0, 0.2)
    for i in range(3):
        axs[i, j].grid(alpha=.3)
    axs[2, j].set_xlabel('alpha [deg]')
axs[0, 0].legend(fontsize=7)
axs[0, 0].set_ylabel('cl'); axs[1, 0].set_ylabel('cd'); axs[2, 0].set_ylabel('cl/cd')
fig.tight_layout()
fig.savefig(os.path.join(HERE, 'check_polars.png'), dpi=80)

#!/usr/bin/env python3
"""XFOIL (Ncrit 9 and 5) vs measured UIUC LSAT Vol. 3 Clark-Y (B): check of the XFOIL setup used for all polars."""
import sys, glob
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, pandas as pd
cache = Path(sys.argv[1] if len(sys.argv) > 1 else "/tmp/uiuc/build/xfoil_cache")
out = Path(sys.argv[2] if len(sys.argv) > 2 else "clarky_xfoil_vs_lsat.png")
meas = pd.read_csv(Path(__file__).resolve().parents[1] / "airfoils/clarky_uiuc_lsat_polars.csv", comment="#") \
    if len(sys.argv) < 4 else pd.read_csv(sys.argv[3], comment="#")
fig, axs = plt.subplots(1, 2, figsize=(12, 5))
cols = {60e3: "C0", 100e3: "C1", 200e3: "C2"}
for re_, c in cols.items():
    m = meas[meas.re == re_]
    axs[0].plot(m.cd, m.cl, "s", c=c, label="UIUC LSAT meas. Re %gk" % (re_ / 1e3))
    axs[1].plot(m.alpha_deg, m.cl, "s", c=c)
    for nc, ls in ((9, "-"), (5, "--")):
        f = glob.glob(str(cache / ("clarky_*_re%d_n%d.csv" % (re_, nc))))
        if f:
            x = pd.read_csv(f[0])
            axs[0].plot(x.cd, x.cl, ls, c=c, lw=1, label="XFOIL Ncrit %d Re %gk" % (nc, re_ / 1e3))
            axs[1].plot(x.alpha_deg, x.cl, ls, c=c, lw=1)
axs[0].set_xlim(0, 0.06); axs[0].set_xlabel("cd"); axs[0].set_ylabel("cl"); axs[0].legend(fontsize=7); axs[0].grid(alpha=.3)
axs[1].set_xlim(-6, 16); axs[1].set_xlabel("alpha [deg]"); axs[1].set_ylabel("cl"); axs[1].grid(alpha=.3)
fig.suptitle("Clark-Y 11.7 %: XFOIL 6.97 (clarky.dat) vs UIUC LSAT Vol. 3 Clark-Y (B) measurement")
plt.tight_layout(); plt.savefig(out, dpi=80)

"""Sanity checks and overview plot of the ntnu_bt1 case -> extraction/check_ntnu_bt1.png"""
import json
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.dirname(HERE)
meta = json.load(open(os.path.join(CASE, "case.json")))
rd = lambda f: pd.read_csv(os.path.join(CASE, f), comment="#")
g = rd("geometry.csv")
R = meta["R"]
assert g.r_m.is_monotonic_increasing and g.r_m.max() <= R + 1e-9, "radius"
s826 = g[g.airfoil == "s826"]
assert (np.diff(s826.twist_deg) < 0).all(), "twist must decrease outboard"
assert g.chord_m.max() < 0.1 and g.chord_m.min() > 0.01, "chord in m"
assert set(g.airfoil) <= set(meta["airfoils"]) | {"transition"}
print("geometry ok: %d rows, r %.4f..%.4f m, chord %.4f..%.4f m, twist %.2f..%.2f deg" %
      (len(g), g.r_m.min(), g.r_m.max(), g.chord_m.min(), g.chord_m.max(), s826.twist_deg.min(), s826.twist_deg.max()))
P = rd("airfoils/s826_polars.csv")
print("polars: Re", sorted(P.re.unique()), "alpha %.1f..%.1f" % (P.alpha_deg.min(), P.alpha_deg.max()))
assert (P.cd > 0).all()
exps = {e["file"]: rd(e["file"]) for e in meta["experiments"]}
for f, d in exps.items():
    assert (d.cp.dropna() < 16 / 27).all(), f + ": Cp above Betz"
    print("%-62s %3d rows  cp max %.3f  ct max %.3f" % (f, len(d), d.cp.max(), d.get("ct", d.get("ct_incl_tower")).max()))

# cross-check: K&L Table II vs BT1 Fig. 3 at common TSR (both on R = 0.447 m now)
a = exps["experiment/performance_krogstad_lund2012_tableII.csv"]
b = exps["experiment/performance_bt1_krogstad_eriksen2013_fig3.csv"]
x = a.tsr[(a.tsr > b.tsr.min()) & (a.tsr < b.tsr.max())]
dcp = np.interp(x, b.tsr, b.cp) - a.cp[x.index]
dct = np.interp(x, b.tsr, b.ct) - a.ct[x.index]
print("BT1 Fig.3 minus K&L Table II: dcp %.3f..%.3f, dct %.3f..%.3f" % (dcp.min(), dcp.max(), dct.min(), dct.max()))

fig, ax = plt.subplots(2, 3, figsize=(16, 9))
ax[0, 0].plot(g.r_m, g.chord_m, "o-"); ax[0, 0].set_ylabel("chord [m]")
ax[0, 0].axvline(meta["Rhub"], color="k", ls=":", label="Rhub"); ax[0, 0].axvline(R, color="r", ls=":", label="R")
ax[0, 0].legend()
ax2 = ax[0, 0].twinx(); ax2.plot(g.r_m, g.twist_deg, "s-", color="C1"); ax2.set_ylabel("twist [deg]", color="C1")
ax[0, 0].set_xlabel("r [m]"); ax[0, 0].set_title("geometry.csv")
for re_, d in P.groupby("re"):
    ax[0, 1].plot(d.alpha_deg, d.cl, ".-", label="%g" % re_)
    ax[0, 2].plot(d.alpha_deg, d.cd, ".-")
N = rd("airfoils/s826_polars_ntnu_bartl2019.csv")
d = N[N.re == 1e5]
ax[0, 1].plot(d.alpha_deg, d.cl, "k--", lw=2, label="NTNU 1e5 (alt)")
ax[0, 2].plot(d.alpha_deg, d.cd, "k--", lw=2)
X = rd("airfoils/s826_polars_bt5_xfoil_re1e5.csv")
ax[0, 1].plot(X.alpha_deg, X.cl, "m:", lw=2, label="XFOIL BT5 1e5"); ax[0, 2].plot(X.alpha_deg, X.cd, "m:", lw=2)
ax[0, 1].set_xlim(-12, 26); ax[0, 2].set_xlim(-12, 26); ax[0, 2].set_ylim(0, 0.5)
ax[0, 1].legend(fontsize=7, ncol=2); ax[0, 1].set_title("S826 cl: DTU Sarlak 2018 (coloured), NTNU (default, dashed)"); ax[0, 2].set_title("S826 cd")
for f, d in exps.items():
    lab = os.path.basename(f).replace(".csv", "")
    if "yaw_sweep" in f or "table2" in f:
        ax[1, 2].plot(d.yaw_deg, d.cp, "o-", label=lab + " cp"); ax[1, 2].plot(d.yaw_deg, d.ct, "s--", label=lab + " ct")
    elif "yaw_krog" in f:
        for y, dd in d.groupby("yaw_deg"):
            ax[1, 1].plot(dd.tsr, dd.cp, ".-", label="yaw %d" % y)
    elif "re_dep" in f:
        for v, dd in d.groupby("v_ms"):
            ax[1, 0].plot(dd.tsr, dd.ct_incl_tower, ".:", lw=0.7, label="ct incl. tower U=%.2f" % v)
    else:
        ax[1, 0].plot(d.tsr, d.cp, "o-", label=lab + " cp"); ax[1, 0].plot(d.tsr, d.ct, "s-", label=lab + " ct")
ax[1, 0].set_xlabel("tsr"); ax[1, 0].legend(fontsize=6); ax[1, 0].set_title("performance (yaw 0)")
ax[1, 1].set_xlabel("tsr"); ax[1, 1].set_ylabel("cp"); ax[1, 1].legend(fontsize=7); ax[1, 1].set_title("K&A 2012 Fig. 5, cp vs tsr at yaw")
ax[1, 2].set_xlabel("yaw [deg]"); ax[1, 2].legend(fontsize=6); ax[1, 2].set_title("TSR 6 yaw sweeps")
for a_ in ax.flat:
    a_.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "check_ntnu_bt1.png"), dpi=90)
print("wrote check_ntnu_bt1.png")

#!/usr/bin/env python3
"""Sanity checks and check plots (extraction/check_<case>.png) for the UIUC propeller cases."""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
CASES = HERE.parents[2] if HERE.parent.name == "extraction" else Path(
    "/Users/miha/Documents/dev/Python-BEM/validation/cases")


def rd(p):
    return pd.read_csv(p, comment="#", skipinitialspace=True)


def check(case_dir):
    meta = json.load(open(case_dir / "case.json"))
    R = meta["R"]
    g = rd(case_dir / meta["geometry_file"])
    probs = []
    if abs(g.r_m.iloc[-1] - R) > 1e-6:
        probs.append("geometry does not end at R (%.5f vs %.5f)" % (g.r_m.iloc[-1], R))
    if not (g.chord_m.max() < 0.5 * R):
        probs.append("chord suspicious")
    if not g.r_m.is_monotonic_increasing:
        probs.append("r not increasing")
    for a in set(g.airfoil) - {"transition"}:
        if a not in meta["airfoils"]:
            probs.append("airfoil %s missing in case.json" % a)
    exps = []
    for e in meta["experiments"]:
        d = rd(case_dir / e["file"])
        if "eta" in d and (d.eta > 1).any():
            probs.append("%s: eta > 1" % e["file"])
        if "eta" in d:
            bad = (np.abs(d.J * d.ct / d.cp - d.eta) > 0.01) & (d.cp.abs() > 0.01)
            if bad.any():
                probs.append("%s: eta != J ct/cp in %d rows" % (e["file"], bad.sum()))
        exps.append((e, d))
    pols = {}
    for n, a in meta["airfoils"].items():
        p = rd(case_dir / a["polars"])
        pols[n] = p
        for re_, gg in p.groupby("re"):
            if not (gg.alpha_deg.min() < 0 < gg.alpha_deg.max()):
                probs.append("polar %s Re=%g does not span 0" % (n, re_))
    # ---------------- plot
    perf = [x for x in exps if x[0]["x"] == "J"]
    stat = [x for x in exps if x[0]["x"] == "rpm"]
    nrow = 2 + len(perf)
    fig, axs = plt.subplots(nrow, 3, figsize=(16, 3.4 * nrow))
    ax = axs[0, 0]
    ax.plot(g.r_m / R, g.chord_m / R, "k.-", label="geometry.csv")
    for k, f in meta.get("alternative_geometry", {}).items():
        ga = rd(case_dir / f)
        ax.plot(ga.r_m / R, ga.chord_m / R, ".--", label=k)
    ax.set_xlabel("r/R"); ax.set_ylabel("c/R"); ax.legend(fontsize=7); ax.set_xlim(0, 1.02); ax.grid(alpha=.3)
    ax.set_title("%s  R=%.4f m  B=%d" % (meta["id"], R, meta["B"]), fontsize=9)
    ax = axs[0, 1]
    ax.plot(g.r_m / R, g.twist_deg, "k.-", label="geometry.csv")
    for k, f in meta.get("alternative_geometry", {}).items():
        ga = rd(case_dir / f)
        ax.plot(ga.r_m / R, ga.twist_deg, ".--", label=k)
    for i, (r, a) in enumerate(zip(g.r_m / R, g.airfoil)):
        if a != "transition" and (i == 0 or g.airfoil.iloc[i - 1] != a):
            ax.annotate(a, (r, g.twist_deg.iloc[i]), fontsize=6, rotation=60)
    ax.set_xlabel("r/R"); ax.set_ylabel("twist = beta [deg]"); ax.legend(fontsize=7); ax.grid(alpha=.3)
    ax = axs[0, 2]
    if "thickness_ratio" in g:
        ax.plot(g.r_m / R, g.thickness_ratio, "k.-")
    ax.set_xlabel("r/R"); ax.set_ylabel("t/c"); ax.grid(alpha=.3)
    # polars: first used airfoil and last used
    used = [n for n in meta["airfoils"] if meta["airfoils"][n].get("used_in_geometry", True)]
    show = [used[0], used[-1]] if len(used) > 1 else used
    for j, n in enumerate(show[:2]):
        p = pols[n]
        ax = axs[1, j]
        for re_, gg in p.groupby("re"):
            ax.plot(gg.alpha_deg, gg.cl, ".-", ms=2, lw=.8, label="%gk" % (re_ / 1e3))
        ax.set_title("%s cl (XFOIL Ncrit 9)" % n, fontsize=8); ax.set_xlabel("alpha"); ax.grid(alpha=.3)
        ax.legend(fontsize=6, ncol=2)
    ax = axs[1, 2]
    for n in show:
        p = pols[n]
        for re_, gg in p.groupby("re"):
            if re_ in (30e3, 60e3, 100e3):
                ax.plot(gg.cd, gg.cl, ".-", ms=2, lw=.8, label="%s %gk" % (n, re_ / 1e3))
    if "clarky_uiuc_lsat" in pols:
        for re_, gg in pols["clarky_uiuc_lsat"].groupby("re"):
            if re_ in (60e3, 100e3, 200e3):
                ax.plot(gg.cd, gg.cl, "s", ms=3, label="Clark-Y meas %gk" % (re_ / 1e3))
    ax.set_xlim(0, 0.1); ax.set_xlabel("cd"); ax.set_ylabel("cl"); ax.legend(fontsize=6); ax.grid(alpha=.3)
    for i, (e, d) in enumerate(perf):
        for j, y in enumerate(["ct", "cp", "eta"]):
            ax = axs[2 + i, j]
            for rpm, gg in d.groupby("rpm"):
                ax.plot(gg.J, gg[y], ".-", ms=3, lw=.8, label="%d rpm" % rpm)
            ax.set_xlabel("J"); ax.set_ylabel(y); ax.grid(alpha=.3)
            ax.set_title("%s (pitch %+g)" % (e["file"], e["defaults"].get("pitch_deg", 0)), fontsize=8)
            if j == 0:
                ax.legend(fontsize=6)
        for e2, s in stat:
            if e2["defaults"].get("pitch_deg") == e["defaults"].get("pitch_deg"):
                ax = axs[2 + i, 0]
                ax.plot(np.zeros(len(s)), s.ct, "kx", ms=4, label="static")
                axs[2 + i, 1].plot(np.zeros(len(s)), s.cp, "kx", ms=4)
    plt.tight_layout()
    out = case_dir / "extraction" / ("check_%s.png" % meta["id"])
    plt.savefig(out, dpi=75)
    plt.close(fig)
    return probs, out


if __name__ == "__main__":
    ids = sys.argv[1:] or sorted(p.name for p in CASES.iterdir() if p.name.startswith("uiuc_"))
    for cid in ids:
        probs, out = check(CASES / cid)
        print(cid, "OK" if not probs else probs, out.name)

"""Cross-check of PV_r (electrical) against the Jan-1993 Risoe time series in sources/ny_r_*.zip
(10-s block averages of electrical power vs met-mast wind speed, 1 m/s bins). Also prints the
air density of each campaign from the logged pressure and temperature."""
import os, glob, subprocess
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(os.path.dirname(HERE), "sources")
blocks = []
for z in sorted(glob.glob(os.path.join(SRC, "ny_r_*.zip"))):
    # old 'implode' compression -> use the unzip tool (Python's zipfile cannot read it)
    names = subprocess.run(["unzip", "-Z1", z], capture_output=True, text=True).stdout.split()
    name = [n for n in names if "PROF" in n.upper()][0]
    L = subprocess.run(["unzip", "-p", z, name], capture_output=True).stdout.decode("latin-1").splitlines()
    d = np.array([[float(x) for x in l.split()] for l in L[5:] if l.strip()])
    rho = d[:, 1] * 100 / (287.05 * (d[:, 2] + 273.15))
    print(os.path.basename(z), "rho %.3f  V %.2f  rpm %.2f  Pel %.1f kW" % (rho.mean(), d[:, 3].mean(), d[:, 6].mean(), d[:, 7].mean()))
    n = 250
    m = len(d) // n
    b = d[:m * n].reshape(m, n, -1).mean(1)
    blocks.append(b[:, [3, 7]])
b = np.vstack(blocks)
pv = np.loadtxt(os.path.join(HERE, "raw", "PV_r"), skiprows=5)
print(" bin  n  Pel_timeseries  PV_r_elec  diff%")
for v0 in np.arange(6, 17, 1.0):
    s = (b[:, 0] >= v0 - 0.5) & (b[:, 0] < v0 + 0.5)
    if s.sum() > 3:
        ref = np.interp(b[s, 0].mean(), pv[:, 0], pv[:, 1])
        print("%4.0f %3d %8.1f %10.1f %7.1f" % (v0, s.sum(), b[s, 1].mean(), ref, 100 * (b[s, 1].mean() / ref - 1)))

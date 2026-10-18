"""Generate the NACA 63-8xx section coordinates of the Bahaj et al. (2007) 800 mm rotor: one section per Table 1
station, at the station's exact t/c (24.0, 22.5, 20.7, ... 12.6 %).

Why one section per station: the blades were lofted from NACA 63-812/815/818/821/824 coordinates interpolated to the
Table 1 thickness of each of the 17 stations (Bahaj et al. 2007, Sec. 3.1), and the authors' BEM interpolated the
section data "for the correct thickness ratio" (Batten et al. 2008, Sec. 3). The runner assigns each annulus the
airfoil of the nearest tabulated station and its 'transition' rule interpolates by section index, not t/c, so the
thickness variation is represented directly: a 63-series section at the exact t/c for every station.

Generator: naca456 v1.1 (R. L. Carmichael, Public Domain Aeronautical Software, https://www.pdas.com/naca456.html,
AIAA 2001-5235), the modern revision of the NASA Langley Ladson-Brooks programs (NASA TM X-3069, TM 4741); source
copied to extraction/naca456/src. Build: `gfortran -fdefault-real-8 naca456.f90 -o naca456` (gfortran 16.2).
63-series thickness (conformal-mapping epsilon/psi functions of the 63-series basic thickness form, scaled linearly
to the t/c - the NACA/Abbott-von Doenhoff practice for 6-series sections of other thicknesses), mean line NACA
6-series a = 1.0, cl_i = 0.8, 'coarse' output (dencode = 1: the standard NACA ordinate stations, 51 points).
Check: the package's 63-021 sample reproduces the Abbott & von Doenhoff ordinates to 1e-6.
(Finer naca456 spacings put ~0.01 %c points ahead of the LE that make XFOIL's BL march diverge, so the
standard stations are used; XFOIL repanels to 240 panels anyway.)

Mean line choice: by NACA convention (Abbott & von Doenhoff, Theory of Wing Sections, Sec. 6.4; NACA TR 824) a
6-series designation without "a = ..." uses the a = 1.0 mean line; the Southampton papers write "NACA 63-8xx" without
an a-value. The drawing of the NACA 63-815 foil tested by the same group (Molland et al. 2004, Fig. 2) agrees with
a = 1.0 and not with a = 0.8 (camber 0.054 instead of 0.044): see overlay_molland2004_fig2.png.

Writes ../airfoils/naca63-8_tNNN.dat (NNN = 10 x t/c in %) and overlay_molland2004_fig2.png.
"""
import os
import subprocess

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'airfoils')
N456 = os.path.join(HERE, 'naca456')
EXE = os.environ.get('NACA456_EXE', '/tmp/naca456/naca456')

# Table 1 of Bahaj et al. 2007 (= Table 1 of Batten et al. 2008): t/c in % at r/R = 0.20 ... 1.00
TC_TABLE1 = [24.0, 22.5, 20.7, 19.5, 18.7, 18.1, 17.6, 17.1, 16.6, 16.1, 15.6, 15.1, 14.6, 14.1, 13.6, 13.1, 12.6]


def foil_name(tc):
    return f'naca63-8_t{round(tc * 10):03d}'


def build_exe():
    if not os.path.exists(EXE):
        os.makedirs(os.path.dirname(EXE), exist_ok=True)
        subprocess.run(['gfortran', '-fdefault-real-8', '-O2', 'naca456.f90', '-o', EXE],
                       cwd=os.path.join(N456, 'src'), check=True)


def run_naca456(name, tc, a=1.0):
    nml = os.path.join(N456, name + '.nml')
    with open(nml, 'w') as f:
        f.write(f"{name}: NACA 63-8xx at t/c = {tc / 100:.3f}; 63-series thickness, a={a} mean line, cl_i = 0.8\n"
                f"&NACA\n  name = '{name}',\n  profile = '63',\n  toc = {tc / 100:.4f},\n  camber = '6',\n"
                f"  cl = 0.8,\n  a = {a},\n  dencode = 1/\n")
    r = subprocess.run([EXE], input=os.path.basename(nml) + '\n', text=True, capture_output=True, cwd=N456)
    assert 'End of naca456' in r.stdout, r.stdout + r.stderr
    os.replace(os.path.join(N456, 'naca.out'), os.path.join(N456, name + '.out'))
    for junk in ('naca.gnu', 'naca.dbg'):
        p = os.path.join(N456, junk)
        if os.path.exists(p):
            os.remove(p)
    return os.path.join(N456, name + '.out')


def read_out(path):
    """Upper/lower surface (x, y) from the 'xupper yupper xlower ylower' block of a naca456 output."""
    lines = open(path).read().splitlines()
    i0 = next(i for i, L in enumerate(lines) if L.split()[:4] == ['xupper', 'yupper', 'xlower', 'ylower'])
    up, lo = [], []
    for L in lines[i0 + 1:]:
        p = L.split()
        if len(p) != 5:
            break
        xu, yu, xl, yl = map(float, p[1:])
        up.append((xu, yu))
        lo.append((xl, yl))
    return np.array(up), np.array(lo)


def main():
    build_exe()
    os.makedirs(OUT, exist_ok=True)
    for tc in TC_TABLE1:
        name = foil_name(tc)
        up, lo = read_out(run_naca456(name, tc))
        pts = np.vstack([up[::-1], lo[1:]])          # TE -> upper -> LE -> lower -> TE (Selig)
        xx = np.linspace(0, 1, 2001)
        tmax = np.max(np.interp(xx, up[:, 0], up[:, 1]) - np.interp(xx, lo[:, 0], lo[:, 1]))
        with open(os.path.join(OUT, name + '.dat'), 'w') as f:
            f.write(f'NACA 63-8xx t/c={tc / 100:.3f} (63-series thickness, a=1.0 mean line, cl_i=0.8; naca456 v1.1)\n')
            for x, y in pts:
                f.write(f' {x:.6f} {y: .6f}\n')
        print(f'{name}: {len(pts)} points, t/c (vertical) = {tmax:.4f}')

    # ---- overlay on Molland et al. 2004 Fig. 2 (middle panel: their tested 63-815 "experiment profile") ----------
    o10 = read_out(run_naca456('cmp_63-815_a10', 15.0, 1.0))
    o08 = read_out(run_naca456('cmp_63-815_a08', 15.0, 0.8))
    img = plt.imread(os.path.join(HERE, 'molland2004_fig2_sections.png'))
    gray = img[..., :3].mean(axis=2) if img.ndim == 3 else img
    # panel frame = longest dark columns/rows in the middle panel (frame x = -0.1 .. 1.1, y = 0.15 .. -0.1)
    band = gray[420:660, :]
    colsum = (band < 0.6).sum(axis=0)
    xcols = [c for c in range(gray.shape[1]) if colsum[c] > 0.6 * band.shape[0]]
    rowsum = (gray[:, 180:1340] < 0.6).sum(axis=1)
    yrows = [r for r in range(400, 670) if rowsum[r] > 0.5 * 1160]

    def groups(v):
        g, out = [], []
        for c in v:
            if g and c - g[-1] > 2:
                out.append(np.mean(g))
                g = []
            g.append(c)
        if g:
            out.append(np.mean(g))
        return out

    gx, gy = groups(xcols), groups(yrows)
    assert len(gx) == 2 and len(gy) == 2, (gx, gy)
    ax_, bx_ = np.polyfit([-0.1, 1.1], gx, 1)
    ay_, by_ = np.polyfit([0.15, -0.1], gy, 1)
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.imshow(img)
    for (up, lo), col, lab in [(o10, 'r', 'naca456 63-815 a=1.0 (used)'), (o08, 'b', 'naca456 63-815 a=0.8')]:
        ax.plot(ax_ * up[:, 0] + bx_, ay_ * up[:, 1] + by_, '-', color=col, lw=0.8, label=lab)
        ax.plot(ax_ * lo[:, 0] + bx_, ay_ * lo[:, 1] + by_, '-', color=col, lw=0.8)
    ax.set_xlim(150, 1371)
    ax.set_ylim(700, 400)
    ax.legend(loc='lower left', fontsize=8)
    ax.set_title('Molland et al. 2004 Fig. 2 (NACA 63-815 experiment profile) vs generated sections')
    ax.set_axis_off()
    fig.savefig(os.path.join(HERE, 'overlay_molland2004_fig2.png'), dpi=120, bbox_inches='tight')
    print('wrote overlay_molland2004_fig2.png')


if __name__ == '__main__':
    main()

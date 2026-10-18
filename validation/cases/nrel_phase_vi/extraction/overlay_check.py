"""Overlay the FINAL experiment CSV values back onto renders of the source figures (audit of the whole chain:
vector extraction -> calibration -> CSV). For every figure the CSV value is mapped back to page coordinates with
the same per-panel calibration the extractor used, then drawn as a red ring on a 200-dpi render of the page.

Writes overlay_simms2001_p27..p40.png, overlay_sorensen2002_p6..p9.png, overlay_lindenburg2003_p63.png.
"""
import os

import fitz
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import extract_simms2001 as es
import extract_sorensen2002 as eso

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.dirname(HERE)
DPI = 200
Z = DPI / 72.0

pc = pd.read_csv(os.path.join(CASE, "experiment", "power_curve.csv"), comment="#")
ld = pd.read_csv(os.path.join(CASE, "experiment", "loads.csv"), comment="#")
sp = pd.read_csv(os.path.join(CASE, "experiment", "spanwise_cn_ct.csv"), comment="#")


def render(page):
    pix = page.get_pixmap(dpi=DPI)
    img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)
    return img


def to_px(page, x, y):
    p = fitz.Point(x, y) * page.rotation_matrix
    return p.x * Z, p.y * Z


def save(img, pts, name, crop=None):
    fig, ax = plt.subplots(figsize=(img.shape[1] / 100, img.shape[0] / 100), dpi=100)
    ax.imshow(img)
    for (px, py, lab) in pts:
        ax.plot(px, py, "o", ms=14, mfc="none", mec="red", mew=1.6)
        if lab:
            ax.text(px + 12, py - 12, lab, color="red", fontsize=9)
    ax.set_axis_off()
    if crop:
        ax.set_xlim(crop[0], crop[2]); ax.set_ylim(crop[3], crop[1])
    fig.subplots_adjust(0, 0, 1, 1)
    fig.savefig(os.path.join(HERE, name), dpi=100)
    plt.close(fig)


def simms():
    doc = fitz.open(es.PDF)
    sim_series = {
        "lss_torque_nm": pc[pc.campaign.notna()].set_index("v_nominal").torque_nm,
        "root_flap_moment_nm": ld.set_index("v_nominal").flap_moment_gauge_axes_nm,
    }
    for rr, col in (("030", 0.30), ("063", 0.633), ("095", 0.95)):
        s = sp[np.isclose(sp.r_R, col)].set_index("v_nominal")
        sim_series["cn_" + rr] = s.cn_simms2001
        sim_series["ct_" + rr] = s.ct_simms2001
        sim_series["cm_" + rr] = s.cm
        sim_series["q_local_pa_" + rr] = s.q_local_pa
    vpos = pc[pc.campaign.notna()].set_index("v_nominal").v_ms
    for page_no, (fig, qty) in es.FIGURES.items():
        page = doc[page_no - 1]
        dr = page.get_drawings()
        words = page.get_text("words")
        grids = [x for x in dr if abs((x.get("width") or 0) - 0.12) < 0.01 and len(x["items"]) > 10]
        pts = []
        for g in grids:
            cw, cv, _ = es.calibrate(g, words)
            for vn, val in sim_series[qty].items():
                y = (vpos[vn] - cw[0]) / cw[1]
                x = (val - cv[0]) / cv[1]
                px, py = to_px(page, x, y)
                pts.append((px, py, ""))
        save(render(page), pts, f"overlay_simms2001_p{page_no}.png")


def sorensen():
    doc = fitz.open(eso.PDF)
    jobs = {6: [pc[pc.campaign.notna()].set_index("v_nominal").torque_sorensen2002_nm],
            7: [ld.set_index("v_nominal").root_flap_sorensen2002_nm, ld.set_index("v_nominal").edge_moment_gauge_axes_nm]}
    cases = ["S0700000", "S1000000", "S1300000", "S1500000", "S2000000", "S2500000"]
    for page_no in (6, 7, 8, 9):
        page = doc[page_no - 1]
        dr = page.get_drawings()
        words = page.get_text("words")
        L = eso.lines(dr)
        frames = [d["rect"] for d in dr if len(d["items"]) == 1 and d["items"][0][0] == "qu" and d["rect"].width > 100]
        frames += [d["rect"] for d in dr if len(d["items"]) > 20 and d["items"][0][0] == "l"
                   and d["rect"].width > 180 and 70 < d["rect"].height < 76]
        frames = sorted(frames, key=lambda r: r.y0)
        pts = []
        for k, fr in enumerate(frames):
            cx, cy, *_ = eso.calibrate(fr, L, words)
            if page_no in jobs:
                ser = jobs[page_no][k]
                for vn, val in ser.dropna().items():
                    pts.append((*to_px(page, (vn - cx[0]) / cx[1], (val - cy[0]) / cy[1]), ""))
                    if page_no == 7 and k == 1:
                        sd = ld.set_index("v_nominal").edge_moment_gauge_axes_nm_std[vn]
                        for s in (+sd, -sd):
                            pts.append((*to_px(page, (vn - cx[0]) / cx[1], (val + s - cy[0]) / cy[1]), ""))
            else:
                s = sp[sp.campaign == cases[k]]
                col = "cn_chord" if page_no == 8 else "ct_chord"
                for _, r in s.iterrows():
                    xr = {0.466: 0.47, 0.633: 0.63}.get(round(r.r_R, 3), r.r_R)
                    pts.append((*to_px(page, (xr - cx[0]) / cx[1], (r[col] - cy[0]) / cy[1]), ""))
        save(render(page), pts, f"overlay_sorensen2002_p{page_no}.png")


def lindenburg():
    path = os.path.join(CASE, "sources", "Lindenburg2003_ECN-C-03-025.pdf")
    page = fitz.open(path)[62]
    raw = pd.read_csv(os.path.join(HERE, "lindenburg2003_measured_raw.csv"))
    # recompute the calibration exactly as the extractor does: reuse its fitted mapping through the raw table
    # (raw v/value -> page) by fitting the inverse from the tick marks again
    import extract_lindenburg2003 as el  # noqa: F401  (documentation: same calibration code)
    dr = page.get_drawings()
    words = page.get_text("words")
    pts = []
    for fig, qty, ylo, yhi, ser in (("Fig. 6.2", "lss_torque_knm", 0, 280, pc.set_index("v_ms").torque_lindenburg2003_nm / 1000),
                                    ("Fig. 6.3", "root_flap_moment_knm", 300, 560, ld.set_index("v_ms").root_flap_lindenburg2003_nm / 1000)):
        segs = [(it[1].x, it[1].y, it[2].x, it[2].y) for d in dr if ylo <= d["rect"].y0 and d["rect"].y1 <= yhi
                for it in d["items"] if it[0] == "l"]
        xt = sorted({round(a, 2) for a, b, c, d in segs if abs(a - c) < 1e-3 and 3.5 < abs(b - d) < 4.5})
        yt = sorted({round(b, 2) for a, b, c, d in segs if abs(b - d) < 1e-3 and 3.5 < abs(a - c) < 4.5})
        ybot = max(max(b, d) for a, b, c, d in segs if abs(a - c) < 1e-3 and 3.5 < abs(b - d) < 4.5)
        xleft = min(min(a, c) for a, b, c, d in segs if abs(b - d) < 1e-3 and 3.5 < abs(a - c) < 4.5)
        xp, xv, yp, yv = [], [], [], []
        for x in xt:
            c = [w for w in words if el.num(w[4]) is not None and ybot < w[1] < ybot + 6 and abs((w[0] + w[2]) / 2 - x) < 4]
            if len(c) == 1:
                xp.append(x); xv.append(el.num(c[0][4]))
        for y in yt:
            c = [w for w in words if el.num(w[4]) is not None and xleft - 25 < w[0] < xleft and abs((w[1] + w[3]) / 2 - y) < 4]
            if len(c) == 1:
                yp.append(y); yv.append(el.num(c[0][4]))
        cx = np.polyfit(xv, xp, 1)
        cy = np.polyfit(yv, yp, 1)
        for v, val in ser.dropna().items():
            pts.append((*to_px(page, np.polyval(cx, v), np.polyval(cy, val)), ""))
    save(render(page), pts, "overlay_lindenburg2003_p63.png", crop=(250, 0, 1500, 1500))


if __name__ == "__main__":
    simms()
    sorensen()
    lindenburg()
    print("overlays written")

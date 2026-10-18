"""Exact extraction of Krogstad & Adaramola (2012), Wind Energy 15:743-756, doi:10.1002/we.502, Fig. 4
(Reynolds-number dependence: CP and CT vs TSR at U = 7.22 ... 10.30 m/s, yaw 0, pitch 0) from the VECTOR
drawing in the PDF with PyMuPDF. PDF page index 4.

All series are drawn in black; they are told apart by marker shape (legend):
 7.22 filled circle | 7.72 open circle | 8.24 open triangle | 8.76 circle with plus | 9.33 asterisk
 9.75 square with plus | 10.30 open star
Marker paths: filled circle = 'fs' black fill, 4 Bezier curves; open circle = 'fs' white fill, 4 curves;
triangle = 'fs' white, 3 lines; circle-plus = 's' 4 curves (+ 2 separate lines); square-plus = 's' rectangle;
star = 'fs' white, 10 lines; asterisk = separate stroked lines only -> its centre is the midpoint of its two
diagonal strokes (no other marker has diagonal strokes; an asterisk may coincide with another marker).
Axis calibration from the frame (= end ticks, checked against tick labels):
 (a) x: 110.3 pt -> 0, 287.9 pt -> 12; y: 204.5 pt -> 0, 72.3 pt -> 0.5
 (b) x: 326.5 pt -> 0, 504.1 pt -> 12; y: 204.5 pt -> 0, 72.1 pt -> 1.4
Precision ~0.1 pt: 0.007 in TSR, 0.0004 in CP, 0.001 in CT.
"""
import os
import fitz
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "..", "sources", "krogstad_adaramola_2012_we502.pdf")
page = fitz.open(PDF)[4]

PANELS = {
    "cp": dict(x0=110.3, x1=287.9, y0=204.5, y1=72.3, yv=0.5, legend=(175, 140, 199, 203)),
    "ct": dict(x0=326.5, x1=504.1, y0=204.5, y1=72.1, yv=1.4, legend=(440, 140, 469, 203)),
}
U = {"fc": 7.22, "oc": 7.72, "tri": 8.24, "cplus": 8.76, "ast": 9.33, "splus": 9.75, "star": 10.30}


def close(c, ref, tol=0.05):
    return c is not None and all(abs(a - b) < tol for a, b in zip(c, ref))


def shape(d):
    kinds = tuple(i[0] for i in d["items"])
    fill = d.get("fill")
    if d["type"] == "fs" and kinds == ("c",) * 4:
        return "oc" if close(fill, (1, 1, 1)) else "fc"
    if d["type"] == "fs" and kinds == ("l",) * 3:
        return "tri"
    if d["type"] == "fs" and kinds == ("l",) * 10:
        return "star"
    if d["type"] == "s" and kinds == ("c",) * 4:
        return "cplus"
    if d["type"] == "s" and kinds == ("re",):
        return "splus"
    if d["type"] == "s" and kinds == ("l",):
        return "line"
    return None


rows = []
for name, P in PANELS.items():
    found, diag = [], []
    for d in page.get_drawings():
        r = d["rect"]
        cx, cy = (r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2
        if not (P["x0"] - 3 <= cx <= P["x1"] + 3 and P["y1"] - 3 <= cy <= P["y0"] + 3):
            continue
        lx0, ly0, lx1, ly1 = P["legend"]
        if lx0 <= cx <= lx1 and ly0 <= cy <= ly1:
            continue
        s = shape(d)
        if s == "line":
            # asterisk strokes: only the asterisk has DIAGONAL strokes (the plus signs of the circle-plus and
            # square-plus markers are horizontal/vertical), so a ~4 pt diagonal stroke marks an asterisk centre
            # (45 deg: equal width and height; the two diagonals of one asterisk have the same bounding box)
            if 3.0 < r.width < 5.0 and abs(r.width - r.height) < 0.15:
                diag.append((cx, cy))
            continue
        if s is None or not (2.0 < r.width < 5.0 and 2.0 < r.height < 5.0):
            continue
        found.append((s, cx, cy))
    # asterisk centres: merge the (two) diagonal strokes that share a midpoint
    ast = []
    D = np.array(diag) if diag else np.zeros((0, 2))
    for x, y in diag:
        if np.sum(np.hypot(D[:, 0] - x, D[:, 1] - y) < 0.05) < 2:   # need both diagonals of the 'x'
            continue
        if not any(np.hypot(x - a, y - b) < 0.8 for a, b in ast):
            ast.append((x, y))
    found += [("ast", x, y) for x, y in ast]
    for s, cx, cy in found:
        tsr = (cx - P["x0"]) / (P["x1"] - P["x0"]) * 12
        val = (cy - P["y0"]) / (P["y1"] - P["y0"]) * P["yv"]
        rows.append(dict(panel=name, v_ms=U[s], tsr=round(tsr, 3), value=round(val, 4), px=cx, py=cy))

df = pd.DataFrame(rows).drop_duplicates(subset=["panel", "v_ms", "tsr", "value"]).sort_values(["panel", "v_ms", "tsr"])
print(df.groupby(["panel", "v_ms"]).size())
df.drop(columns=["px", "py"]).to_csv(os.path.join(HERE, "krogstad_adaramola2012_fig4_raw.csv"), index=False)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
clip = fitz.Rect(85, 65, 510, 230)
pix = page.get_pixmap(dpi=250, clip=clip)
img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)
k = 250 / 72
fig = plt.figure(figsize=(pix.width / 100, pix.height / 100), dpi=100)
ax = fig.add_axes([0, 0, 1, 1]); ax.imshow(img); ax.axis("off")
cols = dict(zip(sorted(U.values()), ["red", "orange", "gold", "green", "cyan", "blue", "magenta"]))
for v, g in df.groupby("v_ms"):
    ax.plot((g.px - clip.x0) * k, (g.py - clip.y0) * k, "x", color=cols[v], ms=6, mew=1.2, label="%.2f" % v)
ax.legend(loc="upper center", ncol=7, fontsize=8)
fig.savefig(os.path.join(HERE, "overlay_krogstad_adaramola2012_fig4.png"))

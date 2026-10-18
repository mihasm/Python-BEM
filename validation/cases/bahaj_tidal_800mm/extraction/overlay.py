"""Overlay check: extracted marker centres (red) drawn on the 300 dpi render of each figure page, and the
de-shifted Fig. 8/9 data (blue) drawn into the unshifted lower panels, where they must scatter about the
authors' curve fits. Writes overlay_fig7.png, overlay_fig8.png, overlay_fig9.png, overlay_fig10_11.png."""
import json
import os

import fitz
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, 'raw_markers.json')))
doc = fitz.open(os.path.join(HERE, '..', 'sources', 'bahaj2007_power_thrust_mct.pdf'))
DPI = 300
S = DPI / 72.0
EXP = os.path.join(HERE, '..', 'experiment')


def read(name):
    return pd.read_csv(os.path.join(EXP, name), comment='#')


def to_px(panel, x, y):
    c = R[panel]['calib']
    return (np.asarray(x) - c['bx']) / c['ax'] * S, (np.asarray(y) - c['by']) / c['ay'] * S


def page_img(pno):
    pix = doc[pno].get_pixmap(dpi=DPI)
    return np.frombuffer(pix.samples, np.uint8).reshape(pix.height, pix.width, pix.n)


def draw(pno, panels, extra, out, clip):
    img = page_img(pno)
    fig, ax = plt.subplots(figsize=(12, 12 * (clip[3] - clip[1]) / (clip[2] - clip[0])))
    ax.imshow(img)
    for p in panels:
        px, py = to_px(p, [m['x'] for m in R[p]['markers']], [m['y'] for m in R[p]['markers']])
        ax.plot(px, py, 'o', mfc='none', mec='r', ms=7, mew=0.8)
    for p, x, y in extra:
        px, py = to_px(p, x, y)
        ax.plot(px, py, '+', color='b', ms=7, mew=1.0)
    ax.set_xlim(clip[0] * S, clip[2] * S)
    ax.set_ylim(clip[3] * S, clip[1] * S)
    ax.set_axis_off()
    fig.savefig(os.path.join(HERE, out), dpi=110, bbox_inches='tight')
    plt.close(fig)
    print('wrote', out)


draw(10, ['fig7a_cp', 'fig7b_ct'], [], 'overlay_fig7.png', (60, 60, 400, 600))
allsets = pd.read_csv(os.path.join(HERE, 'fig8_fig9_deshifted_all_sets.csv'))
f8 = allsets[allsets.figure.astype(str) == '8']
draw(11, ['fig8_cp_upper', 'fig8_ct_upper'],
     [('fig8_cp_lower', f8.tsr, f8.cp), ('fig8_ct_lower', f8.tsr, f8.ct)], 'overlay_fig8.png', (80, 60, 400, 390))
f9 = allsets[allsets.figure.astype(str) == '9']
draw(12, ['fig9_cp_upper', 'fig9_ct_upper'],
     [('fig9_cp_lower', f9.tsr, f9.cp), ('fig9_ct_lower', f9.tsr, f9.ct)], 'overlay_fig9.png', (80, 60, 400, 385))
draw(13, ['fig10a_cp', 'fig10b_ct', 'fig11a_cp', 'fig11b_ct'], [], 'overlay_fig10_11.png', (30, 60, 430, 490))

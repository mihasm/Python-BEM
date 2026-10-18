"""Find where each digitised curve changes from dashed (low-J, faired/extrapolated) to solid.
Samples the raw (normalised) image along the traced path, skipping samples within 4 px of a gridline;
a 1-px-step path sample is a 'gap' when the max darkness across +-2 px is below 0.15. The solid start is the
J above which no 30-px window contains >= 4 gap samples."""
import json, sys, numpy as np
from scipy import ndimage as ndi
from calib import grids2
from digilib import darkness_map

def scan(name, page, box):
    D = json.load(open('dig_%s.json' % name))
    g, vl, hl = grids2(page, box)
    d = ndi.maximum_filter(darkness_map(g, atten=1.0), size=5)
    vl = np.array(vl); hl = np.array(hl)
    out = {}
    for lab, c in D['curves'].items():
        p = np.array(c['px']); x = np.array(c['x'])
        o = np.argsort(p[:, 0]); p = p[o]; x = x[o]
        seg = np.r_[0, np.cumsum(np.hypot(*np.diff(p, axis=0).T))]
        s = np.arange(0, seg[-1], 1.0)
        xs = np.interp(s, seg, p[:, 0]); ys = np.interp(s, seg, p[:, 1]); js = np.interp(s, seg, x)
        near = (np.abs(xs[:, None] - vl[None]).min(1) < 4) | (np.abs(ys[:, None] - hl[None]).min(1) < 4)
        v = ndi.map_coordinates(d, [ys, xs], order=1)
        gap = (v < 0.15) & ~near
        cnt = np.convolve(gap.astype(int), np.ones(30, int), 'same')
        dashed = cnt >= 4
        jd = js[dashed]
        # only consider the low-J part (dashes are at the left end); ignore isolated hits beyond J=1
        jd = jd[jd < 0.8]
        out[lab] = float(jd.max()) if len(jd) else None
    print(name, {k: (round(v, 3) if v is not None else None) for k, v in out.items()})
    return out

if __name__ == '__main__':
    figs = {'f8_ct_2b': (9, (246, 1995, 2170, 3030)), 'f6_cp_2b': (8, (390, 1640, 2200, 3005)),
            'f9_ct_3b': (10, (285, 310, 2270, 1610)), 'f11_cp_3b': (11, (290, 300, 2210, 2275)),
            'f13_ct_4b': (12, (270, 290, 2260, 1890)), 'f15_cp_4b': (13, (290, 430, 2215, 2950))}
    res = {k: scan(k, *v) for k, v in figs.items()}
    json.dump(res, open('dash_ends.json', 'w'), indent=1)

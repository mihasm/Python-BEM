import sys, numpy as np
from digilib import *
def normalise(g, size=41):
    from scipy import ndimage as ndi
    bg = ndi.maximum_filter(ndi.median_filter(g[::4, ::4], 9), 5)
    bg = ndi.zoom(bg, 4, order=1)[:g.shape[0], :g.shape[1]]
    return np.clip(g / np.maximum(bg, 1) * 250.0, 0, 255)

def grids(page, box, thr=200, minlen=120, norm=False):
    g = load_gray('../pg%d-000.png' % page)
    gl = normalise(g) if norm else g
    H, V = line_maps(gl, thr=thr, minlen=minlen)
    x0, y0, x1, y1 = box
    vl = grid_positions(V, 'v', x0, x1, y0, y1, minfrac=0.4)
    hl = grid_positions(H, 'h', y0, y1, x0, x1, minfrac=0.4)
    return (gl if norm else g), vl, hl
if __name__ == '__main__':
    p = int(sys.argv[1]); box = tuple(map(int, sys.argv[2:6]))
    g, vl, hl = grids(p, box)
    print('v', len(vl), np.round(vl, 1).tolist()); print('dv', np.round(np.diff(vl), 1).tolist())
    print('h', len(hl), np.round(hl, 1).tolist()); print('dh', np.round(np.diff(hl), 1).tolist())


def lattice(g, box, axis, spacing=None):
    """Regular gridline lattice: mean darkness profile across the box, its period by autocorrelation,
    then every lattice line refined to the local profile peak (+-4 px). Returns list of positions."""
    x0, y0, x1, y1 = box
    d = np.clip(250.0 - g[y0:y1, x0:x1], 0, None)
    prof = d.mean(0) if axis == 'v' else d.mean(1)
    off = x0 if axis == 'v' else y0
    prof = prof - np.median(prof)
    if spacing is None:
        ac = np.correlate(prof, prof, 'full')[len(prof) - 1:]
        lo = 30
        spacing = lo + int(np.argmax(ac[lo:200]))
    # phase: maximise summed profile on the lattice
    best = None
    for sp in np.arange(spacing - 1.0, spacing + 1.01, 0.05):
        for ph in np.arange(0, sp, 0.5):
            pos = np.arange(ph, len(prof), sp)
            s = np.interp(pos, np.arange(len(prof)), prof).mean()
            if best is None or s > best[0]:
                best = (s, ph, sp)
    spacing = best[2]
    pos = np.arange(best[1], len(prof), spacing)
    out = []
    for p in pos:
        i0, i1 = int(max(0, p - 4)), int(min(len(prof), p + 5))
        seg = prof[i0:i1]
        j = int(np.argmax(seg)) + i0
        lo, hi = max(0, j - 2), min(len(prof), j + 3)
        w = np.clip(prof[lo:hi], 0, None)
        c = (np.arange(lo, hi) * w).sum() / w.sum() if w.sum() > 0 else j
        if prof[j] > 0.25 * np.percentile(prof[np.round(pos).astype(int).clip(0, len(prof) - 1)], 50):
            out.append(off + float(c))
    return out, spacing


def grids2(page, box, norm=True):
    g = load_gray('../pg%d-000.png' % page)
    gl = normalise(g) if norm else g
    vl, sx = lattice(gl, box, 'v'); hl, sy = lattice(gl, box, 'h')
    return gl, vl, hl

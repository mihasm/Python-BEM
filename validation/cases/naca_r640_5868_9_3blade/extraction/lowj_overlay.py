"""Low-J (J<0.62) zoom of each chart with the digitised curves drawn over it, for checking."""
import json, sys, numpy as np
from PIL import Image, ImageDraw
from calib import normalise
from digilib import load_gray
COLS = [(255,0,0),(0,150,0),(0,0,255),(200,0,200),(255,128,0),(0,170,170),(120,60,0)]
def lowj(name, page, jmax=0.62, out=None):
    D = json.load(open('dig_%s.json' % name))
    g = normalise(load_gray('../pg%d-000.png' % page))
    c0 = list(D['curves'].values())[0]
    px = np.array(c0['px']); x = np.array(c0['x']); k = np.polyfit(x, px[:, 0], 1)
    x0 = int(np.polyval(k, 0)) - 10; x1 = int(np.polyval(k, jmax))
    sel = [np.array(c['px'])[np.array(c['x']) < jmax] for c in D['curves'].values()]
    allp = np.vstack([s for s in sel if len(s)]); y0 = int(allp[:, 1].min()) - 20; y1 = int(allp[:, 1].max()) + 5
    im = Image.fromarray(g[y0:y1, x0:x1].astype(np.uint8)).convert('RGB'); dr = ImageDraw.Draw(im)
    for (lab, c), col in zip(D['curves'].items(), COLS):
        p = np.array(c['px']); xx = np.array(c['x']); m = xx < jmax + 0.02
        q = [(a - x0, b - y0) for a, b in p[m]]
        if len(q) > 1: dr.line(q, fill=col, width=1)
        o = np.argsort(xx); Y = np.interp(0.0, xx[o], p[o, 1]) - y0
        dr.text((12, Y - 14), lab, fill=col)
    im.save(out or 'lowj_%s.png' % name)
    return im
if __name__ == '__main__':
    lowj(sys.argv[1], int(sys.argv[2]))

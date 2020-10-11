import math

import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as mp3d
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


def create_face(p1, p2, p3, p4, *args, **kwargs):
    coords = [p1, p2, p3, p4]
    face = mp3d.art3d.Poly3DCollection(
        [coords], facecolors=facecolors, alpha=.5, linewidth=0.0, *args, **kwargs)
    return face


def rotate(origin, point, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle must be given in radians.
    """
    ox, oy = origin
    px, py = point

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return qx, qy


def rotate_array(x, y, xy, angle):
    angle = math.radians(angle)
    x_out, y_out = [], []
    for _x, _y in zip(x, y):
        _xq, _yq = rotate(xy, (_x, _y), angle)
        x_out.append(_xq)
        y_out.append(_yq)
    return x_out, y_out


def scale_and_normalize(foil_x, foil_y, scale, centroid):
    foil_x, foil_y = np.array(foil_x), np.array(foil_y)
    foil_x = foil_x-centroid[0]
    foil_x, foil_y = foil_x*scale, foil_y*scale
    return foil_x, foil_y


def create_3d_blade(SET_INIT, flip_turning_direction=False, propeller_geom=False):
    theta = SET_INIT["theta"]
    r = SET_INIT["r"]
    c = SET_INIT["c"]
    foil = SET_INIT["foils"]
    airfoils = SET_INIT["airfoils"]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    out_x, out_y, out_z = [], [], []
    data = []
    for i in range(len(r)):
        _r = r[i]
        _c = c[i]
        _airfoil = foil[i]
        _theta = theta[i]  # - because of direction

        if _airfoil != "transition": #the only exception
            if propeller_geom:
                _theta = -_theta

            _airfoil_x, _airfoil_y = airfoils[_airfoil]["x"], airfoils[_airfoil]["y"]

            _centroid_x, _centroid_y = airfoils[_airfoil]["centroid_x"], airfoils[_airfoil]["centroid_y"]

            if flip_turning_direction:
                _airfoil_x = -np.array(_airfoil_x)
                _theta = -_theta
                _centroid_x = -_centroid_x

            _centroid = (_centroid_x, _centroid_y)
            
            _airfoil_x, _airfoil_y = scale_and_normalize(_airfoil_x, _airfoil_y, _c, _centroid)
            _airfoil_x, _airfoil_y = rotate_array(_airfoil_x, _airfoil_y, (0, 0), _theta)

            list_x, list_y = [], []
            for _x, _y in zip(_airfoil_x, _airfoil_y):
                out_x.append(_x)
                out_y.append(_y)
                out_z.append(_r)
                list_x.append(_x)
                list_y.append(_y)

            data.append([_r, np.array(list_x), np.array(list_y)])

    X, Y, Z = np.array(out_x), np.array(out_y), np.array(out_z)
    max_range = np.array(
        [X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max() / 2.0
    mid_x = (X.max()+X.min()) * 0.5
    mid_y = (Y.max()+Y.min()) * 0.5
    mid_z = (Z.max()+Z.min()) * 0.5
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)
    # ax.set_aspect("equal")
    Axes3D.mouse_init(ax, rotate_btn=1, zoom_btn=2)

    ax.scatter(X, Y, Z)

    # plt.show()

    data = np.array(data)
    return {"data": data, "X": X, "Y": Y, "Z": Z}

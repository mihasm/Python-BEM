import math

import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as mp3d
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy import interpolate

from utils import get_transition_foils


def create_face(p1, p2, p3, p4, *args, **kwargs):
    """

    :param p1:
    :param p2:
    :param p3:
    :param p4:
    :param args:
    :param kwargs:
    :return:
    """
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
    """

    :param x:
    :param y:
    :param xy:
    :param angle: in degrees
    :return:
    """
    angle = math.radians(angle)
    x_out, y_out = [], []
    for _x, _y in zip(x, y):
        _xq, _yq = rotate(xy, (_x, _y), angle)
        x_out.append(_xq)
        y_out.append(_yq)
    return x_out, y_out


def scale_and_normalize(foil_x, foil_y, scale, centroid):
    """

    :param foil_x:
    :param foil_y:
    :param scale:
    :param centroid:
    :return:
    """
    foil_x, foil_y = np.array(foil_x), np.array(foil_y)
    foil_x = foil_x - centroid[0]
    foil_x, foil_y = foil_x * scale, foil_y * scale
    return foil_x, foil_y


def closest_point(x_point, y_point, x_list, y_list):
    """
    Find the closest point in zip(x_list,y_list), to the point x_point,y_point.
    """
    dist = [0] * len(x_list)

    min_index = None
    min_distance = None

    for i in range(len(x_list)):
        dist_v = np.sqrt((x_list[i] - x_point) ** 2 + (y_list[i] - y_point) ** 2)
        if min_distance == None or dist_v < min_distance:
            min_index = i
            min_distance = dist_v

    return min_index, x_list[min_index], y_list[min_index]


def point_along_line(x1, y1, x2, y2, t):
    """
    if t < 0.5, new point is closer to x1,y1
    if t = 1.0, new point is == x2,y2
    if t = 0.0, new point is == x1,y1
    """
    x, y = (1 - t) * x1 + t * x2, (1 - t) * y1 + t * y2
    return x, y


def interpolate_foils(foil_1_x, foil_1_y, foil_2_x, foil_2_y, t):
    """

    :param foil_1_x:
    :param foil_1_y:
    :param foil_2_x:
    :param foil_2_y:
    :param t:
    :return:
    """
    out_x, out_y = [], []

    tck, u = interpolate.splprep([foil_1_x, foil_1_y], s=0)
    xspline, yspline = interpolate.splev(np.linspace(0, 1, 500), tck, der=0)

    tck2, u2 = interpolate.splprep([foil_2_x, foil_2_y], s=0)
    xspline2, yspline2 = interpolate.splev(np.linspace(0, 1, 100), tck2, der=0)

    # get list of zeros
    zeros_idx_1 = np.where(np.diff(np.sign(yspline)))[0]
    zeros_idx_2 = np.where(np.diff(np.sign(yspline2)))[0]

    zero_1_index = zeros_idx_1[np.abs(zeros_idx_1 - len(zeros_idx_1) / 2).argmin()]
    zero_2_index = zeros_idx_2[np.abs(zeros_idx_2 - len(zeros_idx_2) / 2).argmin()]

    for i in range(len(xspline2)):
        x_point, y_point = xspline2[i], yspline2[i]
        if i <= zero_2_index:
            xspline_in, yspline_in = xspline[:zero_1_index], yspline[:zero_1_index]
        else:
            xspline_in, yspline_in = xspline[zero_1_index:], yspline[zero_1_index:]
        index, x_found, y_found = closest_point(x_point, y_point, xspline_in, yspline_in)
        x_new, y_new = point_along_line(x_found, y_found, x_point, y_point, t)
        out_x.append(x_new)
        out_y.append(y_new)

    return out_x, out_y


def create_3d_blade(input_data, flip_turning_direction=False, propeller_geom=False):
    """

    :param input_data:
    :param flip_turning_direction:
    :param propeller_geom:
    :return:
    """
    theta = input_data["theta"]
    r = input_data["r"]
    c = input_data["c"]
    foil = input_data["foils"]
    airfoils = input_data["airfoils"]
    R = input_data["R"]
    Rhub = input_data["Rhub"]

    r = [Rhub] + list(r) + [R]
    c = [c[0]] + list(c) + [c[-1]]
    theta = [theta[0]] + list(theta) + [theta[-1]]
    foil = [foil[0]] + list(foil) + [foil[-1]]
    transition_foils = get_transition_foils(foil)

    out_x, out_y, out_z = [], [], []
    data = []
    for i in range(len(r)):
        _r = r[i]
        _c = c[i]
        _airfoil = foil[i]
        _theta = theta[i]  # - because of direction

        if _airfoil == "transition":
            _prev_foil = transition_foils[i][0]
            _next_foil = transition_foils[i][1]
            _transition_coefficient = transition_foils[i][2]
            _airfoil_x_prev, _airfoil_y_prev = airfoils[_prev_foil]["x"], airfoils[_prev_foil]["y"]
            _airfoil_x_next, _airfoil_y_next = airfoils[_next_foil]["x"], airfoils[_next_foil]["y"]
            _airfoil_x, _airfoil_y = interpolate_foils(_airfoil_x_prev, _airfoil_y_prev, _airfoil_x_next,
                                                       _airfoil_y_next, _transition_coefficient)
            _centroid_x_prev, _centroid_y_prev = airfoils[_prev_foil]["centroid_x"], airfoils[_prev_foil]["centroid_y"]
            _centroid_x_next, _centroid_y_next = airfoils[_next_foil]["centroid_x"], airfoils[_next_foil]["centroid_y"]
            _centroid_x, _centroid_y = point_along_line(_centroid_x_prev, _centroid_y_prev, _centroid_x_next,
                                                        _centroid_y_next, _transition_coefficient)
        else:
            _airfoil_x, _airfoil_y = airfoils[_airfoil]["x"], airfoils[_airfoil]["y"]
            _centroid_x, _centroid_y = airfoils[_airfoil]["centroid_x"], airfoils[_airfoil]["centroid_y"]

        if propeller_geom:
            _theta = -_theta

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

    # DRAW ########
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X, Y, Z = np.array(out_x), np.array(out_y), np.array(out_z)
    max_range = np.array(
        [X.max() - X.min(), Y.max() - Y.min(), Z.max() - Z.min()]).max() / 2.0
    mid_x = (X.max() + X.min()) * 0.5
    mid_y = (Y.max() + Y.min()) * 0.5
    mid_z = (Z.max() + Z.min()) * 0.5
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)
    Axes3D.mouse_init(ax, rotate_btn=1, zoom_btn=2)
    ax.scatter(X, Y, Z)
    # plt.show()
    ###############

    # data = np.array(data)
    return {"data": data, "X": X, "Y": Y, "Z": Z}

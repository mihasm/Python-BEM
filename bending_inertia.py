import math

import numpy as np
from scipy.interpolate import interp1d

pi = math.pi


def PointsInCircum(r, n=1000):
    """
    Generates circle x,y coordinates.
    :param r: float: radius
    :param n: int: number of points
    :return: tuple: (list of x points, list of y points)
    """
    x_points = [math.cos(2 * pi / n * x) * r for x in range(0, n + 1)]
    y_points = [math.sin(2 * pi / n * x) * r for x in range(0, n + 1)]
    return x_points, y_points


def interpolate_airfoil(x, y, num_interp=100):
    """
    Interpolates airfoil x,y data. Airfoil has to be non-rotated.
    :param x: list of floats: x coordinates
    :param y: list of floats: y coordinates
    :param num_interp: int: number of interpolated points
    :return: tuple: (list of x points, list of y points)
    """
    cross = np.where(np.diff(np.signbit(np.gradient(x))))[0][0]

    x_up = x[:cross + 2]
    y_up = y[:cross + 2]
    x_down = x[cross + 1:]
    y_down = y[cross + 1:]

    interp_up = interp1d(x_up, y_up, 'linear')
    interp_down = interp1d(x_down, y_down, 'linear')

    _x = np.linspace(np.min(x), np.max(x), num_interp)
    _y_up = interp_up(_x)
    _y_down = interp_down(_x)
    x_out = np.concatenate((np.flip(_x), _x))
    y_out = np.concatenate((np.flip(_y_up), _y_down))
    return x_out, y_out


def calculate_bending_inertia_2(x, y):
    """
    Calculates the bending intertia (second moment of area) for a given set of points.
    Any polygon equations: https://en.wikipedia.org/wiki/Second_moment_of_area
    Area equation: https://en.wikipedia.org/wiki/Polygon#Area
    :param x: list of floats: x coordinates
    :param y: list of floats: y coordinates
    :return: tuple: (Ix,Iy,Ixy,A)
    """
    Iy = 0
    Ix = 0
    Ixy = 0
    A = 0
    for i in range(len(x) - 1):
        Iy += 1 / 12 * (x[i] * y[i + 1] - x[i + 1] * y[i]) * (x[i] ** 2 + x[i] * x[i + 1] + x[i + 1] ** 2)
        Ix += 1 / 12 * (x[i] * y[i + 1] - x[i + 1] * y[i]) * (y[i] ** 2 + y[i] * y[i + 1] + y[i + 1] ** 2)
        Ixy += 1 / 24 * (x[i] * y[i + 1] - x[i + 1] * y[i]) * (
                x[i] * y[i + 1] + 2 * x[i] * y[i] + 2 * x[i + 1] * y[i + 1] + x[i + 1] * y[i])
        A += 0.5 * (x[i] * y[i + 1] - x[i + 1] * y[i])
    return Ix, Iy, Ixy, A


def generate_hollow_foil(x, y, thickness):
    """
    Generates hollow airfoil from x,y coordinates.

    thickness should be given in p.u.

    This operation must be done BEFORE ROTATION.
    (otherwise, criterion should be modified)

    :param x:
    :param y:
    :param thickness:
    :return:
    """
    xout, yout = [], []

    cross = np.where(np.diff(np.signbit(np.gradient(x))))[0][0]

    x_up = x[:cross + 2]
    y_up = y[:cross + 2]
    x_down = x[cross + 1:]
    y_down = y[cross + 1:]

    interp_up = interp1d(x_up, y_up, 'linear', fill_value="extrapolate")
    interp_down = interp1d(x_down, y_down, 'linear', fill_value="extrapolate")

    for i in range(1, len(x)):
        if interp_up(x[i]) - interp_down(x[i]) > 2 * thickness:
            dx = x[i] - x[i - 1]
            dy = y[i] - y[i - 1]
            x_90 = -dy
            y_90 = dx
            vec_len = np.sqrt(x_90 ** 2 + y_90 ** 2)
            a = thickness / vec_len
            pristeto_x = x[i] + x_90 * a
            pristeto_y = y[i] + y_90 * a
            xout.append(pristeto_x)
            yout.append(pristeto_y)

    return xout, yout

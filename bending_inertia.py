import math

import numpy as np
from scipy.interpolate import interp1d

pi = math.pi


def PointsInCircum(r, n=1000):
    """
    Generate circle x,y coordinates...
    """
    x_points = [math.cos(2 * pi / n * x) * r for x in range(0, n + 1)]
    y_points = [math.sin(2 * pi / n * x) * r for x in range(0, n + 1)]
    return x_points, y_points


def interpolate_airfoil(x, y, num_interp=100):
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


def calculate_bending_inertia(x, y, num_interp=100):
    """
    Function that calculates the bending inertia (Second moment of Area)
    for a given airfoil, defined by its top and bottom line.

    Calculation result is given in p.u.^4 -> (per unit of length)^4

    Sources:https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/systems-labs-06/spl10b.pdf
            http://www.wingbike.nl/Wingbike_Hydrofoil/Background_files/Method%20to%20determine%20Section%20Modulus%20and%20Bending%20Inertia%20equations.pdf

    Example data [De Lannoy, Method to determine correction factors ...]:
    #x_up = [0,0.75,1.5,2.63,3.75,5.63,7.5,11.25,15,22.5,30,37.5,45,60,75,90,105,120,135,150]
    #y_up = [5.98,8.38,9.49,10.51,11.28,12.39,13.5,15.13,16.41,18.26,19.42,19.83,20,19.49,17.98,15.64,12.56,8.92,4.79,0.21]
    #x_down=[0,0.75,1.5,2.63,3.75,5.63,7.5,11.25,15,22.5,30,37.5,45,60,75,90,105,120,135,150]
    #y_down = [5.98,4.79,3.93,3.16,2.51,2.05,1.5,0.72,0.26,0.05,0,0,0,0,0,0,0,0,0,0]
    """

    cross = np.where(np.diff(np.signbit(np.gradient(x))))[0][0]

    x_up = x[:cross + 2]
    y_up = y[:cross + 2]
    x_down = x[cross + 1:]
    y_down = y[cross + 1:]

    interp_up = interp1d(x_up, y_up, 'linear', fill_value="extrapolate")
    interp_down = interp1d(x_down, y_down, 'linear', fill_value="extrapolate")

    x = np.linspace(np.min(x), np.max(x), num_interp)
    y_up = interp_up(x)
    y_down = interp_down(x)
    y_up = np.nan_to_num(y_up)
    y_down = np.nan_to_num(y_down)

    A = 0
    zsum = 0
    I = 0

    for i in range(len(x) - 1):
        dx = (x[i + 1] - x[i])
        dy_up = (y_up[i + 1] + y_up[i])
        dy_down = (y_down[i + 1] + y_down[i])

        A_up = dy_up * dx * 0.5
        A_down = dy_down * dx * 0.5
        A += A_up - A_down

        zsum += 0.5 * (y_up[i + 1] ** 2 - y_down[i + 1] ** 2) * dx

    z = zsum / A
    I = 0
    for i in range(len(x) - 1):
        dx = (x[i + 1] - x[i])
        I += 1 / 3 * ((y_up[i + 1] - z) ** 3 - (y_down[i + 1] - z) ** 3) * dx

    return I, A


def calculate_bending_inertia_2(x, y):
    """
    Any polygon equations: https://en.wikipedia.org/wiki/Second_moment_of_area
    Area equation: https://en.wikipedia.org/wiki/Polygon#Area
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

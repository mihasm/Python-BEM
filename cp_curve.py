__author__ = "Miha Smrekar"
__credits__ = ["Miha Smrekar"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Miha Smrekar"
__email__ = "miha.smrekar9@gmail.com"
__status__ = "Development"

"""
To build into exe, use pyinstaller, be sure to include scipy paths and hidden import:

$pyinstaller preracun_vetrnica.py
>>--onefile
>>--paths C:\\Users\\miha\\AppData\\Local\\Programs\\Python\\Python36-32\\lib\\site-packages\\scipy\\extra-dll
>>--hidden-import='scipy._lib.messagestream'
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy
from mpl_toolkits.mplot3d import Axes3D
import numbers
from PyQt5 import QtGui, QtCore, QtWidgets
import sys
from table import Table
from utils import dict_to_ar, transpose

from induction_factors import Calculator
from main import MainWindow

a = Axes3D  # only for passing code inspection -> Axes3D needs to be imported


def calculate_power(speed_wind, rpm, sections_radius, chord_lengths, chord_angles, dr, R, Rhub, B, f_c_L, f_c_D,
                    add_angle=None):
    """
    Returns calculated power using BEM analysis.

    Inputs are wind speed, rotational velocity, blade geometry, number of blades, and
    functions for calculating lift and drag coefficients.

    Output is a dictionary with all results.

    :param speed_wind: wind speed [m/s]
    :param rpm: rotational velocity [rpm]
    :param sections_radius: np array of section sections_radius [m]
    :param chord_lengths: np array of section chord lengths [m]
    :param chord_angles: np array of chord angles [degrees]
    :param dr: np array of section heights [m]
    :param R: outer (tip) radius [m]
    :param Rhub: hub radius [m]
    :param B: number of blades
    :param f_c_L: function for calculating lift coefficient
    :param f_c_D: function for calculating drag coefficient
    :param add_angle: [degrees]
    :return: dict with results
    """
    if add_angle != None:
        chord_angles = chord_angles + add_angle
    results = Calculator(f_c_L, f_c_D).run_array(chord_angle=chord_angles, B=B, c=chord_lengths, r=sections_radius,
                                                 dr=dr, rpm=rpm, v=speed_wind, R=R, Rhub=Rhub)
    return results


def calculate_power_3d(sections_radius, chord_lengths, chord_angles, dr, R, Rhub, B, f_c_L, f_c_D, add_angle=0.0):
    """
    Calculates power for given geometry and data for every windspeed and rpm.

    Returns dictionary with arrays with data for every point.

    :param sections_radius: np array of section sections_radius [m]
    :param chord_lengths: np array of section chord lengths [m]
    :param chord_angles: np array of chord angles [degrees]
    :param dr: np array of section heights [m]
    :param R: outer (tip) radius [m]
    :param Rhub: hub radius [m]
    :param B: number of blades
    :param f_c_L: function for calculating lift coefficient
    :param f_c_D: function for calculating drag coefficient
    :param add_angle: [degrees]
    :return: dictionary with all results stored as numpy arrays
    """

    results_3d = {}

    for v in list(numpy.linspace(start=3, stop=20, num=10)):
        for rpm in list(numpy.linspace(start=100, stop=3000, num=10)):
            print("Calculating power for v", v, "rpm", rpm)
            _results = calculate_power(v, rpm, sections_radius, chord_lengths, chord_angles, dr, add_angle=add_angle,
                                       R=R, B=B, f_c_L=f_c_L, f_c_D=f_c_D, Rhub=Rhub)
            print("    Power:", _results["power"], "Cp:", _results["cp"])
            if _results != None and _results["power"]:
                if 0.0 < _results["cp"] <= 0.6:
                    for key, value in _results.items():
                        if key not in results_3d:
                            results_3d[key] = []
                        results_3d[key].append(value)

    for k, v in results_3d.items():
        results_3d[k] = v

    return results_3d


def max_calculate(X, Y, Z):
    """
    Calculates maximum power for every wind speed.

    Returns only points that provide maximum power for given wind speed.

    :param X: Wind speed
    :param Y: RPM
    :param Z: Power
    :return: X,Y,Z (filtered)
    """
    X_un = numpy.unique(X)

    max_x = []
    max_y = []
    max_z = []

    for i in range(len(X_un)):
        X_max = 0.0
        Y_max = 0.0
        Z_max = 0.0
        for j in numpy.where(X == X_un[i])[0]:
            if Z[j] > Z_max:
                Z_max = Z[j]
                X_max = X[j]
                Y_max = Y[j]
        max_x.append(X_max)
        max_y.append(Y_max)
        max_z.append(Z_max)

    max_x = numpy.array(max_x)
    max_y = numpy.array(max_y)
    max_z = numpy.array(max_z)
    return max_x, max_y, max_z


def sort_xy(array_x, array_y):
    """
    Sorts two arrays by values in the first array.
    Useful for sorting pairs of x,y points.

    Returns sorted arrays.

    Arrays have to be the same length.

    :param array_x: x values
    :param array_y: y values
    :return: x,y (sorted)
    """
    out = []
    if len(array_x) == len(array_y):
        for i in range(len(array_x)):
            out.append((array_x[i], array_y[i]))
        out = sorted(out, key=lambda k: k[0])
        out_x, out_y = [], []
        for n in range(len(out)):
            out_x.append(out[n][0])
            out_y.append(out[n][1])
        return out_x, out_y
    else:
        raise Exception("Cannot create XY pairs with arrays with different num of elements")


def run(sections_radius, chord_lengths, chord_angles, dr, B, R, Rhub, f_c_L, f_c_D):
    """
    Main function of this .py file.

    Runs power calculation for given wind turbine geometry and draws results in three graphs.

    :param sections_radius: np array of section sections_radius [m]
    :param chord_lengths: np array of section chord lengths [m]
    :param chord_angles: np array of chord angles [degrees]
    :param dr: np array of section heights [m]
    :param R: outer (tip) radius [m]
    :param Rhub: hub radius [m]
    :param B: number of blades
    :param f_c_L: function for calculating lift coefficient
    :param f_c_D: function for calculating drag coefficient
    :return dictionary with results
    """

    results_3d = calculate_power_3d(sections_radius, chord_lengths, chord_angles, dr, B=B, R=R, Rhub=Rhub,
                                    f_c_D=f_c_D, f_c_L=f_c_L)

    X, Y, Z = results_3d["v"], results_3d["rpm"], results_3d["power"]
    max_x, max_y, max_z = max_calculate(X, Y, Z)  # dodajanje tock maksimalne poweri

    # RISANJE 3D SCATTER PLOT
    fig = plt.figure()
    ax_3D = fig.gca(projection='3d')
    p0 = ax_3D.plot_trisurf(X, Y, Z, cmap=plt.cm.CMRmap)
    plt.title('Theoretical wind turbine power')
    cbar = plt.colorbar(p0)
    cbar.set_label("Power [W]")
    ax_3D.set_xlabel('Wind speed [m/s]')
    ax_3D.set_ylabel('Rotational velocity [rpm]')
    ax_3D.set_zlabel('Power [W]')
    ax_3D.set_zlim(0, 2e3)

    # RISANJE GRAF P-V
    fig2 = plt.figure()
    ax_pv = fig2.add_subplot(1, 2, 1)

    ax_pv.plot(max_x, max_z, "-")
    plt.title("power vs. wind")
    ax_pv.set_xlabel("Wind speed [m/s]")
    ax_pv.set_ylabel("Power [W]")

    # RISANJE CP curve
    CP, TSR = results_3d["cp"], results_3d["TSR"]
    X1, Y1 = sort_xy(TSR, CP)

    ax_cp = fig2.add_subplot(1, 2, 2)
    ax_cp.plot(X1, Y1, '-')
    ax_cp.set_xlabel('lambda')
    ax_cp.set_ylabel('Cp')
    ax_cp.set_xlim(0, 15)
    ax_cp.set_ylim(0.0, 1.0)

    CT, A = results_3d["Ct"], results_3d["a"]
    # RISANJE CT CURVE
    fig4 = plt.figure()
    ax_ct = fig4.gca()
    ax_ct.scatter(A, CT)
    plt.title('Ct curve')

    # plt.show()
    return results_3d


def run_main(sections_radius, chord_lengths, chord_angles, dr, B, R, Rhub, f_c_L, f_c_D):
    """
    Main function of this .py file.

    Runs power calculation for given wind turbine geometry and draws results in three graphs.

    :param sections_radius: np array of section sections_radius [m]
    :param chord_lengths: np array of section chord lengths [m]
    :param chord_angles: np array of chord angles [degrees]
    :param dr: np array of section heights [m]
    :param R: outer (tip) radius [m]
    :param Rhub: hub radius [m]
    :param B: number of blades
    :param f_c_L: function for calculating lift coefficient
    :param f_c_D: function for calculating drag coefficient
    :return dictionary with results
    """

    results_3d = calculate_power_3d(sections_radius, chord_lengths, chord_angles, dr, B=B, R=R, Rhub=Rhub,
                                    f_c_D=f_c_D, f_c_L=f_c_L)

    X, Y, Z = results_3d["v"], results_3d["rpm"], results_3d["power"]
    max_x, max_y, max_z = max_calculate(X, Y, Z)

    app = QtWidgets.QApplication([])
    screen = app.primaryScreen()
    size = screen.size()
    main = MainWindow(size.width(), size.height())

    geom = [
        ["r"] + list(sections_radius),
        ["c"] + list(chord_lengths),
        ["theta"] + list(chord_angles),
        ["dr"] + list(dr)
    ]

    geom = transpose(geom)
    t_geom = Table(geom)
    main.tab_widget.add_tab_widget(t_geom, "Geom")

    alpha = numpy.linspace(-25, 25, 100)
    cL = f_c_L(alpha)
    cD = f_c_D(alpha)

    f = main.tab_widget.add_tab_figure("Cl/Cd")
    main.tab_widget.add_2d_plot_to_figure(f, alpha, cL, 121, "cL", "alpha", "cL")
    main.tab_widget.add_2d_plot_to_figure(f, alpha, cD, 122, "cD", "alpha", "cD")

    f2 = main.tab_widget.add_tab_figure("Moč in Cp")
    main.tab_widget.add_2d_plot_to_figure(f2, max_x, max_z, 121, "Moč vs veter", "veter [m/s]", "moč [W]")

    TSR, CP = sort_xy(results_3d["TSR"], results_3d["cp"])
    main.tab_widget.add_2d_plot_to_figure(f2, TSR, CP, 122, "Cp krivulja", "lambda", "Cp")

    f3 = main.tab_widget.add_tab_figure("3D moč")
    main.tab_widget.add_surface_plot(f3, X, Y, Z, 111, "Moč (veter,rpm)", "veter[m/s]", "rpm", "moč [W]")

    f4 = main.tab_widget.add_tab_figure("Ct krivulja")
    main.tab_widget.add_2d_plot_to_figure(f4, results_3d["a"], results_3d["Ct"], 111, "ct krivulja", "a", "ct",
                                          look="o")

    t = Table(dict_to_ar(results_3d))
    main.tab_widget.add_tab_widget(t, "Data")

    app.exec_()
    return results_3d

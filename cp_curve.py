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

import tkinter as tk
from tkinter import filedialog

import matplotlib
import matplotlib.pyplot as plt
import numpy
from matplotlib.widgets import Button
from mpl_toolkits.mplot3d import Axes3D

from induction_factors import Calculator

matplotlib.use('TkAgg')
a = Axes3D  # only for passing code inspection -> Axes3D needs to be imported
root = tk.Tk()
root.withdraw()


def save_file_dialog(text2save):
    """
    Opens window dialog for saving file.
    :param text2save: File contents.
    :return:
    """
    f = filedialog.asksaveasfile(
        title="Select where to export",
        defaultextension=".csv",
        filetypes=(("Comma-separated values", "*.csv"), ("all files", "*.*")))
    if f is None:
        return
    f.write(text2save)
    f.close()
    return


def calculate_power(speed_wind, rpm, sections_radius, chord_lengths, chord_angles, dr, R, Rhub, B, f_c_L, f_c_D,
                    add_angle=0.0):
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
    chord_angles = chord_angles + add_angle
    results = Calculator(f_c_L, f_c_D).run_array(chord_angle=chord_angles, B=B, c=chord_lengths, r=sections_radius,
                                                 dr=dr, rpm=rpm, v=speed_wind, R=R, Rhub=Rhub)
    return results


def calculate_power_3d(sections_radius, chord_lengths, chord_angles, dr, R, Rhub, B, f_c_L, f_c_D, add_angle=0.0):
    """
    Calculates power for given geometry and data for every windspeed and rpm.

    Returns arrays with data for every point.

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
    :return: wind_speed,rpm,power,cp,TSR for every point
    """
    x = []
    y = []
    z = []
    cp_ar = []
    tsr_ar = []
    for v in list(numpy.linspace(start=3, stop=20, num=10)):
        for rpm in list(numpy.linspace(start=100, stop=3000, num=10)):
            _results = calculate_power(v, rpm, sections_radius, chord_lengths, chord_angles, dr, add_angle=add_angle,
                                       R=R, B=B,
                                       f_c_L=f_c_L, f_c_D=f_c_D, Rhub=Rhub)
            if _results != None:
                power = _results["power"]
                cp = _results["cp"]
                TSR = _results["TSR"]
                print("v:", "%.2f" % round(v, 2), "m/s", "rpm:", "%.2f" % round(rpm, 2), "RPM", "power:",
                      "%.2f" % round(power, 2), "W", "TSR:", "%.2f" % round(TSR, 2), "Cp:", "%.2f" % round(cp, 2))
                if power != None:
                    if 0.0 <= cp <= 1.0:
                        x.append(v)
                        y.append(rpm)
                        z.append(power)
                        cp_ar.append(cp)
                        tsr_ar.append(TSR)

    X = numpy.array(x)
    Y = numpy.array(y)
    Z = numpy.array(z)
    CP_AR = numpy.array(cp_ar)
    TSR_AR = numpy.array(tsr_ar)

    return X, Y, Z, CP_AR, TSR_AR


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
    """

    X, Y, Z, CP, TSR = calculate_power_3d(sections_radius, chord_lengths, chord_angles, dr, B=B, R=R, Rhub=Rhub,
                                          f_c_D=f_c_D, f_c_L=f_c_L)

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
    axcut = plt.axes([0.9, 0.0, 0.1, 0.075], label="test")

    def save_ax_3D():
        """
        Helper function that handles export button
        """
        out = "wind_speed[m/s],rpm,power[W]\n"
        for i in range(len(X)):
            out += "%s,%s,%s,\n" % (X[i], Y[i], Z[i])
        return save_file_dialog(out)

    button = Button(axcut, 'Export')
    button.on_clicked(save_ax_3D)

    # ax_3D.scatter(X2, Y2, Z2, c='g',marker='^',label='tocke za Cp krivuljo')

    # ax_3D.scatter(max_x, max_y, max_z, c='r',label='power (max per v)')

    # RISANJE GRAF P-V
    fig2 = plt.figure()
    ax_pv = fig2.add_subplot(1, 2, 1)

    ax_pv.plot(max_x, max_z, "-")
    plt.title("power vs. wind")
    ax_pv.set_xlabel("Wind speed [m/s]")
    ax_pv.set_ylabel("Power [W]")
    axcut2 = plt.axes([0.4, 0.0, 0.08, 0.04], label="test2")

    def save_ax_pv():
        """
        Helper function that handles export button
        """
        out = "wind_speed[m/s],power[W],\n"
        for i in range(len(X)):
            out += "%s,%s,\n" % (max_x, max_z)
        return save_file_dialog(out)

    button2 = Button(axcut2, 'Export')
    button2.on_clicked(save_ax_pv)

    # RISANJE CP curve
    # X1,Y1 = cp_xy
    # print_sort_xy(X1,Y1)
    # fig2=plt.figure()
    # ax_cp=fig2.gca()
    X1, Y1 = sort_xy(TSR, CP)
    ax_cp = fig2.add_subplot(1, 2, 2)
    ax_cp.plot(X1, Y1, '-')
    ax_cp.set_xlabel('lambda')
    ax_cp.set_ylabel('Cp')
    ax_cp.set_xlim(0, 15)
    ax_cp.set_ylim(0.0, 1.0)

    axcut3 = plt.axes([0.9, 0.0, 0.08, 0.04], label="test3")

    def save_ax_cp():
        """
        Helper function that handles export button
        """
        out = "lambda,Cp,\n"
        for i in range(len(X)):
            out += "%s,%s,\n" % (X1, Y1)
        return save_file_dialog(out)

    button3 = Button(axcut3, 'Export')
    button3.on_clicked(save_ax_cp)

    # Risanje ct curve
    # fig4=plt.figure()
    # ax_ct=fig4.gca()
    # ax_ct.scatter(TSR_A,CTL_A)
    # plt.title('Ct curve')

    plt.show()


def optimize_runner(target_speed, sections_radius, chord_lengths, chord_angles, dr, R, Rhub, B, f_c_L, f_c_D, rpm,
                    min_add_angle=-30, max_add_angle=30, step=0.5):
    """
    This function calculates the optimum pitch angle of the blade for the given wind speed and rotational velocity.

    :param target_speed: target speed [m/s]
    :param sections_radius: np array of section sections_radius [m]
    :param chord_lengths: np array of section chord lengths [m]
    :param chord_angles: np array of chord angles [degrees]
    :param dr: np array of section heights [m]
    :param R: outer (tip) radius [m]
    :param Rhub: hub radius [m]
    :param B: number of blades
    :param f_c_L: function for calculating lift coefficient
    :param f_c_D: function for calculating drag coefficient
    :param rpm: rotational velocity [RPM]
    :param min_add_angle: lower limit [degrees]. Default: -30
    :param max_add_angle: upper limit [degrees]. Default: +30
    :param step: step for changing angle [degrees]
    :return: best angle [degrees]
    """
    print("\nOptimizing for wind speed of ", target_speed, "m/s...")
    power_orig = \
        calculate_power(target_speed, rpm, sections_radius, chord_lengths, chord_angles, dr, R, B=B, Rhub=Rhub, f_c_D=f_c_D,
                        f_c_L=f_c_L)["power"]
    print("Without turning the blade, the power is:", "%.2f" % round(power_orig), "Watts at wind speed", target_speed,
          "m/s and rotational velocity", rpm, "rpm")
    current_add_angle = min_add_angle
    results = []
    while current_add_angle <= max_add_angle:
        print("Testing pitch:", current_add_angle, "° at rpm", rpm)
        power = \
            calculate_power(target_speed, rpm, sections_radius, chord_lengths, chord_angles, dr, add_angle=current_add_angle,
                            R=R,
                            Rhub=Rhub, B=B, f_c_D=f_c_D, f_c_L=f_c_L)["power"]
        results.append((current_add_angle, power))
        current_add_angle += step
        print("---")
    best_angle, best_power = sorted(results, key=lambda x: x[1], reverse=True)[0]
    print("The best angle to turn would be:", "%.3f" % round(best_angle, 2), "°, at this pitch angle the power is:",
          "%.2f" % round(best_power), "Watts")
    print("That is a power gain of ", "%.2f" % round(best_power - power_orig, 2), "Watts, or",
          round(((best_power - power_orig) / power_orig) * 100, 2), "% better than the original version.")
    return best_angle

__author__ = "Miha Smrekar"
__credits__ = ["Miha Smrekar"]
__license__ = "GPL"
__version__ = "0.3.0"
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
import sys, traceback

import numpy
from mpl_toolkits.mplot3d import Axes3D
import time
import datetime
from math import pi

from induction_factors import Calculator
from utils import Printer

a = Axes3D  # only for passing code inspection -> Axes3D needs to be imported


def calculate_power(inp_args):
    """
    Returns calculated power using BEM analysis.

    Inputs are wind speed, rotational velocity, blade geometry, number of blades, and
    functions for calculating lift and drag coefficients.

    Output is a dictionary with all results.

    :return: dict with results
    """
    p = Printer(inp_args["return_print"])

    for f in inp_args["foils_in"]:
        if not f in inp_args["curves"].keys():
            p.print("Section foil %s does not exist in airfoil list." % f)
            raise Exception("Section foil not matching airfoil list error")

    if "add_angle" in inp_args:
        if inp_args["add_angle"] != None:
            inp_args["theta"] = inp_args["theta"] + inp_args["add_angle"]
    try:
        results = Calculator(inp_args["curves"]).run_array(**inp_args)
        return results
    except:
        p.print("!!!!EOF!!!!")
        traceback.print_exc(file=sys.stdout)
        raise


def calculate_power_3d(inp_args, print_eof=True, prepend="", print_out=True):
    """
    Calculates power for given geometry and data for every windspeed and rpm.

    Returns dictionary with arrays with data for every point.
    :return: dictionary with all results stored as numpy arrays
    """

    p = Printer(inp_args["return_print"])
    return_results = inp_args["return_results"]
    results_3d = {}
    speeds = list(
        numpy.linspace(
            start=inp_args["v_min"],
            stop=inp_args["v_max"],
            num=inp_args["v_num"]
        ))
    rpms = list(
        numpy.linspace(
            start=inp_args["rpm_min"],
            stop=inp_args["rpm_max"],
            num=inp_args["rpm_num"],
        )
    )
    total_iterations = int(len(speeds) * len(rpms))
    i = 0

    time_start = time.time()
    for v in speeds:
        for rpm in rpms:

            if print_out:
                p.print(prepend + "Calculating power for v:", v, "[m/s], rpm:", rpm, "[RPM], lambda:",
                        rpm / 60 * 2 * pi * inp_args["R"] / v, "J:", v / (rpm / 60 * inp_args["R"] * 2))
            _inp_args = {**inp_args}
            _inp_args["v"] = v
            _inp_args["rpm"] = rpm
            try:
                _results = calculate_power(_inp_args)
            except:
                p.print("Error in running calculate_power")
                p.print("!!!!EOF!!!!")
                raise
            if _results == "!!!!EOF!!!!":
                return None
            if _results != None and _results["power"]:
                if print_out:
                    p.print(prepend + "    TSR:", _results["TSR"], "J:", _results["J"], "cp:", _results["cp_w"],
                            "ct_p:", _results["ct_p"])
                for key, value in _results.items():
                    if key not in results_3d:
                        results_3d[key] = []
                for key, value in _results.items():
                    results_3d[key].append(value)
                return_results.append(results_3d)

            i += 1
            t_now = int(time.time() - time_start)
            t_left = int((total_iterations / i - 1) * t_now)
            t_left_str = str(datetime.timedelta(seconds=t_left))
            eta_seconds = datetime.datetime.now() + datetime.timedelta(seconds=t_left)
            eta = str(eta_seconds).split(".")[0]
            p.print("    ### Time left:", t_left_str, "ETA:", eta, "###")

    return_results.append(results_3d)
    if print_eof:
        p.print("!!!!EOF!!!!")
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

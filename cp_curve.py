__author__ = "Miha Smrekar"
__credits__ = ["Miha Smrekar"]
__license__ = "GPL"
__version__ = "0.2.6"
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

import numpy
from mpl_toolkits.mplot3d import Axes3D

from induction_factors import Calculator

a = Axes3D  # only for passing code inspection -> Axes3D needs to be imported


def calculate_power(
    inp_args
):
    """
    Returns calculated power using BEM analysis.

    Inputs are wind speed, rotational velocity, blade geometry, number of blades, and
    functions for calculating lift and drag coefficients.

    Output is a dictionary with all results.

    :param return_results: lst, used for returning results to main class
    :param return_print: lst, used for printing using main class
    :param print_all: prints every iteration
    :param relaxation_factor: relaxation factor
    :param method: method of calculating induction factors
    :param rho: air density [kg/m^3]
    :param convergence_limit: convergence criterion
    :param max_iterations: maximum number of iterations
    :param cascade_correction: uses cascade correction
    :param new_hub_loss: uses new tip loss correction
    :param new_tip_loss: uses new tip loss correction
    :param hub_loss: uses Prandtl tip loss correction
    :param tip_loss: uses Prandtl tip loss correction
    :param print_out: prints detailed iteration data
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
    if "add_angle" in inp_args:
        if inp_args["add_angle"] != None:
            inp_args["chord_angles"] = inp_args["chord_angles"] + inp_args["add_angle"]
    results = Calculator(inp_args["f_c_L"], inp_args["f_c_D"], inp_args["inverse_f_c_L"]).run_array(**inp_args)
    return results


def calculate_power_3d(inp_args):
    """
    Calculates power for given geometry and data for every windspeed and rpm.

    Returns dictionary with arrays with data for every point.
    :param return_results: lst, used for returning results to main class
    :param return_print: lst, used for printing using main class
    :param print_all: prints every iteration
    :param rpm_num: number of rpm points
    :param rpm_max: maximum rpm [RPM]
    :param rpm_min: minimum rpm [RPM]
    :param v_num: number of speed points
    :param v_max: maximum speed [m/s]
    :param v_min: minimum speed [m/s]
    :param relaxation_factor: relaxation factor
    :param method: method of calculating induction factors
    :param rho: air density [kg/m^3]
    :param convergence_limit: convergence criterion
    :param max_iterations: maximum number of iterations
    :param cascade_correction: uses cascade correction
    :param new_hub_loss: uses new tip loss correction
    :param new_tip_loss: uses new tip loss correction
    :param hub_loss: uses Prandtl tip loss correction
    :param tip_loss: uses Prandtl tip loss correction
    :param print_out: prints detailed iteration data
    :param r: np array of section sections_radius [m]
    :param c: np array of section chord lengths [m]
    :param theta: np array of chord angles [degrees]
    :param dr: np array of section heights [m]
    :param R: outer (tip) radius [m]
    :param Rhub: hub radius [m]
    :param B: number of blades
    :param f_c_L: function for calculating lift coefficient
    :param f_c_D: function for calculating drag coefficient
    :param add_angle: [degrees]
    :return: dictionary with all results stored as numpy arrays
    """
    
    return_print = inp_args["return_print"]
    return_results = inp_args["return_results"]
    results_3d = {}

    for v in list(numpy.linspace(start=inp_args["v_min"], stop=inp_args["v_max"], num=inp_args["v_num"])):
        for rpm in list(numpy.linspace(start=inp_args["rpm_min"], stop=inp_args["rpm_max"], num=inp_args["rpm_num"])):
            _p = (
                "\nCalculating power for v: "
                + str(v)
                + "[m/s] rpm: "
                + str(rpm)
                + "[RPM]\n"
            )
            inp_args["v"] = v
            inp_args["rpm"] = rpm
            return_print.append(_p)
            _results = calculate_power(
                inp_args
            )

            _p = (
                "    Power: "
                + str(_results["power"])
                + "[W] Cp: "
                + str(_results["cp"])
                + "\n"
            )
            return_print.append(_p)

            if _results != None and _results["power"]:
                if 0.0 < _results["cp"] <= 0.6:
                    for key, value in _results.items():
                            if key not in results_3d:
                                results_3d[key] = []
                            results_3d[key].append(value)    

    for k, v in results_3d.items():
        results_3d[k] = v

    return_results.append(results_3d)
    return_print.append("!!!!EOF!!!!")
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

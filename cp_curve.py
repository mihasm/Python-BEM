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

import numpy
from mpl_toolkits.mplot3d import Axes3D

from induction_factors import Calculator

a = Axes3D  # only for passing code inspection -> Axes3D needs to be imported


def calculate_power(speed_wind, rpm, sections_radius, chord_lengths, chord_angles, dr, R, Rhub, B, f_c_L, f_c_D,
                    add_angle=None, print_out=False, tip_loss=False, hub_loss=False, new_tip_loss=False,
                    new_hub_loss=False, cascade_correction=False, max_iterations=100, convergence_limit=0.001,
                    rho=1.225, method=10, relaxation_factor=0.3,print_all=False):
    """
    Returns calculated power using BEM analysis.

    Inputs are wind speed, rotational velocity, blade geometry, number of blades, and
    functions for calculating lift and drag coefficients.

    Output is a dictionary with all results.

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
    if add_angle != None:
        chord_angles = chord_angles + add_angle
    results = Calculator(f_c_L, f_c_D).run_array(chord_angle=chord_angles, B=B, c=chord_lengths, r=sections_radius,
                                                 dr=dr, rpm=rpm, v=speed_wind, R=R, Rhub=Rhub, print_out=print_out,
                                                 tip_loss=tip_loss, hub_loss=hub_loss, new_tip_loss=new_tip_loss,
                                                 new_hub_loss=new_hub_loss, cascade_correction=cascade_correction,
                                                 max_iterations=max_iterations, convergence_limit=convergence_limit,
                                                 rho=rho, method=method, relaxation_factor=relaxation_factor, print_all=print_all)
    return results


def calculate_power_3d(r, c, theta, dr, R, Rhub, B, f_c_L, f_c_D, add_angle=0.0,
                       print_out=False, tip_loss=False, hub_loss=False, new_tip_loss=False,
                       new_hub_loss=False, cascade_correction=False, max_iterations=100, convergence_limit=0.001,
                       rho=1.225, method=10, v_min=3, v_max=20, v_num=10, rpm_min=100, rpm_max=3000,
                       rpm_num=10, relaxation_factor=0.3, print_all=False, *args, **kwargs):
    """
    Calculates power for given geometry and data for every windspeed and rpm.

    Returns dictionary with arrays with data for every point.

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
    print("--------- RUNNING CALCULATION FOR FOLLOWING PARAMETERS --------")
    print("r", r)
    print("c", c)
    print("theta", theta)
    print("dr", dr)
    print("R", R)
    print("Rhub", Rhub)
    print("B", B)
    print("print_out", print_out)
    print("tip_loss", tip_loss)
    print("hub_loss", hub_loss)
    print("new_tip_loss", new_tip_loss)
    print("new_hub_loss", new_hub_loss)
    print("cascade_correction", cascade_correction)
    print("max_iterations", max_iterations)
    print("convergence_limit", convergence_limit)
    print("rho", rho)
    print("method", method)
    print("v_min", v_min)
    print("v_max", v_max)
    print("v_num", v_num)
    print("rpm_min", rpm_min)
    print("rpm_max", rpm_max)
    print("rpm_num", rpm_num)
    print("---------------------------------------------------------------")
    results_3d = {}

    for v in list(numpy.linspace(start=v_min, stop=v_max, num=v_num)):
        for rpm in list(numpy.linspace(start=rpm_min, stop=rpm_max, num=rpm_num)):
            print("Calculating power for v", v, "rpm", rpm)
            _results = calculate_power(v, rpm, r, c, theta, dr, add_angle=add_angle,
                                       R=R, B=B, f_c_L=f_c_L, f_c_D=f_c_D, Rhub=Rhub, print_out=print_out,
                                       tip_loss=tip_loss, hub_loss=hub_loss, new_tip_loss=new_tip_loss,
                                       new_hub_loss=new_hub_loss, cascade_correction=cascade_correction,
                                       max_iterations=max_iterations, convergence_limit=convergence_limit,
                                       rho=rho, method=method, relaxation_factor=relaxation_factor,print_all=print_all)
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

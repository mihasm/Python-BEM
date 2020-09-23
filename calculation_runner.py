import datetime
from math import pi
import sys
import time
import traceback

from calculation import Calculator
from mpl_toolkits.mplot3d import Axes3D
import numpy
from utils import Printer, generate_v_and_rpm_from_tsr, generate_v_and_rpm_from_J


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
        if f not in inp_args["airfoils"].keys():
            if f != "transition":
                p.print("Section foil %s does not exist in airfoil list." % f)
                raise Exception("Section foil not matching airfoil list error")

    try:
        c = Calculator(inp_args)
        results = c.run_array(**inp_args)
        return results
    except:
        inp_args["EOF"].value = True
        traceback.print_exc(file=sys.stdout)
        raise


def calculate_power_3d(inp_args, print_eof=False, prepend="", print_progress=True):
    """
    Calculates power for given geometry and data for every windspeed and rpm.

    Returns dictionary with arrays with data for every point.
    :return: dictionary with all results stored as numpy arrays
    """

    p = Printer(inp_args["return_print"])
    inp_args["print_progress"] = print_progress
    return_results = inp_args["return_results"]
    results_3d = {}

    # get parameters
    pitches = list(numpy.linspace(start=inp_args["pitch_min"], stop=inp_args["pitch_max"], num=int(inp_args["pitch_num"])))
    tsr_list = list(numpy.linspace(start=inp_args["tsr_min"], stop=inp_args["tsr_max"], num=int(inp_args["tsr_num"])))
    j_list = list(numpy.linspace(start=inp_args["J_min"], stop=inp_args["J_max"], num=int(inp_args["J_num"])))
    constant_speed, constant_rpm, constant_pitch = inp_args["constant_speed"], inp_args["constant_rpm"], inp_args["pitch"]
    variable_selection = inp_args["variable_selection"]
    constant_selection = inp_args["constant_selection"]
    R = inp_args["R"]

    if variable_selection == 0:
        speeds = list(numpy.linspace(start=inp_args["v_min"], stop=inp_args["v_max"], num=int(inp_args["v_num"])))
        rpms = list(numpy.linspace(start=inp_args["rpm_min"], stop=inp_args["rpm_max"], num=int(inp_args["rpm_num"])))
        pitches = [constant_pitch]

    elif variable_selection == 1:
        #TSR mode
        if constant_selection == 0:
            #constant speed, so change constant rpm to None
            constant_rpm = None
        else:
            constant_speed=None
        speeds, rpms = generate_v_and_rpm_from_tsr(tsr_list=tsr_list,R=R,v=constant_speed,rpm=constant_rpm)
        pitches = [constant_pitch]

    elif variable_selection == 2:
        #J mode
        if constant_selection == 0:
            #constant speed, so change constant rpm to None
            constant_rpm = None
        else:
            constant_speed=None
        speeds, rpms = generate_v_and_rpm_from_J(J_list=j_list,R=R,v=constant_speed,rpm=constant_rpm)
        pitches = [constant_pitch]

    elif variable_selection == 3:
        #pitches mode
        speeds, rpms = [constant_speed],[constant_rpm]

    total_iterations = int(len(speeds) * len(rpms))
    i = 0

    time_start = time.time()
    for v in speeds:
        for rpm in rpms:
            for pitch in pitches:
                if print_progress:
                    _lambda = rpm / 60 * 2 * pi * inp_args["R"] / v
                    _advance_ratio = v / (rpm / 60 * inp_args["R"] * 2)
                    p.print(prepend + "v=%.1f m/s, n=%.0f RPM, λ=%.2f, J=%.2f" % (v,rpm,_lambda,_advance_ratio))
                _inp_args = {**inp_args}
                _inp_args["v"] = v
                _inp_args["rpm"] = rpm
                _inp_args["pitch"] = pitch
                try:
                    _results = calculate_power(_inp_args)
                except Exception as e:
                    var = traceback.format_exc()
                    p.print("Error in running optimizer: %s \n %s" % (str(e),var))
                    inp_args["EOF"].value = True
                    raise

                if _results != None and _results["power"]:
                    if print_progress:
                        p.print(prepend + "    TSR:", _results["TSR"], "J:", _results["J"], "cp:", _results["cp"],
                                "ct:", _results["ct"])
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
                #p.print("    ### Time left:", t_left_str, "ETA:", eta, "###")
                if print_progress:
                    p.print("")

    return_results.append(results_3d)
    if print_eof:
        inp_args["EOF"].value = True
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

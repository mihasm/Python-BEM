import datetime
import time
import traceback
from math import pi

import numpy
from mpl_toolkits.mplot3d import Axes3D

from calculation import Calculator
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
    except Exception as e:
        var = traceback.format_exc()
        p.print("Error in running analysis: %s \n %s" % (str(e), var))
        inp_args["EOF"].value = True
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

    try:
        # get parameters
        pitches = list(
            numpy.linspace(start=inp_args["pitch_min"], stop=inp_args["pitch_max"], num=int(inp_args["pitch_num"])))
        tsr_list = list(
            numpy.linspace(start=inp_args["tsr_min"], stop=inp_args["tsr_max"], num=int(inp_args["tsr_num"])))
        j_list = list(numpy.linspace(start=inp_args["J_min"], stop=inp_args["J_max"], num=int(inp_args["J_num"])))

        constant_speed, constant_rpm, constant_pitch = inp_args["constant_speed"], inp_args["constant_rpm"], inp_args[
            "pitch"]
        variable_selection = inp_args["variable_selection"]
        constant_selection = inp_args["constant_selection"]
        geometry_scale = inp_args["geometry_scale"]
        R = inp_args["R"]

        if variable_selection == 0:
            speeds = list(numpy.linspace(start=inp_args["v_min"], stop=inp_args["v_max"], num=int(inp_args["v_num"])))
            rpms = list(
                numpy.linspace(start=inp_args["rpm_min"], stop=inp_args["rpm_max"], num=int(inp_args["rpm_num"])))
            pitches = [constant_pitch]

        elif variable_selection == 1:
            # TSR mode
            if constant_selection == 0:
                # constant speed, so change constant rpm to None
                constant_rpm = None
            else:
                constant_speed = None
            speeds, rpms = generate_v_and_rpm_from_tsr(tsr_list=tsr_list, R=R, geometry_scale=geometry_scale,
                                                       v=constant_speed, rpm=constant_rpm)
            pitches = [constant_pitch]

        elif variable_selection == 2:
            # J mode
            if constant_selection == 0:
                # constant speed, so change constant rpm to None
                constant_rpm = None
            else:
                constant_speed = None
            speeds, rpms = generate_v_and_rpm_from_J(J_list=j_list, R=R, geometry_scale=geometry_scale,
                                                     v=constant_speed, rpm=constant_rpm, printer=p)
            pitches = [constant_pitch]

        elif variable_selection == 3:
            # pitches mode
            speeds, rpms = [constant_speed], [constant_rpm]

        elif variable_selection == 4:
            # pitch + TSR
            if constant_selection == 0:
                # constant speed, so change constant rpm to None
                constant_rpm = None
            else:
                constant_speed = None
            speeds, rpms = generate_v_and_rpm_from_tsr(tsr_list=tsr_list, R=R, geometry_scale=geometry_scale,
                                                       v=constant_speed, rpm=constant_rpm)

        elif variable_selection == 5:
            # pitch + J
            if constant_selection == 0:
                # constant speed, so change constant rpm to None
                constant_rpm = None
            else:
                constant_speed = None
            speeds, rpms = generate_v_and_rpm_from_J(J_list=j_list, R=R, geometry_scale=geometry_scale,
                                                     v=constant_speed, rpm=constant_rpm, printer=p)

        total_iterations = int(len(speeds) * len(rpms))

        i = 0

        pitch_change_list = []

        time_start = time.time()
        for pitch in pitches:
            p.print("Pitch:", pitch)
            pitch_change_list.append(i)
            for v in speeds:
                for rpm in rpms:
                    print_progress_message(v, rpm, inp_args, p, prepend, print_progress)

                    _inp_args = {**inp_args, "v": v, "rpm": rpm, "pitch": pitch}
                    _results = calculate_power(_inp_args)

                    # if results are valid, add them to results list
                    if _results != None and _results["power"]:
                        print_result_message(print_progress, p, prepend, _results)

                        # append the value of the _results to the results_3d list
                        for key, value in _results.items():
                            if key not in results_3d:
                                results_3d[key] = []
                            results_3d[key].append(value)

                    i += 1
                    eta = process_time(time_start, i, total_iterations)
                    # p.print("    ### Time left:", t_left_str, "ETA:", eta, "###")
                    if print_progress:
                        p.print("")

        results_3d["pitch_change_list"] = pitch_change_list
        return_results.append(results_3d)

        if print_eof:
            inp_args["EOF"].value = True

        return results_3d

    except Exception as e:
        var = traceback.format_exc()
        p.print("Error in running analysis: %s \n %s" % (str(e), var))
        inp_args["EOF"].value = True
        raise


def process_time(time_start, i, total_iterations):
    """

    :param time_start:
    :param i:
    :param total_iterations:
    """
    t_now = int(time.time() - time_start)
    t_left = int((total_iterations / i - 1) * t_now)
    t_left_str = str(datetime.timedelta(seconds=t_left))
    eta_seconds = datetime.datetime.now() + datetime.timedelta(seconds=t_left)
    eta = str(eta_seconds).split(".")[0]


def print_progress_message(v, rpm, inp_args, p, prepend, print_progress):
    """

    :param v:
    :param rpm:
    :param inp_args:
    :param p:
    :param prepend:
    :param print_progress:
    """
    if print_progress:
        if v > 0:
            _lambda = rpm / 60 * 2 * pi * inp_args["R"] * inp_args["geometry_scale"] / v
        else:
            _lambda = 0.0
        _advance_ratio = v / (rpm / 60 * 2 * inp_args["R"] * inp_args["geometry_scale"])
        # pitch = inp_args["pitch"]
        p.print(prepend + "v=%.1f m/s, n=%.0f RPM, λ=%.2f, J=%.2f" % (v, rpm, _lambda, _advance_ratio))


def print_result_message(print_progress, p, prepend, _results):
    """

    :param print_progress:
    :param p:
    :param prepend:
    :param _results:
    """
    if print_progress:
        p.print(prepend + "    cp:", _results["cp"], "ct:", _results["ct"], "eff:", _results["eff"])


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

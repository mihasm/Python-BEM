from turbine_data import SET_INIT
from calculation import Calculator
from utils import Printer
from numpy import radians, degrees
import numpy as np
import numpy
from math import pi, exp
import traceback
from matplotlib import pyplot as plt
from scipy.signal import argrelextrema


def optimize_angles_genetic(inp_args):

    # brute force method

    p = Printer(inp_args["return_print"])
    try:
        p.print("Optimizing angles for target variable:",
                inp_args["optimization_variable"])
        p.print("Using genetic algorithm")
        return_results = inp_args["return_results"]

        v = inp_args["target_speed"]
        rpm = inp_args["target_rpm"]
        # print(v,pi,rpm)
        omega = 2 * pi * rpm / 60
        # optimization_variable = "dT"
        optimization_variable = inp_args["optimization_variable"]
        #optimization_variable = "dQ"
        p.print("Optimization variable is", optimization_variable)
        p.print("Propeller mode:", inp_args["propeller_mode"])

        output_angles = []
        output_alphas = []

        inp_args["theta_in"] = np.array([120]*len(inp_args["theta_in"]))

        C = Calculator(inp_args["airfoils"])

        p.print("Input section radiuses:")
        for _r in inp_args["r_in"]:
            p.print(_r)

        p.print("Starting calculation...")
        for section_number in range(len(inp_args["r_in"])):
            p.print("  Section_number is", section_number)

            _r = inp_args["r_in"][section_number]
            _c = inp_args["c_in"][section_number]
            _theta = inp_args["theta_in"][section_number]
            _dr = inp_args["dr"][section_number]
            _airfoil = inp_args["foils_in"][section_number]
            max_thickness = inp_args["airfoils"][_airfoil]["max_thickness"] * _c
            _airfoil_dat = _airfoil + ".dat"

            def fobj(x):
                global dT_max
                global dQ_max
                d = C.calculate_section(v=v, omega=omega, _airfoil=_airfoil, _airfoil_dat=_airfoil_dat,
                                        max_thickness=max_thickness, _r=_r, _c=_c, _dr=_dr, _theta=radians(
                                            x),
                                        printer=p, **inp_args)
                if d == None or d == False:
                    p.print("none")
                    return -1e10

                if optimization_variable == "max dT min dQ":
                    return d["dQ"]/d["dT"]

                return d[optimization_variable]

            it = list(de2(fobj, bounds=[(-10, 45)], printer=p))

            d_final = C.calculate_section(v=v, omega=omega, _airfoil=_airfoil, _airfoil_dat=_airfoil_dat,
                                          max_thickness=max_thickness, _r=_r, _c=_c, _dr=_dr, _theta=radians(
                                              it[-1]),
                                          printer=p, **inp_args)

            output_angles.append(it[-1])
            output_alphas.append(d_final["alpha"])

            p.print("    final theta is", it[-1])
        p.print("Final angles:")
        for a in output_angles:
            p.print(a)

        p.print("!!!!EOF!!!!")
    except Exception as e:
        p.print("Error in running optimizer: %s" % str(e))
        p.print("!!!!EOF!!!!")
        raise


#https://pablormier.github.io/2017/09/05/a-tutorial-on-differential-evolution-with-python/#

def de2(function, bounds, M=0.8, num_individuals=30, iterations=50, printer=None):
    """
    Function that uses the DE genetic optimisation algorithm to maximize
    the fitness function "function".

    Inputs:
    function:function: Python function that is used to determine
        the fitness score of each indicidual. Can be one or N-dimensional.
    bounds:list: List of tuple pairs that correspond to the upper and
        lower boundary for each dimension information.
    M:float: Mutation coefficient.
    num_individuals:int: Number of individuals in the calculation.
    iterations:int: Maximum number of iterations per calculation.

    Output:
    return:float:Individual with the highest fitness.
    """
    p = printer
    dimensions = len(bounds)
    min_bound, max_bound = np.asarray(bounds).T
    population = np.random.uniform(
        min_bound, max_bound, (num_individuals, dimensions))
    fitness = np.asarray([function(p) for p in population])
    best_i = np.argmax(fitness)
    best = population[best_i]
    for i in range(iterations):
        for j in range(num_individuals):
            other_i = list(set(range(num_individuals))-set([j]))
            a, b, c = population[np.random.choice(other_i, 3,
                                                  replace=False)]
            mutation_vector = a+M*(b-c)
            for k in range(dimensions):
                if mutation_vector[k] < min_bound[k]:
                    mutation_vector[k] = min_bound[k]
                if mutation_vector[k] > max_bound[k]:
                    mutation_vector[k] = max_bound[k]
            random_locations = np.random.choice(a=[False, True],
                                                size=(1, dimensions))[0]
            trial = np.where(random_locations, mutation_vector, population[j])
            f = function(trial)
            if f > fitness[j]:
                fitness[j] = f
                population[j] = trial
                if f > fitness[best_i]:
                    best_i = j
                    best = trial
        p.print(best)
    return best


def optimal_pitch(inp_args):
    p = Printer(inp_args["return_print"])
    try:
        p.print("Optimizing angles for target variable:",
                inp_args["optimization_variable"])
        return_results = inp_args["return_results"]

        v = inp_args["target_speed"]
        rpm = inp_args["target_rpm"]
        omega = 2 * pi * rpm / 60
        optimization_variable = inp_args["optimization_variable"]

        output_angles = []

        C = Calculator(inp_args["airfoils"])

        done_pitches = {}
        _pitch = 0  # is in radians

        args = {**inp_args}

        for dpitch in [45, 30, 20, 10, 5, 2, 1, 0.1]:
            while True:
                # middle angle
                if not _pitch in done_pitches:
                    p.print("   calculating out (pitch:%s)" % _pitch)
                    try:
                        args["pitch"] = _pitch
                        out = C.run_array(**args, rpm=rpm, v=v)
                    except Exception as e:
                        p.print(e)
                        p.print(traceback.format_exc())
                        out = None
                    if out == False or out == None:
                        break
                    done_pitches[_pitch] = out
                else:
                    p.print(_pitch, "already calculated, reusing...")
                    out = done_pitches[_pitch]

                _pitch_up = _pitch + dpitch
                # upper angle
                if not _pitch_up in done_pitches:
                    p.print("   calculating out_up (pitch:%s)" % _pitch_up)
                    try:
                        args["pitch"] = _pitch_up
                        out_up = C.run_array(**args, rpm=rpm, v=v)
                    except Exception as e:
                        p.print()
                        p.print(traceback.format_exc())
                        out_up = None
                    if out_up == False or out_up == None:
                        break
                    done_pitches[_pitch_up] = out_up
                else:
                    p.print(_pitch_up, "already calculated, reusing...")
                    out_up = done_pitches[_pitch_up]

                _pitch_down = _pitch - dpitch
                # lower angle
                if not _pitch_down in done_pitches:
                    p.print("   calculating out_down (pitch:%s)" % _pitch_down)
                    try:
                        args["pitch"] = _pitch_down
                        out_down = C.run_array(**args, rpm=rpm, v=v)
                    except Exception as e:
                        p.print(e)
                        p.print(traceback.format_exc())
                        out_down = None
                    if out_down == False or out_down == None:
                        break
                    done_pitches[_pitch_down] = out_down
                else:
                    p.print(_pitch_down, "already calculated, reusing...")
                    out_down = done_pitches[_pitch_down]

                if out_up == False or out_down == False or out == False:
                    p.print("   one is False, breaking...")
                    break
                if out_up == None or out_down == None or out == None:
                    p.print("   one is None, breaking...")
                    break

                var = np.sum(out[optimization_variable])
                var_up = np.sum(out_up[optimization_variable])
                var_down = np.sum(out_down[optimization_variable])
                p.print("   %s" % optimization_variable, var, "%s_up" % optimization_variable, var_up,
                        "%s_down" % optimization_variable, var_down)

                if var_up <= var and var_down <= var:
                    p.print("   none is bigger, breaking...")
                    break

                if var_up > var > var_down:
                    p.print("   going up")
                    _pitch = _pitch_up
                    var = var_up
                    out = out_up

                if var_down > var > var_up:
                    p.print("   going down")
                    _pitch = _pitch_down
                    var = var_down
                    out = out_down

                if var_down > var and var_up > var:
                    if var_up > var_down:
                        p.print("   both up and down are bigger, going up")
                        _pitch = _pitch_up
                        var = var_up
                        out = out_up
                    else:
                        p.print("   both up and down are smaller, going down")
                        _pitch = _pitch_down
                        var = var_down
                        out = out_down
                if var_up == var and var_down == var:
                    p.print("   both are equal, breaking")
                    break
        p.print("Final pitch:", _pitch)
        p.print("Angles:")
        for t in args["theta"]:
            p.print(t+_pitch)
        p.print("!!!!EOF!!!!")
    except:
        p.print("Error in running optimizer")
        p.print("!!!!EOF!!!!")
        raise
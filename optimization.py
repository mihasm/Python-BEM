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
import traceback

#for display
from pyqtgraph.Qt import QtGui, QtCore
import numpy as np
import pyqtgraph as pg
from multiprocessing import Process, Manager, Queue
import sched, time, threading
import time



def optimize_angles_genetic(inp_args,queue_pyqtgraph):

    # brute force method

    p = Printer(inp_args["return_print"])
    try:
        return_results = inp_args["return_results"]

        v = inp_args["target_speed"]
        rpm = inp_args["target_rpm"]
        omega = 2 * pi * rpm / 60

        optimization_variable = inp_args["optimization_variable"]
        pitch_optimization = inp_args["pitch_optimization"]
        min_bound = inp_args["min_bound"]
        max_bound = inp_args["max_bound"]
        mut_coeff = inp_args["mut_coeff"]
        population_size = int(inp_args["population"])
        num_iter = int(inp_args["num_iter"])
        weight_dt = inp_args["weight_dt"]
        weight_dq = inp_args["weight_dq"]

        p.print("Optimization variable is", optimization_variable)
        p.print("Pitch optimization:",pitch_optimization)
        p.print("Propeller mode:", inp_args["propeller_mode"])

        output_theta = []
        output_chord = []

        inp_args["theta_in"] = np.array([120]*len(inp_args["theta_in"]))

        C = Calculator(inp_args["airfoils"])

        p.print("Input section radiuses:")
        for _r in inp_args["r_in"]:
            p.print(_r)

        if pitch_optimization:
            if optimization_variable == "dQ":
                optimization_variable = "Msum"
            elif optimization_variable == "dT":
                optimization_variable = "thrust"
            else:
                p.print("This is not implemented yet....")
                inp_args["EOF"].value = True
                return

            p.print("Optimization variable is (pitch mode)", optimization_variable)

            del inp_args["pitch"]


            def fobj(x):
                out = C.run_array(**inp_args, rpm=rpm, v=v, pitch=x)
                if out == None or out == False:
                    return 0.0
                return out[optimization_variable]

            results = de2(fobj, bounds=[(min_bound, max_bound)], M=mut_coeff, num_individuals=population_size,iterations=num_iter,printer=p, queue=queue_pyqtgraph)
            it = list(results)
            best = it[-1]
        else:
            p.print("Starting calculation...")
            for section_number in range(len(inp_args["r"])):
                p.print("  Section_number is", section_number)

                _r = inp_args["r"][section_number]
                _c = inp_args["c"][section_number]
                _theta = inp_args["theta"][section_number]
                _dr = inp_args["dr"][section_number]
                _airfoil = inp_args["foils"][section_number]
                max_thickness = inp_args["airfoils"][_airfoil]["max_thickness"] * _c
                _airfoil_dat = _airfoil + ".dat"

                #worst_value = 0.0

                def fobj(x):
                    global dT_max
                    global dQ_max
                    d = C.calculate_section(v=v, omega=omega, _airfoil=_airfoil, _airfoil_dat=_airfoil_dat,
                                            max_thickness=max_thickness, _r=_r, _c=_c, _dr=_dr, _theta=radians(
                                                x),
                                            printer=p, **inp_args)
                    if d == None or d == False:
                        return -1e50

                    if optimization_variable == "dQ-dT":
                        return d["dQ"]*weight_dq-d["dT"]*weight_dt
                    if optimization_variable == "dT-dQ":
                        return d["dT"]*weight_dt-d["dQ"]*weight_dq
                    
                    return d[optimization_variable]

                it = list(de2(fobj, bounds=[(min_bound, max_bound)], M=mut_coeff, num_individuals=population_size,iterations=num_iter,printer=p, queue=queue_pyqtgraph))
                
                d_final = C.calculate_section(v=v, omega=omega, _airfoil=_airfoil, _airfoil_dat=_airfoil_dat,
                                              max_thickness=max_thickness, _r=_r, _c=_c, _dr=_dr, _theta=radians(
                                                  it[0]),
                                              printer=p, **inp_args)

                output_theta.append(it[0])
            
            p.print("Number of final angles:",len(output_theta))
            p.print("Final angles:")
            for a in output_theta:
                p.print(a)

            #p.print("Final chords:")
            #for c in output_chord:
            #    p.print(c)
        p.print("Done!")
        time.sleep(0.5)
        inp_args["EOF"].value = True
    except Exception as e:
        var = traceback.format_exc()
        p.print("Error in running optimizer: %s \n %s" % (str(e),var))
        inp_args["EOF"].value = True
        raise


#https://pablormier.github.io/2017/09/05/a-tutorial-on-differential-evolution-with-python/#

def de2(function, bounds, M=0.8, num_individuals=30, iterations=50, printer=None, queue=None):
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
            queue[0]=[population.flatten(),fitness,population[best_i][0],fitness[best_i]]
        if len(set(population.flatten())) == 1:
            break
        p.print(best,fitness[best_i])


    return best
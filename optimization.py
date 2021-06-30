import time
import traceback
from math import pi

import numpy as np
from numpy import radians

from calculation import Calculator
from utils import Printer

import scipy.optimize


def optimize(inp_args, queue_pyqtgraph):
    """

    :param inp_args:
    :param queue_pyqtgraph:
    :return:
    """
    p = Printer(inp_args["return_print"])
    try:
        return_results = inp_args["return_results"]

        inp_args["v"] = inp_args["target_speed"]
        inp_args["rpm"] = inp_args["target_rpm"]
        inp_args["omega"] = 2 * pi * inp_args["rpm"] / 60

        num_sections = len(inp_args["theta"])

        input_variables = inp_args["optimization_inputs"]
        output_variables = inp_args["optimization_outputs"]
        target_variables = inp_args["optimization_targets"]
        
        mut_coeff = inp_args["mut_coeff"]
        population_size = int(inp_args["population"])
        num_iter = int(inp_args["num_iter"])

        p.print("Input variables:", input_variables)
        p.print("Output variables:", output_variables)
        p.print("Target_variables:", target_variables)

        p.print("Propeller mode:", inp_args["propeller_mode"])

        C = Calculator(inp_args)

        output_list = []

        mode = 0

        if mode == 0:

            for n in range(len(inp_args["r"])):
                p.print("  Section_number is", n)

                list_queue_internal_x = []
                list_queue_internal_y = []

                section_inp_args = {}

                section_inp_args["_r"] = inp_args["r"][n]
                section_inp_args["_c"] = inp_args["c"][n]
                section_inp_args["_theta"] = inp_args["theta"][n]
                section_inp_args["_dr"] = inp_args["dr"][n]

                section_inp_args["transition"] = C.transition_array[n]
                section_inp_args["_airfoil"] = C.airfoils_list[n]
                section_inp_args["_airfoil_prev"] = C.transition_foils[n][0]
                section_inp_args["_airfoil_next"] = C.transition_foils[n][1]
                section_inp_args["transition_coefficient"] = C.transition_foils[n][2]
                section_inp_args["max_thickness"] = C.max_thickness_array[n]

                # worst_value = 0.0

                inputs_list = [i[0] for i in input_variables]
                bounds_list = [(b[1],b[2]) for b in input_variables]

                def fobj(input_numbers):
                    for i in range(len(input_numbers)):
                        section_inp_args[inputs_list[i]]=input_numbers[i]

                    args = {**inp_args,**section_inp_args}

                    d = C.calculate_section(**args,printer=p)
                    if d == None or d == False:
                        return 1e10

                    fitness = 0

                    for var, coeff in output_variables:
                        value = d[var]*coeff
                        fitness -= value
                    for var, target_value, coeff in target_variables:
                        comparison = abs(d[var]-target_value)*coeff
                        fitness += comparison

                    #p.print(queue_pyqtgraph)

                    list_queue_internal_x.append(args["_theta"])
                    list_queue_internal_y.append(fitness)

                    queue_pyqtgraph[0] = [list_queue_internal_x,list_queue_internal_y,args["_theta"],fitness]

                    #p.print(args["_theta"],fitness)
                    #p.print(queue_pyqtgraph)

                    return fitness

                decreasing = {"_theta"}

                if n > 0:
                    # use previous iteration to set new bound
                    for _vname in decreasing:
                        index = inputs_list.index(_vname)
                        bounds_list[index] = (bounds_list[index][0],output_list[n-1][index]) #construct new tuple
                
                p.print("Bounds:",bounds_list)

                initial_guess = [top for bottom,top in bounds_list] # guess the top guess

                #it = list(de2(fobj, bounds=bounds_list, iterations=num_iter, M=mut_coeff, num_individuals=population_size, printer=p, queue=queue_pyqtgraph))
                it = list(scipy.optimize.minimize(fobj,initial_guess,method="powell",bounds=bounds_list).x)

                p.print("best combination",it)
                output_list.append(it)

            p.print("Final output:")
            p.print([v[0] for v in output_variables])
            for i in output_list:
                if len(i) > 1:
                    p.print(i)
                else:
                    p.print(i[0])
            
            p.print("Done!")
            time.sleep(0.5)
            inp_args["EOF"].value = True

        elif mode == 1:
            inputs_list = [i[0] for i in input_variables]
            bounds_list = []
            for b in input_variables:
                for i in range(num_sections):
                    bounds_list.append((b[1],b[2]))

            p.print(bounds_list)

            inputs_list = [["theta",-45,45]]
            output_variables = [["cp",1.0]]

            def fobj(*input_numbers):
                input_numbers = input_numbers[0]
                k = 0
                for i in range(len(inputs_list)):
                    for j in range(num_sections):
                        inp_args[inputs_list[i][0]][j]=input_numbers[k]
                        k+=1


                d = C.run_array(**inp_args,printer=p)

                if d == None or d == False:
                    p.print("d is None or False")
                    return -1e50

                fitness = 0
                for var, coeff in output_variables:
                    fitness += d[var]*coeff
                return fitness

            #it = list(de2(fobj, bounds=bounds_list, iterations=num_iter, M=mut_coeff, num_individuals=population_size, printer=p, queue=queue_pyqtgraph))
            it = list(scipy.optimize.differential_evolution(fobj,bounds=bounds_list,maxiter=num_iter,popsize=num_individuals))

    except Exception as e:
        var = traceback.format_exc()
        p.print("Error in running optimizer: %s \n %s" % (str(e), var))
        inp_args["EOF"].value = True
        raise
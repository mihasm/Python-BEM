import time
import traceback
from math import pi

import numpy as np
from numpy import radians

from calculation import Calculator
from utils import Printer


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
        
        mut_coeff = inp_args["mut_coeff"]
        population_size = int(inp_args["population"])
        num_iter = int(inp_args["num_iter"])

        p.print("Input variables:", input_variables)
        p.print("Output variables:", output_variables)
        p.print("Propeller mode:", inp_args["propeller_mode"])

        C = Calculator(inp_args)

        output_list = []

        mode = 0

        if mode == 0:

            for n in range(len(inp_args["r"])):
                p.print("  Section_number is", n)

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
                        p.print("d is None or False")
                        return -1e50

                    fitness = 0
                    for var, coeff in output_variables:
                        fitness += d[var]*coeff
                    return fitness

                decreasing = {"_theta"}

                if n > 0:
                    # use previous iteration to set new bound
                    for _vname in decreasing:
                        index = inputs_list.index(_vname)
                        p.print("boundslist",bounds_list[index][1])
                        p.print("output_list",output_list[n-1][index])
                        bounds_list[index] = (bounds_list[index][0],output_list[n-1][index]) #construct new tuple
                
                p.print("Bounds:",bounds_list)

                it = list(de2(fobj, bounds=bounds_list, iterations=num_iter, M=mut_coeff, num_individuals=population_size, printer=p, queue=queue_pyqtgraph))

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

            it = list(de2(fobj, bounds=bounds_list, iterations=num_iter, M=mut_coeff, num_individuals=population_size, printer=p, queue=queue_pyqtgraph))


    except Exception as e:
        var = traceback.format_exc()
        p.print("Error in running optimizer: %s \n %s" % (str(e), var))
        inp_args["EOF"].value = True
        raise


# https://pablormier.github.io/2017/09/05/a-tutorial-on-differential-evolution-with-python/#

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
    fitness = np.asarray([function(po) for po in population])
    best_i = np.argmax(fitness)
    best = population[best_i]
    for i in range(iterations):
        for j in range(num_individuals):
            other_i = list(set(range(num_individuals)) - {j})
            a, b, c = population[np.random.choice(other_i, 3, replace=False)]
            mutation_vector = a + M * (b - c)
            for k in range(dimensions):
                if mutation_vector[k] < min_bound[k]:
                    mutation_vector[k] = min_bound[k]
                if mutation_vector[k] > max_bound[k]:
                    mutation_vector[k] = max_bound[k]
            random_locations = np.random.choice(a=[False, True], size=(1, dimensions))[0]
            trial = np.where(random_locations, mutation_vector, population[j])
            f = function(trial)
            if f > fitness[j]:
                fitness[j] = f
                population[j] = trial
                if f > fitness[best_i]:
                    best_i = j
                    best = trial
            if dimensions == 1:
                queue[0] = [population.flatten(), fitness.flatten(), population[best_i][0], fitness[best_i]]
            else:
                queue[0] = [np.arange(0,len(population)),fitness.flatten(),best_i,fitness[best_i]]
        if len(set(population.flatten())) == 1:
            break
        p.print(best, fitness[best_i])

    return best

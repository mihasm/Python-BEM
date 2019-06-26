import numpy as np
from matplotlib import pyplot as plt
from turbine_data import SET_INIT
from utils import Printer
from numpy import pi, radians
from calculation import Calculator




# noinspection PyBroadException
def optimize_angles_genetic(inp_args):

    #brute force method

    p = Printer(inp_args["return_print"])
    try:
        p.print("Optimizing angles for target variable:", inp_args["optimization_variable"])
        return_results = inp_args["return_results"]

        v = inp_args["target_speed"]
        rpm = inp_args["target_rpm"]
        # print(v,pi,rpm)
        omega = 2 * pi * rpm / 60
        # optimization_variable = "dT"
        #optimization_variable = inp_args["optimization_variable"]
        optimization_variable = "dQ"
        p.print("Optimization variable is",optimization_variable)
        p.print("Propeller mode:",inp_args["propeller_mode"])

        output_angles = []
        output_alphas = []

        inp_args["theta_in"] = np.array([120]*len(inp_args["theta_in"]))

        C = Calculator(inp_args["airfoils"])

        p.print("Input section radiuses:")
        for _r in inp_args["r_in"]:
            p.print(_r)

        p.print("Starting calculation...")
        for section_number in range(len(inp_args["r_in"])):
            p.print("section_number is", section_number)

            _r = inp_args["r_in"][section_number]
            _c = inp_args["c_in"][section_number]
            _theta = inp_args["theta_in"][section_number]
            _dr = inp_args["dr"][section_number]
            _airfoil = inp_args["foils_in"][section_number]
            max_thickness = inp_args["airfoils"][_airfoil]["max_thickness"] * _c
            _airfoil_dat = _airfoil + ".dat"

            
            #out = C.calculate_section(v=v, omega=omega, _airfoil=_airfoil, _airfoil_dat=_airfoil_dat,
            #      max_thickness=max_thickness, _r=_r, _c=_c, _dr=_dr, _theta=radians(_theta),
            #      printer=p, **inp_args)

            
            def fobj(x):
                value = 0
                for i in range(len(x)):
                    d =  C.calculate_section(v=v, omega=omega, _airfoil=_airfoil, _airfoil_dat=_airfoil_dat,
                        max_thickness=max_thickness, _r=_r, _c=_c, _dr=_dr, _theta=radians(x[i]),
                        printer=p, **inp_args)
                    var_opt = d[optimization_variable]
                    value += var_opt
                return value / len(x)

            it = list(de(fobj, bounds=[(-45, 45)]))
            #print(it)

            d_final = C.calculate_section(v=v, omega=omega, _airfoil=_airfoil, _airfoil_dat=_airfoil_dat,
                        max_thickness=max_thickness, _r=_r, _c=_c, _dr=_dr, _theta=radians(it[-1][0][0]),
                        printer=p, **inp_args)
            #print(it[-1])
            plt.plot(done_angles,done_thrusts,"b.")
            plt.show()

            output_angles.append(it[-1])
            output_alphas.append(d["alpha"])

            p.print("final theta is", it[-1])
            p.print("*******************************")

        p.print("Final angles:")
        for a in output_angles:
            print(a)

        p.print("!!!!EOF!!!!")
    except Exception as e:
        p.print("Error in running optimizer: %s" % str(e))
        p.print("!!!!EOF!!!!")
        raise


#https://pablormier.github.io/2017/09/05/a-tutorial-on-differential-evolution-with-python/#

def de(fobj, bounds, mut=0.8, crossp=0.3, num_individuals=100, its=100):
    #stevilo dimenzij
    dimensions = len(bounds)

    #populacija
    pop = np.random.rand(num_individuals, dimensions)

    #normalizacija
    min_bound, max_bound = np.asarray(bounds).T
    diff = np.fabs(min_bound - max_bound)
    pop_denorm = min_bound + pop * diff

    #zacetni fintes
    fitness = np.asarray([fobj(ind) for ind in pop_denorm])

    #najmocnejsi posameznik
    best_idx = np.argmax(fitness)

    #vrednost najmocnejsega posameznika
    best = pop_denorm[best_idx]

    for i in range(its):
        #print("i",i)
        for j in range(num_individuals):
            #print("\tj",j)
            #print("\tj_value",pop[j])

            #vsi posamezniki, ki niso trenutni posameznik
            idxs = [idx for idx in range(num_individuals) if idx != j]
            #print("\t\tidxs",idxs)

            #tri nakljucni ostali posamezniki
            a, b, c = pop[np.random.choice(idxs, 3, replace = False)]
            #print("\t\ta,b,c",a,b,c)

            #mutacija
            mutant = np.clip(a + mut * (b - c), 0, 1)
            #print("\t\tmutant",mutant)

            #izbira zamenjave
            cross_points = np.random.rand(dimensions) < crossp
            if not np.any(cross_points):
                cross_points[np.random.randint(0, dimensions)] = True
            cross_points = [False]
            #print("\t\tcross_points",cross_points)

            #poskusni posameznik
            trial = np.where(cross_points, mutant, pop[j])
            #print("\t\ttrial",trial)

            #normalizacija poskusnega posameznika
            trial_denorm = min_bound + trial * diff

            #preverjanje ciljne funkcije
            f = fobj(trial_denorm)
            if f > fitness[j]:
                #ce je novi fitnes boljsi, izberi novega posameznika
                fitness[j] = f
                pop[j] = trial
                if f > fitness[best_idx]:
                    best_idx = j
                    best = trial_denorm
        #yield best, fitness[best_idx]
    return best


#SET_INIT["return_print"] = []
#SET_INIT["return_results"] = []
#optimize_angles_genetic(SET_INIT)


def fobj(x):
    value = 0
    for i in range(len(x)):
        value += -x[i]**2
    return value / len(x)

#fobj = lambda x: sum(x**2)/len(x)

it = list(de(fobj, bounds=[(-100, 100)]))

print(it)
"""

#plt.plot(range(len(it)),np.array(it)[:,0])
#plt.show()
"""
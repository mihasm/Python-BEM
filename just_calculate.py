from calculation import Calculator
from numpy import linspace
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from calculation_runner import calculate_power_3d



def simple_calculate():
	from turbine_data import SET_INIT
	SET_INIT["convergence_limit"] = 1e-3
	SET_INIT["relaxation_factor"] = 0.1
	SET_INIT["cascade_correction"] = False
	SET_INIT["max_iterations"] = 1000
	SET_INIT["yaw_angle"] = 0.0
	SET_INIT["skewed_wake_correction"] = False
	SET_INIT["return_print"] = []
	SET_INIT["return_results"] = []
	results = calculate_power_3d(SET_INIT,print_progress=False)
	return results

simple_calculate()
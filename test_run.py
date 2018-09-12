__author__ = "Miha Smrekar"
__credits__ = ["Miha Smrekar"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Miha Smrekar"
__email__ = "miha.smrekar9@gmail.com"
__status__ = "Development"

import numpy
from scipy import interpolate
import numbers
from openpyxl import Workbook
import csv
from cp_curve import run, calculate_power, run_main
from parse_data import parse_podatki, parse_sections
from optimisation import Optimizer
import multiprocessing as mp

podatki = parse_podatki()

f_c_L = interpolate.interp1d(podatki['cL_x'], podatki['cL_y'], fill_value=(podatki['cL_y'][0], podatki['cL_y'][-1]),
                             bounds_error=False)
f_c_D = interpolate.interp1d(podatki['cD_x'], podatki['cD_y'], fill_value=(podatki['cD_y'][0], podatki['cD_y'][-1]),
                             bounds_error=False)

sections_radius, chord_lengths, chord_angles, dr = parse_sections()

run_main(sections_radius, chord_lengths, chord_angles, dr, B=5, R=0.776, Rhub=0.1, f_c_L=f_c_L, f_c_D=f_c_D)

# o =  Optimizer(sections_radius, chord_lengths, chord_angles, dr, B=5, R=0.776, Rhub=0.1, f_c_L=f_c_L, f_c_D=f_c_D)

# angles_new_pitch = o.optimize_pitch(7)
# run_main(sections_radius, chord_lengths, angles_new_pitch, dr, B=5, R=0.776, Rhub=0.1, f_c_L=f_c_L, f_c_D=f_c_D)

# angles_optimized = o.optimize_angles(7,300)
# run_main(sections_radius, chord_lengths, angles_optimized, dr, B=5, R=0.776, Rhub=0.1, f_c_L=f_c_L, f_c_D=f_c_D)


# PODATKI NORVEŽAN
r = numpy.array(
    [0.055, 0.068, 0.083, 0.098, 0.113, 0.128, 0.143, 0.158, 0.173, 0.188, 0.203, 0.218, 0.233, 0.248, 0.278, 0.293,
     0.308, 0.323, 0.338, 0.368, 0.383, 0.398, 0.413, 0.428, 0.443, 0.450])
theta = numpy.array(
    [38.000, 37.055, 32.544, 28.677, 25.262, 22.430, 19.988, 18.034, 16.349, 14.663, 13.067, 11.829, 10.753, 8.883,
     7.988, 7.253, 6.565, 5.919, 4.719, 4.132, 3.544, 2.943, 2.219, 0.110, -0.717, -0.717])
delta_r = numpy.array(
    [0.013, 0.013, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.030, 0.015,
     0.015, 0.015, 0.015, 0.030, 0.015, 0.015, 0.015, 0.015, 0.015, 0.008])
c = numpy.array(
    [0.050, 0.081, 0.080, 0.077, 0.073, 0.069, 0.065, 0.061, 0.058, 0.054, 0.051, 0.048, 0.046, 0.044, 0.040, 0.038,
     0.036, 0.035, 0.032, 0.031, 0.030, 0.029, 0.028, 0.027, 0.026, 0.026])

alfa = [-20, -19.5, -19, -18.5, -18, -17.5, -17, -16.5, -16, -15.5, -15, -14.5, -14, -13.5, -13, -12.5, -12, -11.5, -11,
        -10.5, -10, -9.5, -9, -8.5, -8, -7.5, -7, -6.5, -6, -5.5, -5, -4.5, -4, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1,
        1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14,
        14.5, 15, 15.5, 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20]
cl = [-0.4959, -0.486, -0.4759, -0.4655, -0.4551, -0.4447, -0.4343, -0.424, -0.4142, -0.4049, -0.3961, -0.3878, -0.38,
      -0.3727, -0.366, -0.3599, -0.355, -0.3515, -0.3499, -0.3387, -0.3177, -0.2991, -0.2796, -0.2703, -0.2658, -0.2754,
      -0.2462, -0.1767, -0.1239, -0.0597, 0.0149, 0.0897, 0.2015, 0.339, 0.3936, 0.449, 0.5041, 0.5583, 0.6144, 0.6706,
      0.7259, 0.7835, 0.8396, 0.895, 0.9505, 1.0061, 1.0595, 1.1126, 1.1645, 1.217, 1.2675, 1.3139, 1.3569, 1.3818,
      1.4185, 1.4427, 1.4633, 1.4841, 1.4969, 1.5173, 1.5236, 1.5117, 1.5225, 1.5278, 1.529, 1.5212, 1.5075, 1.4911,
      1.47, 1.4423, 1.4096, 1.3761, 1.3448, 1.3169, 1.2924, 1.2704, 1.2491, 1.2251, 1.1614]
cd = [0.23497, 0.229, 0.22313, 0.21734, 0.21158, 0.20593, 0.20044, 0.19509, 0.18932, 0.18328, 0.17716, 0.17097, 0.1647,
      0.15833, 0.15186, 0.14532, 0.13874, 0.13222, 0.12575, 0.11772, 0.10808, 0.09812, 0.0869, 0.07361, 0.05991,
      0.04021, 0.03044, 0.02194, 0.01743, 0.0145, 0.01331, 0.01065, 0.00767, 0.00828, 0.0088, 0.00864, 0.00879, 0.00902,
      0.00917, 0.00924, 0.0094, 0.00968, 0.00977, 0.00991, 0.00995, 0.01026, 0.01029, 0.0104, 0.01065, 0.01101, 0.01111,
      0.01148, 0.01189, 0.01344, 0.01744, 0.01905, 0.0209, 0.0228, 0.02531, 0.02743, 0.02985, 0.03532, 0.03921, 0.04397,
      0.04948, 0.05624, 0.06421, 0.07323, 0.08338, 0.09493, 0.10775, 0.1209, 0.1356, 0.15183, 0.1694, 0.18803, 0.20837,
      0.2317, 0.27961]

# f_c_L = interpolate.interp1d(alfa, cl, fill_value=(cl[0], cl[-1]), bounds_error=False)
# f_c_D = interpolate.interp1d(alfa, cd, fill_value=(cd[0], cd[-1]), bounds_error=False)
# res = run_main(r, c, theta, delta_r, B=3, R=0.45, Rhub=0.05, f_c_L=f_c_L, f_c_D=f_c_D)

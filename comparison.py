from turbine_data import SET_INIT
from induction_factors import Calculator
from cp_curve import calculate_power_3d
from utils import sort_xy
import numpy
from main import METHODS_STRINGS
import matplotlib.pyplot as plt

c = Calculator(SET_INIT["curves"])
SET_INIT["return_print"] = []
SET_INIT["return_results"] = []
#SET_INIT["tip_loss"] = True
#SET_INIT["hub_loss"] = True
#SET_INIT["new_tip_loss"] = True
#SET_INIT["new_hub_loss"] = True
SET_INIT["cascade_correction"] = True
SET_INIT["rotational_augmentation_correction"] = True

TSR_exp = [0.50732460, 0.91750193, 1.30609098, 1.45720894, 1.62991519, 2.01850424, 2.32074017, 2.70932922, 2.77409406, 3.01156515, 3.14109483, 3.48650732, 3.65921357, 3.76715497, 3.89668466, 4.26368543, 4.56592136, 4.73862760, 4.93292213, 5.17039322, 5.32151118, 5.75327679, 6.18504241, 6.63839630, 7.24286816, 7.60986893, 8.02004626, 8.36545875, 8.75404780, 8.99151889, 9.22898998, 9.67154973, 10.12490362, 10.50269854, 10.93446415, 11.47417116]
CP_exp = [0.00754310, 0.01508621, 0.02262931, 0.02715517, 0.03469828, 0.05129310, 0.06864224, 0.09579741, 0.10107759, 0.15614224, 0.19008621, 0.28211207, 0.33491379, 0.36961207, 0.38017241, 0.40581897, 0.43221983, 0.44730603, 0.45786638, 0.47144397, 0.47672414, 0.48502155, 0.48577586, 0.47974138, 0.46088362, 0.44278017, 0.41637931, 0.38922414, 0.35377155, 0.32963362, 0.30398707, 0.24967672, 0.19461207, 0.14633621, 0.08599138, 0.00150862]
#Experimental data from Karlsen, NTNU, p.31

plt.plot(TSR_exp,CP_exp,label="Experimental data",color="black")

for i in [0,1,2,4,6,7,8,9,10]: #+1 for counting max number in
    SET_INIT["method"] = i
    res = calculate_power_3d(SET_INIT)
    TSR, CP = sort_xy(res["TSR"], res["cp"])
    plt.plot(TSR,CP,label=METHODS_STRINGS[str(i)])
    plt.xlabel("TSR")
    plt.ylabel('Cp')

plt.legend()
plt.show()
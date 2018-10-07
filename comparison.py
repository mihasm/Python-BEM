from turbine_data import SET_INIT
from induction_factors import Calculator
from cp_curve import calculate_power_3d
from utils import sort_xy
import numpy
from main import METHODS_STRINGS
import matplotlib.pyplot as plt
from cycler import cycler
from utils import Printer


#cycler_map = (cycler('color', ['k']) * cycler('linestyle', ['-.','--', ':']) * cycler('marker', ['^',',','.','x'])) #Monochrome
cycler_map = (cycler('linestyle', ['-','-.','--', ':']) * cycler('color', ['r','g','b','c']) * cycler('marker', [',']))
c = Calculator(SET_INIT["curves"])
SET_INIT["return_print"] = []
SET_INIT["return_results"] = []

#SET_INIT["tip_loss"] = True
#SET_INIT["hub_loss"] = True
#SET_INIT["new_tip_loss"] = True
#SET_INIT["new_hub_loss"] = True
#SET_INIT["cascade_correction"] = True
#SET_INIT["rotational_augmentation_correction"] = True

def title_gen(SET_INIT):
    title = "Comparison of BEM methods"
    add=""
    i = 0
    if SET_INIT["tip_loss"]:
        if i>0:
            add+=", "
        add+="Prandtl's tip loss"
        i+=1
    if SET_INIT["hub_loss"]:
        if i>0:
            add+=", "
        add+="Prandtl's hub loss"
        i+=1
    if SET_INIT["new_tip_loss"]:
        if i>0:
            add+=", "
        add+="New tip loss"
        i+=1
    if SET_INIT["new_hub_loss"]:
        if i>0:
            add+=", "
        add+="New hub loss"
        i+=1
    if SET_INIT["cascade_correction"]:
        if i>0:
            add+=", "
        add+="Cascade effect"
        i+=1
    if SET_INIT["rotational_augmentation_correction"]:
        if i>0:
            add+=", "
        add+="Rotational augmentation"
        i+=1
    if i > 0:
        add=" (with "+add
    if i > 1:
        add+=" corrections)"
    elif i==1:
        add+=" correction)"
    return title+add

#Experimental data from Karlsen, NTNU, p.31
TSR_exp = [0.50732460, 0.91750193, 1.30609098, 1.45720894, 1.62991519, 2.01850424, 2.32074017, 2.70932922, 2.77409406, 3.01156515, 3.14109483, 3.48650732, 3.65921357, 3.76715497, 3.89668466, 4.26368543, 4.56592136, 4.73862760, 4.93292213, 5.17039322, 5.32151118, 5.75327679, 6.18504241, 6.63839630, 7.24286816, 7.60986893, 8.02004626, 8.36545875, 8.75404780, 8.99151889, 9.22898998, 9.67154973, 10.12490362, 10.50269854, 10.93446415, 11.47417116]
CP_exp = [0.00754310, 0.01508621, 0.02262931, 0.02715517, 0.03469828, 0.05129310, 0.06864224, 0.09579741, 0.10107759, 0.15614224, 0.19008621, 0.28211207, 0.33491379, 0.36961207, 0.38017241, 0.40581897, 0.43221983, 0.44730603, 0.45786638, 0.47144397, 0.47672414, 0.48502155, 0.48577586, 0.47974138, 0.46088362, 0.44278017, 0.41637931, 0.38922414, 0.35377155, 0.32963362, 0.30398707, 0.24967672, 0.19461207, 0.14633621, 0.08599138, 0.00150862]
error_x=[0.463678516,0.927357032,1.391035549,1.808346213,2.225656878,2.819165379,3.301391036,3.774343122,4.256568779,4.738794436,5.230293663,5.693972179,6.194744977,6.658423493,7.122102009,7.622874807,8.068006182,8.550231839,9.041731066,9.505409583,9.98763524,10.49768161,10.93353941,11.41576507]
error_y=[0.007567568,0.015135135,0.024216216,0.040864865,0.062054054,0.104432432,0.225513514,0.373837838,0.405621622,0.447243243,0.47372973,0.483567568,0.485081081,0.479027027,0.466162162,0.444972973,0.413189189,0.371567568,0.324648649,0.271675676,0.211891892,0.146810811,0.084756757,0.009837838]
error_size=numpy.array([0.004540541,0.007567568,0.015135135,0.024216216,0.015135135,0.027243243,0.028756757,0.056,0.485837838,0.354162162,0.060540541,0.059027027,0.045405405,0.042378378,0.042378378,0.049945946,0.051459459,0.059027027,0.056,0.051459459,0.048432432,0.049945946,0.046918919,0.052972973])/2.

            #[tiploss,hubloss,newtiploss,newhubloss,cascadeffects]
variants = [[False,False,False,False,False],
            [True,False,False,False,False],
            [True,True,False,False,False],
            [False,False,True,False,False],
            [False,False,True,True,False],
            [False,False,False,False,True],
            [True,False,False,False,True],
            [True,True,False,False,True],
            [False,False,True,False,True],
            [False,False,True,True,True]]

def comparison_runner(inp_settings,TSR_exp,CP_exp,error_x,error_y,error_size):
    p = Printer(inp_settings["return_print"])

    for i in range(len(variants)):

        plt.figure(i)

        p.print("Running variant",i,"/",len(variants))
        inp_settings["tip_loss"],inp_settings["hub_loss"],inp_settings["new_tip_loss"],inp_settings["new_hub_loss"],inp_settings["cascade_correction"] = variants[i]
        
        title =title_gen(inp_settings)
        fig, ax = plt.subplots(1,1,figsize=(10, 5))
        plt.plot(TSR_exp,CP_exp,label="Experimental data",color="black")
        plt.errorbar(error_x,error_y,error_size,capsize=2,color="black")
        ax.set_prop_cycle(cycler_map)

        for j in [0,1,2,4,6,7,8,9,10,13]: #+1 for counting max number in
            inp_settings["method"] = j
            p.print("    Calculating using method",j,"(%s)"%METHODS_STRINGS[str(j)])
            res = calculate_power_3d(inp_settings,False,"    ",False)
            TSR, CP = sort_xy(res["TSR"], res["cp"])
            plt.plot(TSR,CP,label=METHODS_STRINGS[str(j)])

        plt.xlabel(r"$\lambda_r$")
        plt.ylabel(r"$C_p$")
        plt.title(title,fontsize="x-small")
        plt.grid(which="both")
        plt.legend(loc=2,fontsize="x-small") # 1="upper right" 2="upper left"
        out_title = './out/Fig%s.png' % str(i)
        print(out_title)
        fig.savefig(out_title, bbox_inches='tight')
    plt.show()

comparison_runner(SET_INIT,TSR_exp,CP_exp,error_x,error_y,error_size)
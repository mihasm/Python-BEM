from turbine_data import SET_INIT
from calculation import Calculator
from utils import Printer
from numpy import radians, degrees
import numpy as np,numpy
from math import pi
import traceback
from matplotlib import pyplot as plt


# noinspection PyBroadException
def optimize_angles(inp_args):
    p = Printer(inp_args["return_print"])
    try:
        p.print("Optimizing angles for target variable:", inp_args["optimization_variable"])
        return_results = inp_args["return_results"]

        v = inp_args["target_speed"]
        rpm = inp_args["target_rpm"]
        # print(v,pi,rpm)
        omega = 2 * pi * rpm / 60
        # optimization_variable = "dT"
        optimization_variable = inp_args["optimization_variable"]

        output_angles = []

        inp_args["theta_in"] = np.array([120]*len(inp_args["theta_in"]))

        C = Calculator(inp_args["airfoils"])

        for section_number in range(len(inp_args["r_in"])):
            p.print("section_number is", section_number)

            _r = inp_args["r_in"][section_number]
            _c = inp_args["c_in"][section_number]
            _theta = inp_args["theta_in"][section_number]
            _dr = inp_args["dr"][section_number]
            _airfoil = inp_args["foils_in"][section_number]
            max_thickness = inp_args["airfoils"][_airfoil]["max_thickness"] * _c
            _airfoil_dat = _airfoil + ".dat"

            done_angles = {}  # key is theta, value is out

            p.print("initial theta is", _theta)
            got_through = False
            for dtheta in [90,60,45,30,10, 5, 1, 0.5, 0.1]:
                p.print("dtheta", dtheta)
                while True:
                    initial = False
                    if not initial:
                        while _theta > -90.0:
                            # initial angle
                            try:
                                out = C.calculate_section(v=v, omega=omega, _airfoil=_airfoil, _airfoil_dat=_airfoil_dat,
                                                          max_thickness=max_thickness, _r=_r, _c=_c, _dr=_dr, _theta=radians(_theta),
                                                          printer=p, **inp_args)
                            except Exception as e:
                                p.print(e)
                                p.print(traceback.format_exc())
                                out = None
                            if out == False or out == None:
                                _theta = _theta - 10
                            else:
                                done_angles[_theta] = out
                                initial = True
                                break

                    # middle angle
                    if not _theta in done_angles:
                        p.print("   calculating out")
                        try:
                            out = C.calculate_section(v=v, omega=omega, _airfoil=_airfoil, _airfoil_dat=_airfoil_dat,
                                                      max_thickness=max_thickness, _r=_r, _c=_c, _dr=_dr, _theta=radians(_theta),
                                                      printer=p, **inp_args)
                        except Exception as e:
                            p.print(e)
                            p.print(traceback.format_exc())
                            out = None
                        if out == False or out == None:
                            break
                        done_angles[_theta] = out
                    else:
                        p.print(_theta, "already calculated, reusing...")
                        out = done_angles[_theta]

                    # upper angle
                    if not _theta + dtheta in done_angles:
                        p.print("   calculating out_up")
                        try:
                            out_up = C.calculate_section(v=v, omega=omega, _airfoil=_airfoil, _airfoil_dat=_airfoil_dat,
                                                         max_thickness=max_thickness, _r=_r, _c=_c, _dr=_dr,
                                                         _theta=radians(_theta + dtheta), printer=p, **inp_args)
                        except Exception as e:
                            p.print(e)
                            p.print(traceback.format_exc())
                            out_up = None
                        if out_up == False or out_up == None:
                            break
                        done_angles[_theta + dtheta] = out_up
                    else:
                        p.print(_theta + dtheta, "already calculated, reusing...")
                        out_up = done_angles[_theta + dtheta]

                    # lower angle
                    if not _theta - dtheta in done_angles:
                        p.print("   calculating out_down")
                        try:
                            out_down = C.calculate_section(v=v, omega=omega, _airfoil=_airfoil,
                                                           _airfoil_dat=_airfoil_dat, max_thickness=max_thickness,
                                                           _r=_r, _c=_c, _dr=_dr, _theta=radians(_theta - dtheta),
                                                           printer=p, **inp_args)
                        except Exception as e:
                            p.print(e)
                            p.print(traceback.format_exc())
                            out_down = None
                        if out_down == False or out_down == None:
                            break
                        done_angles[_theta - dtheta] = out_down
                    else:
                        p.print(_theta - dtheta, "already calculated, reusing...")
                        out_down = done_angles[_theta - dtheta]

                    if out_up == False or out_down == False or out == False:
                        p.print("   one is False, breaking...")
                        break
                    if out_up == None or out_down == None or out == None:
                        p.print("   one is None, breaking...")
                        break
                    got_through = True

                    var = out[optimization_variable]
                    var_up = out_up[optimization_variable]
                    var_down = out_down[optimization_variable]
                    p.print("   %s" % optimization_variable, var, "%s_up" % optimization_variable, var_up,
                            "%s_down" % optimization_variable, var_down)

                    if var_up <= var and var_down <= var:
                        p.print("   none is bigger, breaking...")
                        break

                    if var_up > var > var_down:
                        p.print("   going up")
                        _theta = _theta + dtheta
                        var = var_up
                        out = out_up

                    if var_down > var > var_up:
                        p.print("   going down")
                        _theta = _theta - dtheta
                        var = var_down
                        out = out_down

                    if var_down > var and var_up > var:
                        if var_up > var_down:
                            p.print("   both up and down are bigger, going up")
                            _theta = _theta + dtheta
                            var = var_up
                            out = out_up
                        else:
                            p.print("   both up and down are bigger, going down")
                            _theta = _theta - dtheta
                            var = var_down
                            out = out_down
                    if var_up == var and var_down == var:
                        p.print("   both are equal, breaking")
                        break
                    p.print("   ***")
                p.print("***")
            if not got_through:
                p.print('optimization failed for section',section_number)
                p.print("!!!!EOF!!!!")
                return

            output_angles.append(_theta)
            p.print("final theta is", _theta)
            p.print("*******************************")
        p.print("angles:")
        p.print(output_angles)
        for a in output_angles:
            p.print(a)
        p.print("!!!!EOF!!!!")
    except Exception as e:
        p.print("Error in running optimizer: %s" % str(e))
        p.print("!!!!EOF!!!!")
        #raise

def maximize_for_both(inp_args):
    p = Printer(inp_args["return_print"])
    try:
        p.print("Optimizing angles for both propeller and generator operation...")
        return_results = inp_args["return_results"]

        v = inp_args["target_speed"]
        rpm = inp_args["target_rpm"]
        rpm_prop = inp_args["target_rpm_propeller"]
        omega = 2 * pi * rpm / 60
        omega_prop = 2*pi*rpm / 60
        optimization_variable = inp_args["optimization_variable"]

        output_angles = []

        C = Calculator(inp_args["airfoils"])

        for section_number in range(len(inp_args["r_in"])):
            p.print("section_number is", section_number)

            _r = inp_args["r_in"][section_number]
            _c = inp_args["c_in"][section_number]
            _theta = inp_args["theta_in"][section_number]
            _dr = inp_args["dr"][section_number]
            _airfoil = inp_args["foils_in"][section_number]
            max_thickness = inp_args["airfoils"][_airfoil]["max_thickness"] * _c
            _airfoil_dat = _airfoil + ".dat"

            #done_angles = {}  # key is theta, value is out

            #p.print("initial theta is", degrees(_theta))
            theta_array,dT_array,dQ_array = [],[],[]
            for _theta in np.linspace(0,90,100):
                dT,dQ = None, None
                
                #Test for wind turbine mode
                inp_args["propeller_mode"] = False
                out = C.calculate_section(v=v, omega=omega, _airfoil=_airfoil,_airfoil_dat=_airfoil_dat, max_thickness=max_thickness,_r=_r, _c=_c, _dr=_dr, _theta=radians(_theta),printer=p, **inp_args)
                if out != None and out != False:
                    dQ = out["dQ"]
                
                #Test for propeller
                inp_args["propeller_mode"] = True
                out_prop = C.calculate_section(v=0.01, omega=omega_prop, _airfoil=_airfoil,_airfoil_dat=_airfoil_dat, max_thickness=max_thickness,_r=_r, _c=_c, _dr=_dr, _theta=radians(_theta),printer=p, **inp_args)
                if out_prop != None and out_prop != False:
                    dT = out_prop["dT"]

                if dT != None and dQ != None:
                    theta_array.append(_theta)
                    dT_array.append(dT)
                    dQ_array.append(dQ)

            theta_array = np.array(theta_array)
            dT_array = np.array(dT_array)
            dQ_array = np.array(dQ_array)

            #plt.plot(theta_array,dT_array,"g-",label="dT"+str(section_number))
            #plt.plot(theta_array,dQ_array,"r-",label="dQ"+str(section_number))

            #indexes_positive_values_dT = numpy.where(dT_array > 0)
            #indexes_positive_values_dQ = numpy.where(dQ_array > 0)
            
            positive_indexes = numpy.logical_and(dT_array > 0, dQ_array > 0)
            positive_dT_array = dT_array[positive_indexes]
            positive_dQ_array = dQ_array[positive_indexes]
            positive_theta_array = theta_array[positive_indexes]
            #plt.plot(positive_theta_array,positive_dT_array,"g",label="dT")
            #plt.plot(positive_theta_array,positive_dQ_array,"r",label="dQ")
            normalized_dT_array = (positive_dT_array - positive_dT_array.min())/positive_dT_array.max()
            normalized_dQ_array = (positive_dQ_array - positive_dQ_array.min())/positive_dQ_array.max()
            crossings = get_crossings(positive_theta_array,normalized_dT_array,normalized_dQ_array)
            _max_i = np.where(crossings[:,1] == crossings[:,1].max())[0][0]
            #print("crossings",crossings)
            #print('_max_i',_max_i)

            dT_only_rising = np.all(np.diff(normalized_dT_array) > 0)
            dQ_only_rising = np.all(np.diff(normalized_dQ_array) > 0)
            dT_only_falling = np.all(np.diff(normalized_dT_array) < 0)
            dQ_only_falling = np.all(np.diff(normalized_dQ_array) < 0)

            if dT_only_rising and dQ_only_rising:
                _max_theta = positive_theta_array[-1]
            elif dT_only_falling and dQ_only_falling:
                _max_theta = positive_theta_array[0]
            else:
                _max_theta = crossings[:,0][_max_i]

            p.print("max_theta",_max_theta)
            output_angles.append(_max_theta)
            #plt.plot(positive_theta_array,normalized_dT_array,"g",label="dT"+str(section_number))
            #plt.plot(positive_theta_array,normalized_dQ_array,"b",label="dQ"+str(section_number))
            #plt.plot(crossings[:,0],crossings[:,1],'r*')
            #plt.axvline(_max_theta)

            #plt.legend()
            #plt.show()
        p.print("Angles:")
        for a in output_angles:
            p.print(a)

        p.print("!!!!EOF!!!!")
        #plt.legend()
        #plt.show()
    except Exception as e:
        p.print(str(e))
        p.print(traceback.format_exc())
        p.print("!!!!EOF!!!!")

"""
x = np.linspace(0,2*np.pi,40)
y = np.sin(x)
y_2 = np.cos(x)
plt.plot(x,y,'b.')
plt.plot(x,y_2,'g-')
"""


def get_crossings(x,y_1,y_2):
    """
    This function fetches x-es and y-s where y_1 and y_2 intersect.
    """
    #initial status
    top_indexes = []
    for i in range(len(x)):
        if i == 0:
            status = y_1[i] > y_2[i]
        else:
            new_status = y_1[i] > y_2[i]
            if new_status != status:
                top_indexes.append(i)
                status = new_status
    top_indexes = np.array(top_indexes)
    bottom_indexes = top_indexes-1
    out = []
    for i in range(len(top_indexes)):
        i1 = bottom_indexes[i]
        i2 = top_indexes[i]
        x1 = x[i1]
        y1 = y_1[i1]
        x2 = x[i2]
        y2 = y_1[i2]
        x3 = x1
        y3 = y_2[i1]
        x4 = x2
        y4 = y_2[i2]
        x_out,y_out = findIntersection(x1,y1,x2,y2,x3,y3,x4,y4)
        out.append([x_out,y_out])
    out = np.array(out)
    return out

def findIntersection(x1,y1,x2,y2,x3,y3,x4,y4):
    #https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
    px= ( (x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4) ) / ( (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4) ) 
    py= ( (x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4) ) / ( (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4) )
    return [px, py]

# noinspection PyBroadException
def optimal_pitch(inp_args):
    p = Printer(inp_args["return_print"])
    try:
        p.print("Optimizing angles for target variable:", inp_args["optimization_variable"])
        return_results = inp_args["return_results"]

        v = inp_args["target_speed"]
        rpm = inp_args["target_rpm"]
        omega = 2 * pi * rpm / 60
        optimization_variable = inp_args["optimization_variable"]

        output_angles = []

        C = Calculator(inp_args["airfoils"])

        done_pitches = {}
        _pitch = 0 # is in radians

        args = {**inp_args}

        for dpitch in [45,30,20,10,5,2,1,0.1]:
            while True:
                # middle angle
                if not _pitch in done_pitches:
                    p.print("   calculating out (pitch:%s)" % _pitch)
                    try:
                        args["pitch"] = _pitch
                        out = C.run_array(**args,rpm=rpm,v=v)
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
                    p.print("   calculating out_up (pitch:%s)" %_pitch_up)
                    try:
                        args["pitch"] = _pitch_up
                        out_up = C.run_array(**args,rpm=rpm,v=v)
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
                        out_down = C.run_array(**args,rpm=rpm,v=v)
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
        p.print("Final pitch:",_pitch)
        p.print("Angles:")
        for t in args["theta"]:
            p.print(t+_pitch)
        p.print("!!!!EOF!!!!")
    except:
        p.print("Error in running optimizer")
        p.print("!!!!EOF!!!!")
        raise


"""
from polars import scrape_data,get_extrapolated_data
from utils import sort_data
data = scrape_data("http://airfoiltools.com/airfoil/details?airfoil=s826-nr")
data = get_extrapolated_data(data)
data = sort_data(data)
settings = SET_INIT
#SET_INIT["fix_reynolds"] = True
#SET_INIT["reynolds"] = 200000
SET_INIT["method"] = SET_INIT["method"]+1
SET_INIT["return_print"] = []
SET_INIT["return_results"] = []
SET_INIT["airfoils"]["s826"]["gathered_curves"] = data
#SET_INIT["airfoils"]["s826"]["interp_function_cd"] = f_cd
maximize_for_both(SET_INIT)
"""


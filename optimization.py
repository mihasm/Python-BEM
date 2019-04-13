from turbine_data import SET_INIT
from induction_factors import Calculator
from utils import Printer
from numpy import radians, degrees
from math import pi
import traceback



# noinspection PyBroadException
def optimize_angles(inp_args):
    p = Printer(inp_args["return_print"])
    try:
        p.print("Optimizing angles for target variable:",inp_args["optimization_variable"])
        return_results = inp_args["return_results"]

        v = inp_args["target_speed"]
        rpm = inp_args["target_rpm"]
        #print(v,pi,rpm)
        omega = 2 * pi * rpm / 60
        # optimization_variable = "dT"
        optimization_variable = inp_args["optimization_variable"]

        output_angles = []

        C = Calculator(inp_args["airfoils"])

        for section_number in range(len(inp_args["r_in"])):
            p.print("section_number is", section_number)

            _r = inp_args["r_in"][section_number]
            _c = inp_args["c_in"][section_number]
            _theta = radians(inp_args["theta_in"][section_number])
            _dr = inp_args["dr"][section_number]
            _airfoil = inp_args["foils_in"][section_number]
            max_thickness = inp_args["airfoils"][_airfoil]["max_thickness"] * _c
            _airfoil_dat = _airfoil + ".dat"

            done_angles = {}  # key is theta, value is out

            p.print("initial theta is", degrees(_theta))
            for dtheta in [10, 5, 1, 0.5, 0.1]:
                p.print("dtheta", dtheta)

                while True:

                    # middle angle
                    if not _theta in done_angles:
                        p.print("   calculating out")
                        try:
                            out = C.calculate_section(v=v,omega=omega,_airfoil=_airfoil,_airfoil_dat=_airfoil_dat, max_thickness=max_thickness,
                                                    _r=_r, _c=_c, _dr=_dr, _theta=_theta, printer=p, **inp_args)
                        except Exception as e:
                            print(e)
                            print(traceback.format_exc())
                            out = None
                        if out == False or out == None:
                            break
                        done_angles[_theta] = out
                    else:
                        p.print(_theta, "already calculated, reusing...")
                        out = done_angles[_theta]

                    # upper angle
                    if not _theta + radians(dtheta) in done_angles:
                        p.print("   calculating out_up")
                        try:
                            out_up = C.calculate_section(v=v,omega=omega,_airfoil=_airfoil,_airfoil_dat=_airfoil_dat, max_thickness=max_thickness,
                                                    _r=_r, _c=_c, _dr=_dr, _theta=_theta+radians(dtheta), printer=p, **inp_args)
                        except Exception as e:
                            print()
                            print(traceback.format_exc())
                            out_up = None
                        if out_up == False or out_up == None:
                            break
                        done_angles[_theta + radians(dtheta)] = out_up
                    else:
                        p.print(_theta + radians(dtheta), "already calculated, reusing...")
                        out_up = done_angles[_theta + radians(dtheta)]

                    # lower angle
                    if not _theta - radians(dtheta) in done_angles:
                        p.print("   calculating out_down")
                        try:
                            out_down = C.calculate_section(v=v,omega=omega,_airfoil=_airfoil,_airfoil_dat=_airfoil_dat, max_thickness=max_thickness,
                                                    _r=_r, _c=_c, _dr=_dr, _theta=_theta-radians(dtheta), printer=p, **inp_args)
                        except Exception as e:
                            print(e)
                            print(traceback.format_exc())
                            out_down = None
                        if out_down == False or out_down == None:
                            break
                        done_angles[_theta - radians(dtheta)] = out_down
                    else:
                        p.print(_theta - radians(dtheta), "already calculated, reusing...")
                        out_down = done_angles[_theta - radians(dtheta)]

                    if out_up == False or out_down == False or out == False:
                        p.print("   one is False, breaking...")
                        break
                    if out_up == None or out_down == None or out == None:
                        p.print("   one is None, breaking...")
                        break

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
                        _theta = _theta + radians(dtheta)
                        var = var_up
                        out = out_up

                    if var_down > var > var_up:
                        p.print("   going down")
                        _theta = _theta - radians(dtheta)
                        var = var_down
                        out = out_down

                    if var_down > var and var_up > var:
                        if var_up > var_down:
                            p.print("   both up and down are bigger, going up")
                            _theta = _theta + radians(dtheta)
                            var = var_up
                            out = out_up
                        else:
                            p.print("   both up and down are bigger, going down")
                            _theta = _theta - radians(dtheta)
                            var = var_down
                            out = out_down
                    if var_up == var and var_down == var:
                        p.print("   both are equal, breaking")
                        break
                    p.print("   ***")
                p.print("***")
            output_angles.append(_theta)
            p.print("final theta is", degrees(_theta))
            p.print("*******************************")
        p.print("angles:")
        p.print(degrees(output_angles))
        for a in degrees(output_angles):
            p.print(a)
        p.print("!!!!EOF!!!!")
    except:
        p.print("Error in running optimizer")
        p.print("!!!!EOF!!!!")
        raise

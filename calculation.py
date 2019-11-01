import numbers
from math import sin, cos, atan, acos, pi, exp, sqrt, radians, atan2, degrees, tan
import os

import numpy
import numpy as np
from scipy import interpolate

from utils import Printer, generate_dat, sort_data, normalize_angle
from popravki import *
from xfoil import run_xfoil_analysis, xfoil_runner
from interpolator import interp

numpy.seterr(all="raise")
numpy.seterr(invalid="raise")

OUTPUT_VARIABLES_LIST = {
    "a":{"type":"array","name":"Axial induction factor","symbol":"a","unit":""},
    "a'":{"type":"array","name":"Tangential induction factor","symbol":"a'","unit":""},
    "cL":{"type":"array","name":"Lift coefficient","symbol":r"$C_L$","unit":""},
    "cD":{"type":"array","name":"Drag coefficient","symbol":r"$C_L$","unit":""},
    "alpha":{"type":"array","name":"Angle of attack","symbol":r"$\alpha$","unit":"°"},
    "phi":{"type":"array","name":"Relative wind angle","symbol":r"$\phi$","unit":"°"},
    "F":{"type":"array","name":"Tip loss correction factor","symbol":"F","unit":""},
    "dFt":{"type":"array","name":"Incremental tangential force","symbol":r"$dF_t$","unit":"N"},
    "M":{"type":"array","name":"Torque","symbol":"M","unit":"Nm"},
    "lambda_r":{"type":"array","name":"Local tip speed ratio","symbol":r"$\lambda_r$","unit":""},
    "Ct":{"type":"array","name":"Tangential coefficient","symbol":r"$C_t$","unit":""},
    "dFn":{"type":"array","name":"Incremental normal force","symbol":r"$dF_n$","unit":"N"},
    #"foils":{"type":"string_array","name":"Airfoil name","symbol":"airfoil_name","unit":""},
    "dT":{"type":"array","name":"Incremental thrust","symbol":"dT","unit":"N"},
    "dQ":{"type":"array","name":"Incremental torque","symbol":"dM","unit":"N"},
    "Re":{"type":"array","name":"Reynolds number","symbol":"Re","unit":""},
    "U1":{"type":"array","name":"Far-upwind speed","symbol":"U1","unit":"m/s"},
    "U2":{"type":"array","name":"Near-upwind speed","symbol":"U2","unit":"m/s"},
    "U3":{"type":"array","name":"Near-downwind speed","symbol":"U3","unit":"m/s"},
    "U4":{"type":"array","name":"Far-downwind speed","symbol":"U4","unit":"m/s"},


    "r":{"type":"array","name":"Section radius","symbol":"r","unit":"m"},
    "M":{"type":"array","name":"Section torque","symbol":"M","unit":"Nm"},
    "dr":{"type":"array","name":"Section height","symbol":"dr","unit":"m"},
    "c":{"type":"array","name":"Section chord length","symbol":"c","unit":"m"},
    "theta":{"type":"array","name":"Section twist angle","symbol":r"$\Theta$","unit":"°"},

    "R":{"type":"float","name":"Turbine radius","symbol":"R","unit":"m"},
    "rpm":{"type":"float","name":"Turbine rotational velocity","symbol":r"$\Omega$","unit":"RPM"},
    "v":{"type":"float","name":"Wind speed","symbol":"v","unit":"m/s"},
    "cp":{"type":"float","name":"Power coefficient","symbol":r"$C_P$","unit":""},
    "ct":{"type":"float","name":"Thrust coefficient","symbol":r"$C_T$","unit":""},
    "TSR":{"type":"float","name":"Tip speed ratio","symbol":r"$\lambda$","unit":""},
    "Ft":{"type":"float","name":"Tangential force","symbol":r"$F_t$","unit":"N"},
    "omega":{"type":"float","name":"Rotational velocity","symbol":r"$\Omega$","unit":r"$rad^{-1}$"},
    "Msum":{"type":"float","name":"Torque sum","symbol":r"$M_Sum$","unit":"Nm"},
    "power":{"type":"float","name":"Power","symbol":"P","unit":"W"},
    "thrust":{"type":"float","name":"Thrust","symbol":"T","unit":"N"},
    "Rhub":{"type":"float","name":"Hub radius","symbol":r"$R_hub$","unit":"m"},
    "B":{"type":"float","name":"Number of blades","symbol":"B","unit":""},
    "J":{"type":"float","name":"Advance ratio","symbol":"J","unit":""},
    "eff":{"type":"float","name":"Efficiency","symbol":r"$\eta_p$","unit":""},
}


class Calculator:
    """
    Class for calculation of induction factors using BEM theory.
    """

    def __init__(self, airfoils):
        self.airfoils = airfoils
        for blade_name in self.airfoils:
            self.airfoils[blade_name]["alpha_zero"] = 0.0  # TODO FIX
            generate_dat(
                blade_name, self.airfoils[blade_name]["x"], self.airfoils[blade_name]["y"])

            ncrit_selected = self.airfoils[blade_name]["ncrit_selected"]

            data = self.airfoils[blade_name]["gathered_curves"]
            data = data[np.in1d(data[:,1],ncrit_selected)]
            data = sort_data(data)
            print(data)

            re = data[:, 0].flatten()
            alpha = data[:, 2].flatten()
            cl = data[:, 3].flatten()
            cd = data[:, 4].flatten()

            def interpolation_function_cl(x, y, re=re, alpha=alpha, cl=cl):
                return interp(x, y, re, alpha, cl)

            def interpolation_function_cd(x, y, re=re, alpha=alpha, cd=cd):
                return interp(x, y, re, alpha, cd)

            self.airfoils[blade_name]["interp_function_cl"] = interpolation_function_cl
            self.airfoils[blade_name]["interp_function_cd"] = interpolation_function_cd

    def printer(self, _locals, p):
        p.print("----Running induction calculation for following parameters----")
        for k, v in _locals.items():
            if isinstance(v, dict):
                for k2, v2 in v.items():
                    _p2 = "    " + k2 + ":" + str(v2)
                    p.print(_p2)
            elif isinstance(v, list):
                for l in v:
                    _l = "    " + l
                    p.print(_l)
            else:
                _p = k + ":" + str(v)
                p.print(_p)
        p.print("--------------------------------------------------------------")
        return

    def convert_to_array(self, theta, c, r):
        """
        Converts integers or floats into numpy arrays.
        :param theta: int or float
        :param c: int or float
        :param r: int or float
        :return: np.array(theta),np.array(c),np.array(r)
        """
        if isinstance(theta, numpy.ndarray) and isinstance(c, numpy.ndarray) and isinstance(r, numpy.ndarray):
            return theta, c, r
        else:
            if isinstance(theta, numbers.Real) and isinstance(c, numbers.Real) and isinstance(r, numbers.Real):
                return numpy.array([theta]), numpy.array([c]), numpy.array([r])
            return None

    # noinspection PyUnusedLocal,PyUnusedLocal
    def run_array(self, theta, B, c, r, foils, dr, R, Rhub, rpm, v, pitch, method, propeller_mode, print_out, tip_loss,mach_number_correction, 
                  hub_loss, new_tip_loss, new_hub_loss, cascade_correction, max_iterations, convergence_limit, rho,
                  relaxation_factor, print_all, rotational_augmentation_correction,rotational_augmentation_correction_method,
                   fix_reynolds, reynolds, yaw_angle, skewed_wake_correction, print_progress=False, return_print=[], return_results=[], 
                   *args,**kwargs):
        """
        Calculates induction factors using standard iteration methods.

        Different methods are available as different fInductionCoefficients functions.

        ANGLES REPRESENTATION SHOWN IN
        https://cmm2017.sciencesconf.org/129068/document
        alpha - angle of attack
        phi - angle of relative wind
        beta - theta

        :param reynolds: Reynolds number (when forced) [float]
        :param fix_reynolds: Force Reynolds number [bool]
        :param mach_number_correction: use only for propeller [bool]
        :param foils: list of airfoils [str]
        :param propeller_mode: if calculating propeller thrust [bool]
        :param rotational_augmentation_correction_method:
        :param rotational_augmentation_correction:
        :param return_results: lst, used for returning results to main class
        :param return_print: lst, used for printing using main class
        :param print_all: prints every iteration
        :param relaxation_factor: relaxation factor
        :param method: method of calculating induction factors
        :param rho: air density [kg/m^3]
        :param convergence_limit: convergence criterion
        :param max_iterations: maximum number of iterations
        :param cascade_correction: uses cascade correction
        :param new_hub_loss: uses new tip loss correction
        :param new_tip_loss: uses new tip loss correction
        :param hub_loss: uses Prandtl tip loss correction
        :param tip_loss: uses Prandtl tip loss correction
        :param print_out: bool; if true, prints iteration data, default: False
        :param v: wind speed [m]
        :param r: sections radiuses [m]
        :param c: sections chord lengths [m]
        :param pitch: blade pitch (twist) [degrees]
        :param theta: twist - theta [deg]
        :param rpm: rotational velocity [rpm]
        :param dr: np array of section heights [m]
        :param R: outer (tip) radius [m]
        :param Rhub: hub radius [m]
        :param B: number of blades
        :return: dicitonary with results
        """

        p = Printer(return_print)

        if print_all:
            self.printer(locals(), p=p)

        theta, c, r = self.convert_to_array(theta, c, r)

        # create results array placeholders
        results = {}
        arrays = ["a", "a'", "cL", "alpha", "phi", "F", "dFt", "M", "lambda_r",
                  "Ct", "dFn", "foils", "dT", "dQ", "Re", "U1", "U2", "U3", "U4","cD"]
        for array in arrays:
            results[array] = numpy.array([])

        # set constants that are section-independent
        omega = rpm * 2 * pi / 60
        TSR = omega * R / v  # tip speed ratio
        J = v / (rpm / 60 * R * 2)
        kin_viscosity = 1.4207E-5  # Kinematic viscosity

        section_number = 0

        for n in range(len(theta)):
            section_number += 1

            _r = r[n]
            _c = c[n]
            _theta = radians(theta[n])
            _airfoil = foils[n]
            _airfoil_dat = _airfoil + ".dat"
            _dr = dr[n]

            if print_out:
                p.print("    r", _r, "(" + str(section_number) + ")")

            # get max thickness
            max_thickness = self.airfoils[_airfoil]["max_thickness"] * _c

            # Coning angle (PROPX: Definitions,Derivations, Data Flow, p.22)
            psi = 0.0

            _locals = locals()
            del _locals["self"]

            out_results = self.calculate_section(**_locals, printer=p)

            if print_progress:
                p.print("*", add_newline=False)

            if out_results == None:
                return None

            results["a"] = numpy.append(results["a"], out_results["a"])
            results["a'"] = numpy.append(results["a'"], out_results["aprime"])
            results["cL"] = numpy.append(results["cL"], out_results["Cl"])
            results["cD"] = numpy.append(results["cD"], out_results["Cd"])
            results["alpha"] = numpy.append(
                results["alpha"], out_results["alpha"])
            results["phi"] = numpy.append(results["phi"], out_results["phi"])
            results["F"] = numpy.append(results["F"], out_results["F"])
            results["dFt"] = numpy.append(results["dFt"], out_results["dFt"])
            results["Ct"] = numpy.append(results["Ct"], out_results["Ct"])
            results["dFn"] = numpy.append(results["dFn"], out_results["dFn"])
            results["foils"] = numpy.append(
                results["foils"], out_results["_airfoil"])
            results["dT"] = numpy.append(results["dT"], out_results["dT"])
            results["dQ"] = numpy.append(results["dQ"], out_results["dQ"])
            results["Re"] = numpy.append(results["Re"], out_results["Re"])
            results["U1"] = numpy.append(results["U1"], out_results["U1"])
            results["U2"] = numpy.append(results["U2"], out_results["U2"])
            results["U3"] = numpy.append(results["U3"], out_results["U3"])
            results["U4"] = numpy.append(results["U4"], out_results["U4"])
            results["lambda_r"] = numpy.append(results["lambda_r"], out_results["lambda_r"])

        dFt = results["dFt"]
        Ft = numpy.sum(dFt)
        M = B * dFt * r  # momenti po prerezih
        dQ = results["dQ"]
        Q = numpy.sum(dQ)  # moment for propeller
        power_p = Q * omega
        Msum = numpy.sum(M)
        power = numpy.sum(M) * omega
        Pmax = 0.5 * rho * v ** 3 * pi * R ** 2
        cp_w = power / Pmax
        cp_p = power_p / (rho * (rpm / 60) ** 3 * (2 * R) ** 5)

        dFn = results["dFn"]
        Fn = numpy.sum(dFn)
        dT = results["dT"]
        T = numpy.sum(dT)
        ct_w = T / (0.5 * rho * v ** 2 * pi * R ** 2)
        ct_p = T / (rho * (2 * R) ** 4 * (rpm / 60) ** 2)

        cq_p = Q / (rho * (2*R)**5 * (rpm/60)**2)
        eff = J/2/pi*ct_p/cq_p

        #floats
        results["R"] = R
        results["rpm"] = rpm
        results["v"] = v
        if propeller_mode:
            results["cp"] = cp_p
            results["ct"] = ct_p
        else:
            results["cp"] = cp_w
            results["ct"] = ct_w
        
        results["TSR"] = TSR
        results["Ft"] = Ft
        results["omega"] = omega
        results["Msum"] = Msum
        results["power"] = power
        results["thrust"] = T
        results["Rhub"] = Rhub
        results["B"] = B        
        results["J"] = J
        results["eff"] = eff

        #arrays
        results["r"] = r
        results["M"] = M
        results["dr"] = dr
        results["c"] = c
        results["theta"] = theta
        return results

    def calculate_section(self, v, omega, _r, _c, _theta, _dr, B, R, _airfoil_dat, _airfoil, max_thickness, Rhub,
                          propeller_mode, pitch=0.0, psi=0.0, fix_reynolds=False, reynolds=1e6, tip_loss=False, new_tip_loss=False,
                          hub_loss=False, new_hub_loss=False, cascade_correction=False, rotational_augmentation_correction=False,
                          rotational_augmentation_correction_method=0, mach_number_correction=False, method=5,
                          kin_viscosity=1.4207E-5, rho=1.225, convergence_limit=0.001, max_iterations=100, relaxation_factor=0.3,
                          printer=None, print_all=False, print_out=False, yaw_angle=0.0, tilt_angle=0.0, skewed_wake_correction=False, *args, **kwargs):
        """
        Function that calculates each section of the blade.

        Inputs:
        v:float: Wind speed [m/s].
        omega:float: Rotational velocity [s^-1].
        _r:float: Section radius [m].
        _c:float: Section chord length [m].
        _theta:float: Section angle [rad].
        _dr:float: Section height [m].
        B:int: Number of blades.
        R:float: Wind turbine radius [m].
        _airfoil_dat:str: Airfoil file name.
        _airfoil:str: Airfoil name.
        max_thickness:float: Foil thickness in percentage of chord length.
        propeller_mode:bool: Boolean that turns on propeller calculation mode if True.
        pitch:float: Pitch in degrees that adds fixed angle to _theta.
        psi:float: Coning angle.
        fix_reynolds:bool: True if only data at one Reynolds number is to be used.
        reynolds:float: Reynolds number if fix_reynolds is True.
        tip_loss:bool: Prandtl tip loss.
        hub_loss:bool: Prandtl hub loss.
        new_tip_loss:bool: Shen tip loss.
        new_hub_loss:bool: Shen hub loss.
        cascade_correction:bool: Cascade correction.
        rotational_augmentation_correction:bool: Rotational augmentation correction.
        rotational_augmentation_correction_method:int: Rotational augmentation correction method.
        mach_number_correction:bool: Mach number correction (for propellers).
        method:int: Method for calculating axial and tangential induction.
        kin_viscosity:float: Kinematic viscosity.
        rho:float: Air density [kgm^-3].
        convergence_limit:float: Convergence limit.
        max_iterations:int: Maximum iterations.
        relaxation_factor: Relaxation factor.
        yaw_angle:float: Yaw angle in degrees.
        tilt_angle:float: Tilting angle in degrees.
        """

        p = printer

        # local speed ratio
        lambda_r = omega * _r / v

        # solidity
        sigma = _c * B / (2 * pi * _r)

        # initial guess
        a = 1/3
        aprime = 0.01

        # iterations counter
        i = 0

        # tip mach number
        M = omega * _r / 343

        # convert pitch to radians
        _pitch = radians(pitch)

        # convert yaw to radians
        yaw_angle = radians(yaw_angle) # [Radians]
        psi = radians(psi) #Coning angle [Radians]
        tilt_angle = radians(tilt_angle) # [Radians]

        ############ START ITERATION ############
        while True:
            # update counter
            i = i + 1

            # for pretty-printing only
            prepend = ""
            #Equations for Vx and Vy from https://pdfs.semanticscholar.org/5e7d/9c6408b7dd8841692d950d08bce90c676dc1.pdf
            Vx = v*((cos(yaw_angle)*sin(tilt_angle)+sin(yaw_angle))*sin(psi)+cos(yaw_angle)*cos(psi)*cos(tilt_angle))
            Vy = omega*_r*cos(psi)+v*(cos(yaw_angle)*sin(tilt_angle)-sin(yaw_angle))

            # wind components
            if propeller_mode:
                Un = Vx * (1 + a)
                Ut = Vy * (1 - aprime)
            else:
                Un = Vx * (1 - a)
                Ut = Vy * (1 + aprime)

            Vrel_norm = sqrt(Un ** 2 + Ut ** 2)

            if fix_reynolds:
                Re = reynolds
            else:
                Re_next = Vrel_norm * _c / kin_viscosity
                if Re_next > 1e7:
                    Re_next = 2e5
                Re = int(Re_next)

            # relative wind
            phi = atan2(Un, Ut)

            F = 1
            # Prandtl tip loss
            if tip_loss:
                F = F * fTipLoss(B, _r, R, phi)

            # New tip loss
            elif new_tip_loss:
                F = F * newTipLoss(B, _r, R, phi, lambda_r)

            # Prandtl hub loss
            if hub_loss:
                F = F * fHubLoss(B, _r, Rhub, phi)

            # New hub loss
            elif new_hub_loss:
                F = F * newHubLoss(B, _r, R, phi, lambda_r)

            # angle of attack
            if propeller_mode:
                alpha = (_theta + _pitch) - phi
            else:
                alpha = phi - (_theta + _pitch)

            alpha = radians(normalize_angle(degrees(alpha)))

            # cascade correction
            if cascade_correction:
                alpha = cascadeEffectsCorrection(alpha=alpha, v=v, omega=omega, r=_r, R=R, c=_c, B=B, a=a,
                                                 aprime=aprime, max_thickness=max_thickness)

            """
            # For xFoil cL,cD
            xfoil_return = xfoil_runner(airfoil=_airfoil_dat, reynolds=Re, alpha=alpha, printer=p, print_all=print_all)

            if xfoil_return == False:
                p.print("        Xfoil failed")
                return None

            #Cl, Cd = xfoil_return["CL"], xfoil_return["CD"] #direct xfoil calculation - no interpolation
            """

            Cl, Cd = self.airfoils[_airfoil]["interp_function_cl"](Re, degrees(alpha)), self.airfoils[_airfoil][
                "interp_function_cd"](Re, degrees(alpha))

            if print_all:
                p.print("        CL:", Cl, "Cd:", Cd)

            if rotational_augmentation_correction:
                if print_all:
                    p.print("--")
                    p.print("  Cl:", Cl, "Cd:", Cd)
                Cl, Cd = calc_rotational_augmentation_correction(alpha=alpha, Cl=Cl, Cd=Cd, omega=omega, r=_r, R=R,
                                                                 c=_c, theta=_theta, v=v, Vrel=Vrel_norm,
                                                                 method=rotational_augmentation_correction_method,
                                                                 alpha_zero=radians(-5))

                if print_all:
                    p.print("  Cl_cor:", Cl, "Cd_cor:", Cd)
                    p.print("--")

            if mach_number_correction:
                Cl = machNumberCorrection(Cl, M)

            # normal and tangential coefficients
            C_norm = Cl * cos(phi) + Cd * sin(phi)
            C_tang = Cl * sin(phi) - Cd * cos(phi)

            input_arguments = {"F": F, "lambda_r": lambda_r, "phi": phi, "sigma": sigma, "C_norm": C_norm,
                               "C_tang": C_tang, "Cl": Cl, "Cd": Cd, "B": B, "c": _c, "r": _r, "R": R, "psi": 0.0,
                               "aprime_last": aprime, "omega": omega, "v": v, "a_last": a,
                               # "alpha_zero": airfoils[_airfoil]["alpha_zero"],
                               "method": method, "alpha": alpha, "alpha_deg": degrees(alpha)}

            if print_all:
                args_to_print = sorted(
                    [key for key, value in input_arguments.items()])
                p.print("            i", i)
                for argument in args_to_print:
                    p.print("            ", argument, input_arguments[argument])
                p.print("             --------")

            # calculate new induction coefficients
            coeffs = calculate_coefficients(method, input_arguments)
            if coeffs == None:
                return None

            # save old values
            a_last = a
            aprime_last = aprime

            #set new values
            a, aprime, Ct = coeffs

            if skewed_wake_correction:
                a_skewed = skewed_wake_correction_calculate(yaw_angle,a,_r,R)
                a_no_skew = a
                a = a_skewed

            # force calculation
            dFL = Cl * 0.5 * rho * Vrel_norm ** 2 * _c * _dr  # lift force
            dFD = Cd * 0.5 * rho * Vrel_norm ** 2 * _c * _dr  # drag force
            dFt = dFL * sin(phi) - dFD * cos(phi)  # tangential force
            dFn = dFL * cos(phi) + dFD * sin(phi)  # normal force

            # thrust and torque - Wiley, WE 2nd, p.124
            dT_MT = F * 4 * pi * _r * rho * v ** 2 * a * (1 - a) * _dr
            dT_BET = 0.5 * rho * B * _c * Vrel_norm ** 2 * \
                (Cl * cos(phi) + Cd * sin(phi)) * _dr
            dQ_MT = F * 4 * aprime * (1 - a) * rho * \
                v * pi * _r ** 3 * omega * _dr
            dQ_BET = B * 0.5 * rho * Vrel_norm ** 2 * \
                (Cl * sin(phi) - Cd * cos(phi)) * _c * _dr * _r

            # thrust-propeller
            dT_MT_p = 4 * pi * _r * rho * v ** 2 * (1 + a) * a * _dr
            dQ_MT_p = 4 * pi * _r ** 3 * rho * \
                v * omega * (1 + a) * aprime * _dr
            dT_BET_p = 0.5 * rho * v ** 2 * _c * B * (1 + a) ** 2 / (sin(phi) ** 2) * (
                Cl * cos(phi) - Cd * sin(phi)) * _dr
            dQ_BET_p = 0.5 * rho * v * _c * B * omega * _r ** 2 * (1 + a) * (1 - aprime) / (sin(phi) * cos(phi)) * (
                Cl * sin(phi) + Cd * cos(phi)) * _dr

            # from http://www.aerodynamics4students.com/propel.m
            dT_BET_p_2 = 0.5*rho*Vrel_norm**2*B*_c * \
                (Cl * cos(phi) - Cd * sin(phi))*_dr
            dQ_BET_p_2 = 0.5*rho*Vrel_norm**2*B*_c * \
                _r*(Cl * sin(phi) + Cd * cos(phi)) * _dr
            dT_MT_p_2 = 4*pi*_r*rho*v**2*(1+a)
            dQ_MT_p_2 = 4*pi*_r**3*rho*v*(1+a)*omega

            if propeller_mode:
                dT = dT_BET_p
                dQ = dQ_BET_p
            else:
                dT = dT_BET
                dQ = dQ_BET

            # wind after
            if propeller_mode:
                U1 = v
                U2 = None
                U3 = U1*(1+a)
                U4 = U1*(1+2*a)
            else:
                U1 = v
                U2 = U1*(1-a)
                U3 = None
                U4 = U1 * (1 - 2 * a)

            if skewed_wake_correction:
                # replace a with no skew
                a = a_no_skew

            # check convergence
            if abs(a - a_last) < convergence_limit:
                break

            # check iterations limit
            if i >= max_iterations:
                if print_out:
                    p.print("-*-*-*-*-*-*-*-*-*-*-*-*-*-\n", "|max iterations exceeded\n", "|------>a:", a, " aprime",
                            aprime, )
                    prepend = "|"
                return None

            # relaxation
            # aprime=aprime_last+relaxation_factor*(aprime-aprime_last)
            a = a_last + relaxation_factor * (a - a_last)

        ############ END ITERATION ############

        if print_out:
            p.print(prepend, "        iters: ", i)
            p.print(prepend, "        alpha: ", degrees(alpha)), "Cl", str(Cl)
            p.print(prepend, "        a: ", a, "a'", str(aprime))
            p.print(prepend, "        LSR: ", lambda_r)
            p.print(prepend, "        Vrel: ", Vrel_norm)
            p.print(prepend, "        Re:", Re)
            p.print(prepend, "        foil:", _airfoil_dat)
            p.print(prepend, "        Cl:", Cl)
            p.print(prepend, "        Cd:", Cd)
            p.print(prepend, "        U1:", U1)
            p.print(prepend, "        U2:", U2)
            p.print(prepend, "        U3:", U3)
            p.print(prepend, "        U4:", U4)
            p.print(prepend, "        dT_MT %.2f dT_BET %.2f" %
                    (dT_MT, dT_BET))
            p.print(prepend, "        dQ_MT %.2f dQ_BET %.2f" %
                    (dQ_MT, dQ_BET))
            p.print(prepend, "        dT_MT_p %.5f dT_BET_p %.5f" %
                    (dT_MT_p, dT_BET_p))
            p.print(prepend, "        dQ_MT_p %.5f dQ_BET_p %.5f" %
                    (dQ_MT_p, dQ_BET_p))
            p.print(prepend, "        dT_p %.2f dQ_p %.2f" % (dT_p, dQ_p))
            p.print(prepend, "    ----------------------------")

        out = {"a": a, "aprime": aprime, "Cl": Cl, "Cd":Cd, "alpha": alpha, "phi": phi, "F": F, "dFt": dFt, "Ct": Ct, "dFn": dFn,
               "_airfoil": _airfoil_dat, "dT": dT, "dQ": dQ, "Re": Re, 'U1': U1, 'U2': U2, 'U3': U3, 'U4': U4, "lambda_r":lambda_r}
        return out

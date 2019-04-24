__author__ = "Miha Smrekar"
__credits__ = ["Miha Smrekar"]
__license__ = "GPL"
__version__ = "0.3.0"
__maintainer__ = "Miha Smrekar"
__email__ = "miha.smrekar9@gmail.com"
__status__ = "Development"

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


class Calculator:
    """
    Class for calculation of induction factors using BEM theory.
    """

    def __init__(self, airfoils):
        self.airfoils = airfoils
        for blade_name in self.airfoils:
            self.airfoils[blade_name]["alpha_zero"] = 0.0  # TODO FIX
            generate_dat(blade_name, self.airfoils[blade_name]["x"], self.airfoils[blade_name]["y"])

            data = self.airfoils[blade_name]["gathered_curves"]
            data = sort_data(data)

            re = data[:,0].flatten()
            alpha = data[:,2].flatten()
            cl = data[:,3].flatten()
            cd = data[:,4].flatten()

            def interpolation_function_cl(x,y,re=re,alpha=alpha,cl=cl):
                return interp(x,y,re,alpha,cl)
            def interpolation_function_cd(x,y,re=re,alpha=alpha,cd=cd):
                return interp(x,y,re,alpha,cd)
            
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
        p.print("-------------------------------------------------------------")
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
    def run_array(self, theta, B, c, r, foils, dr, R, Rhub, rpm, v, pitch, method, propeller_mode, print_out, tip_loss,
            hub_loss, new_tip_loss, new_hub_loss, cascade_correction, max_iterations, convergence_limit, rho,
            relaxation_factor, print_all, return_print, return_results, rotational_augmentation_correction,
            rotational_augmentation_correction_method, mach_number_correction, fix_reynolds, reynolds, *args,
            **kwargs, ):
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
        arrays = ["a", "a'", "cL", "alpha", "phi", "F", "dFt", "M", "TSR", "Ct", "dFn", "foils", "dT", "dQ", "Re", "U1", "U2", "U3", "U4"]
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

            if not print_all and not print_out:
                p.print("*", add_newline=False)

            if out_results == None:
                return None

            results["a"] = numpy.append(results["a"], out_results["a"])
            results["a'"] = numpy.append(results["a'"], out_results["aprime"])
            results["cL"] = numpy.append(results["cL"], out_results["Cl"])
            results["alpha"] = numpy.append(results["alpha"], out_results["alpha"])
            results["phi"] = numpy.append(results["phi"], out_results["phi"])
            results["F"] = numpy.append(results["F"], out_results["F"])
            results["dFt"] = numpy.append(results["dFt"], out_results["dFt"])
            results["Ct"] = numpy.append(results["Ct"], out_results["Ct"])
            results["dFn"] = numpy.append(results["dFn"], out_results["dFn"])
            results["foils"] = numpy.append(results["foils"], out_results["_airfoil"])
            results["dT"] = numpy.append(results["dT"], out_results["dT"])
            results["dQ"] = numpy.append(results["dQ"], out_results["dQ"])
            results["Re"] = numpy.append(results["Re"], out_results["Re"])
            results["U1"] = numpy.append(results["U1"], out_results["U1"])
            results["U2"] = numpy.append(results["U2"], out_results["U2"])
            results["U3"] = numpy.append(results["U3"], out_results["U3"])
            results["U4"] = numpy.append(results["U4"], out_results["U4"])

        if not print_all and not print_out:
            p.print("")

        dFt = results["dFt"]
        Ft = numpy.sum(dFt)
        M = B * dFt * r  # momenti po prerezih
        dQ = results["dQ"]
        Q = numpy.sum(dQ)  # moment for propeller
        power_p = Q * omega
        Msum = numpy.sum(M)
        power = numpy.sum(M) * omega
        p.print(power_p,power)
        Pmax = 0.5 * rho * v ** 3 * pi * R ** 2
        cp_w = power / Pmax
        cp_p = power_p / (rho * (rpm / 60) ** 3 * (2 * R) ** 5)

        dFn = results["dFn"]
        Fn = numpy.sum(dFn)
        dT = results["dT"]
        T = numpy.sum(dT)
        ct_w = T / (0.5 * rho * v ** 2 * pi * R ** 2)
        ct_p = T / (rho * (2 * R) ** 4 * (rpm / 60) ** 2)

        results["R"] = R
        results["rpm"] = rpm
        results["v"] = v
        results["cp_w"] = cp_w
        results["cp_p"] = cp_p
        results["ct_w"] = ct_w
        results["ct_p"] = ct_p
        results["TSR"] = TSR
        results["Ft"] = Ft
        results["r"] = r
        results["omega"] = omega
        results["M"] = M
        results["Msum"] = Msum
        results["power"] = power
        results["thrust"] = T
        results["dFt"] = dFt
        results["Rhub"] = Rhub
        results["B"] = B
        results["dr"] = dr
        results["c"] = c
        results["theta"] = theta
        results["J"] = J
        return results

    def calculate_section(self, v, omega, _r, _c, _theta, _dr, B, R, _airfoil_dat, _airfoil, max_thickness,
            propeller_mode, pitch=0.0, psi=0.0, fix_reynolds=False, reynolds=1e6, tip_loss=False, new_tip_loss=False,
            hub_loss=False, new_hub_loss=False, cascade_correction=False, rotational_augmentation_correction=False,
            rotational_augmentation_correction_method=0, mach_number_correction=False, method=5,
            kin_viscosity=1.4207E-5, rho=1.225, convergence_limit=0.001, max_iterations=100, relaxation_factor=0.3,
            printer=None, print_all=False, print_out=False, *args, **kwargs):
        # print(v)

        p = printer

        # local speed ratio
        lambda_r = omega * _r / v

        # solidity
        # sigma=_c*B/(2*pi*_r)

        # solidity
        sigma = _c * B / (2 * pi * _r)  # * abs(cos(_theta)) implemented in QBlade/XBEM/BDATA.cpp

        ## initial guess
        a = 0.01
        aprime = 0.01
        # a,aprime = guessInductionFactors(lambda_r,sigma,_theta)

        # iterations counter
        i = 0

        # tip mach number
        M = omega * _r / 343

        #convert pitch to radians
        _pitch = radians(pitch)

        ############ START ITERATION ############
        while True:
            # update counter
            i = i + 1

            # for pretty-printing only
            prepend = ""

            # wind components
            if propeller_mode:
                Un = v * (1 + a)
                Ut = omega * _r * (1 - aprime)
            else:
                Un = v * (1 - a)
                Ut = omega * _r * (1 + aprime)

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

            # if print_all:
            #    p.print("        Running xfoil for %s,Re=%s,alpha=%s" % (_airfoil_dat, Re, alpha))

            """
            xfoil_return = xfoil_runner(airfoil=_airfoil_dat, reynolds=Re, alpha=alpha, printer=p, print_all=print_all)

            if xfoil_return == False:
                p.print("        Xfoil failed")
                return None
            """

            # Cl, Cd = xfoil_return["CL"], xfoil_return["CD"] #direct xfoil calculation - no interpolation
            Cl, Cd = self.airfoils[_airfoil]["interp_function_cl"](Re, degrees(alpha)), self.airfoils[_airfoil][
                "interp_function_cd"](Re, degrees(alpha))
            #p.print(Re,Cl,Cd)

            if print_all:
                p.print("        CL:", Cl, "Cd:", Cd)

            if rotational_augmentation_correction:
                if print_all:
                    p.print("--")
                    p.print("  Cl:", Cl, "Cd:", Cd)
                Cl, Cd = calc_rotational_augmentation_correction(alpha=alpha, Cl=Cl, Cd=Cd, omega=omega, r=_r, R=R,
                                                                 c=_c, theta=_theta, v=v, Vrel=Vrel_norm,
                                                                 method=rotational_augmentation_correction_method, )

                if print_all:
                    p.print("  Cl_cor:", Cl, "Cd_cor:", Cd)
                    p.print("--")

            if mach_number_correction:
                Cl = machNumberCorrection(Cl, M)

            # normal and tangential coefficients
            C_norm = Cl * cos(phi) + Cd * sin(phi)
            C_tang = Cl * sin(phi) - Cd * cos(phi)

            # save old values, calculate new values of induction factors
            a_last = a
            aprime_last = aprime

            input_arguments = {"F": F, "lambda_r": lambda_r, "phi": phi, "sigma": sigma, "C_norm": C_norm,
                               "C_tang": C_tang, "Cl": Cl, "Cd": Cd, "B": B, "c": _c, "r": _r, "R": R, "psi": 0.0,
                               "aprime_last": aprime, "omega": omega, "v": v, "a_last": a_last,
                               # "alpha_zero": airfoils[_airfoil]["alpha_zero"],
                               "method": method, "alpha": alpha, "alpha_deg": degrees(alpha)}

            if print_all:
                args_to_print = sorted([key for key, value in input_arguments.items()])
                p.print("            i", i)
                for a in args_to_print:
                    p.print("            ", a, input_arguments[a])
                p.print("             --------")

            # calculate induction coefficients
            coeffs = calculate_coefficients(method, input_arguments)
            if coeffs == None:
                return None
            else:
                a, aprime, Ct = coeffs

            # force calculation
            dFL = Cl * 0.5 * rho * Vrel_norm ** 2 * _c * _dr  # lift force
            dFD = Cd * 0.5 * rho * Vrel_norm ** 2 * _c * _dr  # drag force
            dFt = dFL * sin(phi) - dFD * cos(phi)  # tangential force
            dFn = dFL * cos(phi) + dFD * sin(phi)  # normal force

            # thrust and torque - Wiley, WE 2nd, p.124
            dT_MT = F * 4 * pi * _r * rho * v ** 2 * a * (1 - a) * _dr
            dT_BET = 0.5 * rho * B * _c * Vrel_norm ** 2 * (Cl * cos(phi) + Cd * sin(phi)) * _dr
            dQ_MT = F * 4 * aprime * (1 - a) * rho * v * pi * _r ** 3 * omega * _dr
            dQ_BET = B * 0.5 * rho * Vrel_norm ** 2 * (Cl * sin(phi) - Cd * cos(phi)) * _c * _dr * _r

            # thrust and torque from https://apps.dtic.mil/dtic/tr/fulltext/u2/1013408.pdf
            dT_p = B * rho * (omega * _r / cos(phi) * cos(_theta)) ** 2 * _c * _dr * (Cl * cos(phi) - Cd * sin(phi))
            dQ_p = B * rho * (omega * _r / cos(phi) * cos(_theta)) ** 2 * _c * _r * _dr * (
                    Cl * sin(phi) + Cd * cos(phi))

            # thrust and torque from http://www.icas.org/ICAS_ARCHIVE/ICAS2010/PAPERS/434.PDF
            # dT_p = sigma*pi*rho*v**2*(1+a)**2/(sin(phi)**2)*(Cl*cos(phi)-Cd*sin(phi))*_r*_dr
            # dQ_p = sigma*pi*rho*v**2*(1+a)**2/(sin(phi)*+2)*(Cl*sin(phi)+Cd*cos(phi))*_r**2*_dr

            # thrust-propeller
            dT_MT_p = 4 * pi * _r * rho * v ** 2 * (1 + a) * a * _dr
            dQ_MT_p = 4 * pi * _r ** 3 * rho * v * omega * (1 + a) * aprime * _dr
            dT_BET_p = 0.5 * rho * v ** 2 * _c * B * (1 + a) ** 2 / (sin(phi) ** 2) * (
                    Cl * cos(phi) - Cd * sin(phi)) * _dr
            dQ_BET_p = 0.5 * rho * v * _c * B * omega * _r ** 2 * (1 + a) * (1 - aprime) / (sin(phi) * cos(phi)) * (
                    Cl * sin(phi) + Cd * cos(phi)) * _dr

            if propeller_mode:
                dT = dT_BET_p
                dQ = dQ_BET_p
            else:
                dT = dT_BET
                dQ = dQ_BET

            # wind after
            if propeller_mode:
                #p.print('propeller_mode_U')
                U1 = v
                U2 = None
                U3 = U1*(1+a)
                U4 = U1*(1+2*a)
            else:
                #p.print('wind turbine mode_a')
                U1 = v
                U2 = U1*(1-a)
                U3 = None
                U4 = U1 * (1 - 2 * a)

            # check convergence
            if abs(a - a_last) < convergence_limit:
                break

            # p.print("dT_MT %.2f dT_BET %.2f" % (dT_MT,dT_BET))

            # check iterations limit
            if i >= max_iterations:
                if print_out:
                    p.print("-*-*-*-*-*-*-*-*-*-*-*-*-*-\n", "|max iterations exceeded\n", "|------>a:", a, " aprime",
                            aprime, )
                    prepend = "|"
                return None

            # relaxation
            a = a_last + relaxation_factor * (a - a_last)  # aprime=aprime_last+relaxation_factor*(aprime-aprime_last)

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
            p.print(prepend, "        dT_MT %.2f dT_BET %.2f" % (dT_MT, dT_BET))
            p.print(prepend, "        dQ_MT %.2f dQ_BET %.2f" % (dQ_MT, dQ_BET))
            p.print(prepend, "        dT_MT_p %.5f dT_BET_p %.5f" % (dT_MT_p, dT_BET_p))
            p.print(prepend, "        dQ_MT_p %.5f dQ_BET_p %.5f" % (dQ_MT_p, dQ_BET_p))
            p.print(prepend, "        dT_p %.2f dQ_p %.2f" % (dT_p, dQ_p))
            p.print(prepend, "    ----------------------------")

        out = {"a": a, "aprime": aprime, "Cl": Cl, "alpha": alpha, "phi": phi, "F": F, "dFt": dFt, "Ct": Ct, "dFn": dFn,
               "_airfoil": _airfoil_dat, "dT": dT, "dQ": dQ, "Re": Re, 'U1': U1, 'U2':U2, 'U3':U3, 'U4':U4}
        return out

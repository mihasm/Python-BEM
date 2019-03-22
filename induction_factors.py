__author__ = "Miha Smrekar"
__credits__ = ["Miha Smrekar"]
__license__ = "GPL"
__version__ = "0.3.0"
__maintainer__ = "Miha Smrekar"
__email__ = "miha.smrekar9@gmail.com"
__status__ = "Development"

import numbers
from math import sin, cos, atan, acos, pi, exp, sqrt, radians, atan2, degrees, tan
from scipy import interpolate

import numpy
from utils import Printer
from popravki import *
from xfoil import run_xfoil_analysis

numpy.seterr(all="raise")
numpy.seterr(invalid="raise")

def generate_dat(name,x,y):
    out = ""
    out += name+"\n"
    #out += "   x/c        y/c\n"
    #array_dat = self.table_dat.get_values()
    for i in range(len(x)):
            _x = float(x[i])
            _y = float(y[i])
            if _y >= 0:
                out += "%.6f   %.6f\n" % (_x,_y)
            else:
                out += "%.6f  %.6f\n" % (_x,_y)

    print(out)
    f = open(name+".dat","w")
    f.write(out)
    f.close()
    return out

def xfoil_runner(airfoil,reynolds,alpha,printer):
    #printer.print("xfoil runner")
    out = run_xfoil_analysis(airfoil,reynolds,alpha)
    if out == False:
        printer.print("    Convergence failed, modifying parameters...")
        original_alpha = alpha
        original_reynolds = reynolds
        while True:
            razmak_x = 0.0

            while True:
                if razmak_x >= 5:
                    break
                razmak_x += 0.01
                alpha_spodnji = original_alpha-razmak_x
                alpha_zgornji = original_alpha+razmak_x
                out_spodnji = run_xfoil_analysis(airfoil,reynolds,alpha_spodnji)
                out_zgornji = run_xfoil_analysis(airfoil,reynolds,alpha_zgornji)
                if out_spodnji != False and out_zgornji != False:
                    CL_spodnji = out_spodnji["CL"]
                    CL_zgornji = out_zgornji["CL"]
                    dy = CL_zgornji-CL_spodnji
                    dx = alpha_zgornji-alpha_spodnji
                    k = dy/dx
                    n = CL_zgornji-k*alpha_zgornji
                    CL_interpoliran = k*original_alpha+n

                    CD_spodnji = out_spodnji["CD"]
                    CD_zgornji = out_zgornji["CD"]
                    dy_2 = CD_zgornji - CD_spodnji
                    k_2 = dy_2/dx
                    n = CD_zgornji-k_2*alpha_zgornji
                    CD_interpoliran = k*original_alpha+n
                    return {"CD":CD_interpoliran,"CL":CL_interpoliran,"out":out_spodnji["out"]+out_zgornji["out"]}

            reynolds = reynolds + 1e3


    else:
        return out


class Calculator:
    """
    Class for calculation of induction factors using BEM theory.
    """

    def __init__(self, curves):
        self.curves = curves
        for blade_name in self.curves:
            self.curves[blade_name]["alpha_zero"] = 0.0 #TODO FIX
            generate_dat(blade_name,self.curves[blade_name]["x"],self.curves[blade_name]["y"])


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
        if (
                isinstance(theta, numpy.ndarray)
                and isinstance(c, numpy.ndarray)
                and isinstance(r, numpy.ndarray)
        ):
            return theta, c, r
        else:
            if (
                    isinstance(theta, numbers.Real)
                    and isinstance(c, numbers.Real)
                    and isinstance(r, numbers.Real)
            ):
                return numpy.array([theta]), numpy.array([c]), numpy.array([r])
            return None

    # noinspection PyUnusedLocal,PyUnusedLocal
    def run_array(
            self,
            theta,
            B,
            c,
            r,
            foils,
            dr,
            R,
            Rhub,
            rpm,
            v,
            method,
            print_out,
            tip_loss,
            hub_loss,
            new_tip_loss,
            new_hub_loss,
            cascade_correction,
            max_iterations,
            convergence_limit,
            rho,
            relaxation_factor,
            print_all,
            return_print,
            return_results,
            rotational_augmentation_correction,
            rotational_augmentation_correction_method,
            *args,
            **kwargs,
    ):
        """
        Calculates induction factors using standard iteration methods.

        Different methods are available as different fInductionCoefficients functions.

        ANGLES REPRESENTATION SHOWN IN
        https://cmm2017.sciencesconf.org/129068/document
        alpha - angle of attack
        phi - angle of relative wind
        beta - theta

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
        arrays = ["a", "a'", "cL", "alpha",
                  "phi", "F", "dFt", "M", "TSR", "Ct", "dFn", "foils"]
        for array in arrays:
            results[array] = numpy.array([])

        # set constants that are section-independent
        omega = rpm * 2 * pi / 60
        TSR = omega * R / v  # tip speed ratio
        kin_viscosity = 1.4207E-5 # Kinematic viscosity

        for n in range(len(theta)):
            # grab local radius, chord length and twist angle
            _r = r[n]
            _c = c[n]
            _theta = radians(theta[n])
            _airfoil = foils[n]

            #get max thickness
            max_thickness = self.curves[_airfoil]["max_thickness"]*_c

            # local speed ratio
            lambda_r = omega * _r / v

            # solidity
            # sigma=_c*B/(2*pi*_r)

            # solidity - implemented in QBlade/XBEM/BDATA.cpp
            sigma = _c * B / (2 * pi * _r) * abs(cos(_theta))

            # Coning angle (PROPX: Definitions,Derivations, Data Flow, p.22)
            psi = 0.0

            # initial guess
            a = 0.1
            aprime = 0.01
            #a,aprime = guessInductionFactors(lambda_r,sigma,_theta,self.f_c_L)

            # iterations counter
            i = 0

            ############ START ITERATION ############
            while True:
                # update counter
                i = i + 1

                # for pretty-printing only
                prepend = ""

                # wind components
                Ut = omega * _r * (1 + aprime)
                Un = v * (1 - a)
                Vrel_norm = sqrt(Un ** 2 + Ut ** 2)
                Re = int(Vrel_norm*_c/kin_viscosity)

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
                alpha = phi - _theta

                if cascade_correction:
                    alpha = cascadeEffectsCorrection(
                        alpha=alpha,
                        v=v,
                        omega=omega,
                        r=_r,
                        R=R,
                        c=_c,
                        B=B,
                        a=a,
                        aprime=aprime,
                        tmax=max_thickness
                    )

                # lift and drag coefficients
                #Cl, Cd = self.f_c_L(degrees(alpha)), self.f_c_D(degrees(alpha))
                #Cl,Cd = self.curves[_airfoil]["f_c_L"](degrees(alpha)), self.curves[_airfoil]["f_c_D"](degrees(alpha))
                if print_out or print_all:
                    p.print("    Running xfoil for %s,Re=%s,alpha=%s" % (_airfoil,Re,alpha))

                xfoil_return = xfoil_runner(airfoil=_airfoil+".dat",reynolds=Re,alpha=alpha,printer=p)

                if print_all and xfoil_return != False:
                    p.print(xfoil_return["out"])
                if xfoil_return == False:
                    p.print("    Xfoil failed")
                    return None
                    break
                Cl,Cd = xfoil_return["CL"],xfoil_return["CD"]
                #p.print(xfoil_return["out"])

                if rotational_augmentation_correction:
                    if print_all:
                        p.print("--")
                        p.print("  Cl:", Cl, "Cd:", Cd)
                    Cl, Cd = calc_rotational_augmentation_correction(
                        alpha=alpha,
                        alpha_zero=self.curves[_airfoil]["alpha_zero"],
                        Cl=Cl,
                        Cd=Cd,
                        omega=omega,
                        r=_r,
                        R=R,
                        c=_c,
                        theta=_theta,
                        v=v,
                        Vrel=Vrel_norm,
                        method=rotational_augmentation_correction_method,
                    )

                    if print_all:
                        p.print("  Cl_cor:", Cl, "Cd_cor:", Cd)
                        p.print("--")

                # normal and tangential coefficients
                cn = Cl * cos(phi) + Cd * sin(phi)
                ct = Cl * sin(phi) - Cd * cos(phi)


                # save old values, calculate new values of induction factors
                a_last = a
                aprime_last = aprime

                input_arguments = {
                    "F": F,
                    "lambda_r": lambda_r,
                    "phi": phi,
                    "sigma": sigma,
                    "cn": cn,
                    "ct": ct,
                    "Cl": Cl,
                    "B": B,
                    "c": _c,
                    "r": _r,
                    "R": R,
                    "psi": 0.0,
                    "aprime_last": aprime,
                    "omega": omega,
                    "v": v,
                    "a_last": a_last,
                    "alpha_zero": self.curves[_airfoil]["alpha_zero"],
                    "method": method,
                }

                if print_all:
                    args_to_print = sorted(
                        [key for key, value in input_arguments.items()]
                    )
                    p.print("            i", i)
                    for a in args_to_print:
                        p.print("            ", a, input_arguments[a])
                    p.print("             --------")

                # calculate induction coefficients
                a, aprime, Ct = calculate_coefficients(method, input_arguments)

                # force calculation
                dFL = Cl * 0.5 * rho * Vrel_norm ** 2 * \
                      _c * dr[n]  # lift force
                dFD = Cd * 0.5 * rho * Vrel_norm ** 2 * \
                      _c * dr[n]  # drag force
                dFt = dFL * sin(phi) - dFD * cos(phi)  # tangential force
                dFn = dFL * cos(phi) + dFD * sin(phi)  # normal force

                # check convergence
                if abs(a - a_last) < convergence_limit:
                    break

                # check iterations limit
                if i >= max_iterations:
                    if print_out:
                        p.print("-*-*-*-*-*-*-*-*-*-*-*-*-*-\n",
                                "|max iterations exceeded\n",
                                "|------>a:", a, " aprime", aprime,
                                )
                        prepend = "|"
                    return None
                    break

                # relaxation
                a = a_last + relaxation_factor * (a - a_last)
                #aprime=aprime_last+relaxation_factor*(aprime-aprime_last)

            ############ END ITERATION ############
            

            if print_out:
                p.print(prepend, "    r", _r)
                p.print(prepend, "        iters: ", i)
                p.print(prepend, "        phi: ", degrees(phi))
                p.print(prepend, "        _theta: ", degrees(_theta))
                p.print(prepend, "        alpha: ", degrees(alpha)), "Cl", str(Cl)
                p.print(prepend, "        a: ", a, "a'", str(aprime))
                p.print(prepend, "        dFt: ", dFt)
                p.print(prepend, "        LSR: ", lambda_r)
                p.print(prepend, "        Ct: ", Ct)
                p.print(prepend, "        Vrel: ", Vrel_norm)
                p.print(prepend, "        foil:",_airfoil)
                p.print(prepend, "    ----------------------------")

            results["a"] = numpy.append(results["a"], a)
            results["a'"] = numpy.append(results["a'"], aprime)
            results["cL"] = numpy.append(results["cL"], Cl)
            results["alpha"] = numpy.append(results["alpha"], alpha)
            results["phi"] = numpy.append(results["phi"], phi)
            results["F"] = numpy.append(results["F"], F)
            results["dFt"] = numpy.append(results["dFt"], dFt)
            results["Ct"] = numpy.append(results["Ct"], Ct)
            results["dFn"] = numpy.append(results["dFn"], dFn)
            results["foils"] = numpy.append(results["foils"], _airfoil)

        dFt = results["dFt"]
        Ft = numpy.sum(dFt)
        M = B * dFt * r  # momenti po prerezih
        Msum = numpy.sum(M)
        power = numpy.sum(M) * omega
        Pmax = 0.5 * rho * v ** 3 * pi * R ** 2
        cp = power / Pmax

        dFn = results["dFn"]
        Fn = numpy.sum(dFn)
        T=Fn*B
        ct = T/(0.5*rho*v**2*pi*R**2)

        results["R"] = R
        results["rpm"] = rpm
        results["v"] = v
        results["cp"] = cp
        results["ct"] = ct
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
        return results

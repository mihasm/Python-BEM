__author__ = "Miha Smrekar"
__credits__ = ["Miha Smrekar"]
__license__ = "GPL"
__version__ = "0.2.9"
__maintainer__ = "Miha Smrekar"
__email__ = "miha.smrekar9@gmail.com"
__status__ = "Development"

import numbers
from math import sin, cos, atan, acos, pi, exp, sqrt, radians, atan2, degrees, tan

import numpy
from utils import Printer

numpy.seterr(all="raise")
numpy.seterr(invalid="raise")


def fTipLoss(B, r, R, phi):
    """
    Prandtl tip loss.
    :param B: number of blades
    :param r: section radius [m]
    :param R: tip radius [m]
    :param phi: angle of relative wind [rad]
    :return: returns tip loss factor [float]
    """
    F = 1
    f = sin(phi)
    g = (R - r) / r
    Ft = 2 / pi * acos(exp(-B / 2 * abs(g / f)))
    F = F * Ft
    return F


def fHubLoss(B, r, Rhub, phi):
    """
    Prandtl hub loss.
    :param B: number of blades
    :param r: section radius [m]
    :param Rhub: hub radius [m]
    :param phi: angle of relative wind [rad]
    :return: returns hub loss factor [float]
    """
    F = 1
    f = sin(phi / 360 * 2 * pi)
    g = (r - Rhub) / r
    Fr = 2 / pi * acos(exp(-B / 2 * abs(g / f)))
    F = F * Fr
    return F


def newTipLoss(B, r, R, phi, lambda_r):
    """
    Prandtl tip loss with correction.
    :param B: number of blades
    :param r: section radius [m]
    :param R: tip radius [m]
    :param phi: angle of relative wind [rad]
    :param lambda_r: local speed ratio [float]
    :return: returns tip loss factor [float]
    """
    F = 1
    f = sin(phi)
    g = (R - r) / r
    Flt = (
            2
            / pi
            * acos(exp(-B / 2 * abs(g / f) * (exp(-0.15 * (B * lambda_r - 21)) + 0.1)))
    )
    F = F * Flt
    return F


def newHubLoss(B, r, Rhub, phi, lambda_r):
    """
    Prandtl hub loss.
    :param lambda_r: local speed ratio [float]
    :param B: number of blades [int]
    :param r: section radius [m]
    :param Rhub: hub radius [m]
    :param phi: angle of relative wind [rad]
    :return: returns hub loss factor [float]
    """
    F = 1
    f = sin(phi)
    g = (Rhub - r) / r
    Flt = (
            2
            / pi
            * acos(exp(-B / 2 * abs(g / f) * (exp(-0.15 * (B * lambda_r - 21)) + 0.1)))
    )
    F = F * Flt
    return F


def newLosses(cn, ct, B, r, R, phi, lambda_r, Rhub=None):
    """
    Combines Prandtl tip and hub losses with corrections.
    :param cn: normal coefficient [float]
    :param ct: thrust coefficient [float]
    :param B: number of blades [int]
    :param r: section radius [m]
    :param R: tip radius [m]
    :param phi: angle of relative wind [rad]
    :param lambda_r: local speed ratio [float]
    :param Rhub: hub radius [m]
    :return: cn,ct with included losses
    """
    tiploss = newTipLoss(B, r, R, phi, lambda_r)
    if Rhub:
        hubloss = newHubLoss(B, r, Rhub, phi, lambda_r)
    else:
        hubloss = 1.0
    Fl = tiploss * hubloss
    cn = cn * Fl
    ct = ct * Fl
    return cn, ct


# noinspection PyUnusedLocal,PyUnusedLocal
def fInductionCoefficients0(F, phi, sigma, cn, ct, *args, **kwargs):
    """
    METHOD FROM http://orbit.dtu.dk/files/86307371/A_Detailed_Study_of_the_Rotational.pdf
    Basically, original method without any corrections.

    :param F: loss factors
    :param phi: relative wind [rad]
    :param sigma: solidity
    :param cn: normal coefficient
    :param ct: thrust coefficient
    :return: axial induction factor, tangential induction factor
    """
    a = (sigma * cn) / (4 * F * sin(phi) ** 2 + sigma * cn)
    aprime = (sigma * ct) / (4 * F * sin(phi) * cos(phi) - sigma * ct)
    Ct = 4*a*(1-a)
    return a, aprime, Ct


# noinspection PyUnusedLocal,PyUnusedLocal
def fInductionCoefficients1(F, phi, sigma, cn, ct, *args, **kwargs):
    """
    Calculates induction coefficients using method using Spera correction only
    https://cmm2017.sciencesconf.org/129068/document

    :param F: loss factors
    :param phi: relative wind [rad]
    :param sigma: solidity
    :param cn: normal coefficient
    :param ct: thrust coefficient
    :return: axial induction factor, tangential induction factor
    """
    a = (sigma * cn) / (4 * F * sin(phi) ** 2 + sigma * cn)
    Ct = 4*a*(1-a)*F

    # Spera's correction
    if a >= 0.2:
        ac = 0.2
        K = (4 * F * sin(phi) ** 2) / (sigma * cn)
        to_sqrt = abs((K * (1 - 2 * ac) + 2) ** 2 + 4 * (K * ac ** 2 - 1))
        a = 1 + 0.5 * K * (1 - 2 * ac) - 0.5 * sqrt(to_sqrt)
        Ct = 4*(ac**2+(1-2*ac)*a)*F

    aprime = (sigma * ct) / (4 * F * sin(phi) * cos(phi) - sigma * ct)
    return a, aprime, Ct


# noinspection PyUnusedLocal,PyUnusedLocal,PyUnusedLocal
def fInductionCoefficients2(F, phi, sigma, cn, Cl, *args, **kwargs):
    """
    Calculates induction coefficients using method from
    Wiley,Wind Energy, p.126 (Method 2) == p.128

    :param Cl: lift coefficient
    :param F: loss factors
    :param phi: relative wind [rad]
    :param sigma: solidity
    :param cn: normal coefficient
    :return: axial induction factor, tangential induction factor
    """
    a = 1 / (1 + 4 * F * sin(phi) ** 2 / (sigma * Cl * cos(phi)))
    aprime = 1 / (4 * cos(phi) / (sigma * Cl) - 1)
    Ct = 4*a*(1-a)*F
    return a, aprime, Ct


# noinspection PyUnusedLocal,PyUnusedLocal,PyUnusedLocal
def fInductionCoefficients3(lambda_r, phi, sigma, Cl, *args, **kwargs):
    """
    Calculates induction coefficients using method from
    Wind Turbine Blade Analysis, Grant Ingram, 2011

    :param lambda_r: local speed ratio
    :param Cl: lift coefficient
    :param phi: relative wind [rad]
    :param sigma: solidity
    :return: axial induction factor, tangential induction factor
    """

    a = (1 + (4 * cos(phi) ** 2) / (sigma * Cl * sin(phi))) ** -1
    aprime = ((sigma * Cl) / (4 * lambda_r * cos(phi))) * (1 - a)
    Ct = 4*a*(1-a)
    return a, aprime, Ct


# noinspection PyUnusedLocal,PyUnusedLocal,PyUnusedLocal
def fInductionCoefficients4(a_last, F, phi, sigma, cn, Cl, *args, **kwargs):
    """
    Calculates induction coefficients using method from
    Wind Energy Explained, Wiley, p.136
    Turbulent Wake State Modeling

    :param cn: normal coefficient
    :param F: loss factors
    :param a_last: axial induction factor
    :param Cl: lift coefficient
    :param phi: relative wind [rad]
    :param sigma: solidity
    :return: axial induction factor, tangential induction factor
    """

    Ct = sigma * (1 - a_last) ** 2 * cn / (sin(phi) ** 2)
    if Ct <= 0.96 * F:
        a = 1 / (1 + (4 * F * sin(phi) ** 2 / (sigma * Cl * cos(phi))))
    else:
        a = (1 / F) * (0.143 + sqrt(abs(0.0203 - 0.6427 * (0.889 - Ct))))
    aprime = 1 / ((4 * F * cos(phi) / (sigma * Cl)) - 1)
    return a, aprime, Ct


# noinspection PyUnusedLocal,PyUnusedLocal,PyUnusedLocal
def fInductionCoefficients5(
        aprime_last, F, lambda_r, phi, sigma, Cl, B, c, psi, r, R, v, omega, *args, **kwargs
):
    """
    Calculates induction coefficients using method from
    PROPX: Definitions, Derivations and Data Flow, C. Harman, 1994

    :param omega: rotational velocity [rad/s]
    :param v: wind speed [m/s]
    :param R: tip radius [m]
    :param r: section radius [m]
    :param psi: coning angle [rad]
    :param c: chord length [m]
    :param B: number of blades
    :param lambda_r: local speed ratio
    :param aprime_last: tangential induction factor
    :param F: loss factors
    :param Cl: lift coefficient
    :param phi: relative wind [rad]
    :param sigma: solidity
    :return: axial induction factor, tangential induction factor
    """

    K = B * c * Cl * cos(phi) * cos(psi) / (8 * pi * r * sin(phi) ** 2)
    a = K / (1 + K)
    CTL = 4 * a * F * (1 - a)
    
    if a > 0.2:
        ac = 0.2
        acprime = 0.0
        Fc = F
        if Fc == 0.0:
            Fc = 0.1
        X = lambda_r
        CTL = B * c * Cl * X * (1 + aprime_last) * cos(psi) / (2 * pi * R) * sqrt(
            (1 - a) ** 2 + (1 + aprime_last) ** 2 * X ** 2)
        CTL_ac = 4 * ac * Fc * (1 - ac)
        cpsi = cos(psi)
        #clen_4 = (1 + (acprime * (1 - 2 * ac)) / ((1 + 2 * acprime) * ac * cos(psi) ** 2))
        #clen_2 = 4 * ac * B / pi * 1 / tan(pi * Fc / 2) 
        #clen_3 = ((R * cpsi - r * cpsi) / (R * cpsi))
        #clen_5 = (1 - ac)
        dCTL_da = 4 * Fc * (1 - 2 * ac) + \
                  4 * ac * B / pi * 1 / tan(pi * Fc / 2) * \
                  ((R * cpsi - r * cpsi) / (R * cpsi)) * \
                  (cos(phi) ** 2 / sin(phi)) * \
                  (1 + (acprime * (1 - 2 * ac)) / ((1 + 2 * acprime) * ac * cos(psi) ** 2)) * \
                  (1 - ac)
        a = ac - (CTL_ac - CTL) / dCTL_da
    
    XL = (r * cos(psi) * omega) / v

    to_sqrt = abs(1 + (4 * a * (1 - a)) / (XL ** 2))
    aprime = 0.5 * (sqrt(to_sqrt) - 1)
    # aprime=0
    # aprime = 1 / ((4 * F * cos(phi) / (sigma * Cl)) - 1)
    return a, aprime, CTL


# noinspection PyUnusedLocal,PyUnusedLocal,PyUnusedLocal
def fInductionCoefficients6(a_last, F, phi, sigma, cn, ct, Cl, *args, **kwargs):
    """
    Calculates induction coefficients using method from
    AeroDyn manual - theory.

    :param cn: normal coefficient
    :param a_last: axial induction factor
    :param F: loss factors
    :param Cl: lift coefficient
    :param phi: relative wind [rad]
    :param sigma: solidity
    :return: axial induction factor, tangential induction factor
    """

    CT=(1+(sigma*(1-a_last)**2*cn)/(sin(phi)**2))
    if CT > 0.96 * F:
        # Modified Glauert correction
        a = (
                    18 * F - 20 - 3 * sqrt(CT * (50 - 36 * F) +
                                           12 * F * (3 * F - 4))
            ) / (36 * F - 50)
        CT = 8/9+(4*F-40/90)*a+(50/9-4*F)*a**2
    else:
        a = (1 + 4 * F * sin(phi) ** 2 / (sigma * cn)) ** -1
    aprime = (-1 + 4 * F * sin(phi) * cos(phi) / (sigma * ct)) ** -1
    return a, aprime, CT


# noinspection PyUnusedLocal,PyUnusedLocal,PyUnusedLocal
def fInductionCoefficients7(a_last,F, lambda_r, phi, sigma, cn, *args, **kwargs):
    """
    Calculates induction coefficients using method from
    QBlade/src/XBEM/BData.cpp

    :param lambda_r: local speed ratio
    :param Ct: local thrust coefficient
    :param cn: normal coefficient
    :param F: loss factors
    :param phi: relative wind [rad]
    :param sigma: solidity
    :return: axial induction factor, tangential induction factor
    """

    Ct = sigma * (1 - a_last) ** 2 * cn / (sin(phi) ** 2)  # Qblade

    if Ct <= 0.96 * F:
        a = 1 / (4 * F * sin(phi) ** 2 / (sigma * cn) + 1)
    else:
        a = (
                    18 * F - 20 - 3 * abs(Ct * (50 - 36 * F) +
                                          12 * F * (3 * F - 4)) ** 0.5
            ) / (36 * F - 50)

    aprime = 0.5 * (abs(1 + 4 / (lambda_r ** 2) * a * (1 - a)) ** 0.5 - 1)

    return a, aprime, Ct


# noinspection PyUnusedLocal,PyUnusedLocal,PyUnusedLocal
def fInductionCoefficients8(F, phi, sigma, cn, Cl, ct, *args, **kwargs):
    """
    Calculates induction coefficients using method from
    http://www.i-asem.org/publication_conf/asem17/1.SM/W4A.2.SM2143_4195F1.pdf

    :param Cl: lift coefficient
    :param F: loss factors
    :param phi: relative wind [rad]
    :param sigma: solidity
    :param cn: normal coefficient
    :return: axial induction factor, tangential induction factor
    """
    a = (1 + 4 * F * sin(phi) ** 2 / (sigma * cn)) ** -1
    Ct = 4*a*(1-a*F)*F

    ac = 1/3

    # Shen's correction
    if a > ac:
        Y1=a**-1
        a = (2+Y1-sqrt(4*Y1*(1-F)+Y1**2))/(2*(1+F*Y1))
        Ct = 4*(ac**2*F+(1-2*ac*F)*a)*F

    aprime = 1 / (4 * F * sin(phi) * cos(phi)/ (sigma * ct) + 1)
    return a, aprime, Ct

def fInductionCoefficients9(a_last,F,phi,Cl,ct,sigma,*args, **kwargs):
    """
    Method from Pratumnopharat,2010

    Glauert method
    """
    a = a_last
    if a >= 1/3:
        Ct = 4*a*F*(1-a)
    else:
        Ct = 4*a*F*(1-0.25*(5-3*a)*a)

    if Ct <= 0.888*F:
        a=(1-sqrt(1-Ct/F))/2
    else:
        Y=(sqrt(1/36*(Ct/F)**2-145/2187*Ct/F+92/2187)+Ct/(6*F)-145/729)**(1/3)
        a = Y-11/(81*Y)+5/9
    aprime = (4*F*sin(phi)*cos(phi)/(sigma*ct)-1)**-1
    return a,aprime,Ct

def fInductionCoefficients10(a_last,F,phi,Cl,ct,sigma,*args, **kwargs):
    """
    Method from Pratumnopharat,2010

    Wilson and Walker method
    """
    a=a_last

    ac = 0.2

    if a <= ac:
        Ct = 4*a*F*(1-a)
    else:
        Ct = 4*F*(ac**2+(1-2*ac)*a)

    if Ct <= 0.64*F:
        a = (1-sqrt(1-Ct/F))/2
    else:
        a = (Ct-4*F*ac**2)/(4*F*(1-2*ac))

    aprime = (4*F*sin(phi)*cos(phi)/(sigma*ct)-1)**-1

    return a,aprime,Ct


def fInductionCoefficients11(a_last,F,phi,Cl,ct,cn,sigma,*args, **kwargs):
    """
    Method from Pratumnopharat,2010

    Classical momentum brake state model
    """
    a=a_last
    Sw = sigma/(8*sin(phi)**2)*cn
    Ct = Sw*(1-a)**2
    a=(2*Sw+F-sqrt(F**2+4*Sw*F*(1-F)))/(2*(Sw+F**2))
    aprime = (4*F*sin(phi)*cos(phi)/(sigma*ct)-1)**-1
    return a,aprime,Ct

def fInductionCoefficients12(a_last,F,phi,Cl,ct,cn,sigma,lambda_r,*args, **kwargs):
    """
    Method from Pratumnopharat,2010

    Advanced brake state model

    Should be the same result as in Aerodyn or QBlade
    """
    a=a_last
    if a <= 0.4:
        Ct = 4*a*F*(1-a)
    else:
        Ct=8/9+(4*F-40/9)*a+(50/9-4*F)*a**2

    a = (18*F-20-3*sqrt(abs(ct*(50-36*F)+12*F*(3*F-4))))/(36*F-50)
    aprime = (4*F*sin(phi)*cos(phi)/(sigma*ct)-1)**-1
    return a, aprime, Ct


def guessInductionFactors(lambda_r, sigma, theta, f_c_L):
    """
    Provides initial guess to the function coefficients.

    Method from Wiley, Wind energy explained.

    :param lambda_r: local speed ratio
    :param sigma: solidity
    :param theta: twist [rad]
    :param f_c_L: function for lift calculation
    :return: axial induction factor, tangential induction factor
    """

    phi = (2 / 3) * atan(1 / lambda_r)
    Cl = f_c_L(degrees(theta))
    a = 1 / (1 + (4 * sin(phi) ** 2) / (sigma * Cl * cos(phi)))
    return a, aprime


def cascadeEffectsCorrection(alpha, v, omega, r, R, c, B, a, aprime):
    """
    Calculates cascade effects and corresponding change in angle of attack.
    Method from PROPX: Definitions, Derivations and Data Flow, C. Harman, 1994
    :param alpha: angle of attack [rad]
    :param v: wind speed [m/s]
    :param omega: rotational velocity [rad/s]
    :param r: section radius [m]
    :param R: tip radius [m]
    :param c: section chord length [m]
    :param B: number of blades [int]
    :param a: axial induction
    :param aprime: tangential induction
    :return: new angle of attack [rad]
    """

    tmax = 0.02
    delta_alpha_1 = (
            1
            / 4
            * (
                    atan((1 - a) * v / ((1 + a * aprime) * r * omega))
                    - atan(((1 - a) * v) / (r * omega))
            )
    )
    delta_alpha_2 = (
            0.109
            * (B * c * tmax * R * omega / v)
            / (R * c * sqrt((1 - a) ** 2 + (r * omega / v) ** 2))
    )
    out = alpha + delta_alpha_1 + delta_alpha_2

    return out


def calc_rotational_augmentation_correction(
        alpha, alpha_zero, Cl, Cd, omega, r, R, c, theta, v, Vrel, method=0
):
    """
    METHODS FROM http://orbit.dtu.dk/files/86307371/A_Detailed_Study_of_the_Rotational.pdf
    """

    fl = 0
    fd = 0

    if method == 1:
        # Snel et al.
        a_s = 3
        h = 2
        fl = a_s * (c / r) ** h
        fd = 0
    if method == 2:
        # Du & Selig
        gama = omega * R / sqrt(abs(v ** 2 - (omega * R) ** 2))
        ad, dd, bd = 1, 1, 1
        fl = (
                1
                / (2 * pi)
                * (
                        1.6
                        * (c / r)
                        / 0.1267
                        * (ad - (c / r) ** (dd * R / gama / r))
                        / (bd + (c / r) ** (dd * R / gama / r))
                        - 1
                )
        )
        fd = (
                -1
                / (2 * pi)
                * (
                        1.6
                        * (c / r)
                        / 0.1267
                        * (ad - (c / r) ** (dd * R / gama / r / 2))
                        / (bd + (c / r) ** (dd * R / gama / r / 2))
                        - 1
                )
        )
    if method == 3:
        # Chaviaropoulos and Hansen
        ah = 2.2
        h = 1.3
        n = 4
        fl = ah * (c / r) ** h * cos(theta) ** n
        fd = fl
    if method == 4:
        # Lindenburg
        al = 3.1
        h = 2
        fl = al * (omega * R / Vrel) ** 2 * (c / r) ** h
        fd = 0
    if method == 5:
        # Dumitrescu and Cardos
        gd = 1.25
        fl = 1 - exp(-gd / (r / c - 1))
        fd = 0
    # Cl_3D = Cl + fl*(2*pi*sin(alpha-alpha_zero)-Cl)
    Cl_3D = Cl + fl * Cl
    Cd_3D = Cd + fd * Cd
    return Cl_3D, Cd_3D


def calculate_coefficients(method, input_arguments):
    if method == 0:
        return fInductionCoefficients0(**input_arguments)
    if method == 1:
        return fInductionCoefficients1(**input_arguments)
    if method == 2:
        return fInductionCoefficients2(**input_arguments)
    if method == 3:
        return fInductionCoefficients3(**input_arguments)
    if method == 4:
        return fInductionCoefficients4(**input_arguments)
    if method == 5:
        return fInductionCoefficients5(**input_arguments)
    if method == 6:
        return fInductionCoefficients6(**input_arguments)
    if method == 7:
        return fInductionCoefficients7(**input_arguments)
    if method == 8:
        return fInductionCoefficients8(**input_arguments)
    if method == 9:
        return fInductionCoefficients9(**input_arguments)
    if method == 10:
        return fInductionCoefficients10(**input_arguments)
    if method == 11:
        return fInductionCoefficients11(**input_arguments)
    if method == 12:
        return fInductionCoefficients12(**input_arguments)
    raise Exception("Method "+str(method)+" does not exist.")


class Calculator:
    """
    Class for calculation of induction factors using BEM theory.
    """

    def __init__(self, curves):
        self.curves = curves
        #self.alpha_zero = self.inverse_f_c_L(0.0) #TODO CHECK

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

        for n in range(len(theta)):
            # grab local radius, chord length and twist angle
            _r = r[n]
            _c = c[n]
            _theta = radians(theta[n])
            _airfoil = foils[n]

            # local speed ratio
            lambda_r = omega * _r / v

            # solidity
            # sigma=_c*B/(2*pi*_r)

            # solidity - implemented in QBlade/XBEM/BDATA.cpp
            sigma = _c * B / (2 * pi * _r) * abs(cos(_theta))

            # Coning angle (PROPX: Definitions,Derivations, Data Flow, p.22)
            psi = 0.0

            # initial guess
            a = 0.0
            aprime = 0.0
            # a,aprime = guessInductionFactors(lambda_r,sigma,_theta,self.f_c_L)

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

                # relative wind
                phi = atan2(Un, Ut)

                F = 1
                # Prandtl tip loss
                if tip_loss:
                    F = F * fTipLoss(B, _r, R, phi)

                # Prandtl hub loss
                if hub_loss:
                    F = F * fHubLoss(B, _r, Rhub, phi)

                # New tip loss
                if new_tip_loss:
                    F = F * newTipLoss(B, _r, R, phi, lambda_r)

                # New hub loss
                if new_hub_loss:
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
                    )

                # lift and drag coefficients
                #Cl, Cd = self.f_c_L(degrees(alpha)), self.f_c_D(degrees(alpha))
                Cl,Cd = self.curves[_airfoil]["f_c_L"](degrees(alpha)), self.curves[_airfoil]["f_c_D"](degrees(alpha))

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
                    break

                # relaxation
                a = a_last + relaxation_factor * (a - a_last)
                # aprime=aprime_last+0.3*(aprime-aprime_last)

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
                p.print(prepend, "----------------------------")

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
        power = numpy.sum(M) * omega
        Pmax = 0.5 * rho * v ** 3 * pi * R ** 2
        cp = power / Pmax

        results["R"] = R
        results["rpm"] = rpm
        results["v"] = v
        results["cp"] = cp
        results["TSR"] = TSR
        results["Ft"] = Ft
        results["r"] = r
        results["omega"] = omega
        results["M"] = M
        results["power"] = power
        results["dFt"] = dFt
        results["Rhub"] = Rhub
        results["B"] = B
        results["dr"] = dr
        results["c"] = c
        results["theta"] = theta
        return results

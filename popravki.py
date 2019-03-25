import numbers
from math import sin, cos, atan, acos, pi, exp, sqrt, radians, atan2, degrees, tan
from scipy import interpolate

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
    Ct = sigma * (1 - a) ** 2 * cn / (sin(phi) ** 2)
    return a, aprime, Ct


# noinspection PyUnusedLocal,PyUnusedLocal
def fInductionCoefficients1(F, phi, sigma, cn, ct, *args, **kwargs):
    """
    Calculates induction coefficients using method using Spera correction only
    https://cmm2017.sciencesconf.org/129068/document

    Spera 1994


    :param F: loss factors
    :param phi: relative wind [rad]
    :param sigma: solidity
    :param cn: normal coefficient
    :param ct: thrust coefficient
    :return: axial induction factor, tangential induction factor
    """
    a = (sigma * cn) / (4 * F * sin(phi) ** 2 + sigma * cn)
    Ct = sigma * (1 - a) ** 2 * cn / (sin(phi) ** 2)

    # Spera's correction
    if a >= 0.2:
        ac = 0.2
        K = (4 * F * sin(phi) ** 2) / (sigma * cn)
        to_sqrt = abs((K * (1 - 2 * ac) + 2) ** 2 + 4 * (K * ac ** 2 - 1))
        a = 1 + 0.5 * K * (1 - 2 * ac) - 0.5 * sqrt(to_sqrt)
        #Ct = 4*(ac**2+(1-2*ac)*a)*F

    aprime = (sigma * ct) / (4 * F * sin(phi) * cos(phi) - sigma * ct)
    return a, aprime, Ct


# noinspection PyUnusedLocal,PyUnusedLocal,PyUnusedLocal
def fInductionCoefficients2(F, phi, sigma, cn, Cl, *args, **kwargs):
    """
    Calculates induction coefficients using method from
    Wiley,Wind Energy, p.126 (Method 2) == p.128

    #the same as original...

    :param Cl: lift coefficient
    :param F: loss factors
    :param phi: relative wind [rad]
    :param sigma: solidity
    :param cn: normal coefficient
    :return: axial induction factor, tangential induction factor
    """
    a = 1 / (1 + 4 * F * sin(phi) ** 2 / (sigma * Cl * cos(phi)))
    aprime = 1 / (4 * cos(phi) / (sigma * Cl) - 1)
    #Ct = 4*a*(1-a)*F
    Ct = sigma * (1 - a) ** 2 * cn / (sin(phi) ** 2)  # Qblade
    return a, aprime, Ct


# noinspection PyUnusedLocal,PyUnusedLocal,PyUnusedLocal
def fInductionCoefficients3(lambda_r, phi, sigma, Cl, cn, *args, **kwargs):
    """
    Calculates induction coefficients using method from
    Wind Turbine Blade Analysis, Grant Ingram, 2011
    
    #the same as original

    :param lambda_r: local speed ratio
    :param Cl: lift coefficient
    :param phi: relative wind [rad]
    :param sigma: solidity
    :return: axial induction factor, tangential induction factor
    """

    a = (1 + (4 * cos(phi) ** 2) / (sigma * Cl * sin(phi))) ** -1
    aprime = ((sigma * Cl) / (4 * lambda_r * cos(phi))) * (1 - a)
    #Ct = 4*a*(1-a)
    Ct = sigma * (1 - a) ** 2 * cn / (sin(phi) ** 2)  # Qblade
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
    if F == 0:
        F = 1e-6
    Ct = sigma * (1 - a_last) ** 2 * cn / (sin(phi) ** 2)
    if Ct <= 0.96 * F:
        a = 1 / (1 + (4 * F * sin(phi) ** 2 / (sigma * Cl * cos(phi))))
    else:
        a = (1 / F) * (0.143 + sqrt(abs(0.0203 - 0.6427 * (0.889 - Ct))))
    aprime = 1 / ((4 * F * cos(phi) / (sigma * Cl)) - 1)
    return a, aprime, Ct


# noinspection PyUnusedLocal,PyUnusedLocal,PyUnusedLocal
def fInductionCoefficients5(
        aprime_last, F, lambda_r, phi, sigma, Cl, cn, B, c, psi, r, R, v, omega, *args, **kwargs
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

    CTL = sigma * (1 - a) ** 2 * cn / (sin(phi) ** 2)  # Qblade
    
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

    This is equal to Advanced brake state model.


    :param cn: normal coefficient
    :param a_last: axial induction factor
    :param F: loss factors
    :param Cl: lift coefficient
    :param phi: relative wind [rad]
    :param sigma: solidity
    :return: axial induction factor, tangential induction factor
    """

    #CT=(1+(sigma*(1-a_last)**2*cn)/(sin(phi)**2))
    CT=((sigma*(1-a_last)**2*cn)/(sin(phi)**2))
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
def fInductionCoefficients8(a_last,F, phi, sigma, lambda_r, B, r, R, cn, Cl, ct, *args, **kwargs):
    """
    Method from Pratumnopharat,2010

    Shen's correction

    """
    if F == 0:
        F = 1e-6
    a=a_last
    Ct = sigma * (1 - a) ** 2 * cn / (sin(phi) ** 2)  # Qblade
    g=exp(-0.125*(B*lambda_r-21))+0.1
    try:
        F1 = (2/pi)*acos(exp(-g*B*(R-r)/(2*r*sin(phi))))
    except:
        F1 = 1.0
    Y1 = 4*F*sin(phi)**2/(sigma*F1*cn)
    Y2 = 4*F*sin(phi)*cos(phi)/(sigma*F1*ct)
    if Ct <= 0.888:
        a=(1-sqrt(1-Ct))/(2*F)
    else:
        a=(2+Y1-sqrt(4*Y1*(1-F)+Y1**2))/(2*(1+F*Y1))
    aprime = 1/((1-a*F)*Y2/(1-a)-1)

    return a, aprime, Ct

def fInductionCoefficients9(a_last,F,phi,Cl,cn, ct,sigma,*args, **kwargs):
    """
    Method from Pratumnopharat,2010

    Glauert method

    """
    
    a = a_last
    Ct = sigma * (1 - a_last) ** 2 * cn / (sin(phi) ** 2)
    if F == 0:
        F = 1e-6
    if Ct <= 0.888*F:
        a=(1-sqrt(1-Ct/F))/2
    else:
        Y=(sqrt(1/36*(Ct/F)**2-145/2187*Ct/F+92/2187)+Ct/(6*F)-145/729)**(1/3)
        a = Y-11/(81*Y)+5/9
    aprime = (4*F*sin(phi)*cos(phi)/(sigma*ct)-1)**-1
    return a,aprime,Ct

def fInductionCoefficients10(a_last,F,phi,Cl,ct, cn,sigma,*args, **kwargs):
    """
    Method from Pratumnopharat,2010

    Wilson and Walker method

    """
    a=a_last

    ac = 0.2

    Ct = sigma * (1 - a_last) ** 2 * cn / (sin(phi) ** 2)
    if F == 0:
        F = 1e-6
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
    Ct = sigma * (1 - a_last) ** 2 * cn / (sin(phi) ** 2)  # Qblade
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
    #if a <= 0.4:
    #    Ct = 4*a*F*(1-a)
    #else:
    #    Ct=8/9+(4*F-40/9)*a+(50/9-4*F)*a**2
    Ct = sigma * (1 - a_last) ** 2 * cn / (sin(phi) ** 2)  # Qblade

    a = (18*F-20-3*sqrt(abs(ct*(50-36*F)+12*F*(3*F-4))))/(36*F-50)
    aprime = (4*F*sin(phi)*cos(phi)/(sigma*ct)-1)**-1
    return a, aprime, Ct

def fInductionCoefficients13(a_last,F,phi,Cl,ct,cn,sigma,lambda_r,*args, **kwargs):
    """
    Method from Pratumnopharat,2010

    Modified advanced brake state model

    """
    a=a_last
    #if a <= 0.4:
    #    Ct = 4*a*F*(1-a)
    #else:
    #    Ct=8/9+(4*F-40/9)*a+(50/9-4*F)*a**2
    Ct = sigma * (1 - a_last) ** 2 * cn / (sin(phi) ** 2)  # Qblade
    if Ct < 0.96*F:
        a = (1-sqrt(1-Ct/F))/2
    else:
        a=0.1432+sqrt(-0.55106+0.6427*Ct/F)
    aprimeprime = (4*a*F*(1-a))/(lambda_r**2)
    if (1+aprimeprime) < 0:
        aprime = 0
    else:
        aprime = (sqrt(1+aprimeprime)-1)/2

    return a, aprime, Ct


def guessInductionFactors(lambda_r, sigma, theta):
    """
    Provides initial guess to the function coefficients.

    Method from Wiley, Wind energy explained (2nd ed.). (3.124,3.125,3.126)

    :param lambda_r: local speed ratio
    :param sigma: solidity
    :param theta: twist [rad]
    :return: axial induction factor, tangential induction factor
    """

    phi = (2 / 3) * atan(1 / lambda_r)
    #Cl = f_c_L(degrees(theta))
    Cl = 0.4
    a = 1 / (1 + (4 * sin(phi) ** 2) / (sigma * Cl * cos(phi)))
    aprime = (1-3*a)/(4*a-1)
    return a, aprime


def cascadeEffectsCorrection(alpha, v, omega, r, R, c, B, a, aprime, tmax):
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

    #tmax = 0.02
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
    if method == 13:
        return fInductionCoefficients13(**input_arguments)
    raise Exception("Method "+str(method)+" does not exist.")


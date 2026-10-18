# Python BEM - Blade Element Momentum Theory Software.

# Copyright (C) 2022 M. Smrekar

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from math import sin, cos, acos, pi, exp, sqrt, atan2, tan, degrees, radians, tanh

METHODS_STRINGS = {"0": "Classic BEM",
                   "1": "Spera Correction",
                   "2": "Buhl Correction (Aerodyn)",
                   "3": "Buhl Correction (QBlade)",
                   "4": "Wilson and Walker model",
                   "5": "Modified ABS model"}


def calculate_coefficients(method, input_arguments):
    """

    :param method:
    :param input_arguments:
    :return:
    """
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
    raise Exception("Method " + str(method) + " does not exist.")


def machNumberCorrection(Cl, Cd, M):
    """
    Prandtl-Glauert compressibility correction. It scales the pressure-driven lift only;
    the (mostly viscous) drag coefficient is left unchanged.

    :param Cl: Lift coefficient
    :param Cd: Drag coefficient
    :param M: local (relative) Mach number, must be < 1
    :return: Lift coefficient, drag coefficient
    """
    Cl = Cl / sqrt(1 - M ** 2)
    return Cl, Cd


def fTipLoss(B, r, R, phi):
    """
    Prandtl tip loss.
    :param B: number of blades
    :param r: section radius [m]
    :param R: tip radius [m]
    :param phi: angle of relative wind [rad]
    :return: returns tip loss factor [float]
    """
    # F = 1
    F = 2 / pi * acos(exp(-B / 2 * abs((R - r) / r / sin(phi))))
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
    f = sin(phi)
    g = (r - Rhub) / r
    Fr = 2 / pi * acos(exp(-B / 2 * abs(g / f)))
    F = F * Fr
    return F


def newTipLoss(B, r, R, phi, tsr):
    """
    Prandtl tip loss with the correction of Shen et al. (2005):
    g = exp(-0.125 * (B * lambda - 21)) + 0.1, where lambda is the rotor tip speed ratio.
    :param B: number of blades
    :param r: section radius [m]
    :param R: tip radius [m]
    :param phi: angle of relative wind [rad]
    :param tsr: tip speed ratio of the rotor [float]
    :return: returns tip loss factor [float]
    """
    F = 1
    f = sin(phi)
    g = (R - r) / r
    Flt = (
            2
            / pi
            * acos(exp(-B / 2 * abs(g / f) * (exp(-0.125 * (B * tsr - 21)) + 0.1)))
    )
    F = F * Flt
    return F


def newHubLoss(B, r, Rhub, phi, tsr):
    """
    Prandtl hub loss with the correction of Shen et al. (2005).
    :param tsr: tip speed ratio of the rotor [float]
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
            * acos(exp(-B / 2 * abs(g / f) * (exp(-0.125 * (B * tsr - 21)) + 0.1)))
    )
    F = F * Flt
    return F


def fAdkinsTipLoss(B, r, R, phi):
    """
    Tip loss of Adkins and Liebeck (1994): f = B/2 * (1 - xi) / sin(phi_t),
    with the flow angle at the tip tan(phi_t) = xi * tan(phi) and xi = r/R.
    :param B: number of blades [int]
    :param r: section radius [m]
    :param R: tip radius [m]
    :param phi: angle of relative wind [rad]
    :return: returns tip loss factor [float]
    """
    xi = r / R
    phi_t = atan2(xi * sin(phi), cos(phi))
    F = 2 / pi * acos(exp(-B / 2 * abs((1 - xi) / sin(phi_t))))
    return F


# noinspection PyUnusedLocal,PyUnusedLocal
def fInductionCoefficients0(F, phi, sigma, C_norm, C_tang, prop_coeff,  *args, **kwargs):
    """
    Calculates induction coefficients using no corrections.

    NAME: Original
    SOURCE: http://orbit.dtu.dk/files/86307371/A_Detailed_Study_of_the_Rotational.pdf

    :param F: loss factors
    :param phi: relative wind [rad]
    :param sigma: solidity
    :param C_norm: normal coefficient
    :param C_tang: tangential coefficient
    :return: axial induction factor, tangential induction factor
    """

    a = (sigma * C_norm) / (4 * F * sin(phi) ** 2 + sigma * C_norm * prop_coeff)
    aprime = (sigma * C_tang) / (4 * F * sin(phi) * cos(phi) - sigma * C_tang * prop_coeff)
    return a, aprime


# noinspection PyUnusedLocal,PyUnusedLocal
def fInductionCoefficients1(F, phi, sigma, C_norm, C_tang, *args, **kwargs):
    """
    Calculates induction coefficients using method using Spera correction.

    NAME: Spera
    SOURCE: https://cmm2017.sciencesconf.org/129068/document
    AUTHOR: Spera 1994

    :param F: loss factors
    :param phi: relative wind [rad]
    :param sigma: solidity
    :param C_tang: tangential coefficient [float]
    :param C_norm: normal coefficient [float]
    :return: axial induction factor, tangential induction factor
    """
    a = (sigma * C_norm) / (4 * F * sin(phi) ** 2 + sigma * C_norm)
    # Spera's correction
    if a >= 0.2:
        ac = 0.2
        K = (4 * F * sin(phi) ** 2) / (sigma * C_norm)
        to_sqrt = abs((K * (1 - 2 * ac) + 2) ** 2 + 4 * (K * ac ** 2 - 1))
        a = 1 + 0.5 * K * (1 - 2 * ac) - 0.5 * sqrt(to_sqrt)

    aprime = (sigma * C_tang) / (4 * F * sin(phi) * cos(phi) - sigma * C_tang)
    return a, aprime


# noinspection PyUnusedLocal,PyUnusedLocal,PyUnusedLocal
def fInductionCoefficients2(a_last, F, phi, sigma, C_norm, C_tang, Cl, Ct_r, *args, **kwargs):
    """
    Calculates induction coefficients using method used in Aerodyn software (Buhl method).

    NAME: AERODYN (BUHL)
    SOURCE: AeroDyn manual - theory.

    This method is equal to Advanced brake state model method.
    """

    if Ct_r > 0.96 * F:
        # Modified Glauert correction
        a = (
                    18 * F - 20 - 3 * sqrt(Ct_r * (50 - 36 * F) +
                                           12 * F * (3 * F - 4))
            ) / (36 * F - 50)
    else:
        a = (1 + 4 * F * sin(phi) ** 2 / (sigma * C_norm)) ** -1
    aprime = (-1 + 4 * F * sin(phi) * cos(phi) / (sigma * C_tang)) ** -1
    return a, aprime


# noinspection PyUnusedLocal,PyUnusedLocal,PyUnusedLocal
def fInductionCoefficients3(a_last, F, lambda_r, phi, sigma, C_norm, C_tang, Ct_r, *args, **kwargs):
    """
    Calculates induction coefficients using Buhl correction (QBlade implementation).

    NAME: QBLADE (Buhl)
    SOURCE: QBlade/src/XBEM/BData.cpp
    AUTHOR: Buhl

    The tangential induction follows from the tangential momentum balance (with drag), not from the
    drag-free orthogonality relation a'(1+a')lambda_r^2 = a(1-a) QBlade uses; the high-thrust correction
    only replaces the axial momentum balance. With this, the method equals method 2.
    """

    if Ct_r <= 0.96 * F:
        a = 1 / (4 * F * sin(phi) ** 2 / (sigma * C_norm) + 1)
    else:
        a = (
                    18 * F - 20 - 3 * abs(Ct_r * (50 - 36 * F) +
                                          12 * F * (3 * F - 4)) ** 0.5
            ) / (36 * F - 50)

    aprime = (sigma * C_tang) / (4 * F * sin(phi) * cos(phi) - sigma * C_tang)

    return a, aprime


def fInductionCoefficients4(a_last, F, phi, Cl, C_tang, C_norm, sigma, Ct_r, *args, **kwargs):
    """
    NAME: Wilson and walker

    Method from Pratumnopharat,2010

    Wilson and Walker method
    """
    a = a_last
    ac = 0.2

    if F == 0:
        F = 1e-6
    if Ct_r <= 0.64 * F:
        a = (1 - sqrt(1 - Ct_r / F)) / 2
    else:
        a = (Ct_r - 4 * F * ac ** 2) / (4 * F * (1 - 2 * ac))

    aprime = (4 * F * sin(phi) * cos(phi) / (sigma * C_tang) - 1) ** -1

    return a, aprime


def fInductionCoefficients5(a_last, F, phi, Cl, C_norm, C_tang, sigma, lambda_r, Ct_r, *args, **kwargs):
    """
    NAME: Modified ABS model
    Method from Pratumnopharat,2010

    Modified advanced brake state model

    The tangential induction follows from the tangential momentum balance (with drag); the ABS
    correction only replaces the axial momentum balance.

    :param lambda_r: local speed ratio
    :param C_norm: normal coefficient
    :param a_last: axial induction factor
    :param F: loss factors
    :param Cl: lift coefficient
    :param phi: relative wind [rad]
    :param sigma: solidity
    :return: axial induction factor, tangential induction factor
    """
    a = a_last
    if Ct_r < 0.96 * F:
        a = (1 - sqrt(1 - Ct_r / F)) / 2
    else:
        a = 0.1432 + sqrt(-0.55106 + 0.6427 * Ct_r / F)
    aprime = (sigma * C_tang) / (4 * F * sin(phi) * cos(phi) - sigma * C_tang)

    return a, aprime


def cascadeEffectsCorrection(alpha, v, omega, r, R, c, B, a, aprime, max_thickness):
    """
    Calculates cascade effects and corresponding change in angle of attack.
    Method from PROPX: Definitions, Derivations and Data Flow, C. Harman, 1994
    :param max_thickness: maximum airfoil thickness [m]
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

    delta_alpha_1 = (
            1 / 4
            * (
                    atan2((1 - a) * v, ((1 + 2 * aprime) * r * omega))
                    - atan2(((1 - a) * v), (r * omega))
            )
    )
    delta_alpha_2 = (
            0.109
            * (B * c * max_thickness * R * omega / v)
            / (R * c * sqrt((1 - a) ** 2 + (r * R * omega / v / R) ** 2))
    )

    out = alpha + delta_alpha_1 + delta_alpha_2

    return out


def calc_rotational_augmentation_correction(
        alpha, alpha_zero, Cl, Cd, omega, r, R, c, theta, v, Vrel_norm, method,
        printer, print_all, Cl_alpha=2 * pi):
    """
    METHODS FROM http://orbit.dtu.dk/files/86307371/A_Detailed_Study_of_the_Rotational.pdf

    Cl_3D = Cl + fl * (Cl_lin - Cl) * w, Cl_lin = Cl_alpha * (alpha - alpha_zero) is the linear attached-flow lift
    of the published models, with the airfoil's own lift slope Cl_alpha [1/rad].

    Drag that goes with the added lift (Eggers et al. 2003, as in AirfoilPrep / WISDEM):
    dCd = dCl * (sin(alpha) - 0.12 cos(alpha)) / (cos(alpha) + 0.12 sin(alpha)), which turns the added force
    almost normal to the chord: delayed separation raises the normal force, not the leading-edge suction, as the
    pressure measurements on rotating blades show (NREL Phase VI, MEXICO). The own drag factors of Du & Selig
    (fd < 0, lower drag) and Chaviaropoulos & Hansen (fd = fl, Cd - Cd0) both gave a far too large torque and edge
    moment on the Phase VI rotor in stall (validation/run_validation.py, 2026-10-18), so they are not used.

    Limits:
    - rotational augmentation is delayed separation: it recovers at most the lift lost to separation, so it acts
      only where the 2D lift lies below the linear lift (same sign), with 0 <= fl <= 1;
    - w blends the correction out in deep stall: 1 for |alpha| <= 25 deg, linearly to 0 at 45 deg
      (AirfoilPrep / WISDEM "linear_25_45"); there is no attached flow left to delay;
    - sections without a lift slope (cylinder roots, Cl_alpha < 1/rad) are not corrected.
    The Gur & Rosen propeller model keeps its own sign convention.
    """
    p = printer

    fl = 0

    if method == 0:
        # Snel et al.
        a_s = 3
        h = 2
        fl = a_s * (c / r) ** h

    if method == 1:
        # Du & Selig
        gama = omega * R / sqrt(v ** 2 + (omega * R) ** 2)
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

    if method == 2:
        # Chaviaropoulos and Hansen
        ah = 2.2
        h = 1.0
        n = 4
        fl = ah * (c / r) ** h * cos(theta) ** n

    if method == 3:
        # Lindenburg
        al = 3.1
        h = 2
        fl = al * (omega * R / Vrel_norm) ** 2 * (c / r) ** h

    if method == 4:
        # Dumitrescu and Cardos
        gd = 1.25
        fl = 1 - exp(-gd / (r / c - 1))

    if method == 5:
        # Snel et al. for propellers
        # method for propellers from http://acoustics.ae.illinois.edu/pdfs/AIAA-Paper-2015-3296.pdf
        fl = (omega * r / Vrel_norm) ** 2 * (c / r) ** 2 * 1.5

    if method == 6:
        # Gur & Rosen
        # method for propellers from Aviv (2005), Propeller Performance at Low Advance Ratio - JoA vol. 42 No. 2
        if degrees(alpha) >= 8:
            fl = tanh(10.73 * (c / r))
        elif degrees(alpha) > degrees(alpha_zero):
            fl = 0
        else:
            fl = -tanh(10.73 * (c / r))

    Cl_lin = Cl_alpha * (alpha - alpha_zero)
    w = rotational_augmentation_blending(alpha)

    if Cl_alpha < 1.0:
        fl = 0.0
    if method != 6:
        separated = (Cl_lin - Cl) * Cl_lin > 0
        if not separated:
            fl = 0.0
        fl = min(max(fl, 0.0), 1.0)

    dCl = fl * (Cl_lin - Cl) * w
    dCd = dCl * (sin(alpha) - 0.12 * cos(alpha)) / (cos(alpha) + 0.12 * sin(alpha))

    if print_all:
        p.print("             Rotational augmentation:")
        p.print("               Cl_lin", Cl_lin)
        p.print("               fl", fl)
        p.print("               w", w)
        p.print("               c/r", c / r)
    return Cl + dCl, Cd + dCd


def rotational_augmentation_blending(alpha):
    """
    Weight of the rotational augmentation: 1 for |alpha| <= 25 deg, linearly to 0 at |alpha| = 45 deg
    (AirfoilPrep / WISDEM Polar.correction3D, blending "linear_25_45").

    :param alpha: angle of attack [rad]
    :return: weight 0..1
    """
    a = abs(alpha)
    if a <= radians(25):
        return 1.0
    if a >= radians(45):
        return 0.0
    return (radians(45) - a) / radians(20)


def skewed_wake_correction_calculate(yaw_angle, a, r, R, azimuth=0.0):
    """
    Pitt and Peters skewed-wake correction (as in AeroDyn), applied to the converged induction:
    a_skew = a * (1 + 15*pi/32 * r/R * tan(chi/2) * cos(psi_dw)), chi = (0.6*a + 1) * yaw (Burton et al.),
    where psi_dw is the azimuth measured from the most downwind blade position. Averaged over one
    revolution the correction is zero; it redistributes the induction over the rotor disc.

    With the solver's azimuth convention (0 = blade up, the in-plane wind V*sin(yaw) blows radially outward
    along the blade at azimuth 90 deg), cos(psi_dw) = sin(azimuth).

    :param yaw_angle: yaw angle [rad]
    :param a: axial induction factor
    :param r: section radius [m]
    :param R: tip radius [m]
    :param azimuth: blade azimuth [rad]
    :return: corrected axial induction factor
    """
    chi = (0.6 * a + 1) * yaw_angle
    a_skew = a * (1 + 15 * pi / 32 * r / R * tan(chi / 2) * sin(azimuth))
    return a_skew

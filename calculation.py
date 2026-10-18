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

import numpy
import numpy as np
import scipy.optimize as optimize
from numpy import radians, degrees

from bending_inertia import generate_hollow_foil, calculate_bending_inertia_2, calculate_centroid
from popravki import *
from utils import Printer, get_curves_functions
from visualization import interpolate_foils

numpy.seterr(all="raise")
numpy.seterr(invalid="raise")

OUTPUT_VARIABLES_LIST = {
    "a": {"type": "array", "name": "Axial induction factor", "symbol": "a", "unit": ""},
    "a'": {"type": "array", "name": "Tangential induction factor", "symbol": "a'", "unit": ""},
    "Cl": {"type": "array", "name": "Lift coefficient", "symbol": r"$C_l$", "unit": ""},
    "Cd": {"type": "array", "name": "Drag coefficient", "symbol": r"$C_d$", "unit": ""},
    "Cn": {"type": "array", "name": "Normal coefficient", "symbol": r"$C_n$", "unit": ""},
    "Ct": {"type": "array", "name": "Tangential coefficient", "symbol": r"$C_t$", "unit": ""},
    "alpha": {"type": "array", "name": "Angle of attack", "symbol": r"$\alpha$", "unit": "°"},
    "phi": {"type": "array", "name": "Relative wind angle", "symbol": r"$\phi$", "unit": "°"},
    "F": {"type": "array", "name": "Tip loss correction factor", "symbol": "F", "unit": ""},
    "lambda_r": {"type": "array", "name": "Local tip speed ratio", "symbol": r"$\lambda_r$", "unit": ""},
    "Ct_r": {"type": "array", "name": "Local thrust coefficient (def.)", "symbol": r"$C_{T_r}$", "unit": ""},
    "Vrel_norm": {"type": "array", "name": "Relative wind speed", "symbol": "W", "unit": "m/s"},

    "dFn": {"type": "array", "name": "Incremental normal force", "symbol": r"$dF_n$", "unit": "N"},
    "dFt": {"type": "array", "name": "Incremental tangential force", "symbol": r"$dF_t$", "unit": "N"},
    "dFc": {"post_processed": True, "type": "array", "name": "Incremental centrifugal force", "symbol": r"$dF_c$",
            "unit": "N"},
    "dFn/n": {"type": "array", "name": "Incremental normal force (per unit length)", "symbol": r"$dF_n/n$",
              "unit": "N/m"},
    "dFt/n": {"type": "array", "name": "Incremental tangential force (per unit length)", "symbol": r"$dF_t/n$",
              "unit": "N/m"},

    "foils": {"post_processed": True, "type": "string_array", "name": "Airfoil name", "symbol": "airfoil_name",
              "unit": ""},
    "dT": {"type": "array", "name": "Incremental thrust", "symbol": "dT", "unit": "N"},
    "dQ": {"type": "array", "name": "Incremental torque", "symbol": "dM", "unit": "Nm"},
    "Re": {"type": "array", "name": "Reynolds number", "symbol": "Re", "unit": ""},
    "U1": {"type": "array", "name": "Far-upwind speed", "symbol": "U1", "unit": "m/s"},
    "U2": {"type": "array", "name": "Near-upwind speed", "symbol": "U2", "unit": "m/s"},
    "U3": {"type": "array", "name": "Near-downwind speed", "symbol": "U3", "unit": "m/s"},
    "U4": {"type": "array", "name": "Far-downwind speed", "symbol": "U4", "unit": "m/s"},

    "Ix": {"post_processed": True, "type": "array", "name": "Bending inertia (normal)", "symbol": r"$I_x$",
           "unit": r"$mm^4$"},
    "Iy": {"post_processed": True, "type": "array", "name": "Bending inertia (tangential)", "symbol": r"$I_y$",
           "unit": r"$mm^4$"},
    "Ixy": {"post_processed": True, "type": "array", "name": "Bending inertia (xy)", "symbol": r"$I_xy$",
            "unit": r"$mm^4$"},
    "A": {"post_processed": True, "type": "array", "name": "Airfoil area", "symbol": "A", "unit": r"$mm^2$"},
    "Ms_t": {"post_processed": True, "type": "array", "name": "Bending moment (tangential)",
             "symbol": r"$M_{bend,tang.}$", "unit": "Nm"},
    "Ms_n": {"post_processed": True, "type": "array", "name": "Bending moment (normal)", "symbol": r"$M_{bend,norm.}$",
             "unit": "Nm"},
    "stress_norm": {"post_processed": True, "type": "array", "name": "Bending stress (normal)", "symbol": r"$\sigma_n$",
                    "unit": "MPa"},
    "stress_tang": {"post_processed": True, "type": "array", "name": "Bending stress (tangential)",
                    "symbol": r"$\sigma_t$", "unit": "MPa"},
    "stress_cent": {"post_processed": True, "type": "array", "name": "Axial stress (centrifugal)",
                    "symbol": r"$\sigma_c$", "unit": "MPa"},
    "stress_von_mises": {"post_processed": True, "type": "array", "name": "Von Mises stress", "symbol": r"$\sigma_y$",
                         "unit": "MPa"},

    "r": {"post_processed": True, "type": "array", "name": "Section radius", "symbol": "r", "unit": "m"},
    "dM": {"post_processed": True, "type": "array", "name": "Section torque", "symbol": "dM", "unit": "Nm"},
    "dr": {"post_processed": True, "type": "array", "name": "Section height", "symbol": "dr", "unit": "m"},
    "c": {"post_processed": True, "type": "array", "name": "Section chord length", "symbol": "c", "unit": "m"},
    "theta": {"post_processed": True, "type": "array", "name": "Section twist angle", "symbol": r"$\Theta$",
              "unit": "°"},

    "stall": {"type": "array", "name": "Stall boolean", "symbol": "", "unit": ""},

    "R": {"type": "float", "name": "Turbine radius", "symbol": "R", "unit": "m"},
    "rpm": {"type": "float", "name": "Turbine rotational velocity", "symbol": r"$\Omega$", "unit": "RPM"},
    "v": {"type": "float", "name": "Wind speed", "symbol": "v", "unit": "m/s"},
    "cp": {"type": "float", "name": "Power coefficient", "symbol": r"$C_P$", "unit": ""},
    "ct": {"type": "float", "name": "Thrust coefficient", "symbol": r"$C_T$", "unit": ""},
    "TSR": {"type": "float", "name": "Tip speed ratio", "symbol": r"$\lambda$", "unit": ""},
    "Ft": {"type": "float", "name": "Tangential force", "symbol": r"$F_t$", "unit": "N"},
    "omega": {"type": "float", "name": "Rotational velocity", "symbol": r"$\Omega$", "unit": "rad/s"},
    "M": {"type": "float", "name": "Torque sum", "symbol": r"$M_Sum$", "unit": "Nm"},
    "Q": {"type": "float", "name": "Aerodynamic torque", "symbol": "Q", "unit": "Nm"},
    "power": {"type": "float", "name": "Power", "symbol": "P", "unit": "W"},
    "thrust": {"type": "float", "name": "Thrust", "symbol": "T", "unit": "N"},
    "Rhub": {"type": "float", "name": "Hub radius", "symbol": r"$R_hub$", "unit": "m"},
    "B": {"type": "float", "name": "Number of blades", "symbol": "B", "unit": ""},

    "pitch": {"type": "float", "name": "Pitch", "symbol": "p", "unit": "°"},
    "blade_stall_percentage": {"type": "float", "name": "Blade stall percentage", "symbol": r"$s_p$", "unit": ""},

    "J": {"type": "float", "name": "Advance ratio", "symbol": "J", "unit": ""},
    "eff": {"type": "float", "name": "Propeller efficiency", "symbol": r"$\eta_p$", "unit": ""},
    "cq": {"type": "float", "name": "Torque coefficient", "symbol": r"$C_q$", "unit": ""},

    "iterations": {"type": "array", "name": "Iterations", "symbol": "i", "unit": ""},
}


class Calculator:
    """
    Class for calculation of induction factors using BEM theory.
    """

    def __init__(self, input_arguments):
        p = Printer(input_arguments["return_print"])
        airfoils, airfoils_list, transition_foils, transition_array, max_thickness_array = get_curves_functions(
            input_arguments)
        self.airfoils, self.airfoils_list, self.transition_foils, self.transition_array, self.max_thickness_array = airfoils, airfoils_list, transition_foils, transition_array, max_thickness_array

    def printer(self, _locals, p):
        """

        :param _locals:
        :param p:
        :return:
        """
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

    @staticmethod
    def convert_to_array(theta, c, r):
        """
        Converts integers or floats into numpy arrays.
        :param theta: int or float
        :param c: int or float
        :param r: int or float
        :return: np.array(theta),np.array(c),np.array(r)
        """
        if isinstance(theta, numpy.ndarray) and \
                isinstance(c, numpy.ndarray) and \
                isinstance(r, numpy.ndarray):
            return theta, c, r
        else:
            if isinstance(theta, numpy.numbers.Real) and \
                    isinstance(c, numpy.numbers.Real) and \
                    isinstance(r, numpy.numbers.Real):
                return numpy.array([theta]), numpy.array([c]), numpy.array([r])
            return None

    # noinspection PyUnusedLocal,PyUnusedLocal
    def run_array(self, theta, B, c, r, foils, dr, R, Rhub, rpm, v, pitch, method, turbine_type, print_out,
                  mach_number_correction, tip_loss_mode, hub_loss_mode,
                  cascade_correction, max_iterations, convergence_limit, rho, kin_viscosity,
                  relaxation_factor, print_all, rotational_augmentation_correction,
                  rotational_augmentation_correction_method,
                  fix_reynolds, reynolds, yaw_angle, skewed_wake_correction, blade_design, blade_thickness,
                  mass_density, use_minimization_solver, invert_alpha,
                  a_initial, aprime_initial,
                  print_progress=False, return_print=None, return_results=None, num_azimuth=16,
                  speed_of_sound=343.0, *args, **kwargs):
        """
        Calculates induction factors using standard iteration methods.

        Different methods are available as different fInductionCoefficients functions.

        ANGLES REPRESENTATION SHOWN IN
        https://cmm2017.sciencesconf.org/129068/document
        alpha - angle of attack
        phi - angle of relative wind

        :param print_progress:
        :param mass_density:
        :param blade_thickness:
        :param blade_design:
        :param skewed_wake_correction:
        :param yaw_angle:
        :param reynolds: Reynolds number (when forced) [float]
        :param fix_reynolds: Force Reynolds number [bool]
        :param mach_number_correction: use only for propeller [bool]
        :param foils: list of airfoils [str]
        :param turbine_type: int, 0 = wind turbine, 1 = propeller
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
        :param num_azimuth: blade positions per revolution over which a yawed rotor is averaged
        :param speed_of_sound: speed of sound in the fluid [m/s]
        :return: dictionary with results
        """

        if return_print is None:
            return_print = []
        if return_results is None:
            return_results = []
        p = Printer(return_print)

        # create results array placeholders
        results = {}
        arrays = [k for k, v in OUTPUT_VARIABLES_LIST.items() if v["type"] == "array" or v["type"] == "string_array"]
        arrays_no_statics = [k for k, v in OUTPUT_VARIABLES_LIST.items() if
                             (v["type"] == "array" or v["type"] == "string_array") and "post_processed" not in v]
        for array in arrays:
            results[array] = numpy.array([])

        theta, c, r = np.array(theta), np.array(c), np.array(r)
        num_sections = len(theta)

        # set constants that are section-independent
        omega = rpm * 2 * pi / 60
        TSR = omega * R / v
        J = v / (rpm / 60 * R * 2)

        # BEM CALCULATION FOR EVERY SECTION
        for n in range(num_sections):

            _r = r[n]
            _c = c[n]
            _theta = theta[n]
            _dr = dr[n]

            if print_out:
                p.print("    r", _r, "(" + str(n) + ")")

            # Coning angle (PROPX: Definitions,Derivations, Data Flow, p.22)
            psi = 0.0

            transition = self.transition_array[n]
            _airfoil = foils[n]
            _airfoil_prev = self.transition_foils[n][0]
            _airfoil_next = self.transition_foils[n][1]
            transition_coefficient = self.transition_foils[n][2]
            max_thickness = self.max_thickness_array[n]

            _locals = locals()
            del _locals["self"]

            out_results = self.calculate_section_averaged(**_locals, printer=p)

            randomize_parameters = False

            if randomize_parameters:
                try:
                    # Method to try and set a,a' and relaxation factor using DE so convergence is reached.
                    # Doesn't work that well. For now, I left it in.
                    if out_results["finished"] == False:
                        def func(inp):
                            """

                            :param inp:
                            :return:
                            """
                            p.print(inp)
                            a_initial, aprime_initial, relaxation_factor = inp
                            _locals_in = dict(_locals)
                            del _locals_in["a_initial"]
                            del _locals_in["aprime_initial"]
                            del _locals_in["relaxation_factor"]
                            res = self.calculate_section(a_initial=a_initial,
                                                         aprime_initial=aprime_initial,
                                                         relaxation_factor=relaxation_factor,
                                                         **_locals_in,
                                                         printer=p)
                            if res["finished"] == False:
                                return res["criterion_value"]
                            else:
                                global optimizer_results
                                optimizer_results = inp
                                raise Exception("Finished")

                        bounds = [(-100, 100.0), (-1, 1.0), (0.001, 1.0)]
                        initial_guess = [a_initial, aprime_initial, relaxation_factor]
                        result = optimize.differential_evolution(func, bounds)

                        # result = optimize.minimize(func,
                        #    initial_guess,
                        #    bounds=bounds,
                        #    method="powell",
                        #    options={
                        #    'ftol': convergence_limit,
                        #    "xtol":convergence_limit,
                        #    'maxiter':max_iterations}
                        # )

                        a_initial, aprime_initial, relaxation_factor = list(result.x)
                except Exception as e:
                    if "Finished" in str(e):
                        a_initial, aprime_initial, relaxation_factor = list(optimizer_results)
                        p.print("\na:", a_initial, "a':", aprime_initial, "RF:", relaxation_factor, )
                        pass
                    else:
                        raise

                out_results = self.calculate_section(**_locals, printer=p)

            if print_progress:
                p.print("*", add_newline=False)

            if out_results == None:
                return None

            for a in arrays_no_statics:
                results[a] = numpy.append(results[a], out_results[a])

        self.statical_analysis(blade_design, blade_thickness, c, dr, foils, mass_density, num_sections, omega, r,
                               results, theta, pitch=pitch, turbine_type=turbine_type)

        dFt = results["dFt"]
        Ft = numpy.sum(dFt)

        dM = B * dFt * r
        M = numpy.sum(dM)

        dQ = results["dQ"]
        Q = numpy.sum(dQ)

        power_p = Q * omega
        power = M * omega

        Pmax = 0.5 * rho * v ** 3 * pi * R ** 2

        cp_w = power / Pmax
        cp_p = power_p / (rho * (rpm / 60) ** 3 * (2 * R) ** 5)

        dFn = results["dFn"]
        Fn = numpy.sum(dFn)

        dT = results["dT"]
        T = numpy.sum(dT)

        ct_w = T / (0.5 * rho * v ** 2 * pi * R ** 2)
        ct_p = T / (rho * (2 * R) ** 4 * (rpm / 60) ** 2)

        cq_p = Q / (rho * (2 * R) ** 5 * (rpm / 60) ** 2)
        cq_w = Q / (0.5 * rho * v ** 2 * pi * R ** 3)

        eff = J / 2 / pi * ct_p / cq_p

        blade_stall_percentage = np.sum(results["stall"]) / len(results["stall"])

        # floats
        results["R"] = R
        results["rpm"] = rpm
        results["v"] = v
        if turbine_type == 1:  # if propeller
            results["cp"] = cp_p
            results["ct"] = ct_p
            results["cq"] = cq_p
            results["power"] = power_p
            if results["ct"] < 0.0:
                p.print("    ct < 0, excluding...")
                return None
            if eff > 1.0:
                p.print("    eff > 1, excluding...")
                return None
        else:
            results["cp"] = cp_w
            results["ct"] = ct_w
            results["cq"] = cq_w
            results["power"] = power

        results["TSR"] = TSR
        results["Ft"] = Ft
        results["omega"] = omega
        results["M"] = M
        results["Q"] = Q
        results["thrust"] = T
        results["Rhub"] = Rhub
        results["B"] = B
        results["J"] = J
        results["eff"] = eff
        results["pitch"] = pitch

        # arrays
        results["r"] = r
        results["dM"] = dM
        results["dr"] = dr
        results["c"] = c
        results["theta"] = theta

        results["blade_stall_percentage"] = blade_stall_percentage
        return results

    def statical_analysis(self, blade_design, blade_thickness, c, dr, foils, mass_density, num_sections, omega, r,
                          results, theta, pitch=0.0, turbine_type=0):
        """
        Beam-theory stresses of one blade (Euler-Bernoulli, no shear, no torsion).

        The section is oriented in the rotor frame by its blade angle (twist + pitch): X points in the direction
        of rotation, Y downwind for a wind turbine and in the thrust direction (upstream) for a propeller. The
        distributed aerodynamic loads q_X, q_Y outboard of a section give the bending moments
        m_X(r) = int_r^R q_X (rho - r) d rho, m_Y(r) = int_r^R q_Y (rho - r) d rho, and the normal stress of
        unsymmetric bending about the centroid is sigma = a*X + b*Y with
        a*Iyy + b*Ixy = -m_X and a*Ixy + b*Ixx = -m_Y (Ixx = int Y^2 dA, Iyy = int X^2 dA, Ixy = int XY dA).
        The centrifugal force of the outboard blade adds the axial stress F_c/A. All three act along the span, so
        the equivalent (von Mises) stress is the largest |sigma| over the section outline.

        Each section carries its strip of width dr; loads and mass of a strip act at the section radius, and the
        outboard half of the section's own strip is included.

        :param pitch: blade pitch [deg]
        :param turbine_type: 0 = wind turbine, 1 = propeller
        """
        dFn = numpy.asarray(results["dFn"], dtype=float)
        dFt = numpy.asarray(results["dFt"], dtype=float)
        # forces in the rotor frame: a propeller's tangential force opposes the rotation
        q_X = dFt if turbine_type == 0 else -dFt
        q_Y = dFn

        sections = []
        for i in range(num_sections):
            sections.append(self.get_section_airfoil_data(i, foils[i], c[i], theta[i] + pitch, blade_design,
                                                          blade_thickness, turbine_type))

        masses = numpy.array([s["A"] * dr[i] * mass_density for i, s in enumerate(sections)])

        for i in range(num_sections):
            s = sections[i]
            Ixx, Iyy, Ixy, A = s["Ix"], s["Iy"], s["Ixy"], s["A"]

            # bending moments from the loads outboard of r[i], including the outer half of the own strip
            m_X = numpy.sum(q_X[i + 1:] * (r[i + 1:] - r[i])) + q_X[i] * dr[i] / 8
            m_Y = numpy.sum(q_Y[i + 1:] * (r[i + 1:] - r[i])) + q_Y[i] * dr[i] / 8

            # centrifugal force of the blade outboard of r[i]
            F_centrifugal = numpy.sum(masses[i + 1:] * omega ** 2 * r[i + 1:]) + \
                            masses[i] / 2 * omega ** 2 * (r[i] + dr[i] / 4)
            stress_cent = F_centrifugal / A / 1e6  # MPa

            D = Ixx * Iyy - Ixy ** 2
            X, Y = s["X"], s["Y"]

            def bending_stress(mX, mY):
                a_coef = (-mX * Ixx + mY * Ixy) / D
                b_coef = (-mY * Iyy + mX * Ixy) / D
                return (a_coef * X + b_coef * Y) / 1e6  # MPa

            stress_norm = numpy.max(numpy.abs(bending_stress(0.0, m_Y)))
            stress_tang = numpy.max(numpy.abs(bending_stress(m_X, 0.0)))
            stress_von_mises = numpy.max(numpy.abs(bending_stress(m_X, m_Y) + stress_cent))

            results["Ms_n"] = numpy.append(results["Ms_n"], m_Y)
            results["Ms_t"] = numpy.append(results["Ms_t"], m_X if turbine_type == 0 else -m_X)
            results["Ix"] = numpy.append(results["Ix"], Ixx * 1e12)  # to mm4
            results["Iy"] = numpy.append(results["Iy"], Iyy * 1e12)  # to mm4
            results["Ixy"] = numpy.append(results["Ixy"], Ixy * 1e12)  # to mm4
            results["A"] = numpy.append(results["A"], A * 1e6)  # to mm2
            results["stress_norm"] = numpy.append(results["stress_norm"], stress_norm)
            results["stress_tang"] = numpy.append(results["stress_tang"], stress_tang)
            results["dFc"] = numpy.append(results["dFc"], F_centrifugal)
            results["stress_cent"] = numpy.append(results["stress_cent"], stress_cent)
            results["stress_von_mises"] = numpy.append(results["stress_von_mises"], stress_von_mises)

    def get_section_airfoil_data(self, i, _airfoil, _c, _angle, blade_design, blade_thickness, turbine_type):
        """
        Cross-section of section i in the rotor frame (cached: it depends only on the geometry).
        A transition section is the interpolated airfoil shape between the previous and the next airfoil.

        :param _angle: blade angle (twist + pitch) [deg]
        :return: dict with Ix, Iy, Ixy [m^4], A [m^2] and the outline X, Y [m] about the centroid
        """
        key = (i, _airfoil, float(_c), float(_angle), blade_design, float(blade_thickness), turbine_type)
        if not hasattr(self, "_section_cache"):
            self._section_cache = {}
        if key in self._section_cache:
            return self._section_cache[key]

        if _airfoil != "transition":
            foil_x, foil_y = self.airfoils[_airfoil]["x"], self.airfoils[_airfoil]["y"]
        else:
            _airfoil_prev, _airfoil_next, transition_coefficient = self.transition_foils[i]
            # transition_coefficient is the relative position between the previous (0) and the next (1) airfoil
            foil_x, foil_y = interpolate_foils(self.airfoils[_airfoil_prev]["x"], self.airfoils[_airfoil_prev]["y"],
                                               self.airfoils[_airfoil_next]["x"], self.airfoils[_airfoil_next]["y"],
                                               transition_coefficient)

        data = self.get_crossection_data(_c, _angle, foil_x, foil_y, blade_design, blade_thickness, turbine_type)
        self._section_cache[key] = data
        return data

    def calculate_section(self, v, omega, _r, _c, _theta, _dr, B, R, _airfoil, max_thickness, Rhub,
                          turbine_type=0, pitch=0.0, psi=0.0, fix_reynolds=False, reynolds=1e6, tip_loss_mode=0,
                          hub_loss_mode=0, cascade_correction=False,
                          rotational_augmentation_correction=False,
                          rotational_augmentation_correction_method=0, mach_number_correction=False, method=5,
                          kin_viscosity=1.4207E-5, rho=1.225, convergence_limit=0.001, max_iterations=100,
                          relaxation_factor=0.3,
                          printer=None, print_all=False, print_out=False, yaw_angle=0.0, tilt_angle=0.0,
                          skewed_wake_correction=False,
                          lambda_r_array=None, invert_alpha=False,
                          transition=False, _airfoil_prev=None, _airfoil_next=None, transition_coefficient=1.0,
                          num_sections=0, use_minimization_solver=False,
                          a_initial=0.3, aprime_initial=0.01, azimuth=0.0, speed_of_sound=343.0,
                          *args, **kwargs):
        """
        Function that calculates each section of the blade using the optimization function way seen in
        https://www.tandfonline.com/doi/pdf/10.1080/20464177.2011.11020241
        https://github.com/kegiljarhus/pyBEMT/blob/master/pybemt/rotor.py

        :param azimuth: blade azimuth [deg], 0 = blade pointing up; only matters in yawed or tilted inflow
        :param speed_of_sound: speed of sound in the fluid [m/s]
        """

        if lambda_r_array is None:
            lambda_r_array = []
        p = printer

        # local speed ratio
        lambda_r = omega * _r / v

        # tip speed ratio of the rotor
        TSR = omega * R / v

        # solidity
        sigma = _c * B / (2 * pi * _r)

        # initial guess
        a = a_initial
        aprime = aprime_initial

        # iterations counter
        i = 0

        # tip mach number
        Mach_number = omega * R / speed_of_sound
        if Mach_number >= 1.0:
            p.print("\nTip mach 1.0 exceeded")

        # convert pitch to radians
        _pitch = radians(pitch)

        # convert to radians
        yaw_angle = radians(yaw_angle)  # [Radians]
        psi = radians(psi)  # Coning angle [Radians]
        tilt_angle = radians(tilt_angle)  # [Radians]
        lambda_r_array = np.array(lambda_r_array)
        _theta = radians(_theta)
        _azimuth = radians(azimuth)

        # loads vanish where the loss factor is zero (tip with tip loss, hub with hub loss): the BEM equations
        # are singular there, so the section is not iterated
        zero_load_section = (tip_loss_mode in (1, 2, 3) and _r >= R) or (hub_loss_mode in (1, 2) and _r <= Rhub)

        # Wind components in the blade coordinate system (CCBlade, Ning 2014,
        # https://pdfs.semanticscholar.org/5e7d/9c6408b7dd8841692d950d08bce90c676dc1.pdf):
        # yaw, tilt, precone (psi) and blade azimuth; Vy includes the rotation
        sa, ca = sin(_azimuth), cos(_azimuth)
        Vx = v * ((cos(yaw_angle) * sin(tilt_angle) * ca + sin(yaw_angle) * sa) * sin(psi) +
                  cos(yaw_angle) * cos(tilt_angle) * cos(psi))
        Vy = omega * _r * cos(psi) + v * (cos(yaw_angle) * sin(tilt_angle) * sa - sin(yaw_angle) * ca)

        def func(inp, last_iteration=False, evaluate_only=False, phi_given=None, re_given=None):
            """

            :param inp: (a, a')
            :param last_iteration: minimization solver: return the results instead of the residual
            :param evaluate_only: compute the section state for the given (a, a') without updating them
            :param phi_given: inflow angle [rad] to use instead of the one from (a, a') (inflow-angle solver)
            :param re_given: Reynolds number to use with phi_given
            :return:
            """
            if print_all:
                p.print("             i", i)
            ############ START ITERATION ############
            # update counter
            # i = i + 1
            a, aprime = inp
            a_last,aprime_last = None, None
            # p.print(inp)

            # for pretty-printing only
            prepend = ""

            # wind components
            if turbine_type == 1:  # propeller
                Un = Vx * (1 + a)
                Ut = Vy * (1 - aprime)
            else:
                Un = Vx * (1 - a)
                Ut = Vy * (1 + aprime)

            # relative wind
            phi = atan2(Un, Ut)

            Vrel_norm = sqrt(Un ** 2 + Ut ** 2)

            if phi_given is not None:
                # inflow-angle solver: the relative speed follows from the Reynolds number of the outer loop
                phi = phi_given
                if re_given is not None:
                    Vrel_norm = re_given * kin_viscosity / _c

            if fix_reynolds:
                Re = reynolds
            elif re_given is not None:
                Re = re_given
            else:
                Re = Vrel_norm * _c / kin_viscosity

            # F: loss factor of the momentum balance; F1: Shen et al. (2005) factor on the blade-element
            # force coefficients (Cn, Ct), used together with the Prandtl factor in momentum
            F = 1
            F1 = 1

            if tip_loss_mode == 1:
                # Prandtl tip loss
                F = F * fTipLoss(B, _r, R, phi)
            elif tip_loss_mode == 2:
                # Shen tip loss
                F = F * fTipLoss(B, _r, R, phi)
                F1 = F1 * newTipLoss(B, _r, R, phi, TSR)
            elif tip_loss_mode == 3:
                # Adkins tip loss
                F = F * fAdkinsTipLoss(B, _r, R, phi)

            if hub_loss_mode == 1:
                # Prandtl hub loss
                F = F * fHubLoss(B, _r, Rhub, phi)
            elif hub_loss_mode == 2:
                # Shen hub loss
                F = F * fHubLoss(B, _r, Rhub, phi)
                F1 = F1 * newHubLoss(B, _r, Rhub, phi, TSR)

            # angle of attack
            if turbine_type == 1:  # propeller
                alpha = (_theta + _pitch) - phi  # radians
            else:
                alpha = phi - (_theta + _pitch)  # radians

            # cascade correction
            if cascade_correction:
                alpha = cascadeEffectsCorrection(alpha=alpha, v=v, omega=omega, r=_r, R=R, c=_c, B=B, a=a,
                                                 aprime=aprime, max_thickness=max_thickness)

            if invert_alpha:
                alpha = -alpha

            # polars are defined on [-180, 180] degrees
            alpha_deg = (degrees(alpha) + 180.0) % 360.0 - 180.0

            if transition:
                Cl1, Cd1 = self.airfoils[_airfoil_prev]["interp_function_cl"](Re, alpha_deg), \
                           self.airfoils[_airfoil_prev][
                               "interp_function_cd"](Re, alpha_deg)
                Cl2, Cd2 = self.airfoils[_airfoil_next]["interp_function_cl"](Re, alpha_deg), \
                           self.airfoils[_airfoil_next][
                               "interp_function_cd"](Re, alpha_deg)

                if (Cl1 == False and Cd1 == False) or (Cl2 == False and Cd2 == False):
                    p.print("\na:", a, "a':", aprime)
                    p.print("No Cl or CD")
                    p.print("Re:", Re)
                    p.print("alpha:", degrees(alpha))
                    return None

                if invert_alpha:
                    Cl1, Cl2 = -Cl1, -Cl2

                # transition_coefficient is the relative position between the previous (0) and the next (1) airfoil
                Cl = Cl1 * (1 - transition_coefficient) + Cl2 * transition_coefficient
                Cd = Cd1 * (1 - transition_coefficient) + Cd2 * transition_coefficient

                # determine min and max angle of attack for attached region
                aoa_min_stall_1 = self.airfoils[_airfoil_prev]["interpolation_function_stall_min"](Re)
                aoa_max_stall_1 = self.airfoils[_airfoil_prev]["interpolation_function_stall_max"](Re)

                aoa_min_stall_2 = self.airfoils[_airfoil_next]["interpolation_function_stall_min"](Re)
                aoa_max_stall_2 = self.airfoils[_airfoil_next]["interpolation_function_stall_max"](Re)

                aoa_min_stall = aoa_min_stall_1 * (1 - transition_coefficient) + \
                                aoa_min_stall_2 * transition_coefficient
                aoa_max_stall = aoa_max_stall_1 * (1 - transition_coefficient) + \
                                aoa_max_stall_2 * transition_coefficient

                def zero_finding_function1(alpha):
                    """

                    :param alpha:
                    :return:
                    """
                    return self.airfoils[_airfoil_prev]["interp_function_cl"](Re, alpha)

                def zero_finding_function2(alpha):
                    """

                    :param alpha:
                    :return:
                    """
                    return self.airfoils[_airfoil_next]["interp_function_cl"](Re, alpha)

                # zero-lift angle and lift slope: needed by the rotational augmentation only
                alpha_zero = 0.0
                if rotational_augmentation_correction:
                    alpha_zero_1 = optimize.bisect(zero_finding_function1, -10, 10, xtol=1e-3, rtol=1e-3)
                    alpha_zero_2 = optimize.bisect(zero_finding_function2, -10, 10, xtol=1e-3, rtol=1e-3)
                    alpha_zero = alpha_zero_1 * (1 - transition_coefficient) + alpha_zero_2 * transition_coefficient
                    Cl_alpha_1 = self.attached_lift_slope(_airfoil_prev, Re, alpha_zero_1)
                    Cl_alpha_2 = self.attached_lift_slope(_airfoil_next, Re, alpha_zero_2)
                    Cl_alpha = Cl_alpha_1 * (1 - transition_coefficient) + Cl_alpha_2 * transition_coefficient

                alpha_zero = radians(alpha_zero)

                if print_all:
                    p.print("Transition detected, combining airfoils.")
                    p.print("Previous airfoil is", _airfoil_prev)
                    p.print("Next airfoil is", _airfoil_next)
                    p.print("Cl1", Cl1, "Cd1", Cd1)
                    p.print("Cl2", Cl2, "Cd2", Cd2)
                    p.print("Cl=", transition_coefficient, "*Cl1+", (1 - transition_coefficient), "*Cl2=", Cl)
                    p.print("Cd=", transition_coefficient, "*Cd1+", (1 - transition_coefficient), "*Cd2=", Cd)
            else:
                Cl, Cd = self.airfoils[_airfoil]["interp_function_cl"](Re, alpha_deg), self.airfoils[_airfoil][
                    "interp_function_cd"](Re, alpha_deg)
                if Cl == False and Cd == False:
                    p.print("\na:", a, "a':", aprime)
                    p.print("No Cl or CD")
                    p.print("Re:", Re)
                    p.print("alpha:", degrees(alpha))
                    return None

                if invert_alpha:
                    Cl = -Cl

                # determine min and max angle of attack for attached region
                aoa_max_stall = self.airfoils[_airfoil]["interpolation_function_stall_max"](Re)
                aoa_min_stall = self.airfoils[_airfoil]["interpolation_function_stall_min"](Re)

                def zero_finding_function(alpha):
                    """

                    :param alpha:
                    :return:
                    """
                    return self.airfoils[_airfoil]["interp_function_cl"](Re, alpha)

                # zero-lift angle and lift slope: needed by the rotational augmentation only
                alpha_zero = 0.0
                if rotational_augmentation_correction:
                    alpha_zero = optimize.bisect(zero_finding_function, -10, 10, xtol=1e-2, rtol=1e-3)
                    Cl_alpha = self.attached_lift_slope(_airfoil, Re, alpha_zero)

                alpha_zero = radians(alpha_zero)

            stall = 0

            # stall region determination
            if degrees(alpha) > aoa_max_stall:
                stall = 1
            # inverse stall region determination
            if degrees(alpha) < aoa_min_stall:
                stall = 1

            if print_all:
                p.print("             Cl:", Cl, "Cd:", Cd)

            if rotational_augmentation_correction:
                Cl, Cd = calc_rotational_augmentation_correction(alpha=alpha, Cl=Cl, Cd=Cd, omega=omega, r=_r, R=R,
                                                                 c=_c, theta=_theta, v=v, Vrel_norm=Vrel_norm,
                                                                 method=rotational_augmentation_correction_method,
                                                                 alpha_zero=alpha_zero, printer=printer,
                                                                 print_all=print_all, Cl_alpha=Cl_alpha)

                if print_all:
                    p.print("             Cl_cor:", Cl, "Cd_cor:", Cd)

            if mach_number_correction:
                # local Mach number of the relative flow at the section
                Mach_number_local = Vrel_norm / speed_of_sound
                if Mach_number_local < 1.0:
                    Cl, Cd = machNumberCorrection(Cl, Cd, Mach_number_local)

            # circulation gamma
            Gamma_B = 0.5 * Vrel_norm * _c * Cl

            # normal and tangential coefficients
            if turbine_type == 1:  # propeller
                C_norm = Cl * cos(phi) - Cd * sin(phi)
                C_tang = Cl * sin(phi) + Cd * cos(phi)
            else:
                C_norm = Cl * cos(phi) + Cd * sin(phi)
                C_tang = Cl * sin(phi) - Cd * cos(phi)

            # Shen et al. (2005): the blade-element forces go to zero at the tip
            C_norm = F1 * C_norm
            C_tang = F1 * C_tang

            # force calculation (per blade); C_norm and C_tang already carry the sign of drag
            # for the selected turbine type
            dFn = C_norm * 0.5 * rho * Vrel_norm ** 2 * _c * _dr  # normal force
            dFt = C_tang * 0.5 * rho * Vrel_norm ** 2 * _c * _dr  # tangential force

            if turbine_type == 1:
                Ct_r = (sigma * (1 + a) ** 2 * C_norm) / (sin(phi) ** 2)
            else:
                Ct_r = (sigma * (1 - a) ** 2 * C_norm) / (sin(phi) ** 2)

            dFn_norm = dFn / _dr
            dFt_norm = dFt / _dr

            if invert_alpha:
                alpha = -alpha

            if turbine_type == 0:
                # wind turbine
                prop_coeff = 1.0
                if invert_alpha:
                    prop_coeff = -1.0
            else:
                # propeller
                prop_coeff = -1.0
                if invert_alpha:
                    prop_coeff = 1.0

            input_arguments = {"F": F, "lambda_r": lambda_r, "phi": phi, "sigma": sigma, "C_norm": C_norm,
                               "C_tang": C_tang, "Cl": Cl, "Cd": Cd, "B": B, "c": _c, "r": _r, "R": R, "psi": psi,
                               "aprime_last": aprime, "omega": omega, "v": v, "a_last": a, "Ct_r": Ct_r,
                               "method": method, "alpha": alpha, "alpha_deg": degrees(alpha),
                               "prop_coeff":prop_coeff}

            if print_all:
                args_to_print = ["a_last", "aprime_last", "alpha_deg", "phi"]
                for argument in args_to_print:
                    p.print("            ", argument, input_arguments[argument])
                p.print("             --------")

            if not use_minimization_solver and not evaluate_only:
                # calculate new induction coefficients
                coeffs = calculate_coefficients(method, input_arguments)
                if coeffs == None:
                    return None

                # save old values
                a_last = a
                aprime_last = aprime

                # set new values
                a, aprime = coeffs

            # wake rotation correction
            k = omega * Gamma_B / (np.pi * v ** 2)
            aprime_vct = k / (4 * lambda_r ** 2)
            relevant_radiuses = lambda_r_array[np.nonzero(lambda_r_array >= lambda_r)]

            # thrust and torque - Wiley, WE 2nd, p.124; the momentum balance uses the axial wind component Vx
            # (equal to v without yaw), consistent with Un = Vx(1 -+ a) in the blade-element equations
            dT_MT = F * 4 * pi * _r * rho * Vx ** 2 * a * (1 - a) * _dr
            dT_BET = 0.5 * rho * B * _c * Vrel_norm ** 2 * \
                     C_norm * _dr
            dQ_MT = F * 4 * aprime * (1 - a) * rho * \
                    Vx * pi * _r ** 3 * omega * _dr
            dQ_BET = B * 0.5 * rho * Vrel_norm ** 2 * \
                     C_tang * _c * _dr * _r

            # thrust and torque - propeller
            dT_MT_p = F * 4 * pi * _r * rho * Vx ** 2 * (1 + a) * a * _dr
            dQ_MT_p = F * 4 * pi * _r ** 3 * rho * Vx * omega * (1 + a) * aprime * _dr
            dT_BET_p = 0.5 * rho * Vx ** 2 * _c * B * (1 + a) ** 2 / (sin(phi) ** 2) * C_norm * _dr
            dQ_BET_p = 0.5 * rho * Vx * _c * B * omega * _r ** 2 * (1 + a) * (1 - aprime) / (
                    sin(phi) * cos(phi)) * C_tang * _dr

            # blade element forces (C_norm and C_tang are defined per turbine type)
            dT = dT_BET
            dQ = dQ_BET

            # wind after (actuator disc: the speed is continuous through the rotor, U2 == U3)
            # (axial components)
            if turbine_type == 1:  # propeller
                U1 = Vx
                U2 = U1 * (1 + a)
                U3 = U2
                U4 = U1 * (1 + 2 * a)
            else:
                U1 = Vx
                U2 = U1 * (1 - a)
                U3 = U2
                U4 = U1 * (1 - 2 * a)

            if a_last == None:
                a_last = a
            if aprime_last == None:
                aprime_last = aprime

            out = {
                "a": a,
                "a'": aprime,
                "Cl": Cl,
                "Cd": Cd,
                "alpha": degrees(alpha),
                "phi": degrees(phi),
                "F": F,
                "dFt": dFt,
                "dFn": dFn,
                "dFt/n": dFt_norm,
                "dFn/n": dFn_norm,
                "_airfoil": _airfoil,
                "dT": dT,
                "dQ": dQ,
                "Re": Re,
                'U1': U1,
                'U2': U2,
                'U3': U3,
                'U4': U4,
                "lambda_r": lambda_r,
                "stall": stall,
                "Ct_r": Ct_r,
                "Vrel_norm": Vrel_norm,
                "Cn":C_norm, "Ct":C_tang,
                "iterations": i,
                "criterion_value": abs(a - a_last)
            }

            if evaluate_only:
                out["finished"] = True
                return out

            if use_minimization_solver:
                if turbine_type == 1:  # propeller
                    g = ((dT_BET_p - dT_MT_p) ** 2 + 100 * (dQ_BET_p - dQ_MT_p) ** 2)
                else:
                    g = (dT_BET - dT_MT) ** 2 + 100 * (dQ_BET - dQ_MT) ** 2

                if last_iteration:
                    return out
                else:
                    return g
            else:
                # check convergence
                if abs(a - a_last) < convergence_limit:
                    if abs(aprime - aprime_last) < convergence_limit:
                        out["finished"] = True
                        return out

                # relaxation
                a = a_last * (1 - relaxation_factor) + a * relaxation_factor
                aprime = aprime_last * (1 - relaxation_factor) + aprime * relaxation_factor
                out["a"] = a
                out["a'"] = aprime
                out["finished"] = False
                return out

        if zero_load_section:
            out = func((0.0, 0.0), evaluate_only=True)
            if out is None:
                return None
            for key in ("dFt", "dFn", "dFt/n", "dFn/n", "dT", "dQ", "Ct_r", "a", "a'"):
                out[key] = 0.0
            out["F"] = 0.0
            return out

        def propeller_inflow_angle_solution():
            """
            Propeller section solved for the inflow angle phi (Ning 2014, CCBlade) instead of by fixed-point
            iteration on (a, a'). With k = sigma Cn / (4 F sin^2 phi) and k' = sigma Ct / (4 F sin phi cos phi)
            the propeller momentum equations (method 0) give 1 + a = 1 / (1 - k) and 1 - a' = 1 / (1 + k'), and
            phi = atan2(Vx (1 + a), Vy (1 - a')) becomes the residual
                R(phi) = sin(phi) (1 - k) - Vx / Vy cos(phi) (1 + k'),
            which is continuous where a itself is singular (k = 1). A bracketed root is always found, while the
            fixed-point iteration can settle on the spurious root a = -1, a' = 1 (no relative flow, no load)
            at high loading. The Reynolds number is updated in an outer loop.
            :return: section results, or None when no root is bracketed (the caller falls back to iterating)
            """
            if Vx <= 0.0 or Vy <= 0.0:
                return None
            state = {"re": sqrt(Vx ** 2 + Vy ** 2) * _c / kin_viscosity}

            def k_factors(phi):
                o = func((a_initial, aprime_initial), evaluate_only=True, phi_given=phi, re_given=state["re"])
                if o is None or o["F"] <= 0.0:
                    return None
                k = sigma * o["Cn"] / (4 * o["F"] * sin(phi) ** 2)
                kp = sigma * o["Ct"] / (4 * o["F"] * sin(phi) * cos(phi))
                return k, kp

            def residual(phi):
                kk = k_factors(phi)
                if kk is None:
                    return float("nan")
                k, kp = kk
                return sin(phi) * (1 - k) - Vx / Vy * cos(phi) * (1 + kp)

            eps = 1e-6
            # the normal propeller working range first, then the windmill-brake and the reversed ranges
            intervals = [(eps, pi / 2 - eps), (-pi / 4, -eps), (pi / 2 + eps, pi - eps)]
            out = None
            for _ in range(4):
                root = None
                for lo, hi in intervals:
                    grid = np.linspace(lo, hi, 25)
                    values = [residual(x) for x in grid]
                    for j in range(len(grid) - 1):
                        if np.isfinite(values[j]) and np.isfinite(values[j + 1]) and values[j] * values[j + 1] <= 0:
                            root = optimize.brentq(residual, grid[j], grid[j + 1], xtol=1e-10)
                            break
                    if root is not None:
                        break
                if root is None:
                    return None
                k, kp = k_factors(root)
                if abs(1 - k) < 1e-12:
                    return None
                out = func((k / (1 - k), kp / (1 + kp)), evaluate_only=True)
                if out is None or fix_reynolds or abs(out["Re"] - state["re"]) <= 1e-3 * state["re"]:
                    break
                state["re"] = out["Re"]
            if out is not None:
                out["iterations"] = 0
                out["criterion_value"] = abs(residual(root))
            return out

        if turbine_type == 1 and method == 0 and not use_minimization_solver and not invert_alpha:
            out = propeller_inflow_angle_solution()
            if out is not None:
                return out

        if use_minimization_solver:
            bounds = [(0, 1.0), (0, 1.0)]
            initial_guess = [a_initial, aprime_initial]
            result = optimize.minimize(func, initial_guess, bounds=bounds, method="powell",
                                       options={"xtol": convergence_limit,
                                                'ftol': convergence_limit,
                                                'maxiter': max_iterations})
            out = func(result.x, True)
            if out is not None:
                out["finished"] = bool(result.success)
        else:
            i = 0
            while True:
                i = i + 1
                out = func((a, aprime))

                if out == None:
                    return None

                if out["finished"] == True:
                    break
                else:
                    a, aprime = out["a"], out["a'"]

                    # check iterations limit
                    if i >= max_iterations:
                        if print_out:
                            p.print("-*-*-*-*-*-*-*-*-*-*-*-*-*-\n", "|max iterations exceeded\n", "|------>a:", a,
                                    " aprime",
                                    aprime, )
                            prepend = "|"
                        break

        ############ END ITERATION ############

        if out is None:
            return None

        # Pitt and Peters skewed-wake correction, applied to the converged induction (as in AeroDyn); the loads
        # follow from the corrected inflow
        if skewed_wake_correction and yaw_angle != 0.0:
            a_skew = skewed_wake_correction_calculate(yaw_angle, out["a"], _r, R, _azimuth)
            out_skew = func((a_skew, out["a'"]), evaluate_only=True)
            if out_skew is None:
                return None
            for key in ("iterations", "finished", "criterion_value"):
                out_skew[key] = out[key]
            out = out_skew

        return out

    def calculate_section_averaged(self, num_azimuth=16, yaw_angle=0.0, tilt_angle=0.0, turbine_type=0, method=0,
                                   cascade_correction=False, skewed_wake_correction=False, *args, **kwargs):
        """
        Solves one blade section. In axial inflow (no yaw, no tilt) the flow is the same at every azimuth and the
        section is solved once. In yawed or tilted inflow the inflow seen by the blade changes over one revolution,
        so the section is solved at num_azimuth equally spaced blade positions and the results are averaged over
        the revolution (time average of the rotor loads).

        Takes the same arguments as calculate_section.
        """
        if turbine_type == 1:
            if method != 0:
                raise ValueError("Propellers are calculated with method 0 (Classic BEM) only: the high-thrust "
                                 "corrections of methods 1-5 describe the turbulent wake state of wind turbines.")
            if cascade_correction:
                raise ValueError("The cascade correction is formulated for wind turbines only.")
            if skewed_wake_correction:
                raise ValueError("The skewed-wake correction is formulated for wind turbines only.")

        common = dict(yaw_angle=yaw_angle, tilt_angle=tilt_angle, turbine_type=turbine_type, method=method,
                      cascade_correction=cascade_correction, skewed_wake_correction=skewed_wake_correction)

        if (yaw_angle == 0.0 and tilt_angle == 0.0) or int(num_azimuth) <= 1:
            return self.calculate_section(*args, azimuth=0.0, **common, **kwargs)

        outs = []
        for k in range(int(num_azimuth)):
            out = self.calculate_section(*args, azimuth=360.0 * k / int(num_azimuth), **common, **kwargs)
            if out is None:
                return None
            outs.append(out)

        averaged = {}
        for key, value in outs[0].items():
            if key == "finished":
                averaged[key] = all(o[key] for o in outs)
            elif key == "criterion_value":
                averaged[key] = max(o[key] for o in outs)
            elif isinstance(value, (bool, str)) or value is None:
                averaged[key] = value
            else:
                averaged[key] = float(np.mean([o[key] for o in outs]))
        return averaged

    def attached_lift_slope(self, airfoil, Re, alpha_zero_deg, delta_deg=2.0):
        """
        Lift slope of the attached flow at the zero-lift angle [1/rad], from the airfoil's own polar.
        """
        f_cl = self.airfoils[airfoil]["interp_function_cl"]
        return (f_cl(Re, alpha_zero_deg + delta_deg) - f_cl(Re, alpha_zero_deg - delta_deg)) / radians(2 * delta_deg)

    @staticmethod
    def to_rotor_frame(x, y, angle, turbine_type):
        """
        Airfoil coordinates (x from leading to trailing edge, y towards the suction side) to the rotor frame:
        X in the direction of rotation, Y downwind (wind turbine) or in the thrust direction (propeller).
        At blade angle 0 the chord lies in the rotor plane with the leading edge pointing in the direction of
        rotation; a positive blade angle turns the leading edge towards the incoming flow.

        :param angle: blade angle (twist + pitch) [deg]
        """
        b = numpy.radians(angle)
        x, y = numpy.asarray(x, dtype=float), numpy.asarray(y, dtype=float)
        X = -numpy.cos(b) * x - numpy.sin(b) * y
        if turbine_type == 1:
            Y = -numpy.sin(b) * x + numpy.cos(b) * y
        else:
            Y = numpy.sin(b) * x - numpy.cos(b) * y
        return X, Y

    def get_crossection_data(self, _c, _angle, foil_x, foil_y, blade_design, blade_thickness, turbine_type=0):
        """
        Section properties about the centroid, in the rotor frame.

        :param _c: chord [m]
        :param _angle: blade angle (twist + pitch) [deg]
        :param foil_x: airfoil x coordinates (chord 1)
        :param foil_y: airfoil y coordinates (chord 1)
        :param blade_design: 0 = solid, 1 = hollow with wall thickness blade_thickness [m]
        :param blade_thickness: wall thickness [m]
        :param turbine_type: 0 = wind turbine, 1 = propeller
        :return: dict with Ix (int Y^2 dA), Iy (int X^2 dA), Ixy (int XY dA) [m^4], A [m^2], outline X, Y [m]
        """
        _airfoil_x, _airfoil_y = numpy.asarray(foil_x, dtype=float) * _c, numpy.asarray(foil_y, dtype=float) * _c

        # the formulas assume counter-clockwise points (positive area, inward offset of the hollow foil)
        if calculate_centroid(_airfoil_x, _airfoil_y)[0] < 0:
            _airfoil_x, _airfoil_y = _airfoil_x[::-1], _airfoil_y[::-1]

        polygons = [(_airfoil_x, _airfoil_y)]
        if blade_design == 1:
            _airfoil_x2, _airfoil_y2 = generate_hollow_foil(_airfoil_x, _airfoil_y, blade_thickness)  # inner foil
            if len(_airfoil_x2) >= 3:  # otherwise too thin for the wall thickness: solid
                polygons.append((numpy.array(_airfoil_x2), numpy.array(_airfoil_y2)))

        # rotor frame; a mirroring transform reverses the point order, which is restored
        frame = []
        for px, py in polygons:
            X, Y = self.to_rotor_frame(px, py, _angle, turbine_type)
            if calculate_centroid(X, Y)[0] < 0:
                X, Y = X[::-1], Y[::-1]
            frame.append((X, Y))

        # move the origin to the centroid of the (hollow) section, so that the inertia is centroidal
        A_sum, Sx, Sy = 0.0, 0.0, 0.0
        for k, (X, Y) in enumerate(frame):
            A_k, Cx_k, Cy_k = calculate_centroid(X, Y)
            sign = 1.0 if k == 0 else -1.0
            A_sum += sign * A_k
            Sx += sign * A_k * Cx_k
            Sy += sign * A_k * Cy_k
        Cx, Cy = Sx / A_sum, Sy / A_sum

        Ix, Iy, Ixy, A = 0.0, 0.0, 0.0, 0.0
        for k, (X, Y) in enumerate(frame):
            sign = 1.0 if k == 0 else -1.0
            Ix_k, Iy_k, Ixy_k, A_k = calculate_bending_inertia_2(X - Cx, Y - Cy)
            Ix += sign * Ix_k
            Iy += sign * Iy_k
            Ixy += sign * Ixy_k
            A += sign * A_k

        X_outer, Y_outer = frame[0]
        return {"Ix": Ix, "Iy": Iy, "Ixy": Ixy, "A": A, "X": X_outer - Cx, "Y": Y_outer - Cy}

import numpy
import numpy as np
import scipy.optimize as optimize
from numpy import radians, degrees

from bending_inertia import generate_hollow_foil, calculate_bending_inertia_2
from popravki import *
from utils import Printer, get_curves_functions
from visualization import scale_and_normalize, rotate_array

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
    "dQ": {"type": "array", "name": "Incremental torque", "symbol": "dM", "unit": "N"},
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
    "omega": {"type": "float", "name": "Rotational velocity", "symbol": r"$\Omega$", "unit": r"$rad^{-1}$"},
    "M": {"type": "float", "name": "Torque sum", "symbol": r"$M_Sum$", "unit": "Nm"},
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
                  mass_density, geometry_scale, use_minimization_solver, invert_alpha,
                  a_initial, aprime_initial,
                  print_progress=False, return_print=None, return_results=None, *args, **kwargs):
        """
        Calculates induction factors using standard iteration methods.

        Different methods are available as different fInductionCoefficients functions.

        ANGLES REPRESENTATION SHOWN IN
        https://cmm2017.sciencesconf.org/129068/document
        alpha - angle of attack
        phi - angle of relative wind

        :param print_progress:
        :param geometry_scale:
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
        R = R * geometry_scale
        Rhub = Rhub * geometry_scale
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

            out_results = self.calculate_section(**_locals, printer=p)

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
                               results, theta)

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

        eff = J / 2 / pi * ct_p / cq_p

        blade_stall_percentage = np.sum(results["stall"]) / len(results["stall"])

        # floats
        results["R"] = R
        results["rpm"] = rpm
        results["v"] = v
        if turbine_type == 1:  # if propeller
            results["cp"] = cp_p
            results["ct"] = ct_p
            if results["ct"] < 0.0:
                p.print("    ct < 0, excluding...")
                return None
            if eff > 1.0:
                p.print("    eff > 1, excluding...")
                return None
        else:
            results["cp"] = cp_w
            results["ct"] = ct_w
        results["cq"] = cq_p

        results["TSR"] = TSR
        results["Ft"] = Ft
        results["omega"] = omega
        results["M"] = M
        results["power"] = power
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
                          results, theta):
        """

        :param blade_design:
        :param blade_thickness:
        :param c:
        :param dr:
        :param foils:
        :param mass_density:
        :param num_sections:
        :param omega:
        :param r:
        :param results:
        :param theta:
        """
        blade_mass = 0
        ### STATICAL ANALYSIS
        for i in range(num_sections):

            ### BENDING INTERTIA AND STRESS CALCULATION
            _c = c[i]
            _theta = theta[i]
            _dr = dr[i]

            transition = self.transition_array[i]
            _airfoil = foils[i]
            _airfoil_prev = self.transition_foils[i][0]
            _airfoil_next = self.transition_foils[i][1]
            transition_coefficient = self.transition_foils[i][2]

            A, Ix, Ixy, Iy, norm_dist, tang_dist = self.get_section_airfoil_data(_airfoil, _airfoil_next, _airfoil_prev,
                                                                                 _c, _theta, blade_design,
                                                                                 blade_thickness,
                                                                                 transition_coefficient)

            section_mass = A * _dr * mass_density
            blade_mass += section_mass

            # BENDING MOMENT CALCULATION
            Ms_n = 0
            Ms_t = 0
            for j in range(i, num_sections):
                Ms_n = Ms_n + results["dFn"][j] * \
                       (r[j] - r[i])
                Ms_t = Ms_t + results["dFt"][j] * \
                       (r[j] - r[i])

            results["Ms_t"] = numpy.append(results["Ms_t"], Ms_t)
            results["Ms_n"] = numpy.append(results["Ms_n"], Ms_n)
            results["Ix"] = numpy.append(results["Ix"], Ix * 1e12)  # to mm4
            results["Iy"] = numpy.append(results["Iy"], Iy * 1e12)  # to mm4
            results["Ixy"] = numpy.append(results["Ixy"], Ixy * 1e12)  # to mm4
            results["A"] = numpy.append(results["A"], A * 1e6)  # to mm2

            # STRESS CALCULATION
            max_tang_dist = numpy.max(numpy.abs(tang_dist))
            max_norm_dist = numpy.max(numpy.abs(norm_dist))
            stress_norm = max_norm_dist * Ms_n / Ix / 1e6  # MPa
            stress_tang = max_tang_dist * Ms_t / Iy / 1e6  # MPa

            results["stress_norm"] = numpy.append(results["stress_norm"], stress_norm)
            results["stress_tang"] = numpy.append(results["stress_tang"], stress_tang)

        # CENTRIFUGAL FORCE CALCULATION
        for i in range(num_sections):
            F_centrifugal = 0
            _A_section = results["A"][i]  # mm2 !
            for j in range(i, num_sections):
                _dr = dr[j]
                _r = r[j]
                _A = results["A"][j]  # mm2 !
                v_tan = _r * omega
                section_mass = _A * 1e-6 * _dr * mass_density
                F_centrifugal_section = section_mass * v_tan ** 2 / _r
                F_centrifugal += F_centrifugal_section
            stress_cent = F_centrifugal / _A_section  # MPa
            results["dFc"] = numpy.append(results["dFc"], F_centrifugal)
            results["stress_cent"] = numpy.append(results["stress_cent"], stress_cent)

        # VON MISES STRESS CALCULATION
        for i in range(num_sections):
            sigma_1 = results["stress_norm"][i]
            sigma_2 = results["stress_tang"][i]
            sigma_3 = results["stress_cent"][i]
            stress_von_mises = numpy.sqrt(
                0.5 * (sigma_1 - sigma_2) ** 2 + (sigma_2 - sigma_3) ** 2 + (sigma_3 - sigma_1) ** 2)
            results["stress_von_mises"] = numpy.append(results["stress_von_mises"], stress_von_mises)

    def get_section_airfoil_data(self, _airfoil, _airfoil_next, _airfoil_prev, _c, _theta, blade_design,
                                 blade_thickness, transition_coefficient):
        """

        :param _airfoil:
        :param _airfoil_next:
        :param _airfoil_prev:
        :param _c:
        :param _theta:
        :param blade_design:
        :param blade_thickness:
        :param transition_coefficient:
        :return:
        """

        if _airfoil != "transition":
            Ix, Iy, Ixy, A, tang_dist, norm_dist = self.get_crossection_data(_c, _theta, _airfoil, blade_design,
                                                                             blade_thickness)
            _centroid_x, _centroid_y = self.airfoils[_airfoil]["centroid_x"], self.airfoils[_airfoil]["centroid_y"]
        else:
            Ix1, Iy1, Ixy1, A1, tang_dist1, norm_dist1 = self.get_crossection_data(_c, _theta, _airfoil_prev,
                                                                                   blade_design, blade_thickness)
            Ix2, Iy2, Ixy2, A2, tang_dist2, norm_dist2 = self.get_crossection_data(_c, _theta, _airfoil_next,
                                                                                   blade_design, blade_thickness)
            _centroid_x1, _centroid_y1 = self.airfoils[_airfoil_prev]["centroid_x"], self.airfoils[_airfoil_prev][
                "centroid_y"]
            _centroid_x2, _centroid_y2 = self.airfoils[_airfoil_next]["centroid_x"], self.airfoils[_airfoil_next][
                "centroid_y"]
            _centroid_x = _centroid_x1 * transition_coefficient + _centroid_x2 * (1 - transition_coefficient)
            _centroid_y = _centroid_y1 * transition_coefficient + _centroid_y2 * (1 - transition_coefficient)
            tang_dist = tang_dist1 * transition_coefficient + tang_dist2 * (1 - transition_coefficient)
            norm_dist = norm_dist1 * transition_coefficient + norm_dist2 * (1 - transition_coefficient)

            Ix = Ix1 * transition_coefficient + Ix2 * (1 - transition_coefficient)
            Iy = Iy1 * transition_coefficient + Iy2 * (1 - transition_coefficient)
            Ixy = Ixy1 * transition_coefficient + Ixy2 * (1 - transition_coefficient)
            A = A1 * transition_coefficient + A2 * (1 - transition_coefficient)
        return A, Ix, Ixy, Iy, norm_dist, tang_dist

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
                          a_initial=0.3, aprime_initial=0.01,
                          *args, **kwargs):
        """
        Function that calculates each section of the blade using the optimization function way seen in
        https://www.tandfonline.com/doi/pdf/10.1080/20464177.2011.11020241
        https://github.com/kegiljarhus/pyBEMT/blob/master/pybemt/rotor.py
        """

        if lambda_r_array is None:
            lambda_r_array = []
        p = printer

        # local speed ratio
        lambda_r = omega * _r / v

        # solidity
        sigma = _c * B / (2 * pi * _r)

        # initial guess
        a = a_initial
        aprime = aprime_initial

        # iterations counter
        i = 0

        # tip mach number
        Mach_number = omega * R / 343
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

        def func(inp, last_iteration=False):
            """

            :param inp:
            :param last_iteration:
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
            # Equations for Vx and Vy from https://pdfs.semanticscholar.org/5e7d/9c6408b7dd8841692d950d08bce90c676dc1.pdf
            Vx = v * ((cos(yaw_angle) * sin(tilt_angle) + sin(yaw_angle)) * sin(psi) + cos(yaw_angle) * cos(psi) * cos(
                tilt_angle))
            Vy = omega * _r * cos(psi) + v * (cos(yaw_angle) * sin(tilt_angle) - sin(yaw_angle))

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

            if fix_reynolds:
                Re = reynolds
            else:
                Re = Vrel_norm * _c / kin_viscosity

            F = 1

            if tip_loss_mode == 1:
                # Prandtl tip loss
                F = F * fTipLoss(B, _r, R, phi)
            elif tip_loss_mode == 2:
                # New tip loss
                F = F * newTipLoss(B, _r, R, phi, lambda_r)
            elif tip_loss_mode == 3:
                # Adkins tip loss
                F = F * fAdkinsTipLoss(B, _r, R, phi)

            if hub_loss_mode == 1:
                # Prandtl hub loss
                F = F * fHubLoss(B, _r, Rhub, phi)
            elif hub_loss_mode == 2:
                # New hub loss
                F = F * newHubLoss(B, _r, R, phi, lambda_r)

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

            if transition:
                Cl1, Cd1 = self.airfoils[_airfoil_prev]["interp_function_cl"](Re, degrees(alpha)), \
                           self.airfoils[_airfoil_prev][
                               "interp_function_cd"](Re, degrees(alpha))
                Cl2, Cd2 = self.airfoils[_airfoil_next]["interp_function_cl"](Re, degrees(alpha)), \
                           self.airfoils[_airfoil_next][
                               "interp_function_cd"](Re, degrees(alpha))

                if Cl1 == False and Cd1 == False:
                    return {"finished": False, "iterations": i, "criterion_value": abs(a - a_last)}
                if Cl2 == False and Cd2 == False:
                    return {"finished": False, "iterations": i, "criterion_value": abs(a - a_last)}

                if invert_alpha:
                    Cl1, Cl2 = -Cl1, -Cl2

                Cl = Cl1 * transition_coefficient + Cl2 * (1 - transition_coefficient)
                Cd = Cd1 * transition_coefficient + Cd2 * (1 - transition_coefficient)

                # determine min and max angle of attack for attached region
                aoa_min_stall_1 = self.airfoils[_airfoil_prev]["interpolation_function_stall_min"](Re)
                aoa_max_stall_1 = self.airfoils[_airfoil_prev]["interpolation_function_stall_max"](Re)

                aoa_min_stall_2 = self.airfoils[_airfoil_next]["interpolation_function_stall_min"](Re)
                aoa_max_stall_2 = self.airfoils[_airfoil_next]["interpolation_function_stall_max"](Re)

                aoa_min_stall = aoa_min_stall_1 * transition_coefficient + aoa_min_stall_2 * (
                        1 - transition_coefficient)
                aoa_max_stall = aoa_max_stall_1 * transition_coefficient + aoa_max_stall_2 * (
                        1 - transition_coefficient)

                def zero_finding_function1(alpha):
                    """

                    :param alpha:
                    :return:
                    """
                    return self.airfoils[_airfoil_prev]["interp_function_cl"](Re, alpha)

                alpha_zero_1 = optimize.bisect(zero_finding_function1, -10, 10, xtol=1e-3, rtol=1e-3)

                def zero_finding_function2(alpha):
                    """

                    :param alpha:
                    :return:
                    """
                    return self.airfoils[_airfoil_next]["interp_function_cl"](Re, alpha)

                alpha_zero_2 = optimize.bisect(zero_finding_function2, -10, 10, xtol=1e-3, rtol=1e-3)

                alpha_zero = alpha_zero_1 * transition_coefficient + alpha_zero_2 * (1 - transition_coefficient)
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
                alpha_deg = degrees(alpha)
                if alpha_deg > 180:
                    alpha_deg = 180
                if alpha_deg < -180:
                    alpha_deg = -180

                Cl, Cd = self.airfoils[_airfoil]["interp_function_cl"](Re, degrees(alpha)), self.airfoils[_airfoil][
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

                alpha_zero = optimize.bisect(zero_finding_function, -10, 10, xtol=1e-2, rtol=1e-3)
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
                                                                 print_all=print_all)

                if print_all:
                    p.print("             Cl_cor:", Cl, "Cd_cor:", Cd)

            if mach_number_correction:
                Cl, Cd = machNumberCorrection(Cl, Cd, Mach_number)

            # circulation gamma
            Gamma_B = 0.5 * Vrel_norm * _c * Cl

            # normal and tangential coefficients
            if turbine_type == 1:  # propeller
                C_norm = Cl * cos(phi) - Cd * sin(phi)
                C_tang = Cl * sin(phi) + Cd * cos(phi)
            else:
                C_norm = Cl * cos(phi) + Cd * sin(phi)
                C_tang = Cl * sin(phi) - Cd * cos(phi)

            # force calculation
            dFL = Cl * 0.5 * rho * Vrel_norm ** 2 * _c * _dr  # lift force
            dFD = Cd * 0.5 * rho * Vrel_norm ** 2 * _c * _dr  # drag force
            dFt = dFL * sin(phi) - dFD * cos(phi)  # tangential force
            dFn = dFL * cos(phi) + dFD * sin(phi)  # normal force

            Ct_r = (sigma * (1 - a) ** 2 * C_norm) / (sin(phi) ** 2)

            dFn_norm = dFn * num_sections
            dFt_norm = dFt * num_sections

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

            if not use_minimization_solver:
                # calculate new induction coefficients
                coeffs = calculate_coefficients(method, input_arguments)
                if coeffs == None:
                    return {"finished": False, "iterations": i, "criterion_value": abs(a - a_last)}

                # save old values
                a_last = a
                aprime_last = aprime

                # set new values
                a, aprime = coeffs

            # wake rotation correction
            k = omega * Gamma_B / (np.pi * v ** 2)
            aprime_vct = k / (4 * lambda_r ** 2)
            relevant_radiuses = lambda_r_array[np.nonzero(lambda_r_array >= lambda_r)]

            if skewed_wake_correction:
                a_skewed = skewed_wake_correction_calculate(yaw_angle, a, _r, R)
                a = a_skewed

            # thrust and torque - Wiley, WE 2nd, p.124
            dT_MT = F * 4 * pi * _r * rho * v ** 2 * a * (1 - a) * _dr
            dT_BET = 0.5 * rho * B * _c * Vrel_norm ** 2 * \
                     C_norm * _dr
            dQ_MT = F * 4 * aprime * (1 - a) * rho * \
                    v * pi * _r ** 3 * omega * _dr
            dQ_BET = B * 0.5 * rho * Vrel_norm ** 2 * \
                     C_tang * _c * _dr * _r

            # thrust and torque - propeller
            dT_MT_p = 4 * pi * _r * rho * v ** 2 * (1 + a) * a * _dr
            dQ_MT_p = 4 * pi * _r ** 3 * rho * v * omega * (1 + a) * aprime * _dr
            dT_BET_p = 0.5 * rho * v ** 2 * _c * B * (1 + a) ** 2 / (sin(phi) ** 2) * C_norm * _dr
            dQ_BET_p = 0.5 * rho * v * _c * B * omega * _r ** 2 * (1 + a) * (1 - aprime) / (
                    sin(phi) * cos(phi)) * C_tang * _dr

            if turbine_type == 1:  # propeller
                dT = dT_MT_p
                dQ = dQ_MT_p
            else:
                dT = dT_BET
                dQ = dQ_BET

            # wind after
            if turbine_type == 1:  # propeller
                U1 = v
                U2 = None
                U3 = U1 * (1 + a)
                U4 = U1 * (1 + 2 * a)
            else:
                U1 = v
                U2 = U1 * (1 - a)
                U3 = None
                U4 = U1 * (1 - 2 * a)

            if a_last == None:
                a_last = a
            if aprime_last == None:
                aprime_last = aprime

            out = {"a": a, "a'": aprime, "Cl": Cl, "Cd": Cd, "alpha": degrees(alpha), "phi": degrees(phi), "F": F,
                   "dFt": dFt, "dFn": dFn, "_airfoil": _airfoil, "dT": dT, "dQ": dQ, "Re": Re,
                   'U1': U1, 'U2': U2, 'U3': U3, 'U4': U4,
                   "lambda_r": lambda_r, "dFt/n": dFt_norm, "dFn/n": dFn_norm, "stall": stall, "Ct_r": Ct_r,
                   "Vrel_norm": Vrel_norm,
                   "Cn":C_norm,"Ct":C_tang,
                   "iterations": i, "criterion_value": abs(a - a_last)}

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
                out["finished"] = False
                return out

        if use_minimization_solver:
            bounds = [(0, 1.0), (0, 1.0)]
            initial_guess = [a_initial, aprime_initial]
            result = optimize.minimize(func, initial_guess, bounds=bounds, method="powell",
                                       options={"xtol": convergence_limit,
                                                'ftol': convergence_limit,
                                                'maxiter': max_iterations})
        else:
            i = 0
            while True:
                i = i + 1
                out = func((a, aprime))

                if out == None:
                    return

                if out["finished"] == True:
                    return out
                else:
                    a, aprime = out["a"], out["a'"]

                    # check iterations limit
                    if i >= max_iterations:
                        if print_out:
                            p.print("-*-*-*-*-*-*-*-*-*-*-*-*-*-\n", "|max iterations exceeded\n", "|------>a:", a,
                                    " aprime",
                                    aprime, )
                            prepend = "|"
                        return out

        ############ END ITERATION ############
        out = func(result.x, True)
        return out

    def get_crossection_data(self, _c, _theta, _airfoil, blade_design, blade_thickness):
        """
        :param _c:
        :param _theta: in degrees
        :param _airfoil:
        :param blade_design:
        :param blade_thickness:
        :return:
        """
        _airfoil_x, _airfoil_y = self.airfoils[_airfoil]["x"], self.airfoils[_airfoil]["y"]
        _centroid_x, _centroid_y = self.airfoils[_airfoil]["centroid_x"], self.airfoils[_airfoil]["centroid_y"]
        _centroid = (_centroid_x, _centroid_y)
        _airfoil_x, _airfoil_y = scale_and_normalize(_airfoil_x, _airfoil_y, _c, _centroid)  # outer foil

        if blade_design == 1:
            _airfoil_x2, _airfoil_y2 = generate_hollow_foil(_airfoil_x, _airfoil_y, blade_thickness)  # inner foil

            _airfoil_x, _airfoil_y = rotate_array(_airfoil_x, _airfoil_y, (0, 0), _theta)  # outer foil
            _airfoil_x2, _airfoil_y2 = rotate_array(_airfoil_x2, _airfoil_y2, (0, 0), _theta)  # inner foil

            Ix1, Iy1, Ixy1, A1 = calculate_bending_inertia_2(_airfoil_x, _airfoil_y)  # outer foil
            Ix2, Iy2, Ixy2, A2 = calculate_bending_inertia_2(_airfoil_x2, _airfoil_y2)  # inner foil

            Ix = Ix1 - Ix2
            Iy = Iy1 - Iy2
            Ixy = Ixy1 - Ixy2
            A = A1 - A2

        else:
            _airfoil_x, _airfoil_y = rotate_array(_airfoil_x, _airfoil_y, (0, 0), _theta)  # outer foil
            Ix, Iy, Ixy, A = calculate_bending_inertia_2(_airfoil_x, _airfoil_y)

        min_x = numpy.min(_airfoil_x)
        max_x = numpy.max(_airfoil_x)
        min_y = numpy.min(_airfoil_y)
        max_y = numpy.max(_airfoil_y)

        return Ix, Iy, Ixy, A, numpy.abs(numpy.array((min_x, max_x))), numpy.abs(numpy.array((min_y, max_y)))

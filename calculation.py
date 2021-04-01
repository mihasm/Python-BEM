import numpy
import numpy as np
from numpy import radians, degrees

from bending_inertia import generate_hollow_foil, calculate_bending_inertia_2
from popravki import *
from utils import Printer, generate_dat, sort_data, get_transition_foils, interp
from visualize import scale_and_normalize, rotate_array
from scipy import interpolate

numpy.seterr(all="raise")
numpy.seterr(invalid="raise")

OUTPUT_VARIABLES_LIST = {
    "a": {"type": "array", "name": "Axial induction factor", "symbol": "a", "unit": ""},
    "a'": {"type": "array", "name": "Tangential induction factor", "symbol": "a'", "unit": ""},
    "cL": {"type": "array", "name": "Lift coefficient", "symbol": r"$C_l$", "unit": ""},
    "cD": {"type": "array", "name": "Drag coefficient", "symbol": r"$C_d$", "unit": ""},
    "alpha": {"type": "array", "name": "Angle of attack", "symbol": r"$\alpha$", "unit": "°"},
    "phi": {"type": "array", "name": "Relative wind angle", "symbol": r"$\phi$", "unit": "°"},
    "F": {"type": "array", "name": "Tip loss correction factor", "symbol": "F", "unit": ""},
    "lambda_r": {"type": "array", "name": "Local tip speed ratio", "symbol": r"$\lambda_r$", "unit": ""},
    "Ct": {"type": "array", "name": "Tangential coefficient", "symbol": r"$C_t$", "unit": ""},

    "dFn": {"type": "array", "name": "Incremental normal force", "symbol": r"$dF_n$", "unit": "N"},
    "dFt": {"type": "array", "name": "Incremental tangential force", "symbol": r"$dF_t$", "unit": "N"},
    "dFc": {"type": "array", "name": "Incremental centrifugal force", "symbol": r"$dF_c$", "unit": "N"},
    "dFn/n": {"type": "array", "name": "Incremental normal force (per unit length)", "symbol": r"$dF_n/n$",
              "unit": "N/m"},
    "dFt/n": {"type": "array", "name": "Incremental tangential force (per unit length)", "symbol": r"$dF_t/n$",
              "unit": "N/m"},

    "foils": {"type": "string_array", "name": "Airfoil name", "symbol": "airfoil_name", "unit": ""},
    "dT": {"type": "array", "name": "Incremental thrust", "symbol": "dT", "unit": "N"},
    "dQ": {"type": "array", "name": "Incremental torque", "symbol": "dM", "unit": "N"},
    "Re": {"type": "array", "name": "Reynolds number", "symbol": "Re", "unit": ""},
    "U1": {"type": "array", "name": "Far-upwind speed", "symbol": "U1", "unit": "m/s"},
    "U2": {"type": "array", "name": "Near-upwind speed", "symbol": "U2", "unit": "m/s"},
    "U3": {"type": "array", "name": "Near-downwind speed", "symbol": "U3", "unit": "m/s"},
    "U4": {"type": "array", "name": "Far-downwind speed", "symbol": "U4", "unit": "m/s"},

    "Ix": {"type": "array", "name": "Bending inertia (normal)", "symbol": r"$I_x$", "unit": r"$mm^4$"},
    "Iy": {"type": "array", "name": "Bending inertia (tangential)", "symbol": r"$I_y$", "unit": r"$mm^4$"},
    "Ixy": {"type": "array", "name": "Bending inertia (xy)", "symbol": r"$I_xy$", "unit": r"$mm^4$"},
    "A": {"type": "array", "name": "Airfoil area", "symbol": "A", "unit": r"$mm^2$"},

    "Ms_t": {"type": "array", "name": "Bending moment (tangential)", "symbol": r"$M_{bend,tang.}$", "unit": "Nm"},
    "Ms_n": {"type": "array", "name": "Bending moment (normal)", "symbol": r"$M_{bend,norm.}$", "unit": "Nm"},
    "stress_norm": {"type": "array", "name": "Bending stress (normal)", "symbol": r"$\sigma_n$", "unit": "MPa"},
    "stress_tang": {"type": "array", "name": "Bending stress (tangential)", "symbol": r"$\sigma_t$", "unit": "MPa"},
    "stress_cent": {"type": "array", "name": "Axial stress (centrifugal)", "symbol": r"$\sigma_c$", "unit": "MPa"},
    "stress_von_mises": {"type": "array", "name": "Von Mises stress", "symbol": r"$\sigma_y$", "unit": "MPa"},

    "r": {"type": "array", "name": "Section radius", "symbol": "r", "unit": "m"},
    "dM": {"type": "array", "name": "Section torque", "symbol": "dM", "unit": "Nm"},
    "dr": {"type": "array", "name": "Section height", "symbol": "dr", "unit": "m"},
    "c": {"type": "array", "name": "Section chord length", "symbol": "c", "unit": "m"},
    "theta": {"type": "array", "name": "Section twist angle", "symbol": r"$\Theta$", "unit": "°"},

    "stall": {"type": "array", "name": "Stall boolean", "symbol":"", "unit":""},

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
    "J": {"type": "float", "name": "Advance ratio", "symbol": "J", "unit": ""},
    "eff": {"type": "float", "name": "Propeller efficiency", "symbol": r"$\eta_p$", "unit": ""},
    "pitch": {"type": "float", "name": "Pitch", "symbol": "p", "unit": "°"},
    "blade_stall_percentage" : {"type":"float", "name":"Blade stall percentage", "symbol":r"$s_p$", "unit":""}
}



class Calculator:
    """
    Class for calculation of induction factors using BEM theory.
    """

    def __init__(self, input_arguments):
        p = Printer(input_arguments["return_print"])

        self.airfoils = input_arguments["airfoils"]  # Define airfoil data
        self.airfoils_list = input_arguments["foils"]  # List of airfoils per section

        for blade_name in self.airfoils:
            self.airfoils[blade_name]["alpha_zero"] = 0.0  # TODO FIX
            generate_dat(
                blade_name, self.airfoils[blade_name]["x"], self.airfoils[blade_name]["y"])

            ncrit_selected = self.airfoils[blade_name]["ncrit_selected"]

            data = self.airfoils[blade_name]["gathered_curves"]
            data = data[np.in1d(data[:, 1], ncrit_selected)]
            data = sort_data(data)

            re = data[:, 0].flatten()
            alpha = data[:, 2].flatten()
            cl = data[:, 3].flatten()
            cd = data[:, 4].flatten()

            def interpolation_function_cl(re_in, alpha_in, re=re, alpha=alpha, cl=cl):
                return interp(re_in, alpha_in, re, alpha, cl)

            def interpolation_function_cd(re_in, alpha_in, re=re, alpha=alpha, cd=cd):
                return interp(re_in, alpha_in, re, alpha, cd)

            self.airfoils[blade_name]["interp_function_cl"] = interpolation_function_cl
            self.airfoils[blade_name]["interp_function_cd"] = interpolation_function_cd
            
            re_stall_list,aoa_min_stall_list,aoa_max_stall_list = self.airfoils[blade_name]["stall_angles"]

            if len(re_stall_list) == 1:
                #only one curve
                def interpolation_function_stall_min(re_in):
                    return aoa_min_stall_list[0]
                def interpolation_function_stall_max(re_in):
                    return aoa_max_stall_list[0]
            else:
                interpolation_function_stall_min = interpolate.interp1d(
                                                        re_stall_list,
                                                        aoa_min_stall_list,
                                                        fill_value=(aoa_min_stall_list[0],aoa_min_stall_list[-1]),
                                                        bounds_error=False)
                interpolation_function_stall_max = interpolate.interp1d(
                                                        re_stall_list,
                                                        aoa_max_stall_list,
                                                        fill_value=(aoa_max_stall_list[0],aoa_max_stall_list[-1]),
                                                        bounds_error=False)

            self.airfoils[blade_name]["interpolation_function_stall_min"] = interpolation_function_stall_min
            self.airfoils[blade_name]["interpolation_function_stall_max"] = interpolation_function_stall_max

        self.transition_foils = get_transition_foils(self.airfoils_list)
        self.transition_array = []  # True,False,False, etc.
        self.max_thickness_array = []

        for n in range(len(self.airfoils_list)):
            _c = input_arguments["c"][n]
            if self.airfoils_list[n] == 'transition':
                _airfoil_prev = self.transition_foils[n][0]
                _airfoil_next = self.transition_foils[n][1]
                transition_coefficient = self.transition_foils[n][2]

                max_thickness = self.airfoils[_airfoil_prev]["max_thickness"] * _c * transition_coefficient + \
                                self.airfoils[_airfoil_next]["max_thickness"] * _c * (1 - transition_coefficient)

                self.transition_array.append(True)
            else:
                _airfoil = self.airfoils_list[n]
                _airfoil_prev, _airfoil_next, transition_coefficient = None, None, None

                max_thickness = self.airfoils[_airfoil]["max_thickness"] * _c
                self.transition_array.append(False)

            self.max_thickness_array.append(max_thickness)

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
    def run_array(self, theta, B, c, r, foils, dr, R, Rhub, rpm, v, pitch, method, propeller_mode, print_out, tip_loss,
                  mach_number_correction,
                  hub_loss, new_tip_loss, new_hub_loss, cascade_correction, max_iterations, convergence_limit, rho,
                  relaxation_factor, print_all, rotational_augmentation_correction,
                  rotational_augmentation_correction_method,
                  fix_reynolds, reynolds, yaw_angle, skewed_wake_correction, blade_design, blade_thickness,
                  mass_density, geometry_scale,
                  print_progress=False, return_print=[], return_results=[], *args, **kwargs):
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
        :return: dictionary with results
        """

        p = Printer(return_print)

        # create results array placeholders
        results = {}
        arrays = [k for k, v in OUTPUT_VARIABLES_LIST.items() if v["type"] == "array" or v["type"] == "string_array"]
        for array in arrays:
            results[array] = numpy.array([])

        theta, c, r = self.convert_to_array(theta, c, r)
        R = R * geometry_scale
        Rhub = Rhub * geometry_scale
        num_sections = len(theta)

        # set constants that are section-independent
        omega = rpm * 2 * pi / 60
        TSR = omega * R / v
        J = v / (rpm / 60 * R * 2)
        kin_viscosity = 1.4207E-5  # Kinematic viscosity

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

            if print_progress:
                p.print("*", add_newline=False)

            if out_results == None:
                return None

            results["a"] = numpy.append(results["a"], out_results["a"])
            results["a'"] = numpy.append(results["a'"], out_results["aprime"])
            results["cL"] = numpy.append(results["cL"], out_results["Cl"])
            results["cD"] = numpy.append(results["cD"], out_results["Cd"])
            results["alpha"] = numpy.append(results["alpha"], out_results["alpha"])
            results["phi"] = numpy.append(results["phi"], out_results["phi"])
            results["F"] = numpy.append(results["F"], out_results["F"])
            results["Ct"] = numpy.append(results["Ct"], out_results["Ct"])

            results["dFn"] = numpy.append(results["dFn"], out_results["dFn"])
            results["dFt"] = numpy.append(results["dFt"], out_results["dFt"])
            results["dFn/n"] = numpy.append(results["dFn/n"], out_results["dFn/n"])
            results["dFt/n"] = numpy.append(results["dFt/n"], out_results["dFt/n"])

            results["foils"] = numpy.append(results["foils"], out_results["_airfoil"])
            results["dT"] = numpy.append(results["dT"], out_results["dT"])
            results["dQ"] = numpy.append(results["dQ"], out_results["dQ"])
            results["Re"] = numpy.append(results["Re"], out_results["Re"])
            results["U1"] = numpy.append(results["U1"], out_results["U1"])
            results["U2"] = numpy.append(results["U2"], out_results["U2"])
            results["U3"] = numpy.append(results["U3"], out_results["U3"])
            results["U4"] = numpy.append(results["U4"], out_results["U4"])
            results["lambda_r"] = numpy.append(results["lambda_r"], out_results["lambda_r"])

            results["stall"] = numpy.append(results["stall"], out_results["stall"])

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

        if propeller_mode and eff > 1.0:
            return None

        blade_stall_percentage = np.sum(results["stall"])/len(results["stall"])

        # floats
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
                _A = results["A"][j] # mm2 !
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
                          propeller_mode, pitch=0.0, psi=0.0, fix_reynolds=False, reynolds=1e6, tip_loss=False,
                          new_tip_loss=False,
                          hub_loss=False, new_hub_loss=False, cascade_correction=False,
                          rotational_augmentation_correction=False,
                          rotational_augmentation_correction_method=0, mach_number_correction=False, method=5,
                          kin_viscosity=1.4207E-5, rho=1.225, convergence_limit=0.001, max_iterations=100,
                          relaxation_factor=0.3,
                          printer=None, print_all=False, print_out=False, yaw_angle=0.0, tilt_angle=0.0,
                          skewed_wake_correction=False,
                          lambda_r_array=[],
                          transition=False, _airfoil_prev=None, _airfoil_next=None, transition_coefficient=1.0,
                          num_sections=0,
                          *args, **kwargs):
        """
        Function that calculates each section of the blade.

        Inputs:
        v:float: Wind speed [m/s].
        omega:float: Rotational velocity [s^-1].
        _r:float: Section radius [m].
        _c:float: Section chord length [m].
        _theta:float: Section angle [deg].
        _dr:float: Section height [m].
        B:int: Number of blades.
        R:float: Wind turbine radius [m].
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
        a = 1 / 3
        aprime = 0.01

        # iterations counter
        i = 0

        # tip mach number
        M = omega * _r / 343

        # convert pitch to radians
        _pitch = radians(pitch)

        # convert to radians
        yaw_angle = radians(yaw_angle)  # [Radians]
        psi = radians(psi)  # Coning angle [Radians]
        tilt_angle = radians(tilt_angle)  # [Radians]
        lambda_r_array = np.array(lambda_r_array)
        _theta = radians(_theta)

        ############ START ITERATION ############
        while True:
            # update counter
            i = i + 1

            # for pretty-printing only
            prepend = ""
            # Equations for Vx and Vy from https://pdfs.semanticscholar.org/5e7d/9c6408b7dd8841692d950d08bce90c676dc1.pdf
            Vx = v * ((cos(yaw_angle) * sin(tilt_angle) + sin(yaw_angle)) * sin(psi) + cos(yaw_angle) * cos(psi) * cos(
                tilt_angle))
            Vy = omega * _r * cos(psi) + v * (cos(yaw_angle) * sin(tilt_angle) - sin(yaw_angle))

            # wind components
            if propeller_mode:
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
                Re_next = Vrel_norm * _c / kin_viscosity
                Re = int(Re_next)

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
                alpha = (_theta + _pitch) - phi  # radians
            else:
                alpha = phi - (_theta + _pitch)  # radians

            # cascade correction
            if cascade_correction:
                alpha = cascadeEffectsCorrection(alpha=alpha, v=v, omega=omega, r=_r, R=R, c=_c, B=B, a=a,
                                                 aprime=aprime, max_thickness=max_thickness)

            if transition:
                Cl1, Cd1 = self.airfoils[_airfoil_prev]["interp_function_cl"](Re, degrees(alpha)), \
                           self.airfoils[_airfoil_prev][
                               "interp_function_cd"](Re, degrees(alpha))
                Cl2, Cd2 = self.airfoils[_airfoil_next]["interp_function_cl"](Re, degrees(alpha)), \
                           self.airfoils[_airfoil_next][
                               "interp_function_cd"](Re, degrees(alpha))

                if Cl1 == False and Cd1 == False:
                    return None
                if Cl2 == False and Cd2 == False:
                    return None

                Cl = Cl1 * transition_coefficient + Cl2 * (1 - transition_coefficient)
                Cd = Cd1 * transition_coefficient + Cd2 * (1 - transition_coefficient)
                
                # determine min and max angle of attack for attached region
                aoa_min_stall_1 = self.airfoils[_airfoil_prev]["interpolation_function_stall_min"](Re)
                aoa_max_stall_1 = self.airfoils[_airfoil_prev]["interpolation_function_stall_max"](Re)

                aoa_min_stall_2 = self.airfoils[_airfoil_next]["interpolation_function_stall_min"](Re)
                aoa_max_stall_2 = self.airfoils[_airfoil_next]["interpolation_function_stall_max"](Re)

                aoa_min_stall = aoa_min_stall_1 * transition_coefficient + aoa_min_stall_2 * (1 - transition_coefficient)
                aoa_max_stall = aoa_max_stall_1 * transition_coefficient + aoa_max_stall_2 * (1 - transition_coefficient)

                if print_all:
                    p.print("Transition detected, combining airfoils.")
                    p.print("Previous airfoil is", _airfoil_prev)
                    p.print("Next airfoil is", _airfoil_next)
                    p.print("Cl1", Cl1, "Cd1", Cd1)
                    p.print("Cl2", Cl2, "Cd2", Cd2)
                    p.print("Cl=", transition_coefficient, "*Cl1+", (1 - transition_coefficient), "*Cl2=", Cl)
                    p.print("Cd=", transition_coefficient, "*Cd1+", (1 - transition_coefficient), "*Cd2=", Cd)
            else:
                Cl, Cd = self.airfoils[_airfoil]["interp_function_cl"](Re, degrees(alpha)), self.airfoils[_airfoil][
                    "interp_function_cd"](Re, degrees(alpha))
                if Cl == False and Cd == False:
                    return None

                # determine min and max angle of attack for attached region
                aoa_max_stall = self.airfoils[_airfoil]["interpolation_function_stall_max"](Re)
                aoa_min_stall = self.airfoils[_airfoil]["interpolation_function_stall_min"](Re)

            stall = 0

            # stall region determination
            if degrees(alpha) > aoa_max_stall:
                stall = 1
            # inverse stall region determination
            if degrees(alpha) < aoa_min_stall:
                stall = 1

            if print_all:
                p.print("        Cl:", Cl, "Cd:", Cd)

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

            # circulation gamma
            Gamma_B = 0.5 * Vrel_norm * _c * Cl

            # normal and tangential coefficients
            if propeller_mode:
                C_norm = Cl * cos(phi) - Cd * sin(phi)
                C_tang = Cl * sin(phi) + Cd * cos(phi)
            else:
                C_norm = Cl * cos(phi) + Cd * sin(phi)
                C_tang = Cl * sin(phi) - Cd * cos(phi)

            input_arguments = {"F": F, "lambda_r": lambda_r, "phi": phi, "sigma": sigma, "C_norm": C_norm,
                               "C_tang": C_tang, "Cl": Cl, "Cd": Cd, "B": B, "c": _c, "r": _r, "R": R, "psi": psi,
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

            # set new values
            a, aprime, Ct = coeffs

            # wake rotation correction
            k = omega * Gamma_B / (np.pi * v ** 2)
            aprime_vct = k / (4 * lambda_r ** 2)
            relevant_radiuses = lambda_r_array[np.nonzero(lambda_r_array >= lambda_r)]

            if skewed_wake_correction:
                a_skewed = skewed_wake_correction_calculate(yaw_angle, a, _r, R)
                a = a_skewed

            # force calculation
            dFL = Cl * 0.5 * rho * Vrel_norm ** 2 * _c * _dr  # lift force
            dFD = Cd * 0.5 * rho * Vrel_norm ** 2 * _c * _dr  # drag force
            dFt = dFL * sin(phi) - dFD * cos(phi)  # tangential force
            dFn = dFL * cos(phi) + dFD * sin(phi)  # normal force

            dFn_norm = dFn * num_sections
            dFt_norm = dFt * num_sections

            # thrust and torque - Wiley, WE 2nd, p.124
            dT_MT = F * 4 * pi * _r * rho * v ** 2 * a * (1 - a) * _dr
            dT_BET = 0.5 * rho * B * _c * Vrel_norm ** 2 * \
                     C_norm * _dr
            dQ_MT = F * 4 * aprime * (1 - a) * rho * \
                    v * pi * _r ** 3 * omega * _dr
            dQ_BET = B * 0.5 * rho * Vrel_norm ** 2 * \
                     C_tang * _c * _dr * _r

            # thrust-propeller
            dT_MT_p = 4 * pi * _r * rho * v ** 2 * (1 + a) * a * _dr
            dQ_MT_p = 4 * pi * _r ** 3 * rho * v * omega * (1 + a) * aprime * _dr
            dT_BET_p = 0.5 * rho * v ** 2 * _c * B * (1 + a) ** 2 / (sin(phi) ** 2) * C_norm * _dr
            dQ_BET_p = 0.5 * rho * v * _c * B * omega * _r ** 2 * (1 + a) * (1 - aprime) / (
                    sin(phi) * cos(phi)) * C_tang * _dr

            # from http://www.aerodynamics4students.com/propel.m
            dT_BET_p_2 = 0.5 * rho * Vrel_norm ** 2 * B * _c * C_norm * _dr
            dQ_BET_p_2 = 0.5 * rho * Vrel_norm ** 2 * B * _c * _r * C_tang * _dr
            dT_MT_p_2 = 4 * pi * _r * rho * v ** 2 * (1 + a)
            dQ_MT_p_2 = 4 * pi * _r ** 3 * rho * v * (1 + a) * omega

            # https://apps.dtic.mil/dtic/tr/fulltext/u2/1013408.pdf
            dT_prop = 0.5 * B * rho * Vrel_norm ** 2 * _c * _dr * C_norm
            dQ_prop = 0.5 * B * rho * Vrel_norm ** 2 * _c * _r * _dr * C_tang

            if propeller_mode:
                dT = dT_prop
                dQ = dQ_prop
            else:
                dT = dT_BET
                dQ = dQ_BET

            # wind after
            if propeller_mode:
                U1 = v
                U2 = None
                U3 = U1 * (1 + a)
                U4 = U1 * (1 + 2 * a)
            else:
                U1 = v
                U2 = U1 * (1 - a)
                U3 = None
                U4 = U1 * (1 - 2 * a)

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
            p.print(prepend, "        foil:", _airfoil)
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
            # p.print(prepend, "        dT_p %.2f dQ_p %.2f" % (dT_p, dQ_p))
            p.print(prepend, "    ----------------------------")

        out = {"a": a, "aprime": aprime, "Cl": Cl, "Cd": Cd, "alpha": degrees(alpha), "phi": degrees(phi), "F": F,
               "dFt": dFt, "Ct": Ct, "dFn": dFn, "_airfoil": _airfoil, "dT": dT, "dQ": dQ, "Re": Re,
               'U1': U1, 'U2': U2, 'U3': U3, 'U4': U4,
               "lambda_r": lambda_r, "dFt/n": dFt_norm, "dFn/n": dFn_norm, "stall": stall}
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

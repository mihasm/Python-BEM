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

import copy
import os
import re
import sys

import numpy
import numpy as np
import scipy
from PyQt5.QtGui import QPalette, QColor
from scipy import interpolate

# determine if application is a script file or frozen exe
if getattr(sys, 'frozen', False):
    application_path = os.path.dirname(sys.executable)
elif __file__:
    application_path = os.path.dirname(__file__)


def transpose(a):
    """
    Transposes matrix-like array a.

    Input:
    [[1,2],
     [3,4],
     [5,6]]

    Output:
    [[1,3,5],
     [2,4,6]]

    :param a: input array
    :return: a^T
    """
    o = []
    for i in range(len(a[0])):
        o.append([])
        for _r in a:
            o[i].append(_r[i])
    return o


def dict_to_ar(inp_dict):
    """

    :param inp_dict:
    :return:
    """
    prep = []
    i = 0

    list_items = inp_dict.items()
    for k, v in list_items:
        if not k == "pitch_change_list":
            prep.append([k])
            for j in v:
                if isinstance(j, np.ndarray):
                    j = np.array2string(j, max_line_width=10000000)
                prep[i].append(str(j))
            i += 1
    prep = transpose(prep)
    return prep


def dict_to_csv(inp_dict, delimiter=";"):
    """
    Creates string from input dict that is formatted like a csv file.

    E.g.:
        inp_dict={"a":[1,2,3],"b":4,"c":[5,6,7]}

        out:
        "a,b,c,
         1,4,5,
         2,4,6,
         3,4,7,"

    :param delimiter: delimiter string (default ;)
    :param inp_dict: input dictionary
    :return: csv-like string
    """
    prep = dict_to_ar(inp_dict)
    out = ""
    for _r in prep:
        for e in _r:
            out += str(e) + delimiter
        out += "\n"
    return out

def array_to_csv(in_ar, delimiter="\t"):
    """

    :param in_ar:
    :param delimiter:
    :return:
    """
    out = ""
    for r in range(len(in_ar)):
        for c in in_ar[r]:
            if c == None:
                c = ""
            out += c + delimiter
        out = out[0:-1]
        out += "\n"
    return out


def sort_xy(array_x, array_y):
    """
    Sorts two arrays by values in the first array.
    Useful for sorting pairs of x,y points.

    Returns sorted arrays.

    Arrays have to be the same length.

    :param array_x: x values
    :param array_y: y values
    :return: x,y (sorted)
    """
    out = []
    if len(array_x) == len(array_y):
        for i in range(len(array_x)):
            out.append((array_x[i], array_y[i]))
        out = sorted(out, key=lambda k: k[0])
        out_x, out_y = [], []
        for n in range(len(out)):
            out_x.append(out[n][0])
            out_y.append(out[n][1])
        return out_x, out_y
    else:
        raise Exception(
            "Cannot create XY pairs with arrays with different num of elements"
        )


def interpolate_geom(r, c, theta, foils, R, Rhub, num=None, linspace_interp=False, geometry_scale=1.0):
    """
    interpolates c,r,theta with num elements:
    """
    c_interpolator = interpolate.interp1d(r, c)
    theta_interpolator = interpolate.interp1d(r, theta)
    r_orig = r.copy()
    foils_orig = foils.copy()
    if linspace_interp:
        r = np.linspace(start=r[0], stop=r[-1], num=int(num) + 1)
        c = c_interpolator(r)
        theta = theta_interpolator(r)
        foils = []
        for _r in r:
            closest_index = find_nearest(r_orig, _r)
            foils.append(foils_orig[closest_index])
    else:
        foils = foils_orig

    # calculate dr
    dr = calculate_dr(r, R, Rhub)

    # scaling
    r = geometry_scale * r
    dr = geometry_scale * dr
    c = geometry_scale * c

    return r, c, theta, foils, dr


def calculate_dr(r, R, Rhub):
    """

    :param r:
    :param R:
    :param Rhub:
    :return:
    """
    # calculate dr
    dr = np.zeros(len(r))
    for i in range(len(r)):
        if i == 0:
            r_between = (r[i] + r[i + 1]) / 2
            _dr = r_between - Rhub
        elif i == len(r) - 1:
            r_between = (r[i] + r[i - 1]) / 2
            _dr = R - r_between
        else:
            r_between_up = (r[i] + r[i + 1]) / 2
            r_between_down = (r[i] + r[i - 1]) / 2
            _dr = r_between_up - r_between_down
        dr[i] = _dr
    return dr


def find_nearest(_array, value):
    """

    :param _array:
    :param value:
    :return:
    """
    _array = np.asarray(_array)
    idx = (np.abs(_array - value)).argmin()
    return idx


def to_float(inpt):
    """

    :param inpt:
    :return:
    """
    if isinstance(inpt, str):
        inpt = inpt.replace(",", ".")
    return float(inpt)


class Printer:
    """

    """
    def __init__(self, arr):
        self.out = arr

    def print(self, *args, add_newline=True):
        """

        :param args:
        :param add_newline:
        :return:
        """
        out_str = ""
        i = 0
        for a in args:
            if i > 0:
                out_str += " "
            if isinstance(a, float):
                a = "%.3f" % round(a, 3)
            out_str += str(a)
            i += 1
        print(out_str)
        if add_newline:
            out_str += "\n"
        self.out.append(out_str)
        return out_str

def filter_3d_results(results_3d):
    """
    Input: 3D results from calculation_runner
    Output: Same, with numpy arrays converted to lists,
    to ease JSON serialization.
    """
    for k, v in results_3d.items():
        if isinstance(v, list):
            for i in range(len(v)):
                if isinstance(results_3d[k][i], np.ndarray):
                    results_3d[k][i] = list(results_3d[k][i])
    return results_3d


def generate_dat(name, x, y):
    """

    :param name:
    :param x:
    :param y:
    :return:
    """
    out = ""
    out += name + "\n"
    for i in range(len(x)):
        _x = float(x[i])
        _y = float(y[i])
        if _y >= 0:
            out += "%.6f   %.6f\n" % (_x, _y)
        else:
            out += "%.6f  %.6f\n" % (_x, _y)

    f = open(os.path.join(application_path, "foils", name + ".dat"), "w")
    f.write(out)
    f.close()
    return out


WHITE = QColor(255, 255, 255)
BLACK = QColor(0, 0, 0)
RED = QColor(255, 0, 0)
PRIMARY = QColor(53, 53, 53)
SECONDARY = QColor(35, 35, 35)
TERTIARY = QColor(42, 130, 218)
LIGHT_GRAY = QColor(180, 180, 180)
DARK_GRAY = QColor(50, 50, 50)


def css_rgb(color, a=False):
    """Get a CSS `rgb` or `rgba` string from a `QtGui.QColor`."""
    return ("rgba({}, {}, {}, {})" if a else "rgb({}, {}, {})").format(*color.getRgb())


class QDarkPalette(QPalette):
    """Dark palette for a Qt application meant to be used with the Fusion theme."""

    def __init__(self, *__args):
        super().__init__(*__args)

        # Set all the colors based on the constants in globals
        self.setColor(QPalette.Window, PRIMARY)
        self.setColor(QPalette.WindowText, WHITE)
        self.setColor(QPalette.Base, SECONDARY)
        self.setColor(QPalette.AlternateBase, PRIMARY)
        self.setColor(QPalette.ToolTipBase, WHITE)
        self.setColor(QPalette.ToolTipText, WHITE)
        self.setColor(QPalette.Text, LIGHT_GRAY)
        self.setColor(QPalette.Button, PRIMARY)
        self.setColor(QPalette.ButtonText, WHITE)
        self.setColor(QPalette.BrightText, RED)
        self.setColor(QPalette.Link, TERTIARY)
        self.setColor(QPalette.Highlight, TERTIARY)
        self.setColor(QPalette.HighlightedText, BLACK)

    @staticmethod
    def set_stylesheet(app):
        """Static method to set the tooltip stylesheet to a `QtWidgets.QApplication`."""
        app.setStyleSheet("QToolTip {{"
                          "color: {white};"
                          "background-color: {tertiary};"
                          "border: 1px solid {white};"
                          "}}".format(white=css_rgb(WHITE), tertiary=css_rgb(TERTIARY)))

    def set_app(self, app):
        """Set the Fusion theme and this palette to a `QtWidgets.QApplication`."""
        app.setStyle("Fusion")
        app.setPalette(self)
        self.set_stylesheet(app)


def sort_data(data, columns=(0, 2)):
    """

    :param data:
    :param columns:
    :return:
    """
    if len(columns) == 0:
        raise Exception("Sorting must be done for more than zero columns.")
    first = False
    for i in columns:
        if not first:
            data = data[data[:, i].argsort()]  # sort by reynolds
            first = True
        else:
            data = data[data[:, i].argsort(kind="mergesort")]  # sort by alpha
    return data


def create_folder(name_path):
    """

    :param name_path:
    """
    if not os.path.exists(name_path):
        os.makedirs(name_path)


def get_centroid_coordinates(x, y):
    """
    Calculates area of a polygon.
    :param x:
    :param y:
    :return: 
    """

    A = 0
    for i in range(0, len(x) - 1):
        A = A + (x[i] * y[i + 1] - x[i + 1] * y[i])
    A = A * 1 / 2

    Cx = 0
    Cy = 0

    for i in range(0, len(x) - 1):
        Cx = Cx + (x[i] + x[i + 1]) * (x[i] * y[i + 1] - x[i + 1] * y[i])
        Cy = Cy + (y[i] + y[i + 1]) * (x[i] * y[i + 1] - x[i + 1] * y[i])

    Cx = Cx * 1 / (6 * A)
    Cy = Cy * 1 / (6 * A)
    return Cx, Cy


def generate_v_and_rpm_from_tsr(tsr_list, R, geometry_scale, v=None, rpm=None):
    """
    TSR = omega * R / v

    TSR = 2*pi*rpm/60 * R / v
    """
    out_v = []
    out_rpm = []
    if v == None:
        # rpm is fixed
        out_rpm.append(rpm)
        for tsr in tsr_list:
            _v = 2 * np.pi * rpm / 60 * (R * geometry_scale) / tsr
            out_v.append(_v)
    elif rpm == None:
        # v is fixed
        out_v.append(v)
        for tsr in tsr_list:
            _rpm = tsr * v * 60 / (R * geometry_scale) / 2 / np.pi
            out_rpm.append(_rpm)
    return out_v, out_rpm


def generate_v_and_rpm_from_J(J_list, R, geometry_scale, v=None, rpm=None, printer=None):
    """
    J = v / (rpm/60 * D) ...

    J = v / (rpm/60 * (2*R))

    J = v / (rpm/60) / (2*R)

    (rpm/60)=v/J/(2*R)

    rpm = 60*v/J/(2*R)
    """

    out_v = []
    out_rpm = []
    if v == None:
        # rpm is fixed
        out_rpm.append(rpm)
        for J in J_list:
            _v = J * (rpm / 60) * (2 * R * geometry_scale)
            out_v.append(_v)
    elif rpm == None:
        # v is fixed
        out_v.append(v)
        for J in J_list:
            _rpm = 60 * v / J / (2 * R * geometry_scale)
            out_rpm.append(_rpm)
    return out_v, out_rpm


def import_dat(file_path):
    """

    :param file_path:
    :return:
    """
    f = open(file_path, "r")
    # f.readlines()
    lines = f.readlines()
    f.close()

    x, y = [], []

    for l in lines:
        l = l.strip()
        l = re.sub(r'\s+', ' ', l).strip()
        splitted = l.split(" ")
        if len(splitted) == 2:
            _x = float(splitted[0])
            _y = float(splitted[1])
            x.append(_x)
            y.append(_y)
    return x, y


def import_nrel_dat(file_path):
    """

    :param file_path:
    :return:
    """
    f = open(file_path, "r")
    lines = f.readlines()
    f.close()

    i = 0
    found = False
    Re = 0.0
    ncrit = 0.0
    startline = 0
    while i < len(lines):
        if "Table of aerodynamics coefficients" in lines[i]:
            found = True
            startline = i + 4
        if "Re" in lines[i][14:33]:
            Re = float(lines[i][0:11]) * 1e6
        i += 1

    print("starting line:", startline)
    print("Reynolds:", Re)
    print("found:", found)
    if found:
        # alpha,cl,cd,cm=[],[],[],[]
        data = []
        i = startline
        while i < len(lines):
            if "!" in lines[i]:
                break
            l = lines[i]
            l = l.strip()
            l = re.sub(r'\s+', ' ', l).strip()
            splitted = l.split(" ")
            if len(splitted) == 4:
                _alpha = float(splitted[0])
                _cl = float(splitted[1])
                _cd = float(splitted[2])
                _cm = float(splitted[3])
                data.append([Re, ncrit, _alpha, _cl, _cd])
            i += 1
        # print(data)
        data = np.array(data)
        return data


def get_transition_foils(foils):
    """

    :param foils:
    :return:
    """
    transition_foils = []
    for j in range(len(foils)):
        if foils[j] == "transition":
            k = j
            while k > 0:
                k = k - 1
                prev_foil = foils[k]
                if prev_foil != "transition":
                    break
            l = j
            while l < len(foils):
                l = l + 1
                next_foil = foils[l]
                if next_foil != "transition":
                    break
            number_of_transition_sections = l - k
            relative_position = j - k
            coefficient_lower = relative_position / number_of_transition_sections
            transition_foils.append([prev_foil, next_foil, coefficient_lower])
        else:
            transition_foils.append([None, None, None])
    return transition_foils


def greek_letters_to_string(string):
    """

    :param string:
    :return:
    """
    dict_letters = {"\\alpha": "α",
                    "\\beta": "β",
                    "\\gamma": "γ",
                    "\\delta": "δ",
                    "\\epsilon": "ε",
                    "\\zeta": "ζ",
                    "\\eta": "η",
                    "\\theta": "Θ",
                    "\\iota": "ι",
                    "\\kappa": "κ",
                    "\\lambda": "λ",
                    "\\mu": "μ",
                    "\\nu": "ν",
                    "\\xi": "ξ",
                    "\\omicron ": "ℴ",
                    "\\pi": "π",
                    "\\rho": "ρ",
                    "\\sigma": "σ",
                    "\\tau": "τ",
                    "\\upsilon": "υ",
                    "\\phi": "ϕ",
                    "\\chi": "χ",
                    "\\psi": "ψ",
                    "\\omega": "ω",
                    "\\Alpha": "A",
                    "\\Beta": "B",
                    "\\Gamma": "Γ",
                    "\\Delta": "Δ",
                    "\\Epsilon": "E",
                    "\\Zeta": "Z",
                    "\\Eta": "H",
                    "\\Theta": "Θ",
                    "\\Iota": "I",
                    "\\Kappa": "K",
                    "\\Lambda": "Λ",
                    "\\Mu": "M",
                    "\\Nu": "N",
                    "\\Xi": "Ξ",
                    "\\Omicron": "O",
                    "\\Pi": "Π",
                    "\\Rho": "P",
                    "\\Sigma": "Σ",
                    "\\Tau": "T",
                    "\\Upsilon": "Υ",
                    "\\Phi": "Φ",
                    "\\Chi": "X",
                    "\\Psi": "Ψ",
                    "\\Omega": "Ω"}
    while True:
        found = False
        for k, v in dict_letters.items():
            if k in string:
                string = string.replace(k, v)
                found = True
        if not found:
            break
    return string


def get_curves_functions(input_arguments):
    """

    :param input_arguments:
    :return:
    """
    airfoils = input_arguments["airfoils"]  # Define airfoil data
    airfoils_list = input_arguments["foils"]  # List of airfoils per section

    for blade_name in airfoils:
        generate_dat(
            blade_name, airfoils[blade_name]["x"], airfoils[blade_name]["y"])

        ncrit_selected = airfoils[blade_name]["ncrit_selected"]

        data = airfoils[blade_name]["gathered_curves"]
        data = np.array(data)
        data = data[np.in1d(data[:, 1], ncrit_selected)]
        data = sort_data(data)

        re = data[:, 0].flatten()
        alpha = data[:, 2].flatten()
        cl = data[:, 3].flatten()
        cd = data[:, 4].flatten()

        def interpolation_function_cl(re_in, alpha_in, re=re, alpha=alpha, cl=cl):
            """

            :param re_in:
            :param alpha_in:
            :param re:
            :param alpha:
            :param cl:
            :return:
            """
            return interp(re_in, alpha_in, re, alpha, cl)

        def interpolation_function_cd(re_in, alpha_in, re=re, alpha=alpha, cd=cd):
            """

            :param re_in:
            :param alpha_in:
            :param re:
            :param alpha:
            :param cd:
            :return:
            """
            return interp(re_in, alpha_in, re, alpha, cd)

        airfoils[blade_name]["interp_function_cl"] = interpolation_function_cl
        airfoils[blade_name]["interp_function_cd"] = interpolation_function_cd
        re_stall_list, aoa_min_stall_list, aoa_max_stall_list = airfoils[blade_name]["stall_angles"]

        if len(re_stall_list) == 1:
            # only one curve
            def interpolation_function_stall_min(re_in):
                """

                :param re_in:
                :return:
                """
                return aoa_min_stall_list[0]

            def interpolation_function_stall_max(re_in):
                """

                :param re_in:
                :return:
                """
                return aoa_max_stall_list[0]
        else:
            interpolation_function_stall_min = interpolate.interp1d(
                re_stall_list,
                aoa_min_stall_list,
                fill_value=(aoa_min_stall_list[0], aoa_min_stall_list[-1]),
                bounds_error=False)
            interpolation_function_stall_max = interpolate.interp1d(
                re_stall_list,
                aoa_max_stall_list,
                fill_value=(aoa_max_stall_list[0], aoa_max_stall_list[-1]),
                bounds_error=False)

        airfoils[blade_name]["interpolation_function_stall_min"] = interpolation_function_stall_min
        airfoils[blade_name]["interpolation_function_stall_max"] = interpolation_function_stall_max

    transition_foils = get_transition_foils(airfoils_list)
    transition_array = []  # True,False,False, etc.
    max_thickness_array = []

    for n in range(len(airfoils_list)):
        _c = input_arguments["c"][n]
        if airfoils_list[n] == 'transition':
            _airfoil_prev = transition_foils[n][0]
            _airfoil_next = transition_foils[n][1]
            transition_coefficient = transition_foils[n][2]

            max_thickness = airfoils[_airfoil_prev]["max_thickness"] * _c * transition_coefficient + \
                            airfoils[_airfoil_next]["max_thickness"] * _c * (1 - transition_coefficient)

            transition_array.append(True)
        else:
            _airfoil = airfoils_list[n]
            _airfoil_prev, _airfoil_next, transition_coefficient = None, None, None

            max_thickness = airfoils[_airfoil]["max_thickness"] * _c
            transition_array.append(False)

        max_thickness_array.append(max_thickness)

    return airfoils, airfoils_list, transition_foils, transition_array, max_thickness_array


### CHORD TWIST GENERATORS ###

def generate_chord_lengths_betz(radiuses, R, Cl_max, B, TSR):
    """
    Source: http://wflportal.amcplaza.com/Research/DYNAM/Resource%20Documents/WT_Theory_2009.pdf
    """
    chords = 16 * np.pi * R / (9 * B * Cl_max) * (TSR * np.sqrt(TSR ** 2 * (radiuses / R) ** 2 + 4 / 9)) ** -1
    return chords


def generate_chord_lengths_schmitz(radiuses, R, Cl_max, B, TSR):
    """
    Source: http://wflportal.amcplaza.com/Research/DYNAM/Resource%20Documents/WT_Theory_2009.pdf
    """
    chords = 16 * np.pi * radiuses / (B * Cl_max) * np.sin(1 / 3 * np.arctan(R / (TSR * radiuses))) ** 2
    return chords


def generate_twists_betz(radiuses, R, TSR, alpha_d, ):
    """
    Source: http://wflportal.amcplaza.com/Research/DYNAM/Resource%20Documents/WT_Theory_2009.pdf
    """
    thetas = np.rad2deg(np.arctan(2 * R / (3 * radiuses * TSR))) + alpha_d
    return thetas


def generate_twists_schmitz(radiuses, R, TSR, alpha_d):
    """
    Source: http://wflportal.amcplaza.com/Research/DYNAM/Resource%20Documents/WT_Theory_2009.pdf
    """
    thetas = 2 / 3 * np.rad2deg(np.arctan(R / (radiuses * TSR))) - alpha_d
    return thetas


def generate_propeller_larabee(radiuses, R, B, RPM, drag_lift_ratio, v, T, rho, cl):
    """
    Source: Larabee, 1979
    """
    radiuses = np.array(radiuses)
    xi = radiuses / R
    omega = RPM * 2 * np.pi / 60
    x = omega * radiuses / v
    TSR = v / omega / R
    f = B / 2 * (np.sqrt(TSR ** 2 + 1) / TSR) * (1 - radiuses / R)
    F = 2 / np.pi * np.arccos(np.exp(-f))
    G = F * x ** 2 / (1 + x ** 2)
    _xi = np.insert(xi, 0, 0, axis=0)
    dxi = np.diff(_xi)
    y = G * (1 - drag_lift_ratio / x) * xi
    y2 = G * (1 - drag_lift_ratio / x) * xi / (x ** 2 + 1)
    I1 = 4 * np.trapz(y, x=xi)
    I2 = 2 * np.trapz(y2, x=xi)
    thrust_coeff = 2 * T / (rho * v ** 2 * np.pi * R ** 2)
    zeta = I1 / (2 * I2) * (1 - np.sqrt(1 - (4 * I2 * thrust_coeff) / I1 ** 2))
    vprime = zeta * v
    c_R = 4 * np.pi / B * TSR * G / (np.sqrt(1 + x ** 2)) * zeta / cl
    c = c_R * R
    theta = np.arctan(TSR / xi * (1 + 0.5 * zeta))
    theta = np.rad2deg(theta)
    return c, theta


alpha_last = None


def generate_propeller_adkins(inp):
    # radiuses, R, B, RPM, v, T, rho, cl, airfoil,
    """METHOD FROM ADKINS: https://arc.aiaa.org/doi/pdf/10.2514/3.23779"""
    Rhub = inp["Rhub"]
    R = inp["R"]
    num_gen_sections = inp["num_gen_sections"]
    r = np.linspace(Rhub, R, int(num_gen_sections))
    B = inp["B"]
    RPM = inp["design_RPM"]
    v = inp["design_velocity"]
    T = inp["design_thrust"]
    use_power_constraint = inp["design_use_power_constraint"]
    P = inp["design_power"]
    rho = inp["design_rho"]
    kin_viscosity = inp["design_kin_viscosity"]
    cl_des = inp["design_cl"]
    airfoil = inp["design_airfoil"]
    relaxation_factor = inp["design_relaxation"]
    iters = int(inp["design_iters"])
    convergence_criterion_adkins = inp["convergence_criterion_adkins"]
    minimize_losses = inp["design_minimize_losses"]

    airfoils, airfoils_list, transition_foils, transition_array, max_thickness_array = get_curves_functions(inp)

    zeta = 0.1  # initial guess

    cl_arr = np.array([cl_des] * len(r))

    for count in range(iters):
        print("count", count)
        xi = r / R
        omega = 2 * np.pi * RPM / 60
        J = v / (RPM / 60 * 2 * R)
        print("J", J)
        x = omega * r / v
        _lambda = v / (omega * R)
        print("_lambda", _lambda)
        phi_t = np.arctan(_lambda * (1 + zeta / 2))
        phi = np.arctan(np.tan(phi_t) / xi)  # direct from adkins (21)
        phi2 = np.arctan((1 + zeta / 2) / x)  # adkins (8)
        phi3 = np.arctan((1 + zeta / 2) * _lambda / xi)  # adkins (8)

        print("phi", np.rad2deg(phi))
        print("phi2", np.rad2deg(phi2))
        print("phi3", np.rad2deg(phi3))
        f = B / 2 * (1 - xi) / np.sin(phi_t)
        F = 2 / np.pi * np.arccos(np.exp(-f))
        # F = 2 / np.pi * np.arccos(np.exp(-B / 2 * np.abs((R - r) / r / np.sin(phi)))) # Typical Prandtl implementation
        print("F", F)
        G = F * x ** 2 / (1 + x ** 2)
        print("G", G)
        G = F * x * np.cos(phi) * np.sin(phi)
        print("G2", G)
        print("cl_arr", cl_arr)
        Wc = 4 * np.pi * _lambda * G * v * R * zeta / (cl_arr * B)
        print("Wc_des", Wc)
        Re = Wc / kin_viscosity
        print("Re_des", Re)
        eps_arr = [None] * len(r)
        alpha_arr = [None] * len(r)
        cl_arr = [None] * len(r)

        for i in range(len(r)):
            print("section", i)

            Re_section = Re[i]

            min_cl = None
            max_cl = None
            for a in range(-45, 45):
                _cl = airfoils[airfoil]["interp_function_cl"](Re_section, a)
                if min_cl == None:
                    min_cl = _cl
                if max_cl == None:
                    max_cl = _cl
                if _cl < min_cl:
                    min_cl = _cl
                if _cl > max_cl:
                    max_cl = _cl

            if minimize_losses:

                def eps_minimize_function(cl):
                    """

                    :param cl:
                    :return:
                    """
                    # finds minimal losses by choosing the right Cl

                    def get_dcl(alpha, cl_req=cl):
                        """

                        :param alpha:
                        :param cl_req:
                        :return:
                        """
                        # returns 0 if alpha produces specified cl, else returns > 0.
                        # used as the zero-finding function

                        cl_actual = airfoils[airfoil]["interp_function_cl"](Re_section, alpha)
                        dcl = cl_actual - cl_req
                        return dcl  # so zero finding finds value of cl

                    lower_bound = 45
                    while get_dcl(lower_bound) > 0:
                        lower_bound -= 1
                        if lower_bound < -90:
                            raise Exception("Too low bound, perhaps too low initial Cl?")

                    upper_bound = -45
                    while get_dcl(upper_bound) < 0:
                        upper_bound += 1
                        if upper_bound > 90:
                            raise Exception("Too high bound, perhaps too high initial Cl?")

                    _alpha = scipy.optimize.ridder(get_dcl, lower_bound, upper_bound, xtol=1e-2, rtol=1e-4)

                    alpha_arr[i] = _alpha
                    cl = airfoils[airfoil]["interp_function_cl"](Re_section, _alpha)
                    cd = airfoils[airfoil]["interp_function_cd"](Re_section, _alpha)

                    cl_arr[i] = cl
                    _eps = cd / cl
                    return _eps

                min_eps = scipy.optimize.minimize(eps_minimize_function, 0.1, bounds=[(0, max_cl)], method="powell",
                                                  options={'ftol': 0.001, "xtol": 0.01, 'maxiter': 5, 'maxfev': 5}).fun
                eps_arr[i] = min_eps
            else:
                def get_dcl(alpha, cl_req=cl_des):
                    """

                    :param alpha:
                    :param cl_req:
                    :return:
                    """
                    # returns 0 if alpha produces specified cl, else returns > 0.
                    # used as the zero-finding function
                    cl_actual = airfoils[airfoil]["interp_function_cl"](Re_section, alpha)
                    dcl = cl_actual - cl_req
                    return dcl  # so zero finding finds value of cl

                lower_bound = 45
                while get_dcl(lower_bound) > 0:
                    lower_bound -= 1
                    if lower_bound < -90:
                        raise Exception("Too low bound, perhaps too low initial Cl?")

                upper_bound = -45
                while get_dcl(upper_bound) < 0:
                    upper_bound += 1
                    if upper_bound > 90:
                        raise Exception("Too high bound, perhaps too high initial Cl?")

                _alpha = scipy.optimize.ridder(get_dcl, lower_bound, upper_bound, xtol=1e-2, rtol=1e-4)
                alpha_arr[i] = _alpha
                cl = airfoils[airfoil]["interp_function_cl"](Re_section, _alpha)
                cd = airfoils[airfoil]["interp_function_cd"](Re_section, _alpha)
                _eps = cd / cl
                cl_arr[i] = cl
                eps_arr[i] = _eps

        alpha_arr = np.array(alpha_arr)
        print("alpha", alpha_arr)
        eps_arr = np.array(eps_arr)
        cl_arr = np.array(cl_arr)
        print("cl", cl_arr)

        a = (zeta / 2) * np.cos(phi) ** 2 * (1 - eps_arr * np.tan(phi))
        print("a", a)

        aprime = (zeta / (2 * x)) * np.cos(phi) * np.sin(phi) * (1 + eps_arr / np.tan(phi))
        print("aprime", aprime)

        phi4 = np.arctan(v * (1 + a) / (omega * r * (1 - aprime)))
        print("phi4", np.rad2deg(phi4))

        phi5 = np.arctan((zeta / 2 - a) / (x * aprime))
        print("phi5", np.rad2deg(phi5))

        W = v * (1 + a) / np.sin(phi)
        print("W", W)
        W2 = np.sqrt((v * (1 + a)) ** 2 + (omega * r * (1 - aprime)) ** 2)
        print("W2", W2)
        c = Wc / W
        print("c", c)

        beta = np.deg2rad(alpha_arr) + phi

        I1p = 4 * xi * G * (1 - eps_arr * np.tan(phi))
        I2p = _lambda * (I1p / (2 * xi)) * (1 + eps_arr / np.tan(phi)) * np.sin(phi) * np.cos(phi)
        J1p = 4 * xi * G * (1 + eps_arr / np.tan(phi))
        J2p = (J1p / 2) * (1 - eps_arr * np.tan(phi)) * np.cos(phi) ** 2

        print("I1p", I1p)
        print("I2p", I2p)
        print("J1p", J1p)
        print("J2p", J2p)

        I1 = np.trapz(I1p, x=xi)
        I2 = np.trapz(I2p, x=xi)
        J1 = np.trapz(J1p, x=xi)
        J2 = np.trapz(J2p, x=xi)

        if use_power_constraint:
            print("Power")
            P_c = 2 * P / (rho * v ** 3 * np.pi * R ** 2)
            zeta_new = -(J1 / (2 * J2)) + ((J1 / (2 * J2)) ** 2 + P_c / J2) ** 0.5
        else:
            print("Thrust")
            T_c = 2 * T / (rho * v ** 2 * np.pi * R ** 2)
            zeta_new = (I1 / (2 * I2)) - ((I1 / (2 * I2)) ** 2 - T_c / I2) ** 0.5

        print("zeta_new", zeta_new)
        if abs(abs(zeta - zeta_new) * zeta_new) < convergence_criterion_adkins:
            print("converged")
            print(c, np.rad2deg(beta))
            return c, np.rad2deg(beta)
        else:
            print("not converged")
            zeta = (1 - relaxation_factor) * zeta + relaxation_factor * zeta_new  # relaxation
    print("Not converged")
    return None


### INTERPOLATION ###

def interp(re_in, alpha_in, re, alpha, cl):
    """
    Interpolation function uses input arrays re, alpha and cl, and input re_in and alpha_in,

    to get the interpolated value of a curve (cl or cd, or any 3D function for that matter).

    In cases that re_in is too large, or too small, the function returns the closest available valid Reynolds data.
    """

    re_list, alpha_list, cl_list = np.unique(
        re), np.unique(alpha), np.unique(cl)

    # print("For itnerpolation:",alpha,cl)

    if re_in >= re_list.max():
        # print("Reyonlds is bigger")
        indexes = np.where(re == re_list.max())
        alpha_selected = alpha[indexes]
        cl_selected = cl[indexes]
        return np.interp(alpha_in, alpha_selected, cl_selected, left=False, right=False)

    if re_in <= re_list.min():
        # print("Reynolds is smaller")
        # print("Alpha_in",alpha_in)
        indexes = np.where(re == re_list.min())
        alpha_selected = alpha[indexes]
        cl_selected = cl[indexes]
        # print("alpha_selected:",alpha_selected)
        oo = np.interp(alpha_in, alpha_selected, cl_selected, left=False, right=False)
        # print(oo)
        return oo

    re_bottom_index = np.where(re_list <= re_in)[0][-1]
    re_bottom = re_list[re_bottom_index]
    re_top = re_list[re_bottom_index + 1]

    indexes_bottom = np.where(re == re_bottom)
    indexes_top = np.where(re == re_top)

    alpha_bottom = alpha[indexes_bottom]
    cl_bottom = cl[indexes_bottom]

    alpha_top = alpha[indexes_top]
    cl_top = cl[indexes_top]

    _cl_1 = np.interp(alpha_in, alpha_bottom, cl_bottom, left=False, right=False)
    _cl_2 = np.interp(alpha_in, alpha_top, cl_top, left=False, right=False)

    if _cl_1 == False or _cl_2 == False:
        return False

    cL = (_cl_2 - _cl_1) / (re_top - re_bottom) * (re_in - re_bottom) + _cl_1
    return cL


def interp_at(x, y, v, xp, yp, algorithm='linear', extrapolate=False):
    """
    Interpolate data onto the specified points.

    Parameters:

    * x, y : 1D arrays
        Arrays with the x and y coordinates of the data points.
    * v : 1D array
        Array with the scalar value assigned to the data points.
    * xp, yp : 1D arrays
        Points where the data values will be interpolated
    * algorithm : string
        Interpolation algorithm. Either ``'cubic'``, ``'nearest'``,
        ``'linear'`` (see scipy.interpolate.griddata)
    * extrapolate : True or False
        If True, will extrapolate values outside of the convex hull of the data
        points.

    Returns:

    * v : 1D array
        1D array with the interpolated v values.

    """
    if algorithm not in ['cubic', 'linear', 'nearest']:
        raise ValueError("Invalid interpolation algorithm: " + str(algorithm))
    grid = interpolate.griddata(
        (x, y), v, (xp, yp), method=algorithm).ravel()
    if extrapolate and algorithm != 'nearest' and numpy.any(numpy.isnan(grid)):
        if xp.size > 2:
            grid = extrapolate_nans(xp, yp, grid)
        else:
            return interpolate.griddata((x, y), v, (xp, yp), method="nearest")
    return grid


def extrapolate_nans(x, y, v):
    """
    Extrapolate the NaNs or masked values in a grid INPLACE using nearest
    value.

    .. warning:: Replaces the NaN or masked values of the original array!

    Parameters:

    * x, y : 1D arrays
        Arrays with the x and y coordinates of the data points.
    * v : 1D array
        Array with the scalar value assigned to the data points.

    Returns:

    * v : 1D array
        The array with NaNs or masked values extrapolated.

    """
    if numpy.ma.is_masked(v):
        nans = v.mask
    else:
        nans = numpy.isnan(v)
    notnans = numpy.logical_not(nans)
    v[nans] = scipy.interpolate.griddata((x[notnans], y[notnans]), v[notnans], (x[nans], y[nans]),
                                         method='nearest').ravel()
    return v


def get_interpolation_function(x, y, z, num_x=10, num_y=360):
    """

    :param x:
    :param y:
    :param z:
    :param num_x:
    :param num_y:
    :return:
    """
    x, y, z = np.array(x), np.array(y), np.array(z)
    xi, yi = np.linspace(x.min(), x.max(), num_x), np.linspace(
        y.min(), y.max(), num_y)
    xi, yi = np.meshgrid(xi, yi)
    zi = interp_at(x, y, z, xi.ravel(), yi.ravel(),
                   algorithm="linear", extrapolate=True)
    fun = scipy.interpolate.interp2d(xi, yi, zi, kind='linear')
    return fun


### SOLIDWORKS MACRO BUILDER ###

def create_macro_text(list_of_files, data):
    """

    :param list_of_files:
    :param data:
    :return:
    """

    template1 = """
    Dim swApp As Object
    Dim swModel As SldWorks.ModelDoc2
    Dim part As Object
    Dim boolstatus As Boolean
    Dim longstatus As Long, longwarnings As Long
    Dim nPtData(%s) As Double
    Dim vPtData As Variant
    Dim swSketchSeg As SldWorks.SketchSegment

    Sub main()

    Set swApp = Application.SldWorks

    template = swApp.GetUserPreferenceStringValue(swDefaultTemplatepart)
    Set part = swApp.NewDocument(template, 0, 0, 0)

    Set swModel = swApp.ActiveDoc

    Dim myModelView As Object
    """ % (len(data) * 3 - 1)

    template2 = ""
    for f in list_of_files:
        template2 += 'boolstatus = part.InsertCurveFile("%s")\n' % f
    i = 0
    j = 0
    for row in data:
        template2 += "nPtData(%s) = %s\n" % (j, data[i][0][0])
        template2 += "nPtData(%s) = %s\n" % (j + 1, data[i][1][0])
        template2 += "nPtData(%s) = %s\n" % (j + 2, data[i][2][0])
        i += 1
        j += 3

    template3 = """
    vPtData = nPtData
    swModel.Insert3DSketch2 True
    Set swSketchSeg = swModel.CreateSpline(vPtData)
    Debug.Assert Not swSketchSeg Is Nothing
    swModel.Insert3DSketch2 True
    """

    template4 = ''

    for i in range(len(data)):
        template4 += "Dim myRefPlane%s As Object\n" % i

    for i in range(len(data)):
        _y = data[i][1][0]
        template4 += 'swModel.ClearSelection\n'
        template4 += 'boolstatus = Part.Extension.SelectByID2("Top Plane", "PLANE", 0, 0, 0, True, 0, Nothing, 0)\n'
        template4 += "Set myRefPlane%s = Part.FeatureManager.InsertRefPlane(8, %s, 0, 0, 0, 0)\n" % (i, _y)

    template5 = "End Sub\n"

    return template1 + "\n" + template2 + "\n" + template3 + "\n" + template4 + "\n" + template5 + "\n"

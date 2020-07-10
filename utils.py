import copy
import numpy
from scipy import interpolate
from numpy import array
import os
import sys
import re

from PyQt5.QtGui import QPalette, QColor
from PyQt5.QtWidgets import (QComboBox, QMainWindow, QPushButton, QTextEdit, QWidget, QFormLayout, QLabel, QLineEdit,
                             QGridLayout, QCheckBox, QStyleFactory, QMessageBox, QAction, QFileDialog, QSlider,
                             QTabWidget, QApplication, QScrollArea, QVBoxLayout)
from PyQt5 import QtCore, QtGui
from PyQt5.QtCore import QThread, QTextStream, pyqtSignal, QProcess, QRect, Qt
import traceback

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
    prep = []
    i = 0

    list_items = inp_dict.items()
    for k, v in list_items:
        prep.append([k])
        for j in v:
            if isinstance(j, numpy.ndarray):
                j = numpy.array2string(j, max_line_width=10000000)
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


def interpolate_geom(r, c, theta, foils, num=None, linspace_interp=False,geometry_scale=1.0):
    """
    interpolates c,r,theta with num elements:
    """
    # print("interpolating")
    # print("foils before",foils)
    c_interpolator = interpolate.interp1d(r, c)
    theta_interpolator = interpolate.interp1d(r, theta)
    r_orig = r.copy()
    foils_orig = foils.copy()
    if linspace_interp:
        r = numpy.linspace(start=r[0], stop=r[-1], num=int(num) + 1)
        c = c_interpolator(r)
        theta = theta_interpolator(r)
        foils = []
        for _r in r:
            closest_index = find_nearest(r_orig, _r)
            foils.append(foils_orig[closest_index])
    else:
        foils = foils_orig

    # calculate dr
    r_shifted = [r[0]]
    for _r in r:
        r_shifted.append(_r)
    r_shifted = array(r_shifted[:-1])
    dr = r - r_shifted
    dr[0] = dr[1]  # TODO: what
    # print("foils after",foils)

    #scaling
    r = geometry_scale * r
    dr = geometry_scale * dr
    c = geometry_scale * c
    
    return r, c, theta, foils, dr


def find_nearest(_array, value):
    _array = numpy.asarray(_array)
    idx = (numpy.abs(_array - value)).argmin()
    return idx


def to_float(inp):
    if isinstance(inp, str):
        inp = inp.replace(",", ".")
    return float(inp)


class Printer:
    def __init__(self, arr):
        self.out = arr

    def print(self, *args, add_newline=True):
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


def fltr(node, vals):
    print(node)
    if isinstance(node, dict):
        retVal = {}
        for key, value in node.items():
            if isinstance(value, numpy.ndarray):
                node[key] = value.tolist()
            if isinstance(key, vals) and isinstance(value, vals):
                retVal[key] = copy.deepcopy(node[key])
            elif isinstance(node[key], list) or isinstance(node[key], dict):
                child = fltr(node[key], vals)
                if child:
                    retVal[key] = child
        if retVal:
            return retVal
        else:
            return None
    elif isinstance(node, list):
        retVal = []
        for entry in node:
            child = fltr(entry, vals)
            if child:
                retVal.append(child)
        if retVal:
            return retVal
        else:
            return None


def generate_dat(name, x, y):
    out = ""
    out += name + "\n"
    for i in range(len(x)):
        _x = float(x[i])
        _y = float(y[i])
        if _y >= 0:
            out += "%.6f   %.6f\n" % (_x, _y)
        else:
            out += "%.6f  %.6f\n" % (_x, _y)

    f = open(os.path.join(application_path,"foils", name + ".dat"), "w")
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


def sort_data(data, columns=[0, 2]):
    if len(columns) == 0:
        raise Exception("Sorting must be done for more than zero columns.")
    first = False
    for i in columns:
        if first == False:
            data = data[data[:, i].argsort()]  # sort by reynolds
            first = True
        else:
            data = data[data[:, i].argsort(kind="mergesort")]  # sort by alpha
    return data


def normalize_angle(angle):
    # reduce the angle
    angle = angle % 360
    # force it to be the positive remainder, so that 0 <= angle < 360
    angle = (angle + 360) % 360
    # force into the minimum absolute value residue class, so that -180 < angle <= 180
    if (angle > 180):
        angle -= 360

    return angle


def create_folder(name_path):
    if not os.path.exists(name_path):
        os.makedirs(name_path)


def get_centroid_coordinates(foil_x, foil_y):
    centroid = (numpy.sum(foil_x) / len(foil_x),
                numpy.sum(foil_y) / len(foil_y))
    return centroid

class MyMessageBox(QMessageBox):
    def __init__(self):
        QtGui.QMessageBox.__init__(self)
        self.setSizeGripEnabled(True)

    def event(self, e):
        result = QtGui.QMessageBox.event(self, e)

        self.setMinimumHeight(0)
        self.setMaximumHeight(16777215)
        self.setMinimumWidth(0)
        self.setMaximumWidth(16777215)
        self.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)

        textEdit = self.findChild(QtGui.QTextEdit)
        if textEdit != None:
            textEdit.setMinimumHeight(0)
            textEdit.setMaximumHeight(16777215)
            textEdit.setMinimumWidth(0)
            textEdit.setMaximumWidth(16777215)
            textEdit.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)

        return result

def ErrorMessageBox():
    msg = MyMessageBox()
    msg.setIcon(QMessageBox.Warning)
    msg.setText("Error while getting settings")
    var = traceback.format_exc()
    msg.setDetailedText(str(var))
    msg.exec_()


def generate_v_and_rpm_from_tsr(tsr_list,R,v=None,rpm=None):
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
            _v=2*numpy.pi*rpm/60*R/tsr
            out_v.append(_v)
    elif rpm == None:
        # v is fixed
        out_v.append(v)
        for tsr in tsr_list:
            _rpm = tsr*v*60/R/2/numpy.pi
            out_rpm.append(_rpm)
    return out_v,out_rpm

def generate_v_and_rpm_from_J(J_list,R,v=None,rpm=None):
    """
    J = v / (rpm / 60 * R * 2)

    J = v / rpm * 60 / R / 2
    """
    out_v = []
    out_rpm = []
    if v == None:
        # rpm is fixed
        out_rpm.append(rpm)
        for J in J_list:
            _v=J*(rpm / 60 * R * 2)
            out_v.append(_v)
    elif rpm == None:
        # v is fixed
        out_v.append(v)
        for J in J_list:
            _rpm = v / J * 60 / R / 2
            out_rpm.append(_rpm)
    return out_v,out_rpm

def import_dat(file_path):
    f = open(file_path,"r")
    #f.readlines()
    lines = f.readlines()
    f.close()

    x,y=[],[]

    for l in lines:
        l = l.strip()
        l = re.sub(r'\s+',' ', l).strip()
        splitted = l.split(" ")
        if len(splitted) == 2:
            _x = float(splitted[0])
            _y = float(splitted[1])
            x.append(_x)
            y.append(_y)
    return x,y

#x,y = import_dat("C:\\Users\\Miha\\Google Drive\\faks\\BEM program\\foils\\DU_91_W2_250.dat")

def get_transition_foils(foils):
    transition_foils = []
    for j in range(len(foils)):
        if foils[j] == "transition":
            k=j
            while k>0:
                k=k-1
                prev_foil = foils[k]
                if prev_foil != "transition":
                    break
            l=j
            while l<len(foils):
                l=l+1
                next_foil = foils[l]
                if next_foil != "transition":
                    break
            number_of_transition_sections = l-k
            relative_position = j-k
            coefficient_lower = relative_position/number_of_transition_sections
            transition_foils.append([prev_foil,next_foil,coefficient_lower])
        else:
            transition_foils.append([None,None,None])
    return transition_foils


def greek_letters_to_string(string):
    dict_letters={"\\alpha":"α",
    "\\beta":"β",
    "\\gamma":"γ",
    "\\delta":"δ",
    "\\epsilon":"ε",
    "\\zeta":"ζ",
    "\\eta":"η",
    "\\theta":"Θ",
    "\\iota":"ι",
    "\\kappa":"κ",
    "\\lambda":"λ",
    "\\mu":"μ",
    "\\nu":"ν",
    "\\xi":"ξ",
    "\\omicron ":"ℴ",
    "\\pi":"π",
    "\\rho":"ρ",
    "\\sigma":"σ",
    "\\tau":"τ",
    "\\upsilon":"υ",
    "\\phi":"ϕ",
    "\\chi":"χ",
    "\\psi":"ψ",
    "\\omega":"ω",
    "\\Alpha":"A",
    "\\Beta":"B",
    "\\Gamma":"Γ",
    "\\Delta":"Δ",
    "\\Epsilon":"E",
    "\\Zeta":"Z",
    "\\Eta":"H",
    "\\Theta":"Θ",
    "\\Iota":"I",
    "\\Kappa":"K",
    "\\Lambda":"Λ",
    "\\Mu":"M",
    "\\Nu":"N",
    "\\Xi":"Ξ",
    "\\Omicron":"O",
    "\\Pi":"Π",
    "\\Rho":"P",
    "\\Sigma":"Σ",
    "\\Tau":"T",
    "\\Upsilon":"Υ",
    "\\Phi":"Φ",
    "\\Chi":"X",
    "\\Psi":"Ψ",
    "\\Omega":"Ω"}
    while True:
        found=False
        for k,v in dict_letters.items():
            if k in string:
                string = string.replace(k,v)
                found=True
        if found==False:
            break
    return string

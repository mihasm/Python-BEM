import copy
import csv
import os
import re
import sys
import traceback

import numpy as np
from PyQt5 import QtCore, QtGui
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtGui import QPalette, QColor
from PyQt5.QtWidgets import (QWidget, QMessageBox, QApplication, QVBoxLayout, QTableWidget, QTableWidgetItem, QMenu)
from numpy import array
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
    prep = []
    i = 0

    list_items = inp_dict.items()
    for k, v in list_items:
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
    _array = np.asarray(_array)
    idx = (np.abs(_array - value)).argmin()
    return idx


def to_float(inpt):
    if isinstance(inpt, str):
        inpt = inpt.replace(",", ".")
    return float(inpt)


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
            if isinstance(value, np.ndarray):
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
    centroid = (np.sum(foil_x) / len(foil_x),
                np.sum(foil_y) / len(foil_y))
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
            _v=2*np.pi*rpm/60*R/tsr
            out_v.append(_v)
    elif rpm == None:
        # v is fixed
        out_v.append(v)
        for tsr in tsr_list:
            _rpm = tsr*v*60/R/2/np.pi
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

def import_nrel_dat(file_path):
    f = open(file_path,"r")
    lines = f.readlines()
    f.close()

    i = 0
    found = False
    Re = 0.0
    ncrit = 0.0
    startline=0
    while (i<len(lines)):
        if "Table of aerodynamics coefficients" in lines[i]:
            found=True
            startline=i+4
        if "Re" in lines[i][14:33]:
            Re=float(lines[i][0:11])*1e6
        i+=1

    print("starting line:",startline)
    print("Reynolds:",Re)
    print("found:",found)
    if found:
        #alpha,cl,cd,cm=[],[],[],[]
        data = []
        i=startline
        while i<len(lines):
            if "!" in lines[i]:
                break
            l = lines[i]
            l = l.strip()
            l = re.sub(r'\s+',' ', l).strip()
            splitted = l.split(" ")
            if len(splitted) == 4:
                _alpha = float(splitted[0])
                _cl = float(splitted[1])
                _cd = float(splitted[2])
                _cm = float(splitted[3])
                data.append([Re, ncrit, _alpha, _cl, _cd])
            i+=1
        #print(data)
        data = np.array(data)
        return data

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

### CHORD TWIST GENERATORS ###

def generate_chord_lengths_betz(radiuses,R,Cl_max,B,TSR):
    """
    Source: http://wflportal.amcplaza.com/Research/DYNAM/Resource%20Documents/WT_Theory_2009.pdf
    """
    chords = 16*np.pi*R/(9*B*Cl_max)*(TSR*np.sqrt(TSR**2*(radiuses/R)**2+4/9))**-1
    return chords

def generate_chord_lengths_schmitz(radiuses,R,Cl_max,B,TSR):
    """
    Source: http://wflportal.amcplaza.com/Research/DYNAM/Resource%20Documents/WT_Theory_2009.pdf
    """
    chords = 16*np.pi*radiuses/(B*Cl_max)*np.sin(1/3*np.arctan(R/(TSR*radiuses)))**2
    return chords

def generate_twists_betz(radiuses,R,TSR,alpha_d):
    """
    Source: http://wflportal.amcplaza.com/Research/DYNAM/Resource%20Documents/WT_Theory_2009.pdf
    """
    thetas = np.rad2deg(np.arctan(2*R/(3*radiuses*TSR)))+alpha_d
    return thetas

def generate_twists_schmitz(radiuses,R,TSR,alpha_d):
    """
    Source: http://wflportal.amcplaza.com/Research/DYNAM/Resource%20Documents/WT_Theory_2009.pdf
    """
    thetas = 2/3*np.rad2deg(np.arctan(R/(radiuses*TSR)))-alpha_d
    return thetas

### INTERPOLATION ###

def interp(re_in, alpha_in, re, alpha, cl):
    """
    Interpolation function uses input arrays re, alpha and cl, and input re_in and alpha_in,

    to get the interpolated value of a curve (cl or cd, or any 3D function for that matter).

    In cases that re_in is too large, or too small, the function returns the closest available valid Reynolds data.
    """


    re_list, alpha_list, cl_list = np.unique(
        re), np.unique(alpha), np.unique(cl)

    #print("For itnerpolation:",alpha,cl)

    if re_in >= re_list.max():
        #print("Reyonlds is bigger")
        indexes = np.where(re == re_list.max())
        alpha_selected = alpha[indexes]
        cl_selected = cl[indexes]
        return np.interp(alpha_in, alpha_selected, cl_selected,left=False, right=False)

    if re_in <= re_list.min():
        #print("Reynolds is smaller")
        #print("Alpha_in",alpha_in)
        indexes = np.where(re == re_list.min())
        alpha_selected = alpha[indexes]
        cl_selected = cl[indexes]
        #print("alpha_selected:",alpha_selected)
        oo = np.interp(alpha_in, alpha_selected, cl_selected,left=False, right=False)
        #print(oo)
        return oo

    re_bottom_index = np.where(re_list < re_in)[0][-1]
    re_bottom = re_list[re_bottom_index]
    re_top = re_list[re_bottom_index+1]

    indexes_bottom = np.where(re == re_bottom)
    indexes_top = np.where(re == re_top)

    alpha_bottom = alpha[indexes_bottom]
    cl_bottom = cl[indexes_bottom]

    alpha_top = alpha[indexes_top]
    cl_top = cl[indexes_top]

    _cl_1 = np.interp(alpha_in, alpha_bottom, cl_bottom,left=False, right=False)
    _cl_2 = np.interp(alpha_in, alpha_top, cl_top,left=False, right=False)

    if _cl_1 == False or _cl_2 == False:
        return False

    cL = (_cl_2-_cl_1)/(re_top-re_bottom)*(re_in-re_bottom)+_cl_1
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
    x, y, z = np.array(x), np.array(y), np.array(z)
    xi, yi = np.linspace(x.min(), x.max(), num_x), np.linspace(
        y.min(), y.max(), num_y)
    xi, yi = np.meshgrid(xi, yi)
    zi = interp_at(x, y, z, xi.ravel(), yi.ravel(),
                   algorithm="linear", extrapolate=True)
    fun = scipy.interpolate.interp2d(xi, yi, zi, kind='linear')
    return fun


### TABLE ###

# noinspection PyArgumentList
class Table(QWidget):
    def __init__(self):
        super().__init__()
        self.selected_array = []
        self.tableWidget = QTableWidget()
        self.tableWidget.setTabKeyNavigation(False)
        self.layout = QVBoxLayout()
        self.initUI()
        self.clip = QApplication.clipboard()
        self.set_headers()

    def set_headers(self):
        self.horizontal_headers = self.tableWidget.horizontalHeader()
        self.horizontal_headers.setContextMenuPolicy(
            QtCore.Qt.CustomContextMenu)
        self.horizontal_headers.customContextMenuRequested.connect(
            self.horizontal_header_popup)
        self.vertical_headers = self.tableWidget.verticalHeader()
        self.vertical_headers.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.vertical_headers.customContextMenuRequested.connect(
            self.vertical_header_popup)

    def set_labels(self, arr):
        self.tableWidget.setHorizontalHeaderLabels(arr)

    def initUI(self):
        self.createEmpty(4, 4)

        # Add box layout, add table to box layout and add box layout to widget
        self.layout.addWidget(self.tableWidget)
        self.setLayout(self.layout)

        # Show widget
        self.show()

    def createTable(self, array):
        if len(array) > 0:
            # Create table
            self.tableWidget.setRowCount(len(array))
            self.tableWidget.setColumnCount(len(array[0]))
            i = 0
            for r in array:
                j = 0
                for c in r:
                    if not isinstance(c, str):
                        c = str(c)
                    self.tableWidget.setItem(i, j, QTableWidgetItem(c))
                    j += 1
                i += 1

            self.tableWidget.move(0, 0)

            # table selection change
            self.tableWidget.clicked.connect(self.on_click)

    def createEmpty(self, x, y):
        # Create table
        self.tableWidget.setRowCount(y)
        self.tableWidget.setColumnCount(x)

        self.tableWidget.move(0, 0)

        # table selection change
        self.tableWidget.clicked.connect(self.on_click)

    @pyqtSlot()
    def on_click(self):
        self.get_selected()

    def get_selected(self):
        self.selected_array = []

        rows_added = sorted(set(index.row()
                                for index in self.tableWidget.selectedIndexes()))
        columns_added = sorted(set(index.column()
                                   for index in self.tableWidget.selectedIndexes()))

        delta_r = rows_added[0]
        delta_c = columns_added[0]

        for r in range(rows_added[-1] - rows_added[0] + 1):
            self.selected_array.append([])
            for c in range(columns_added[-1] - columns_added[0] + 1):
                self.selected_array[r].append(None)
        for currentQTableWidgetItem in self.tableWidget.selectedItems():
            row = currentQTableWidgetItem.row() - delta_r
            column = currentQTableWidgetItem.column() - delta_c
            text = currentQTableWidgetItem.text()
            self.selected_array[row][column] = text
        return self.selected_array

    def keyPressEvent(self, e):
        if e.modifiers() & QtCore.Qt.ControlModifier:
            if e.key() == QtCore.Qt.Key_C:  # copy
                s = array_to_csv(self.get_selected())
                self.clip.setText(s)
        if e.key() == QtCore.Qt.Key_Return or e.key() == QtCore.Qt.Key_Enter:
            self.select_next_row()
        if e.modifiers() & QtCore.Qt.ControlModifier:
            if e.key() == QtCore.Qt.Key_V:  # paste
                self.paste()

        if e.modifiers() & QtCore.Qt.ControlModifier:
            if e.key() == QtCore.Qt.Key_S:  # test
                self.get_values()

        if e.key() == QtCore.Qt.Key_Delete:
            self.delete_data()

    def paste(self):
        results = []
        text = self.clip.text()
        text = text.replace("   ", "\t")
        text = text.replace("  ", "\t")
        if len(text) > 0:
            # change contents to floats
            reader = csv.reader(text.splitlines(), delimiter="\t")
            for row in reader:  # each row is a list
                results.append(row)
            numrows = len(results)
            numcolumns = len(results[0])
            selected_row = sorted(
                set(index.row() for index in self.tableWidget.selectedIndexes()))[0]
            selected_column = sorted(
                set(index.column() for index in self.tableWidget.selectedIndexes()))[0]
            if selected_row + numrows >= self.tableWidget.rowCount():
                self.tableWidget.setRowCount(selected_row + numrows)
            if selected_column + numcolumns >= self.tableWidget.columnCount():
                self.tableWidget.setColumnCount(selected_column + numcolumns)
            currow = selected_row
            for r in results:
                curcolumn = selected_column
                for c in r:
                    self.tableWidget.setItem(
                        currow, curcolumn, QTableWidgetItem(c))
                    curcolumn += 1
                currow += 1
        return

    def delete_data(self):
        rows = sorted(set(index.row() for index in self.tableWidget.selectedIndexes()))
        columns = sorted(set(index.column() for index in self.tableWidget.selectedIndexes()))
        for r in rows:
            for c in columns:
                self.tableWidget.setItem(r, c, QTableWidgetItem(""))

    def clear_table(self):
        num_rows = self.tableWidget.rowCount()
        num_columns = self.tableWidget.columnCount()
        row_indexes = range(num_rows)
        column_indexes = range(num_columns)
        for r in row_indexes:
            for c in column_indexes:
                self.tableWidget.setItem(r, c, QTableWidgetItem(""))


    def select_next_row(self):
        rows = sorted(set(index.row()
                          for index in self.tableWidget.selectedIndexes()))
        columns = sorted(set(index.column()
                             for index in self.tableWidget.selectedIndexes()))
        last_selected_row = rows[-1]
        first_selected_column = columns[0]
        num_rows = self.tableWidget.rowCount()
        if last_selected_row + 1 >= num_rows:
            self.tableWidget.insertRow(num_rows)
        self.tableWidget.setCurrentCell(
            last_selected_row + 1, first_selected_column)

    def get_values(self):
        data = []
        for row in range(self.tableWidget.rowCount()):
            data.append([])
            for column in range(self.tableWidget.columnCount()):
                item = self.tableWidget.item(row, column)
                if item == None:
                    item = ""
                else:
                    item = item.text()
                data[row].append(item)
        return data

    def contextMenuEvent(self, event):
        menu = QMenu(self)
        item = self.tableWidget.itemAt(event.pos())
        if item != None:
            delete_row = menu.addAction("delete row(s)")
            delete_column = menu.addAction("delete column(s)")
        else:
            delete_row = False
            delete_column = False

        insert_row = menu.addAction("insert row")
        insert_column = menu.addAction("insert column")

        action = menu.exec_(self.mapToGlobal(event.pos()))

        rows = sorted(set(index.row()
                          for index in self.tableWidget.selectedIndexes()), reverse=True, )
        columns = sorted(set(index.column()
                             for index in self.tableWidget.selectedIndexes()), reverse=True, )

        if action == insert_row:
            if len(rows) == 0:
                rows.append(0)
            self.tableWidget.insertRow(rows[-1])
            for c in range(self.tableWidget.columnCount()):
                self.tableWidget.setItem(rows[-1], c, QTableWidgetItem(""))
        elif action == delete_row:
            for r in rows:
                self.tableWidget.removeRow(r)
        elif action == insert_column:
            if len(columns) == 0:
                columns.append(0)
            self.tableWidget.insertColumn(columns[-1])
            for r in range(self.tableWidget.rowCount()):
                self.tableWidget.setItem(r, columns[-1], QTableWidgetItem(""))
        elif action == delete_column:
            for c in columns:
                self.tableWidget.removeColumn(c)

    def horizontal_header_popup(self, position):
        pass

    def vertical_header_popup(self, position):
        pass


### SOLIDWORKS MACRO BUILDER ###

def create_macro_text(list_of_files):
    template_start = """
    Dim swApp As Object

    Dim part As Object
    Dim boolstatus As Boolean
    Dim longstatus As Long, longwarnings As Long

    Sub main()

    Set swApp = Application.SldWorks

    template = swApp.GetUserPreferenceStringValue(swDefaultTemplatepart)
    Set part = swApp.NewDocument(template, 0, 0, 0)

    Dim myModelView As Object
    """

    template_end = """
    End Sub
    """
    str_between = ""
    for f in list_of_files:
        str_between += 'boolstatus = part.InsertCurveFile("%s")\n' % f
    return template_start+"\n"+str_between+"\n"+template_end

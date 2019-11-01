# coding=utf-8
__author__ = "Miha Smrekar"
__credits__ = ["Miha Smrekar"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Miha Smrekar"
__email__ = "miha.smrekar9@gmail.com"
__status__ = "Development"

from popravki import METHODS_STRINGS
from utils import get_centroid_coordinates
from sw_macro_builder import create_macro_text
from visualize import create_3d_blade
from xfoil import generate_polars_data
from montgomerie import Montgomerie
from polars import scrape_data
from interpolator import interp_at
from utils import interpolate_geom, to_float, fltr
from optimization import optimize_angles_genetic
from utils import interpolate_geom, to_float, fltr, QDarkPalette, create_folder
from turbine_data import SET_INIT
from table import Table
from results import ResultsWindow
from calculation_runner import calculate_power_3d
from matplotlib import cm
import mpl_toolkits.mplot3d as mp3d
from matplotlib.widgets import Slider, Button, RadioButtons
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from scipy.interpolate import interp1d
from scipy import interpolate
import numpy as np
from numpy import array
from PyQt5.QtGui import QPalette, QColor
from PyQt5.QtWidgets import (QComboBox, QMainWindow, QPushButton, QTextEdit, QWidget, QFormLayout, QLabel, QLineEdit,
                             QGridLayout, QCheckBox, QStyleFactory, QMessageBox, QAction, QFileDialog, QSlider,
                             QTabWidget, QApplication, QScrollArea, QVBoxLayout)
from PyQt5 import QtCore, QtGui
from PyQt5.QtCore import QThread, QTextStream, pyqtSignal, QProcess, QRect, Qt

from multiprocessing import Process, Manager, Queue
import multiprocessing
import json
from pprint import pprint
import time
import sys
import ctypes
import os
import pyqtgraph as pg
from functools import partial
import traceback
import signal

np.set_printoptions(threshold=sys.maxsize)

TITLE_STR = "BEM analiza v%s" % __version__


class MainWindow(QMainWindow):
    """
    Main class that sets up the GUI layout.
    All QWidgets within the MainWindow class are displayed using the custom TabWidget.
    The MainWindow class also holds the functions for saving and loading data from files.

    The MainWindow class stores references to all other subclasses used in the program.
    There are four subclasses in MainWindow called AirfoilManager, WindTurbineProperties, Analysis, and ThreadGetter.
    """
    emitter_add = pyqtSignal(str)
    emitter_done = pyqtSignal()

    def __init__(self, width, height):
        super().__init__()

        self.statusBar()
        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu('&File')
        saveFile = QAction("&Save File", self)
        saveFile.setShortcut("Ctrl+S")
        saveFile.setStatusTip('Save File')
        saveFile.triggered.connect(self.file_save)
        loadFile = QAction("&Load File", self)
        loadFile.setShortcut("Ctrl+L")
        loadFile.setStatusTip('Load File')
        loadFile.triggered.connect(self.file_load)
        getSettings = QAction("Get settings", self)
        getSettings.triggered.connect(self.get_all_settings)
        fileMenu.addAction(saveFile)
        fileMenu.addAction(loadFile)
        fileMenu.addAction(getSettings)

        self.screen_width = width
        self.screen_height = height
        self.setGeometry(width * 0.125, height * 0.125, width * 0.75, height * 0.75)
        self.setWindowTitle(TITLE_STR)
        self.tab_widget = TabWidget(self)
        self.setCentralWidget(self.tab_widget)

        self.curve_manager = AirfoilManager(self)
        self.tab_widget.add_tab(self.curve_manager, "Airfoil management")

        self.wind_turbine_properties = WindTurbineProperties(self)
        self.tab_widget.add_tab(self.wind_turbine_properties, "Turbine info")

        self.analysis = Analysis(self)
        self.tab_widget.add_tab(self.analysis, "Analysis")

        self.getter = ThreadGetter(self)

        self.optimization = Optimization(self)
        self.tab_widget.add_tab(self.optimization, "Optimization")

        self.running = False
        self.manager = Manager()
        self.set_all_settings(SET_INIT)

        create_folder("foils")  # Used by XFoil

        self.show()

    def set_title(self):
        """
        Sets the title of the GUI window, based on the name of the wind turbine.
        """
        s = self.wind_turbine_properties.name.text()
        if s == "":
            self.setWindowTitle(TITLE_STR)
        else:
            self.setWindowTitle(TITLE_STR + " - " + s)

    def file_save(self):
        """
        Saves the wind turbine data to a file.
        """
        name = QFileDialog.getSaveFileName(self, 'Save File')[0]
        if name != "":
            file = open(name, 'w')
            d = self.get_all_settings()
            d_to_save = fltr(d, (float, int, list, str, bool, np.ndarray))
            json_d = json.dumps(d_to_save)
            file.write(json_d)
            file.close()

    def file_load(self):
        """
        Loads the wind turbine data from a file. Also clears the calculation text areas and sets the appropriate title.
        """
        file_path = QFileDialog.getOpenFileName(self, "Load File")[0]
        if file_path != "":
            with open(file_path, "r") as fp:
                data = json.load(fp)
            self.set_all_settings(data)
        self.analysis.clear()
        self.optimization.clear()
        self.set_title()

    def get_all_settings(self):
        """
        Used to save the input configuration for the BEM method calculation.

        It fetches the settings from the four subclasses and combines them into a single settings dictionary.

        After all the settings are combined, it calculates the interpolation of the wind turbine geometry as per
        user choice. (e.g. if you want to increase or decrease the number of sections with regard to the original
        geometry).

        :return: dict: Settings
        """
        try:
            properties = self.wind_turbine_properties.get_settings()
            settings = self.analysis.get_settings()
            opt_settings = self.optimization.get_settings()
            curve_manager_settings = self.curve_manager.get_settings()

            out = {**properties, **settings, **opt_settings, **curve_manager_settings}
            _r = out["r"]
            _c = out["c"]
            _theta = out["theta"]
            _foils = out["foils"]
            r, c, theta, foils, dr = interpolate_geom(_r, _c, _theta, _foils, out["num_interp"], out["linspace_interp"])
            out["r"], out["c"], out["theta"], out["foils"], out["dr"] = r, c, theta, foils, dr
            out["r_in"], out["c_in"], out["theta_in"], out["foils_in"] = _r, _c, _theta, _foils

            return out
        except:
            msg = ErrorMessageBox()
            return None

    # noinspection PyBroadException
    def set_all_settings(self, inp_dict):
        """
        Sets the settings from the appropriate dictionary object by sending the input dictionary object to each of the
        subclasses.
        :param inp_dict: dict: Settings dictionary.
        """
        try:
            self.analysis.set_settings(inp_dict)

            self.optimization.set_settings(inp_dict)

            self.wind_turbine_properties.set_settings(inp_dict)

            self.curve_manager.set_settings(inp_dict)
        except:
            msg = ErrorMessageBox()

    def get_input_params(self):
        """
        Used to fetch the BEM calculation settings, create (reset) the multiprocessing queues used for communication
        between processes and for creating the final dictionary object which is then sent to the calculation module.
        :return: dict: Settings object (with multiprocessing communication entries).
        """
        settings = self.get_all_settings()
        if settings == None:
            return None
        self.return_print = self.manager.list([])
        self.return_results = self.manager.list([])
        self.queue_pyqtgraph = self.manager.list([])
        self.end_of_file = self.manager.Value("EOF", False)
        inp_params = {**settings, "return_print": self.return_print, "return_results": self.return_results,
                      "EOF": self.end_of_file}
        return inp_params

    def set_process_running(self):
        """
        Enables/disables the main stop/start buttons and sets the main boolean self.running to True.

        Used at the beginning of the calculation process.
        """
        self.analysis.buttonRun.setEnabled(False)
        self.optimization.buttonOptimization.setEnabled(False)
        self.analysis.buttonStop.setEnabled(True)
        self.optimization.buttonStop.setEnabled(True)
        self.running = True

    def set_process_stopped(self):
        """
        Enables/disables the main start/stop buttons and sets the main boolean self.running to True.

        Used at the end of the calculation process.
        """
        self.analysis.buttonRun.setEnabled(True)
        self.optimization.buttonOptimization.setEnabled(True)
        self.analysis.buttonStop.setEnabled(False)
        self.optimization.buttonStop.setEnabled(False)
        self.running = False


class WindTurbineProperties(QWidget):
    """
    Class used for storing the main wind turbine information, such as its name, number of blades, radiuses, and
    its geometry information (radius, chord length, twist angle and airfoil name for each wind blade section).
    """
    def __init__(self, parent=None):
        super(WindTurbineProperties, self).__init__(parent)

        self.main = self.parent()

        grid = QGridLayout()
        self.setLayout(grid)

        left = QWidget()
        fbox = QFormLayout()
        left.setLayout(fbox)

        self.table_properties = Table()
        self.table_properties.createEmpty(4, 30)
        self.table_properties.set_labels(["r [m]", "c [m]", "theta [deg]", "airfoil"])

        grid.addWidget(left, 1, 1)
        grid.addWidget(self.table_properties, 1, 2)

        _name = QLabel("Turbine Name")
        self.name = QLineEdit()
        fbox.addRow(_name, self.name)
        self.name.textEdited.connect(self.main.set_title)

        _Rhub = QLabel("Hub radius [m]")
        self.Rhub = QLineEdit()
        self.Rhub.setText("0.1")
        fbox.addRow(_Rhub, self.Rhub)

        _R = QLabel("Tip radius [m]")
        self.R = QLineEdit()
        self.R.setText("0.776")
        fbox.addRow(_R, self.R)

        _B = QLabel("Number of blades")
        self.B = QLineEdit()
        self.B.setText("5")
        fbox.addRow(_B, self.B)

        self.export_button = QPushButton("Export curve data")
        self.export_button.clicked.connect(self.export)
        fbox.addRow("Export:", self.export_button)

        self.flip_turning_direction = QCheckBox()
        fbox.addRow("Flip turning direction", self.flip_turning_direction)

        self.propeller_geom = QCheckBox()
        fbox.addRow("Propeller", self.propeller_geom)

    def get_settings(self):
        """
        Used to get the basic wind turbine settings in a dictionary format, namely the:
        -hub radius (Rhub)
        -tip radius (R)
        -number of blades (B)
        -turbine name (turbine_name)
        -geometry (r,c,theta,foils) as four separate list objects
        :return: dict: Settings dictionary (Basic wind turbine information)
        """
        out_properties = {"Rhub": to_float(self.Rhub.text()), "R": to_float(self.R.text()), "B": int(self.B.text()),
                          "turbine_name": self.name.text(), }
        geom_array = self.table_properties.get_values()
        r, c, theta, foils = [], [], [], []
        for row in geom_array:
            if row[0] != "" and row[1] != "" and row[2] != "":
                r.append(to_float(row[0]))
                c.append(to_float(row[1]))
                theta.append(to_float(row[2]))
                foils.append(row[3])
        out_properties["r"] = array(r)
        out_properties["c"] = array(c)
        out_properties["theta"] = array(theta)
        out_properties["foils"] = foils
        return out_properties

    def set_settings(self, dict_settings):
        """
        Reads the settings from the input dictionary and sets the values in the GUI.
        :param dict_settings: dict: Settings dictionary.
        """
        if "Rhub" in dict_settings:
            t = str(dict_settings["Rhub"])
            self.Rhub.setText(t)
        if "R" in dict_settings:
            t = str(dict_settings["R"])
            self.R.setText(t)
        if "B" in dict_settings:
            t = str(dict_settings["B"])
            self.B.setText(t)
        if "r_in" in dict_settings and "c_in" in dict_settings and "theta_in" in dict_settings and "foils_in" in dict_settings:
            _array = []
            for r in range(len(dict_settings["r_in"])):
                _r = dict_settings["r_in"][r]
                _c = dict_settings["c_in"][r]
                _theta = dict_settings["theta_in"][r]
                _f = dict_settings["foils_in"][r]
                _array.append([_r, _c, _theta, _f])
            self.table_properties.createTable(_array)
        if "turbine_name" in dict_settings:
            t = str(dict_settings["turbine_name"])
            self.name.setText(t)
        else:
            self.name.setText("")

    def export(self):
        """
        Exports the wind turbine geometry data as spatial data (XYZ points), and creates a VB macro, which
        you can use to import the geometry into the Solidworks 3D modeller.
        :return:
        """
        print("Getting settings...")
        settings_fetched = self.parent().parent().parent().get_all_settings()
        if settings_fetched == None:
            return
        data = create_3d_blade(settings_fetched, self.flip_turning_direction.isChecked(), self.propeller_geom.isChecked())
        w = MatplotlibWindow()
        w.ax = w.figure.add_subplot(111, projection="3d")
        w.ax.scatter(data["X"], data["Y"], data["Z"])
        X, Y, Z = array(data["X"]), array(data["Y"]), array(data["Z"])
        max_range = array([X.max() - X.min(), Y.max() - Y.min(), Z.max() - Z.min()]).max() / 2.0
        mid_x = (X.max() + X.min()) * 0.5
        mid_y = (Y.max() + Y.min()) * 0.5
        mid_z = (Z.max() + Z.min()) * 0.5
        w.ax.set_xlim(mid_x - max_range, mid_x + max_range)
        w.ax.set_ylim(mid_y - max_range, mid_y + max_range)
        w.ax.set_zlim(mid_z - max_range, mid_z + max_range)
        # w.ax.set_aspect("equal")

        create_folder("export")
        folder_path = os.path.join("export", SET_INIT["turbine_name"])
        create_folder(folder_path)

        filenames = []
        print("Exporting... (and converting m to mm)")
        for z, x_data, y_data in data["data"]:
            print("Exporting z=" + str(z), "[m]")
            z = z * 1e3  # in mm
            file_name = os.path.join(folder_path, "z_%s.txt" % z)
            filenames.append(os.path.join(os.getcwd(), file_name))
            # print(file_name)
            f = open(os.path.join(folder_path, "z_%s.txt" % z), "w")
            for x, y in zip(x_data, y_data):
                x, y = x * 1e3, y * 1e3  # in mm
                f.write("%s\t%s\t%s\n" % (x, y, z))
            f.close()
        print("Filenames:", filenames)
        macro_text = create_macro_text(filenames)

        print("'===============MACRO START==================")
        print(macro_text)
        print("'===============MACRO   END==================")


class AirfoilManager(QWidget):
    emitter = pyqtSignal(str)

    def __init__(self, parent=None):
        super(AirfoilManager, self).__init__(parent)

        self.main = self.parent()

        self.grid = QGridLayout()
        self.setLayout(self.grid)

        self.tab_widget = TabWidget(self)
        self.grid.addWidget(self.tab_widget, 2, 0)

        self.upper_widget = QWidget()
        self.upper_layout = QGridLayout()
        self.upper_widget.setLayout(self.upper_layout)
        self.grid.addWidget(self.upper_widget, 1, 0)

        self.button_add_foil = QPushButton("Add airfoil")
        self.button_add_foil.clicked.connect(self.add_foil_popup)
        self.button_remove_foil = QPushButton("Remove foil")
        self.button_remove_foil.clicked.connect(self.tab_widget.remove_current_tab)
        self.button_rename_foil = QPushButton("Rename foil")
        self.button_rename_foil.clicked.connect(self.rename_foil_popup)

        self.upper_layout.addWidget(self.button_add_foil, 0, 1)
        self.upper_layout.addWidget(self.button_remove_foil, 0, 2)
        self.upper_layout.addWidget(self.button_rename_foil, 0, 3)

    def add_foil_popup(self):
        self.emitter.connect(self.add_foil)
        self.p = PopupText("foil name", "airfoil_name", self.emitter)
        self.p.setGeometry(QRect(100, 100, 400, 200))
        self.p.show()

    def add_foil(self, string):
        c = Airfoils(string, self)
        self.tab_widget.add_tab(c, string)

    def rename_foil_popup(self):
        self.emitter.connect(self.rename_foil)
        self.p = PopupText("foil name", self.tab_widget.current_tab_name(), self.emitter)
        self.p.setGeometry(QRect(100, 100, 400, 200))
        self.p.show()

    def rename_foil(self, string):
        self.tab_widget.rename_current_tab(string)  # self.tab_widget.tabs

    def get_settings(self):
        out = {}
        i = 0

        # TODO Dont rely on name being set correctly in n!
        for w, n in self.tab_widget.tabs:
            out[n] = w.get_settings()
            i += 1

        return {"airfoils": out}

    def set_settings(self, dict_settings):
        self.tab_widget.remove_all_tabs()
        if "airfoils" in dict_settings:
            if len(dict_settings["airfoils"]) > 0:
                for c_name, c_dict in dict_settings["airfoils"].items():
                    curve_widget = Airfoils(c_name, self)
                    curve_widget.set_settings(c_dict)
                    self.tab_widget.add_tab(curve_widget, c_name)


class PopupText(QWidget):
    def __init__(self, message="message", default_str="", emitter=None):
        QWidget.__init__(self)

        self.emitter = emitter

        self.layout = QGridLayout()
        self.setLayout(self.layout)

        self.message = QLabel(message)

        self.layout.addWidget(self.message, 0, 0)

        self.inp = QLineEdit()
        self.inp.setText(default_str)
        self.layout.addWidget(self.inp, 1, 0)

        self.button = QPushButton("OK")
        self.button.clicked.connect(self.send_signal)
        self.layout.addWidget(self.button, 2, 0)

    def send_signal(self):
        if self.emitter != None:
            self.emitter.emit(self.inp.text())
            self.emitter.disconnect()
        self.close()


class Airfoils(QWidget):
    def __init__(self, airfoil_name, parent=None):
        super(Airfoils, self).__init__(parent)

        self.curves = Curves()

        self.viewer = CurveViewer(self)

        self.airfoil_name = airfoil_name

        self.grid = QGridLayout()
        self.setLayout(self.grid)

        self.left = QWidget()
        self.fbox = QFormLayout()
        self.left.setLayout(self.fbox)

        self.scroll_area = QScrollArea()
        self.scroll_widget = QWidget()
        self.scroll_widget_layout = QVBoxLayout()

        self.scroll_widget.setLayout(self.scroll_widget_layout)
        self.scroll_area.setWidget(self.left)
        self.scroll_area.setWidgetResizable(True)

        self.grid.addWidget(self.scroll_area, 1, 1)

        self.interp_function_cl = None
        self.interp_function_cd = None

        self.table_dat = Table()
        self.table_dat.createEmpty(2, 50)
        self.table_dat.set_labels(["x", "y"])
        self.grid.addWidget(self.table_dat, 1, 2)

        self.plt = plt.figure(figsize=(10, 5))
        self.ax = self.plt.add_subplot(111)

        self.canvas = FigureCanvas(self.plt)
        self.grid.addWidget(self.canvas, 1, 3)
        
        toolbar = NavigationToolbar(self.canvas, self)
        self.grid.addWidget(toolbar, 2, 3)
        
        self.buttonRefresh = QPushButton("Refresh curve")
        self.fbox.addRow(self.buttonRefresh)

        self.buttonRefresh.clicked.connect(self.refresh)
        self.link = QLineEdit("link (airfoiltools.com)")
        self.fbox.addRow(self.link)

        self.button_open_viewer = QPushButton("Open Curve Viewer")
        self.button_open_viewer.clicked.connect(self.open_viewer)
        self.fbox.addRow(self.button_open_viewer)

        self.button_generate_curves_xfoil = QPushButton("Generate xfoil curves [debug]")
        self.button_generate_curves_xfoil.clicked.connect(self.generate_curves_xfoil)
        self.fbox.addRow(self.button_generate_curves_xfoil)

        self.button_generate_curves_link = QPushButton("Generate curves from link")
        self.button_generate_curves_link.clicked.connect(self.generate_curves_link)
        self.fbox.addRow(self.button_generate_curves_link)

        self.button_visualize = QPushButton("Create curve visualization")
        self.button_visualize.clicked.connect(self.visualize)
        self.fbox.addRow(self.button_visualize)

        self.get_centroid_button = QPushButton("Calculate centroid")
        self.get_centroid_button.clicked.connect(self.calculate_centroid)
        self.fbox.addRow(self.get_centroid_button)

        self.centroid_widget = QWidget()
        self.centroid_grid = QGridLayout()
        self.centroid_widget.setLayout(self.centroid_grid)
        self.centroid_label = QLabel("Centroid coordinates:")
        self.fbox.addRow(self.centroid_widget)

        self.centroid_x_edit = QLineEdit()
        self.centroid_y_edit = QLineEdit()

        self.centroid_grid.addWidget(self.centroid_label, 1, 1)
        self.centroid_grid.addWidget(self.centroid_x_edit, 1, 2)
        self.centroid_grid.addWidget(self.centroid_y_edit, 1, 3)

        self.ncrit_selection = QComboBox()
        self._ncrit_selection = QLabel("Ncrit")
        self.fbox.addRow(self._ncrit_selection,self.ncrit_selection)


    def visualize(self):
        print("Visualizing")
        data = self.curves.gather_curves()

        re = data[:, 0]
        alpha = data[:, 2]
        cl = data[:, 3]
        cd = data[:, 4]

        re_min, re_max = data[:, 0].min(), data[:, 0].max()
        alpha_min, alpha_max = data[:, 2].min(), data[:, 2].max()

        x, y = np.linspace(re_min, re_max, 10), np.linspace(alpha_min, alpha_max, 180)
        xi, yi = np.meshgrid(x, y)
        xi, yi = xi.flatten(), yi.flatten()
        z_1 = interp_at(re, alpha, cl, xi, yi)
        z_2 = interp_at(re, alpha, cd, xi, yi)
        w = MatplotlibWindow()
        w.ax = w.figure.add_subplot(111, projection="3d")
        p = w.ax.plot_trisurf(xi, yi, z_1, cmap=cm.coolwarm)
        w.ax.set_xlabel("Reynolds", fontsize=15, labelpad=20)
        w.ax.set_ylabel(r'$\alpha$ [°]', fontsize=15, labelpad=20)
        w.ax.set_zlabel("Cl", fontsize=15, labelpad=20)
        w.ax.xaxis.set_tick_params(labelsize=12)
        w.ax.yaxis.set_tick_params(labelsize=12)
        w.ax.zaxis.set_tick_params(labelsize=12)
        bar = w.figure.colorbar(p)
        bar.ax.set_xlabel('Cl', fontsize=15, labelpad=20)

        w2 = MatplotlibWindow()
        w2.ax = w2.figure.add_subplot(111, projection="3d")
        p = w2.ax.plot_trisurf(xi, yi, z_2, cmap=cm.coolwarm)
        w2.ax.set_xlabel("Reynolds", fontsize=15, labelpad=20)
        w2.ax.set_ylabel(r'$\alpha$ [°]', fontsize=15, labelpad=20)
        w2.ax.set_zlabel("Cd", fontsize=15, labelpad=20)
        w2.ax.xaxis.set_tick_params(labelsize=12)
        w2.ax.yaxis.set_tick_params(labelsize=12)
        w2.ax.zaxis.set_tick_params(labelsize=12)
        bar2 = w2.figure.colorbar(p)
        bar2.ax.set_xlabel('Cd', fontsize=15, labelpad=20)

    def open_viewer(self):
        print("opening viewer")
        self.viewer.show()
        self.viewer.generate_views()

    def generate_interp_functions(self):
        data = self.gather_curves()
        x, y = self.get_x_y()
        self.interp_function_cl, self.interp_function_cd = get_cl_cd_interpolation_function(data, x, y)

    def generate_curves_xfoil(self):
        print("Generating xfoil curves")
        data = generate_polars_data(self.airfoil_name + ".dat")
        self.populate_curve_list(data)
        print("Done")

    def generate_curves_link(self):
        print("Scraping from link...")
        data = scrape_data(self.link.text())
        self.populate_curve_list(data)
        print("Done.")
        self.refresh()

    def populate_curve_list(self, data):
        self.curves.curve_list = []
        x, y = self.get_x_y()
        ncrit_list = np.unique(data[:, 1])
        for ncrit_selected in ncrit_list:
            rows_with_ncrit = data[np.in1d(data[:,1],ncrit_selected)]
            Re_list = np.unique(rows_with_ncrit[:, 0])
            for Re in Re_list:
                rows_with_Re = rows_with_ncrit[np.in1d(rows_with_ncrit[:, 0], Re)]
                _alpha = rows_with_Re[:, 2].flatten()
                _cl = rows_with_Re[:, 3].flatten()
                _cd = rows_with_Re[:, 4].flatten()
                c = Curve()
                c.create(x=x, y=y, Re=Re, ncrit=ncrit_selected, alpha=_alpha, cl=_cl, cd=_cd)
                self.curves.add(c)

    def refresh(self):
        self.ax.clear()
        x_values = []
        y_values = []
        array_dat = self.table_dat.get_values()
        for r in array_dat:
            if r[0] != "" and r[1] != "":
                # noinspection PyBroadException
                try:
                    _x = to_float(r[0])
                    _y = to_float(r[1])
                except:
                    print("Error drawing airfoil because _x or _y isn't a float.")
                    return
                x_values.append(_x)
                y_values.append(_y)
        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(-0.5, 0.5)
        self.ax.plot(x_values, y_values)
        try:
            centroid_x = float(self.centroid_x_edit.text())
            centroid_y = float(self.centroid_y_edit.text())
            self.ax.plot(centroid_x, centroid_y, "r+")
        except:
            msg = ErrorMessageBox()
            

        self.plt.canvas.draw()

        self.refresh_ncrits_combobox()



    def get_max_thickness(self):
        x = []
        y = []
        array_dat = self.table_dat.get_values()
        for r in array_dat:
            if r[0] != "" and r[1] != "":
                _x = to_float(r[0])
                _y = to_float(r[1])
                x.append(_x)
                y.append(_y)
        if len(y) > 0:
            y_max = np.max(y)
            y_min = np.min(y)
            thickness = (abs(y_max) + abs(y_min)) / 1
            return thickness
        return None

    def get_x_y(self):
        x = []
        y = []

        array_dat = self.table_dat.get_values()
        for r in array_dat:
            if r[0] != "" and r[1] != "":
                _x = to_float(r[0])
                _y = to_float(r[1])
                x.append(_x)
                y.append(_y)
        return x, y

    def calculate_centroid(self):
        foil_x, foil_y = self.get_x_y()
        x, y = get_centroid_coordinates(foil_x, foil_y)
        self.centroid_x_edit.setText(str(round(x,10)))
        self.centroid_y_edit.setText(str(round(y,10)))
        return x, y

    def get_ncrits(self):
        curves = self.curves.gather_curves()
        ncrit_list = np.unique(curves[:,1])
        return ncrit_list

    def refresh_ncrits_combobox(self):
        self.ncrit_selection.clear()
        self.ncrit_selection.addItems([str(n) for n in list(self.get_ncrits())])

    def get_settings(self):
        x, y = self.get_x_y()
        centroid_x = float(self.centroid_x_edit.text())
        centroid_y = float(self.centroid_y_edit.text())
        out = {"x": x,
               "y": y,
               "max_thickness": self.get_max_thickness(),
               "link": self.link.text(),
               "interp_function_cl": self.interp_function_cl,
               "interp_function_cd": self.interp_function_cd,
               "curves": self.curves.save_curves(),
               "gathered_curves": self.curves.gather_curves(),
               "centroid_x": centroid_x,
               "centroid_y": centroid_y,
               "ncrit_selected": float(self.ncrit_selection.currentText())}
        return out

    def set_settings(self, dict_settings):
        array_dat = []
        if len(dict_settings["x"]) > 0 and len(dict_settings["y"]) > 0:
            for r in range(len(dict_settings["x"])):
                array_dat.append([str(dict_settings["x"][r]), str(dict_settings["y"][r])])
            self.table_dat.createTable(array_dat)
            self.calculate_centroid()
            

        self.link.setText(dict_settings["link"])
        self.curves.load_curves(dict_settings["curves"])
        self.refresh()


class MatplotlibWindow(QWidget):
    def __init__(self):
        super(MatplotlibWindow, self).__init__(None)
        self.layout = QGridLayout()
        self.setLayout(self.layout)
        self.figure = plt.figure(figsize=(10, 5))
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setMinimumSize(500, 500)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.layout.addWidget(self.canvas)
        self.layout.addWidget(self.toolbar)
        self.show()


class Curves:
    def __init__(self):
        self.curve_list = []

    def add(self, curve):
        self.curve_list.append(curve)

    def save_curves(self):
        out_list = []
        for c in self.curve_list:
            data_curve = c.save_curve()
            out_list.append(data_curve)
        return out_list

    def load_curves(self, out):
        self.curve_list = []
        for data_curve in out:
            c = Curve()
            c.load_curve(data_curve)
            self.curve_list.append(c)

    def gather_curves(self):
        out = []
        for curve in self.curve_list:
            alpha, cl, cd = curve.get_combined_curve()
            for i in range(len(alpha)):
                Re = curve.Re
                ncrit = curve.ncrit
                _alpha = alpha[i]
                _cl = cl[i]
                _cd = cd[i]
                out.append([Re, ncrit, _alpha, _cl, _cd])
        out = array(out)
        return out


class Curve:
    def __init__(self):
        self.x = None
        self.y = None
        self.Re = None
        self.ncrit = None
        self.alpha = None
        self.cl = None
        self.cd = None
        self.A = None
        self.B = None
        self.Am = None
        self.Bm = None
        self.m_CD90 = None
        self.slope = None

    def create(self, x, y, Re, ncrit, alpha, cl, cd):
        self.x = x
        self.y = y
        self.Re = Re
        self.ncrit = ncrit
        self.alpha = alpha
        self.cl = cl
        self.cd = cd
        self.A = 5
        self.B = 5
        self.Am = 5
        self.Bm = 5
        self.m_CD90 = 2.0
        self.slope = 0.106

    def get_cl_curve(self):
        return self.alpha_cl, self.cl

    def get_extrapolated_curve(self):
        M = Montgomerie(x=self.x, y=self.y, alpha=self.alpha, Cl=self.cl, Cd=self.cd, Re=self.Re, A=self.A, Am=self.Am,
                        B=self.B, Bm=self.Bm, m_CD90=self.m_CD90, slope=self.slope)
        alpha, cl, cd = M.calculate_extrapolation()
        return alpha, cl, cd

    def get_combined_curve(self):
        M = Montgomerie(x=self.x, y=self.y, alpha=self.alpha, Cl=self.cl, Cd=self.cd, Re=self.Re, A=self.A, Am=self.Am,
                        B=self.B, Bm=self.Bm, m_CD90=self.m_CD90, slope=self.slope)
        _alpha, _cl, _cd = M.calculate_extrapolation()
        cl_out, cd_out = [], []
        f_cl = interp1d(self.alpha, self.cl, bounds_error=True)
        f_cd = interp1d(self.alpha, self.cd, bounds_error=True)
        for i in range(len(_alpha)):
            # x.append(self.Re)
            # y.append(m_Alpha[i])
            a = _alpha[i]
            try:
                cl = f_cl(a)
            except ValueError:
                cl = _cl[i]
            # try:
            #     cd = f_cd(a)
            # except ValueError:
            #     cd = _cd[i]
            # tukaj vzamem  samo Montgomerie interpolacijo cd, za lazjo interpolacijo
            cd = _cd[i]
            cl_out.append(cl)
            cd_out.append(cd)
        return _alpha, cl_out, cd_out

    def save_curve(self):
        out = {"x": list(self.x), "y": list(self.y), "Re": self.Re, "ncrit": self.ncrit, "alpha": list(self.alpha),
               "cl": list(self.cl), "cd": list(self.cd), "A": self.A, "B": self.B, "Am": self.Am, "Bm": self.Bm,
               "m_CD90": self.m_CD90, "slope": self.slope}
        return out

    def load_curve(self, out):
        self.x = out["x"]
        self.y = out["y"]
        self.Re = out["Re"]
        self.ncrit = out["ncrit"]
        self.alpha = out["alpha"]
        self.cl = out["cl"]
        self.cd = out["cd"]
        self.A = out["A"]
        self.B = out["B"]
        self.Am = out["Am"]
        self.Bm = out["Bm"]
        self.m_CD90 = out["m_CD90"]
        self.slope = out["slope"]


class CurveViewer(QWidget):
    def __init__(self, parent=None):
        super(CurveViewer, self).__init__(None)
        self.resize(1600, 768)
        self.parent = parent
        self.grid = QGridLayout()
        self.setLayout(self.grid)
        self.button = QPushButton("Close")
        self.grid.addWidget(self.button, 1, 1)
        self.button.clicked.connect(self.close)
        self.button_refresh = QPushButton("Refresh")
        self.grid.addWidget(self.button_refresh, 1, 2)
        self.button_refresh.clicked.connect(self.generate_views)

        self.bottom = QWidget()
        self.grid_curves = QGridLayout()
        self.bottom.setLayout(self.grid_curves)

        self.scroll_area = QScrollArea()
        self.scroll_widget = QWidget()
        self.scroll_widget_layout = QVBoxLayout()

        self.scroll_widget.setLayout(self.scroll_widget_layout)
        self.scroll_area.setWidget(self.bottom)
        self.scroll_area.setWidgetResizable(True)
        self.grid.addWidget(self.scroll_area, 2, 1, 2, 2)

        # self.generate_views()

    def generate_views(self):

        # delete stuff already here
        for i in reversed(range(self.grid_curves.count())):
            self.grid_curves.itemAt(i).widget().setParent(None)

        # for i in range(10):
        #    control = CurveControl(self,None)
        #    self.grid_curves.addWidget(control)

        for curve in self.parent.curves.curve_list:
            label = QLabel("Re = " + str(round(curve.Re,2)) + ", Ncrit = " + str(round(curve.ncrit,2)) )
            control = CurveControl(self, curve)
            control.update()
            self.grid_curves.addWidget(label)
            self.grid_curves.addWidget(control)


class CurveControl(QWidget):
    def __init__(self, parent=None, curve=None):
        super(CurveControl, self).__init__(parent)
        # self.setMinimumSize(300,400)
        self.parent = parent

        self.layout = QGridLayout()
        self.setLayout(self.layout)

        self.curve = curve

        self.right = QWidget()
        self.right_layout = QFormLayout()
        self.right.setLayout(self.right_layout)

        self.A = QLineEdit(str(self.curve.A))
        self.B = QLineEdit(str(self.curve.B))
        self.Am = QLineEdit(str(self.curve.Am))
        self.Bm = QLineEdit(str(self.curve.Bm))

        self.A = QSlider(Qt.Horizontal)
        self.A.setMinimum(-10)
        self.A.setMaximum(30)
        self.A.setValue(self.curve.A)
        self.A.setTickPosition(QSlider.TicksBelow)
        self.A.setTickInterval(1)
        self.A.valueChanged.connect(self.update)

        self.B = QSlider(Qt.Horizontal)
        self.B.setMinimum(1)
        self.B.setMaximum(100)
        self.B.setValue(self.curve.B)
        self.B.setTickPosition(QSlider.TicksBelow)
        self.B.setTickInterval(1)
        self.B.valueChanged.connect(self.update)

        self.Am = QSlider(Qt.Horizontal)
        self.Am.setMinimum(1)
        self.Am.setMaximum(80)
        self.Am.setValue(self.curve.Am)
        self.Am.setTickPosition(QSlider.TicksBelow)
        self.Am.setTickInterval(1)
        self.Am.valueChanged.connect(self.update)

        self.Bm = QSlider(Qt.Horizontal)
        self.Bm.setMinimum(1)
        self.Bm.setMaximum(70)
        self.Bm.setValue(self.curve.Bm)
        self.Bm.setTickPosition(QSlider.TicksBelow)
        self.Bm.setTickInterval(1)
        self.Bm.valueChanged.connect(self.update)

        self.m_CD90 = QLineEdit(str(self.curve.m_CD90))
        self.m_CD90.textChanged.connect(self.update)
        self.slope = QLineEdit(str(self.curve.slope))
        self.slope.textChanged.connect(self.update)

        self.right_layout.addRow("A", self.A)
        self.right_layout.addRow("B", self.B)
        self.right_layout.addRow("A-", self.Am)
        self.right_layout.addRow("B-", self.Bm)
        self.right_layout.addRow("CD@90°", self.m_CD90)
        self.right_layout.addRow("Slope", self.slope)

        self.layout.addWidget(self.right, 1, 2)

        self.left = QWidget()
        self.left_layout = QGridLayout()
        self.left.setLayout(self.left_layout)

        self.figure = plt.figure(figsize=(10, 5))
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setMinimumSize(500, 500)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.ax = self.figure.add_subplot(111)
        self.left_layout.addWidget(self.canvas)
        self.left_layout.addWidget(self.toolbar)

        self.layout.addWidget(self.left, 1, 1)

        self.button_update = QPushButton("update")
        self.button_update.clicked.connect(self.draw_extrapolation)
        self.left_layout.addWidget(self.button_update)

        # self.draw_base()

        self.show()

    def clear(self):
        self.ax.cla()

    def draw_base(self):
        self.ax.plot(self.curve.alpha, self.curve.cl)
        self.ax.plot(self.curve.alpha, self.curve.cd, "o-")
        self.canvas.draw()

    def draw_extrapolation(self):
        self.clear()

        self.draw_base()

        alpha, cl, cd = self.curve.get_extrapolated_curve()
        self.ax.plot(alpha, cl, "g.")
        self.ax.plot(alpha, cd, "r.")
        self.canvas.draw()

    def update(self):
        self.curve.A = int(self.A.value())
        self.curve.B = int(self.B.value())
        self.curve.Am = int(self.Am.value())
        self.curve.Bm = int(self.Bm.value())
        try:
            self.curve.m_CD90 = float(self.m_CD90.text())
            self.curve.slope = float(self.slope.text())
            if self.curve.slope == 0:
                self.curve.slope = 1.0
        except:
            print("Error in slope or m_CD90")
        self.draw_extrapolation()


class Analysis(QWidget):
    def __init__(self, parent=None):
        super(Analysis, self).__init__(parent)

        self.main = self.parent()

        self.settings = {"propeller_mode": False, "tip_loss": False, "hub_loss": False, "new_tip_loss": False,
                         "new_hub_loss": False, "cascade_correction": False, "skewed_wake_correction": False,
                         "rotational_augmentation_correction": False, "rotational_augmentation_correction_method": 1,
                         "mach_number_correction": False, "max_iterations": 100, "convergence_limit": 0.001,
                         "rho": 1.225, "method": 10, "linspace_interp": False, "num_interp": 25, "v_min": 3,
                         "v_max": 20, "v_num": 10, "rpm_min": 100, "rpm_max": 3000, "rpm_num": 10, "pitch": 0.0,
                         "relaxation_factor": 0.3, "print_all": False, "print_out": False, "reynolds": 50000,
                         "fix_reynolds": False, "yaw_angle": 0}

        self.settings_to_name = {"propeller_mode": "Propeller mode", "print_out": "Print final iteration data",
                                 "tip_loss": "Prandtl tip loss", "hub_loss": "Prandtl hub loss",
                                 "new_tip_loss": "New tip loss", "new_hub_loss": "New hub loss",
                                 "cascade_correction": "Cascade correction", "max_iterations": "Maximum iterations",
                                 "convergence_limit": "Convergence criteria", "rho": "Air density [kg/m^3]",
                                 "method": "Calculation method", "v_min": "Min calc. wind speed [m/s]",
                                 "v_max": "Max calc. wind speed [m/s]", "v_num": "Number of wind speed points",
                                 "rpm_min": "Min calc. RPM [RPM]", "rpm_max": "Max calc. RPM [RPM]",
                                 "rpm_num": "Number of RPM points", "relaxation_factor": "Relaxation factor",
                                 "print_all": "Print every iteration [debug]",
                                 "num_interp": "Number of sections (interp)",
                                 "linspace_interp": "Custom number of sections",
                                 "rotational_augmentation_correction": "Rot. augmentation cor.",
                                 "rotational_augmentation_correction_method": "Rot. augmentation cor. method",
                                 "fix_reynolds": "Fix Reynolds", "reynolds": "Reynolds",
                                 "mach_number_correction": "Mach number correction", "pitch": "Pitch",
                                 "yaw_angle": "Yaw angle [°]", "skewed_wake_correction": "Skewed Wake Correction"}

        self.list_settings_for_updating_tsr = ["v_min", "v_max", "v_num", "rpm_min", "rpm_max", "rpm_num"]

        self.methods_to_names = METHODS_STRINGS

        self.name_to_methods = {v: k for k, v in self.methods_to_names.items()}
        self.name_to_settings = {v: k for k, v in self.settings_to_name.items()}

        self.grid = QGridLayout()
        self.setLayout(self.grid)

        self.left = QWidget()
        self.fbox = QFormLayout()
        self.left.setLayout(self.fbox)

        self.scroll_area = QScrollArea()
        self.scroll_widget = QWidget()
        self.scroll_widget_layout = QVBoxLayout()

        self.scroll_widget.setLayout(self.scroll_widget_layout)
        self.scroll_area.setWidget(self.left)
        self.scroll_area.setWidgetResizable(True)

        self.grid.addWidget(self.scroll_area, 1, 1)

        self.textEdit = QTextEdit()
        self.textEdit.setReadOnly(True)
        self.grid.addWidget(self.textEdit, 1, 2)

        self.buttonRun = QPushButton("Run")

        self.form_list = []
        self.validator = QtGui.QDoubleValidator()

        for key, value in self.settings.items():
            if key == "method":
                form = QComboBox()
                form.addItems([self.methods_to_names[k] for k, v in self.methods_to_names.items()])
                form.setCurrentIndex(7)
            elif key == "rotational_augmentation_correction_method":
                form = QComboBox()
                form.addItems(["1", "2", "3", "4", "5"])
            elif isinstance(value, bool):
                form = QCheckBox()
                form.setTristate(value)
            else:
                form = QLineEdit()
                form.setValidator(self.validator)
                form.textChanged.connect(self.check_state)
                form.textChanged.emit(form.text())
                form.insert(str(value))
                if key in self.list_settings_for_updating_tsr:
                    form.textChanged.connect(self.update_tsr_and_j)
            key_orig = key
            key = self.settings_to_name[key]
            self.fbox.addRow(key, form)
            self.form_list.append([key, form, key_orig])

        self.tsr_string = QLabel("0")
        self.J_string = QLabel("0")
        self.emptyLabel = QLabel(" ")
        self.buttonRun.clicked.connect(self.run)
        self.buttonClear = QPushButton("Clear screen")
        self.buttonClear.clicked.connect(self.clear)
        self.buttonEOF = QCheckBox()
        self.buttonEOF.setChecked(True)
        self.buttonEOFdescription = QLabel("Scroll to end of screen")
        self.buttonStop = QPushButton("Stop")
        self.buttonStop.clicked.connect(self.terminate)

        self.fbox.addRow(self.emptyLabel, self.buttonRun)
        self.fbox.addRow(self.buttonClear, self.buttonStop)
        self.fbox.addRow(self.buttonEOFdescription, self.buttonEOF)
        self.fbox.addRow("TSR:", self.tsr_string)
        self.fbox.addRow("J:", self.J_string)

    def update_tsr_and_j(self):
        try:
            s = self.get_settings()
            R = float(self.main.wind_turbine_properties.R.text())
            tsr_min = 2 * np.pi * float(s["rpm_min"]) * R / 60 / float(s["v_max"])
            tsr_max = 2 * np.pi * float(s["rpm_max"]) * R / 60 / float(s["v_min"])
            self.tsr_string.setText("%.2f - %.2f" % (tsr_min, tsr_max))
            J_min = float(s["v_min"]) / (float(s["rpm_max"]) / 60 * 2 * R)
            J_max = float(s["v_max"]) / (float(s["rpm_min"]) / 60 * 2 * R)
            self.J_string.setText("%.2f - %.2f" % (J_min, J_max))
        except:
            print("couldnt update tsr min/max or J min/max")

    def check_forms(self):
        out = ""
        for n, f, n_short in self.form_list:
            if isinstance(f, QLineEdit):
                state = self.validator.validate(f.text(), 0)[0]
                if state == QtGui.QValidator.Acceptable:
                    pass
                elif state == QtGui.QValidator.Intermediate:
                    out += "Form %s appears not to be valid.\n" % n
                else:
                    out += "Form %s is not of the valid type.\n" % n
        if out == "":
            return True
        return out

    def check_state(self, *args, **kwargs):
        sender = self.sender()
        validator = sender.validator()
        state = validator.validate(sender.text(), 0)[0]
        if state == QtGui.QValidator.Acceptable:
            color = "#edf5e1"  # green
        elif state == QtGui.QValidator.Intermediate:
            color = "#fff79a"  # yellow
        else:
            color = "#f6989d"  # red
        sender.setStyleSheet("QLineEdit { background-color: %s; color: #000000 }" % color)

    def get_settings(self):
        out_settings = {}
        for name, value, name_short in self.form_list:
            name = self.name_to_settings[name]
            if isinstance(value, QCheckBox):
                value = bool(value.checkState())
            if isinstance(value, QLineEdit):
                value = to_float(value.text())
            if isinstance(value, QComboBox):
                value = int(value.currentIndex())
            out_settings[name] = value
        return out_settings

    def set_settings(self, inp_dict):
        for name_long, item, name in self.form_list:
            if name in inp_dict:
                if isinstance(item, QComboBox):
                    _index = inp_dict[name]
                    if _index >= 0:
                        index = _index
                    else:
                        index = 0
                    item.setCurrentIndex(index)
                elif isinstance(item, QLineEdit):
                    item.setText(str(inp_dict[name]))
                elif isinstance(item, QCheckBox):
                    item.setChecked(inp_dict[name])

    def validate_inputs(self):
        check = self.check_forms()
        if check != True:
            msg = MyMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText("Input validation error.")
            msg.setDetailedText(check)
            msg.exec_()
        return check

    def run(self):
        self.clear()

        if not self.validate_inputs():
            return

        self.runner_input = self.main.get_input_params()

        if self.runner_input == None:
            print("No settings fetched... Exiting.")
            return

        self.main.emitter_add.connect(self.add_text)
        self.main.emitter_done.connect(self.terminate)

        if not self.main.running:
            self.main.set_process_running()
            self.main.getter.start()
            self.p = Process(target=calculate_power_3d, args=[self.runner_input])
            self.p.start()

    def add_text(self, string):
        if self.buttonEOF.checkState() == 2:
            self.textEdit.moveCursor(QtGui.QTextCursor.End)
        self.textEdit.insertPlainText(string)

    def clear(self):
        self.textEdit.clear()

    def terminate(self):
        self.p.terminate()

        self.main.set_process_stopped()
        self.main.getter.quit()

        self.main.emitter_add.disconnect()
        self.main.emitter_done.disconnect()

        self.show_results()

    def show_results(self):
        if len(self.main.return_results) > 0:
            results = self.main.return_results[-1]
            if "v" in results:
                if len(results["v"]) > 0:
                    inp_params = self.runner_input
                    self.r = ResultsWindow(None, self.main.screen_width, self.main.screen_width, results, inp_params, )
                else:
                    print("Not enough points to print results...")
            else:
                print("No results to print...")


class Optimization(QWidget):
    def __init__(self, parent=None):
        super(Optimization, self).__init__(parent)
        self.main = self.parent()

        self.validator = QtGui.QDoubleValidator()
        self.left = QWidget()
        self.textEdit = QTextEdit()
        self.textEdit.setReadOnly(True)

        self.form_list = []

        self.target_speed = QLineEdit()
        self.target_speed.setValidator(self.validator)
        self.target_speed.textChanged.connect(self.check_state)
        self.target_speed.textChanged.emit(self.target_speed.text())
        self._target_speed = QLabel("Target speed [m/s]")
        self.form_list.append([self._target_speed, self.target_speed])

        self.target_rpm = QLineEdit()
        self.target_rpm.setValidator(self.validator)
        self.target_rpm.textChanged.connect(self.check_state)
        self.target_rpm.textChanged.emit(self.target_rpm.text())
        self._target_rpm = QLabel("Target rpm [RPM]")
        self.form_list.append([self._target_rpm, self.target_rpm])

        self.min_bound = QLineEdit()
        self.min_bound.setValidator(self.validator)
        self.min_bound.textChanged.connect(self.check_state)
        self.min_bound.textChanged.emit(self.min_bound.text())
        self._min_bound = QLabel("Minimum boundary")
        self.form_list.append([self._min_bound, self.min_bound])

        self.max_bound = QLineEdit()
        self.max_bound.setValidator(self.validator)
        self.max_bound.textChanged.connect(self.check_state)
        self.max_bound.textChanged.emit(self.max_bound.text())
        self._max_bound = QLabel("Maximum boundary")
        self.form_list.append([self._max_bound, self.max_bound])

        self.mut_coeff = QLineEdit()
        self.mut_coeff.setValidator(self.validator)
        self.mut_coeff.textChanged.connect(self.check_state)
        self.mut_coeff.textChanged.emit(self.mut_coeff.text())
        self._mut_coeff = QLabel("Mutation coefficient")
        self.form_list.append([self._mut_coeff, self.mut_coeff])

        self.population = QLineEdit()
        self.population.setValidator(self.validator)
        self.population.textChanged.connect(self.check_state)
        self.population.textChanged.emit(self.population.text())
        self._population = QLabel("Population")
        self.form_list.append([self._population, self.population])

        self.num_iter = QLineEdit()
        self.num_iter.setValidator(self.validator)
        self.num_iter.textChanged.connect(self.check_state)
        self.num_iter.textChanged.emit(self.num_iter.text())
        self._num_iter = QLabel("Number of iterations")
        self.form_list.append([self._num_iter, self.num_iter])

        self._opt_variable = QLabel("Optimization variable")
        self.opt_variable = QComboBox()
        self.opt_variable.addItems(["Thrust (propeller)", "Torque (Turbine)", "max dQ min dT", "Best pitch"])
        self.form_list.append([self._opt_variable, self.opt_variable])

        self.pitch_optimization = QCheckBox()
        self.pitch_optimization.setChecked(False)
        self._pitch_optimization = QLabel("Pitch optimization")
        self.form_list.append([self._pitch_optimization, self.pitch_optimization])

        self.buttonOptimization = QPushButton("Run optimization")
        self.buttonOptimization.clicked.connect(self.run)

        self.buttonStop = QPushButton("Stop")
        self.buttonStop.clicked.connect(self.terminate)

        self.buttonClear = QPushButton("Clear screen")
        self.buttonClear.clicked.connect(self.clear)

        self.buttonEOF = QCheckBox()
        self.buttonEOF.setChecked(True)
        self.buttonEOFdescription = QLabel("Scroll to end of screen")

        self.grid = QGridLayout()
        self.setLayout(self.grid)
        self.grid.addWidget(self.left, 1, 1)
        self.grid.addWidget(self.textEdit, 1, 2)

        self.fbox = QFormLayout()
        self.left.setLayout(self.fbox)

        for a, b in self.form_list:
            self.fbox.addRow(a, b)

        self.fbox.addRow(self.buttonOptimization)
        self.fbox.addRow("", QLabel())

        self.fbox.addRow(self.buttonClear, self.buttonStop)
        self.fbox.addRow(self.buttonEOFdescription, self.buttonEOF)

        self.win = PyQtGraphWindow(self)
        self.manager_pyqtgraph = Manager()
        self.queue_pyqtgraph = self.manager_pyqtgraph.list()

    def check_forms_angles(self):
        out = ""
        _needed_vars = [[self._target_speed, self.target_speed], [self._target_rpm, self.target_rpm], ]
        for n, f in _needed_vars:
            if isinstance(f, QLineEdit):
                state = self.validator.validate(f.text(), 0)[0]
                if state == QtGui.QValidator.Acceptable:
                    pass
                elif state == QtGui.QValidator.Intermediate:
                    out += "Form %s appears not to be valid.\n" % n.text()
                else:
                    out += "Form %s is not of the valid type.\n" % n.text()
        if out == "":
            return True
        return out

    def check_state(self, *args, **kwargs):
        sender = self.sender()
        validator = sender.validator()
        state = validator.validate(sender.text(), 0)[0]
        if state == QtGui.QValidator.Acceptable:
            color = "#edf5e1"  # green
        elif state == QtGui.QValidator.Intermediate:
            color = "#fff79a"  # yellow
        else:
            color = "#f6989d"  # red
        sender.setStyleSheet("QLineEdit { background-color: %s; color: #000000 }" % color)

    def validate_inputs(self):
        check = self.check_forms_angles()
        check_analysis = self.main.analysis.check_forms()
        if check != True or check_analysis != True:
            if check == True:
                check = ""
            if check_analysis == True:
                check_analysis = ""
            check = check + check_analysis
            msg = MyMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText("Input validation error")
            msg.setDetailedText(check)
            msg.exec_()
            return False
        return True

    def clear(self):
        self.textEdit.clear()

    def add_text(self, string):
        self.textEdit.insertPlainText(string)
        if self.buttonEOF.checkState() == 2:
            self.textEdit.moveCursor(QtGui.QTextCursor.End)

    def run(self):
        self.clear()

        if not self.validate_inputs():
            return

        self.main.emitter_add.connect(self.add_text)
        self.main.emitter_done.connect(self.terminate)

        if not self.main.running:
            self.main.set_process_running()
            self.runner_input = self.main.get_input_params()
            self.main.getter.start()
            self.p = Process(target=optimize_angles_genetic, args=[self.runner_input, self.queue_pyqtgraph])
            self.p.start()
            self.win.show()
            self.win.start_update()

    def terminate(self):
        self.main.set_process_stopped()
        self.p.terminate()
        self.main.getter.quit()

        self.main.emitter_add.disconnect()
        self.main.emitter_done.disconnect()

        self.win.stop_update()

    def get_settings(self):
        out = {}
        out["target_rpm"] = self.target_rpm.text()
        out["target_speed"] = self.target_speed.text()
        out["pitch_optimization"] = bool(self.pitch_optimization.checkState())
        out["min_bound"] = self.min_bound.text()
        out["max_bound"] = self.max_bound.text()
        out["mut_coeff"] = self.mut_coeff.text()
        out["population"] = self.population.text()
        out["num_iter"] = self.num_iter.text()

        if int(self.opt_variable.currentIndex()) == 0:
            out["optimization_variable"] = "dT"
        elif int(self.opt_variable.currentIndex()) == 1:
            out["optimization_variable"] = "dQ"
        elif int(self.opt_variable.currentIndex()) == 2:
            out["optimization_variable"] = "max dT min dQ"
        elif int(self.opt_variable.currentIndex()) == 3:
            out["optimization_variable"] = "best_pitch"
        for k, v in out.items():
            if v == "":
                v = None
            elif v == None:
                pass
            else:
                if not k == "optimization_variable":
                    v = to_float(v)
            out[k] = v
        return out

    def set_settings(self, inp_dict):
        self.target_rpm.setText(str(inp_dict["target_rpm"]))
        self.target_speed.setText(str(inp_dict["target_speed"]))
        self.pitch_optimization.setChecked(inp_dict["pitch_optimization"])
        self.min_bound.setText(str(inp_dict["min_bound"]))
        self.max_bound.setText(str(inp_dict["max_bound"]))
        self.mut_coeff.setText(str(inp_dict["mut_coeff"]))
        self.population.setText(str(inp_dict["population"]))
        self.num_iter.setText(str(inp_dict["num_iter"]))


class ThreadGetter(QThread):
    def __init__(self, parent):
        super(ThreadGetter, self).__init__(parent)
        self.dataCollectionTimer = QtCore.QTimer()
        self.dataCollectionTimer.moveToThread(self)
        self.dataCollectionTimer.timeout.connect(self.updateInProc)

    def run(self):
        self.dataCollectionTimer.start(1)  # 0 causes freeze
        self.loop = QtCore.QEventLoop()
        self.loop.exec_()

    def updateInProc(self):
        if len(self.parent().return_print) > 0:
            t = self.parent().return_print.pop(0)
            self.parent().emitter_add.emit(str(t))
        if self.parent().end_of_file.value == True:
            self.parent().emitter_done.emit()


class DataCaptureThread(QThread):

    def __init__(self, parent, *args, **kwargs):
        QThread.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        self.dataCollectionTimer = QtCore.QTimer()
        self.dataCollectionTimer.moveToThread(self)
        self.dataCollectionTimer.timeout.connect(self.updateInProc)

    def run(self):
        self.dataCollectionTimer.start(1)  # 0 causes freeze
        self.loop = QtCore.QEventLoop()
        self.loop.exec_()

    def updateInProc(self):
        if len(self.parent.parent.queue_pyqtgraph) > 0:
            item = self.parent.parent.queue_pyqtgraph.pop()
            x = item[0]
            y = item[1]
            best_x = [item[2]]
            best_y = [item[3]]
            self.parent.curve.setData(x, y)
            self.parent.curve_red.setData(best_x, best_y)


class PyQtGraphWindow(QMainWindow):
    def __init__(self, parent):
        super(PyQtGraphWindow, self).__init__(parent)
        self.obj = pg.PlotWidget()
        self.setCentralWidget(self.obj)
        self.curve = self.obj.plot(pen=None, symbol='o', symbolPen=None, symbolSize=4, symbolBrush=('b'))
        self.curve_red = self.obj.plot(pen=None, symbol='o', symbolPen=None, symbolSize=5, symbolBrush=('r'))
        self.parent = parent
        self.thread = DataCaptureThread(self)

    def start_update(self):
        self.thread.start()

    def stop_update(self):
        self.thread.quit()


class MyMessageBox(QtGui.QMessageBox):
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


class TabWidget(QTabWidget):
    def __init__(self, parent=None):
        super(TabWidget, self).__init__(parent)
        self.tabs = []

    def add_tab(self, widget, tab_name):
        self.tabs.append([widget, tab_name])
        self.addTab(widget, tab_name)
        return

    def remove_tab(self, index):
        self.removeTab(index)
        del self.tabs[index]

    def remove_all_tabs(self):
        while len(self.tabs) > 0:
            self.remove_tab(0)

    def remove_current_tab(self):
        self.remove_tab(self.currentIndex())

    def rename_current_tab(self, string):
        self.setTabText(self.currentIndex(), string)
        self.tabs[self.currentIndex()][1] = string

    def current_tab_name(self):
        return self.tabText(self.currentIndex())


def main(quick_results=False):
    if sys.platform.startswith("win"):
        # On Windows calling this function is necessary for multiprocessing.
        multiprocessing.freeze_support()
        # To show icon in taskbar
        myappid = 'FEUM.BEM_Analiza.%s' % __version__  # arbitrary string
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    app = QApplication([])
    app_icon = QtGui.QIcon("icon_bem.ico")
    app.setWindowIcon(app_icon)
    app.setStyle("Fusion")
    if sys.platform.startswith("darwin") or True:
        # dark theme fix on OSX
        palette = QDarkPalette()
        palette.set_app(app)
        palette.set_stylesheet(app)
    screen = app.primaryScreen()
    size = screen.size()
    main = MainWindow(size.width(), size.height())
    if quick_results:
        main.analysis.run()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()

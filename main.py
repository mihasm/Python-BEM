# coding=utf-8
__author__ = "Miha Smrekar"
__credits__ = ["Miha Smrekar"]
__license__ = "GPL"
__version__ = "0.4.5"
__maintainer__ = "Miha Smrekar"
__email__ = "miha.smrekar9@gmail.com"
__status__ = "Development"

from functools import partial
import json
from multiprocessing import Process, Manager, Queue
import multiprocessing
import os
from pprint import pprint
import signal
import sys
import time
import traceback

from PyQt5 import QtCore, QtGui
import PyQt5
from PyQt5.QtCore import QThread, QTextStream, pyqtSignal, QProcess, QRect, Qt, QLocale
from PyQt5.QtGui import QPalette, QColor
from PyQt5.QtWidgets import (QComboBox, QMainWindow, QPushButton, QTextEdit, QWidget, QFormLayout, QLabel, QLineEdit,
                             QGridLayout, QCheckBox, QStyleFactory, QMessageBox, QAction, QFileDialog, QSlider,
                             QTabWidget, QApplication, QScrollArea, QVBoxLayout, QSplitter)
from calculation_runner import calculate_power_3d
import ctypes
from matplotlib import cm
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from montgomerie import Montgomerie
import mpl_toolkits.mplot3d as mp3d
from numpy import array
import numpy as np
from optimization import optimize_angles_genetic
from polars import get_x_y_from_link, scrape_data
from popravki import METHODS_STRINGS
import pyqtgraph as pg
from scipy.interpolate import interp1d
from turbine_data import SET_INIT
from utils import (generate_chord_lengths_betz, generate_chord_lengths_schmitz,
                   generate_twists_betz, generate_twists_schmitz,
                   get_centroid_coordinates,interpolate_geom, to_float, fltr, QDarkPalette,
                   create_folder, ErrorMessageBox, MyMessageBox, sort_data,
                   interpolate_geom, to_float, fltr, transpose, import_dat, import_nrel_dat, generate_dat, interp_at,
                   Table,create_macro_text)
from visualize import create_3d_blade
from xfoil import generate_polars_data


#from results import ResultsWindow

np.set_printoptions(threshold=sys.maxsize)

TITLE_STR = "BEM analiza v%s" % __version__

# determine if application is a script file or frozen exe
if getattr(sys, 'frozen', False):
    application_path = os.path.dirname(sys.executable)
elif __file__:
    application_path = os.path.dirname(__file__)



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
        self.tab_widget.add_tab(self.curve_manager, "Airfoil management", "Izbira profilov, cL/cD krivulje")

        self.wind_turbine_properties = WindTurbineProperties(self)
        self.tab_widget.add_tab(self.wind_turbine_properties, "Turbine info", "Nastavitev geometrije turbine/lopatic")

        self.analysis = Analysis(self)
        self.tab_widget.add_tab(self.analysis, "Analysis", "Nastavitev parametrov in BEM analiza")

        self.getter = ThreadGetter(self)

        self.optimization = Optimization(self)
        self.tab_widget.add_tab(self.optimization, "Optimization", "Optimizacija lopatice s pomočjo algoritma diferencialne evolucije")

        self.running = False
        self.manager = Manager()
        self.set_all_settings(SET_INIT)

        create_folder(os.path.join(application_path,"foils"))  # Used by XFoil

        self.set_process_stopped()

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
        name = QFileDialog.getSaveFileName(self, 'Save File',"", "BEM (*.bem)")[0]
        if name != "":
            
            d = self.get_all_settings()
            d_to_save = fltr(d, (float, int, list, str, bool, np.ndarray))
            json_d = json.dumps(d_to_save)

            file = open(name, 'w')
            file.write(json_d)
            file.close()

    def file_load(self):
        """
        Loads the wind turbine data from a file. Also clears the calculation text areas and sets the appropriate title.
        """
        file_path = QFileDialog.getOpenFileName(self, "Load File","", "BEM (*.bem)")[0]
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
        except:
            msg = ErrorMessageBox()
            properties = {}

        try:
            settings = self.analysis.get_settings()
        except:
            msg = ErrorMessageBox()
            settings = {}

        try:
            opt_settings = self.optimization.get_settings()
        except:
            msg = ErrorMessageBox()
            opt_settings = {}

        try:
            curve_manager_settings = self.curve_manager.get_settings()
        except:
            msg = ErrorMessageBox()
            curve_manager_settings = {}

        try:
            out = {**properties, **settings, **opt_settings, **curve_manager_settings}
            return out
        except:
            msg = ErrorMessageBox()
            #pprint(out)
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
        except:
            msg = ErrorMessageBox()

        try:
            self.wind_turbine_properties.set_settings(inp_dict)
        except:
            msg = ErrorMessageBox()

        try:
            self.optimization.set_settings(inp_dict)
        except:
            msg = ErrorMessageBox()

        try:
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
        self.Rhub.setToolTip("Radij pesta. Vpliv ima le pri Hub-Loss popravku.")

        _R = QLabel("Tip radius [m]")
        self.R = QLineEdit()
        self.R.setText("0.776")
        fbox.addRow(_R, self.R)
        self.R.setToolTip("Radij turbine. Vpliv ima na izračun moči ter pri Tip-Loss popravku.")

        _B = QLabel("Number of blades")
        self.B = QLineEdit()
        self.B.setText("5")
        fbox.addRow(_B, self.B)
        self.B.setToolTip("Število lopatic.")

        _blade_design = QLabel("Blade design")
        self.blade_design = QComboBox()
        self.blade_design.addItems(["Filled","Hollow","Spar"])
        fbox.addRow(_blade_design, self.blade_design)
        self.B.setToolTip("Pomembno za statični izračun.")

        _blade_thickness = QLabel("Blade thickness [m]")
        self.blade_thickness = QLineEdit()
        self.blade_thickness.setText("0.001")
        fbox.addRow(_blade_thickness, self.blade_thickness)
        self.blade_thickness.setToolTip("Debelina lopatice / Spar Cap-a.")

        fbox.addRow(QLabel("————— Scale and interpolation —————"))

        _geometry_scale = QLabel("Scale factor")
        self.geometry_scale = QLineEdit()
        self.geometry_scale.setText("1.0")
        fbox.addRow(_geometry_scale, self.geometry_scale)
        self.geometry_scale.setToolTip("Scale factor")

        _linspace_interp = QLabel("Interpolate geometry")
        self.linspace_interp = QCheckBox()
        fbox.addRow(_linspace_interp, self.linspace_interp)
        self.linspace_interp.setToolTip("Interpolate_geom")

        _num_interp = QLabel("Number of interpolation points")
        self.num_interp = QLineEdit()
        self.num_interp.setText("25")
        fbox.addRow(_num_interp, self.num_interp)
        self.num_interp.setToolTip("Number of interpolation points")

        fbox.addRow(QLabel("————— Generate geometry —————"))

        _num_gen_sections = QLabel("Number of gen. sections")
        self.num_gen_sections = QLineEdit()
        self.num_gen_sections.setText("10")
        fbox.addRow(_num_gen_sections, self.num_gen_sections)
        self.num_gen_sections.setToolTip("Število odsekov za generiranje.")

        _design_tsr = QLabel("Design TSR")
        self.design_tsr = QLineEdit()
        self.design_tsr.setText("7")
        fbox.addRow(_design_tsr, self.design_tsr)
        self.design_tsr.setToolTip("Željeni TSR za generiranje.")

        _design_aoa = QLabel("Design AoA [°]")
        self.design_aoa = QLineEdit()
        self.design_aoa.setText("7")
        fbox.addRow(_design_aoa, self.design_aoa)
        self.design_aoa.setToolTip("Željeni AoA za generiranje.")

        _design_cl = QLabel("Cl @ Design AoA")
        self.design_cl = QLineEdit()
        self.design_cl.setText("1.4")
        fbox.addRow(_design_cl, self.design_cl)
        self.design_cl.setToolTip("Koeficient vzgona pri željenem AoA.")

        _design_airfoil = QLabel("Airfoil")
        self.design_airfoil = QLineEdit()
        self.design_airfoil.setText("s826")
        fbox.addRow(_design_airfoil, self.design_airfoil)
        self.design_airfoil.setToolTip("Tekst za v stolpec airfoil.")

        _design_method = QLabel("Design method.")
        self.design_method = QComboBox()
        self.design_method.addItems(["Betz","Schmitz"])
        fbox.addRow(_design_method, self.design_method)
        self.design_method.setToolTip("Metoda dizajniranja.")

        _button_generate_geometry = QLabel("Generate geometry.")
        self.button_generate_geometry = QPushButton("Generate")
        fbox.addRow(_button_generate_geometry, self.button_generate_geometry)
        self.button_generate_geometry.clicked.connect(self.generate_geometry)
        self.button_generate_geometry.setToolTip("Generiraj geometrijo (povozi predhodno!).")

        fbox.addRow(QLabel("————— Export to Solidworks —————"))

        self.export_button = QPushButton("Export curve data")
        self.export_button.clicked.connect(self.export)
        fbox.addRow("Export:", self.export_button)
        self.export_button.setToolTip("Krivulje na vseh radijih lopatice se shranijo v posamezne datoteke. Solidworks makro se nato zgenerira v Python konzoli.")

        self.flip_turning_direction = QCheckBox()
        fbox.addRow("Flip turning direction", self.flip_turning_direction)

        self.propeller_geom = QCheckBox()
        fbox.addRow("Propeller", self.propeller_geom)

        fbox.addRow(QLabel("—————————————————————————"))

        _button_create_geometry_graph = QLabel("Create R,C,θ graph.")
        self.button_create_geometry_graph = QPushButton("Create R,C,θ graph.")
        fbox.addRow(_button_create_geometry_graph, self.button_create_geometry_graph)
        self.button_create_geometry_graph.clicked.connect(self.create_geometry_graph)
        self.button_create_geometry_graph.setToolTip("Izris grafa R,C,θ.")

        self.window = None

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
        out = {"Rhub": to_float(self.Rhub.text()), "R": to_float(self.R.text()), "B": int(self.B.text()),
              "turbine_name": self.name.text(), "geometry_scale":to_float(self.geometry_scale.text()),
              "linspace_interp":self.linspace_interp.isChecked(),
              "num_interp":int(self.num_interp.text()),
              "blade_thickness":to_float(self.blade_thickness.text()),
              "blade_design":self.blade_design.currentIndex()}
        geom_array = self.table_properties.get_values()
        r, c, theta, foils = [], [], [], []
        for row in geom_array:
            if row[0] != "" and row[1] != "" and row[2] != "":
                r.append(to_float(row[0]))
                c.append(to_float(row[1]))
                theta.append(to_float(row[2]))
                foils.append(row[3])
        out["r"] = array(r)
        out["c"] = array(c)
        out["theta"] = array(theta)
        out["foils"] = foils
        _r = out["r"]
        _c = out["c"]
        _theta = out["theta"]
        _foils = out["foils"]
        out["R"] = out["R"]*out["geometry_scale"]
        out["Rhub"] = out["Rhub"]*out["geometry_scale"]
        r, c, theta, foils, dr = interpolate_geom(_r, _c, _theta, _foils, out["num_interp"], out["linspace_interp"], out["geometry_scale"])
        out["r"], out["c"], out["theta"], out["foils"], out["dr"] = r, c, theta, foils, dr
        out["r_in"], out["c_in"], out["theta_in"], out["foils_in"] = _r, _c, _theta, _foils
        return out

    def create_geometry_graph(self):
        out = self.get_settings()
        self.gw = MatplotlibWindow()
        self.gw.setWindowTitle("r,c,θ graph")

        self.gw.ax = self.gw.figure.add_subplot(111)
        self.gw.ax.set_title("c(r) and θ(r)")

        self.gw.ax2 = self.gw.ax.twinx()
        self.gw.ax.plot(out["r"],out["c"],color="b")
        self.gw.ax2.plot(out["r"],out["theta"],color="r")

        self.gw.ax.set_xlabel("Radius r [m]")
        self.gw.ax.set_ylabel("Chord c [m]",color="tab:blue")
        self.gw.ax2.set_ylabel("Twist θ [°]",color="tab:red")
        self.gw.ax.tick_params(axis='y', labelcolor="tab:blue")
        self.gw.ax2.tick_params(axis='y', labelcolor="tab:red")



    def generate_geometry(self):
        array = []
        R = float(self.R.text())
        Rhub = float(self.Rhub.text())
        num_gen_sections = int(self.num_gen_sections.text())
        radiuses = np.linspace(Rhub,R,num_gen_sections)
        Cl_max = float(self.design_cl.text())
        B = float(self.B.text())
        TSR = float(self.design_tsr.text())
        method = self.design_method.currentIndex()
        airfoil = self.design_airfoil.text()
        design_aoa = float(self.design_aoa.text())
        
        if method == 0:
            chords = generate_chord_lengths_betz(radiuses=radiuses,R=R,Cl_max=Cl_max,B=B,TSR=TSR)
            thetas = generate_twists_betz(radiuses=radiuses,R=R,TSR=TSR,alpha_d=design_aoa)

        elif method == 1:
            chords = generate_chord_lengths_schmitz(radiuses=radiuses,R=R,Cl_max=Cl_max,B=B,TSR=TSR)
            thetas = generate_twists_schmitz(radiuses=radiuses,R=R,TSR=TSR,alpha_d=design_aoa)

        for r in range(num_gen_sections):
            array.append([round(radiuses[r],4), round(chords[r],4), round(thetas[r],4), airfoil])
        self.table_properties.createTable(array)

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
        if "blade_thickness" in dict_settings:
            t=str(dict_settings["blade_thickness"])
            self.blade_thickness.setText(t)
        if "blade_design" in dict_settings:
            self.blade_design.setCurrentIndex(dict_settings["blade_design"])

    def export(self):
        """
        Exports the wind turbine geometry data as spatial data (XYZ points), and creates a VB macro, which
        you can use to import the geometry into the Solidworks 3D modeller.
        :return:
        """

        if self.window != None:
            self.window.close()
        self.window = PrintoutWindow(self)
        self.window.setWindowTitle("Solidworks Export Macro")


        print("Getting settings...")
        settings_fetched = self.parent().parent().parent().get_all_settings()
        if settings_fetched == None:
            return
        data = create_3d_blade(settings_fetched, self.flip_turning_direction.isChecked(), self.propeller_geom.isChecked())
        self.w = MatplotlibWindow()
        self.w.setWindowTitle("Export 3D preview")
        self.w.ax = self.w.figure.add_subplot(111, projection="3d")
        self.w.ax.scatter(data["X"], data["Y"], data["Z"])
        X, Y, Z = array(data["X"]), array(data["Y"]), array(data["Z"])

        # Create cubic bounding box to simulate equal aspect ratio
        max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max()
        Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X.max()+X.min())
        Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y.max()+Y.min())
        Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Z.max()+Z.min())
        # Comment or uncomment following both lines to test the fake bounding box:
        for xb, yb, zb in zip(Xb, Yb, Zb):
           self.w.ax.plot([xb], [yb], [zb], 'w')
        # self.w.ax.set_aspect("equal")

        create_folder(os.path.join(application_path,"export"))
        folder_path = os.path.join(application_path,"export", SET_INIT["turbine_name"])
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
        self.tab_widget.setCurrentIndex(len(self.tab_widget.tabs)-1)

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


class Airfoils(QWidget):
    def __init__(self, airfoil_name, parent=None):
        super(Airfoils, self).__init__(parent)

        self.parent = parent

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
        self.grid.addWidget(self.buttonRefresh, 3, 3)
        self.buttonRefresh.setToolTip("Osvežitev grafa krivulje profila (na podlagi tabele na levi strani)")
        self.buttonRefresh.clicked.connect(self.refresh)
        
        self.link = QLineEdit("link (airfoiltools.com)")
        self.fbox.addRow(self.link)
        self.link.setToolTip("Tu lahko downloadamo krivulje cL/cD iz airfoiltools.com. Obliko profila moramo vnesti sami v tabelo (copy-paste iz excela).")

        self.button_generate_curves_link = QPushButton("Scrape curves from link")
        self.button_generate_curves_link.clicked.connect(self.generate_curves_link)
        self.fbox.addRow(self.button_generate_curves_link)
        self.button_generate_curves_link.setToolTip("S pomočjo linka do profila, dostopnega na strani airfoiltools.com, lahko program zdownloada vrednosti, zgenerirane z XFOIL, direktno iz spletne strani")


        self.button_curve_editor = QPushButton("Curve Editor")
        self.button_curve_editor.clicked.connect(self.open_curve_editor)
        self.fbox.addRow(self.button_curve_editor)
        self.button_curve_editor.setToolTip("Ročno spreminjanje/uvažanje/izvažanje cL/cD krivulj")


        self.button_open_viewer = QPushButton("Open Curve Extrapolator (Montgomerie)")
        self.button_open_viewer.clicked.connect(self.open_viewer)
        self.fbox.addRow(self.button_open_viewer)
        self.button_open_viewer.setToolTip("S pomočjo tega okna prilagajamo parametre ekstrapolacije z Montgomerie metodo za vsak dani Reynolds za cL in cD (alpha) krivulji")

        self.button_generate_curves_xfoil = QPushButton("Generate xfoil curves [debug]")
        self.button_generate_curves_xfoil.clicked.connect(self.generate_curves_xfoil)
        self.fbox.addRow(self.button_generate_curves_xfoil)

        self.button_visualize = QPushButton("Create curve visualization")
        self.button_visualize.clicked.connect(self.visualize)
        self.fbox.addRow(self.button_visualize)
        self.button_visualize.setToolTip("Prikaz 3D grafa ekstrapoliranih krivulj, za dodatno verifikacijo vhodnih podatkov v analizo")

        self.get_centroid_button = QPushButton("Calculate centroid")
        self.get_centroid_button.clicked.connect(self.calculate_centroid)
        self.fbox.addRow(self.get_centroid_button)
        self.get_centroid_button.setToolTip("S pomočjo tega gumba izračunamo sredino podanih točk v tabeli (težišče ploskve).")

        self.centroid_widget = QWidget()
        self.centroid_grid = QGridLayout()
        self.centroid_widget.setLayout(self.centroid_grid)
        self.centroid_label = QLabel("Centroid coordinates:")
        self.fbox.addRow(self.centroid_widget)
        self.centroid_label.setToolTip("Okoli te točke se zavrtijo točke pri 3D generaciji geometrije (Solidworks Makro). Na samo analizo nima vpliva.")

        self.centroid_x_edit = QLineEdit()
        self.centroid_y_edit = QLineEdit()

        self.centroid_grid.addWidget(self.centroid_label, 1, 1)
        self.centroid_grid.addWidget(self.centroid_x_edit, 1, 2)
        self.centroid_grid.addWidget(self.centroid_y_edit, 1, 3)

        self.button_import_dat_from_file = QPushButton("Import .dat")
        self.button_import_dat_from_file.clicked.connect(self.dat_importer)
        self.fbox.addRow(self.button_import_dat_from_file)
        self.button_import_dat_from_file.setToolTip("Import .dat file")

        self.button_import_nrel_dat_from_file = QPushButton("Import NREL .dat")
        self.button_import_nrel_dat_from_file.clicked.connect(self.nrel_dat_importer)
        self.fbox.addRow(self.button_import_nrel_dat_from_file)
        self.button_import_nrel_dat_from_file.setToolTip("Import NREL .dat file (cl and cd curves)")

        self.grid.setColumnStretch(1, 1)
        self.grid.setColumnStretch(2, 1)
        self.grid.setColumnStretch(3, 2)

        self.ncrit_selection = QComboBox()
        self._ncrit_selection = QLabel("Ncrit")
        self.fbox.addRow(self._ncrit_selection,self.ncrit_selection)
        self.ncrit_selection.setToolTip("Tu nastavimo N vrednost krivulj, ki jih želimo uporabiti. (oblika mejne plasti (e^N) -> XFOIL)")

        navodila = QLabel("Navodila za uporabo:\n"+
            "1. Na strani http://airfoiltools.com/\n"+
            "izberite poljubni aerodinamični profil.\n"+
            "2. Link vnesite zgoraj in pritisnite 'Scrape'.\n"+
            "Sedaj so cL/cD krivulje naložene v program.\n"+
            "(Ročno jih lahko spremenite s Curve Editor)\n"+
            "3. Sedaj je treba nastaviti koef. ekstrapolacije\n"+
            "z orodjem Curve Extrapolator (Montgomerie)\n"+
            "4. Končane krivulje lahko preverite\n"+
            "z uporabo orodja Create curve visualization,\n"+
            "kjer so prikazane v odvisnosti od Re")
        self.fbox.addRow(navodila)

        self.window = None
        self.curve_editor = CurveEditor(self)

    def dat_importer(self):
        """
        Loads the wind turbine data from a file. Also clears the calculation text areas and sets the appropriate title.
        """
        file_path = QFileDialog.getOpenFileName(self, "Import .dat file","", "dat (*.dat)")[0]
        if file_path != "":
            x,y = import_dat(file_path)
            self.set_x_y(x,y)
            self.refresh()

    def nrel_dat_importer(self):
        """
        Loads the wind turbine data from a file. Also clears the calculation text areas and sets the appropriate title.
        """
        file_path = QFileDialog.getOpenFileName(self, "Import .dat file","", "dat (*.dat)")[0]
        if file_path != "":
            data = import_nrel_dat(file_path)
            self.populate_curve_list(data)

            

    def visualize(self):
        print("Visualizing")
        data = self.curves.gather_curves()
        data = data[np.in1d(data[:,1],float(self.ncrit_selection.currentText()))] #current Ncrit

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
        self.w = MatplotlibWindow()
        self.w.setWindowTitle("Cl(alpha,Re)")
        self.w.ax = self.w.figure.add_subplot(111, projection="3d")
        p = self.w.ax.plot_trisurf(xi, yi, z_1, cmap=cm.coolwarm)
        self.w.ax.set_xlabel("Reynolds", fontsize=15, labelpad=20)
        self.w.ax.set_ylabel(r'$\alpha$ [°]', fontsize=15, labelpad=20)
        self.w.ax.set_zlabel("Cl", fontsize=15, labelpad=20)
        self.w.ax.xaxis.set_tick_params(labelsize=12)
        self.w.ax.yaxis.set_tick_params(labelsize=12)
        self.w.ax.zaxis.set_tick_params(labelsize=12)
        bar = self.w.figure.colorbar(p)
        bar.ax.set_xlabel('Cl', fontsize=15, labelpad=20)

        self.w2 = MatplotlibWindow()
        self.w2.setWindowTitle("Cd(alpha,Re)")
        self.w2.ax = self.w2.figure.add_subplot(111, projection="3d")
        p = self.w2.ax.plot_trisurf(xi, yi, z_2, cmap=cm.coolwarm)
        self.w2.ax.set_xlabel("Reynolds", fontsize=15, labelpad=20)
        self.w2.ax.set_ylabel(r'$\alpha$ [°]', fontsize=15, labelpad=20)
        self.w2.ax.set_zlabel("Cd", fontsize=15, labelpad=20)
        self.w2.ax.xaxis.set_tick_params(labelsize=12)
        self.w2.ax.yaxis.set_tick_params(labelsize=12)
        self.w2.ax.zaxis.set_tick_params(labelsize=12)
        bar2 = self.w2.figure.colorbar(p)
        bar2.ax.set_xlabel('Cd', fontsize=15, labelpad=20)


    def open_viewer(self):
        print("opening viewer")
        self.viewer.show()
        self.viewer.generate_views()

    def open_curve_editor(self):
        print("opening curve editor")
        self.curve_editor.show()
        self.curve_editor.load_curves()

    def generate_interp_functions(self):
        data = self.gather_curves()
        x, y = self.get_x_y()
        self.interp_function_cl, self.interp_function_cd = get_cl_cd_interpolation_function(data, x, y)

    def generate_curves_xfoil(self):
        print("Generating xfoil curves")
        x, y = self.get_x_y()
        generate_dat(self.airfoil_name,x,y)

        if self.window != None:
            self.window.close()
        self.window = PrintoutWindow(self)
        self.thread=XFoilThread(self)
        self.thread.set_params(self.airfoil_name + ".dat")
        self.thread.completeSignal.connect(self.xfoil_completion)
        self.thread.start()
        print("Done")

    def xfoil_completion(self,nothing_important):
        self.populate_curve_list(self.xfoil_generated_data)
        self.refresh()

    def generate_curves_link(self):
        print("Scraping from link...")
        if self.window != None:
            self.window.close()
        self.window = PrintoutWindow(self)
        self.thread=ScrapeThread(self)
        self.thread.set_params(self.link)
        self.thread.completeSignal.connect(self.generate_curves_link_completion)
        self.thread.start()
        
    def generate_curves_link_completion(self,nothing_important):
        self.table_dat.clear_table()
        self.populate_curve_list(self.scraping_generated_data[0])
        self.set_x_y(self.scraping_generated_data[1],self.scraping_generated_data[2])
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
        return 0.1 #default value

    def get_x_y(self):
        """
        Gets x and y values from table.
        """
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

    def set_x_y(self,x,y):
        """
        Sets x and y values from input to table.
        """
        array_dat = []
        if len(x) > 0 and len(y) > 0:
            for r in range(len(x)):
                array_dat.append([str(x[r]), str(y[r])])
            self.table_dat.createTable(array_dat)
            self.calculate_centroid()

    def calculate_centroid(self):
        foil_x, foil_y = self.get_x_y()
        x, y = get_centroid_coordinates(foil_x, foil_y)
        self.centroid_x_edit.setText(str(round(x,6)))
        self.centroid_y_edit.setText(str(round(y,6)))
        return x, y

    def get_ncrits(self):
        curves = self.curves.gather_curves()
        if len(curves)>0:
            ncrit_list = np.unique(curves[:,1])
            return ncrit_list

    def refresh_ncrits_combobox(self):
        self.ncrit_selection.clear()
        ncrits = self.get_ncrits()
        if ncrits is not None:
            self.ncrit_selection.addItems([str(n) for n in list(ncrits)])
        

    def get_settings(self):
        out = {}

        x, y = self.get_x_y()
        try:
            centroid_x = float(self.centroid_x_edit.text())
            centroid_y = float(self.centroid_y_edit.text())
        except:
            centroid_x,centroid_y = 0.0, 0.0

        try:
            ncrit_selected = float(self.ncrit_selection.currentText())
        except:
            ncrit_selected = 0.0
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
               "ncrit_selected": ncrit_selected}
        return out

    def set_settings(self, dict_settings):
        
        x,y = dict_settings["x"], dict_settings["y"]

        self.set_x_y(x,y)
        self.link.setText(dict_settings["link"])
        self.curves.load_curves(dict_settings["curves"])
        self.refresh()


class Curves:
    def __init__(self):
        self.curve_list = []

    def add(self, curve):
        self.curve_list.append(curve)

    def get_curves_sorted(self):
        return sorted(self.curve_list, key=lambda c: (c.ncrit, c.Re))

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
        for curve in self.get_curves_sorted():
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

    def get_curve(self,re_in,ncrit_in):
        out = []
        for curve in self.curve_list:
            re = curve.Re
            ncrit = curve.ncrit
            if ncrit == ncrit_in and re==re_in:
                out.append(curve)

        if len(out) == 0:
            return None
        if len(out) == 1:
            return out[0]
        if len(out) > 1:
            for c in out:
                print(c.Re,c.ncrit)
            raise Exception ("DataError: Multiple curves have same Reynolds and Ncrit...")

    def remove_curve(self,re_in,ncrit_in):
        out = []
        i = 0
        for curve in self.curve_list:
            re = curve.Re
            ncrit = curve.ncrit
            if ncrit == ncrit_in and re==re_in:
                j=i
                out.append(curve)
            i+=1

        if len(out) == 0:
            return None
        if len(out) > 1:
            for c in out:
                print(c.Re,c.ncrit)
            raise Exception ("DataError: Multiple curves have same Reynolds and Ncrit...")
        
        del self.curve_list[j]


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
        self.A = -5
        self.B = 5
        self.Am = 8
        self.Bm = 5
        self.m_CD90 = 1.5
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
            a = _alpha[i]
            try:
                cl = f_cl(a)
            except ValueError:
                cl = _cl[i]
            try:
                cd = f_cd(a)
            except ValueError:
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


class CurveEditor(QWidget):
    def __init__(self, parent=None):
        super(CurveEditor, self).__init__(None)
        self.resize(1600, 768)
        self.parent = parent

        self.grid = QGridLayout()
        self.setLayout(self.grid)

        #self.tab_widget = TabWidget(self)
        #self.grid.addWidget(self.tab_widget, 2, 0)

        self.validator = QtGui.QDoubleValidator()

        self.table = Table()
        self.table.createTable([[0,0,0]])
        self.table.set_labels(["alpha [°]", "cL (lift coeff.)", "cD (drag coeff.)"])
        self.grid.addWidget(self.table,2,0)

        self.upper_widget = QWidget()
        self.upper_layout = QGridLayout()
        self.upper_widget.setLayout(self.upper_layout)
        self.grid.addWidget(self.upper_widget, 1, 0)

        self.button_remove_curve = QPushButton("Remove curve")
        self.button_remove_curve.clicked.connect(self.remove_curve)
        self.upper_layout.addWidget(self.button_remove_curve,2,3)
        self.button_save_curve = QPushButton("Update curve")
        self.button_save_curve.clicked.connect(self.save_curve)
        self.upper_layout.addWidget(self.button_save_curve,1,3)
        self.button_add_curve = QPushButton("Add curve")
        self.button_add_curve.clicked.connect(self.add_curve)
        self.upper_layout.addWidget(self.button_add_curve,4,3)

        self.ncrit_edit = QLineEdit("Insert ncrit")
        self.upper_layout.addWidget(self.ncrit_edit,4,1)
        self.ncrit_edit.setValidator(self.validator)
        self.ncrit_edit.textChanged.connect(self.check_state)
        self.re_edit = QLineEdit("Insert reynolds")
        self.upper_layout.addWidget(self.re_edit,4,2)
        self.re_edit.setValidator(self.validator)
        self.re_edit.textChanged.connect(self.check_state)

        #self.picker_mach_label = QLabel("Mach number:")
        #self.picker_mach = QComboBox()
        #self.picker_mach.setEditable(True)
        #self.upper_layout.addWidget(self.picker_mach_label,2,1)
        #self.upper_layout.addWidget(self.picker_mach,3,1)

        self.picker_ncrit_label = QLabel("NCrit:")
        self.picker_ncrit = QComboBox()
        self.picker_ncrit.setEditable(False)
        #
        self.upper_layout.addWidget(self.picker_ncrit_label,1,1)
        self.upper_layout.addWidget(self.picker_ncrit,2,1)

        self.picker_reynolds_label = QLabel("Reynolds:")
        self.picker_reynolds = QComboBox()
        self.picker_reynolds.setEditable(False)
        #
        self.upper_layout.addWidget(self.picker_reynolds_label,1,2)
        self.upper_layout.addWidget(self.picker_reynolds,2,2)


        #self.picker_mach.lineEdit().returnPressed.connect(self.refresh_dropdowns)
        #self.sig1 = self.picker_reynolds.lineEdit().returnPressed.connect(self.save_curve)
        #self.sig2 = self.picker_ncrit.lineEdit().returnPressed.connect(self.save_curve)

        self.connect()

        #self.load_curves()

    def save_curve(self):

        out_chosen_values = self.get_chosen_values_from_dropdowns()
        re_chosen,ncrit_chosen = out_chosen_values

        data_from_table = self.table.get_values()
        alpha_table,cl_table,cd_table = transpose(data_from_table)
        alpha,cl,cd = [],[],[]

        for i in range(len(alpha_table)):
            if alpha_table[i] != "" and cl_table[i] != "" and cd_table[i] != "":
                alpha.append(float(alpha_table[i]))
                cl.append(float(cl_table[i]))
                cd.append(float(cd_table[i]))

        if self.parent.curves.get_curve(re_in=re_chosen,ncrit_in=ncrit_chosen) == None:
            msg = MyMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText("Curve with these Re - ncrit values does not exist yet. Did you mean to add a new curve?")
            msg.setDetailedText("Curve with Re %s and ncrit %s values does not exist yet. Did you mean to add a new curve?" % (re_chosen,ncrit_chosen))
            msg.exec_()
        else:
            self.current_curve = self.parent.curves.get_curve(re_in=re_chosen,ncrit_in=ncrit_chosen)
            self.current_curve.alpha = alpha
            self.current_curve.cl = cl
            self.current_curve.cd = cd
            self.load_curves()

    def add_curve(self):
        self.disconnect()

        re_chosen,ncrit_chosen = self.re_edit.text(), self.ncrit_edit.text()
        try:
            re_chosen,ncrit_chosen = float(re_chosen),float(ncrit_chosen)
        except:
            msg = MyMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText("Values of Re and ncrit for new curve do not seem to be valid.")
            msg.setDetailedText("Values of Re '%s' and ncrit '%s' could not be converted to float. Please double check the numbers." % (re_chosen,ncrit_chosen))
            msg.exec_()
            return

        data_from_table = self.table.get_values()
        alpha_table,cl_table,cd_table = transpose(data_from_table)
        alpha,cl,cd = [],[],[]

        for i in range(len(alpha_table)):
            if alpha_table[i] != "" and cl_table[i] != "" and cd_table[i] != "":
                alpha.append(float(alpha_table[i]))
                cl.append(float(cl_table[i]))
                cd.append(float(cd_table[i]))

        if self.parent.curves.get_curve(re_in=re_chosen,ncrit_in=ncrit_chosen) == None:
            print("item does not exist,creating new...")
            x,y = self.parent.get_x_y()
            self.current_curve = Curve()
            self.current_curve.create(x,y,re_chosen,ncrit_chosen,alpha,cl,cd)
            self.parent.curves.add(self.current_curve)
            print("self.parent.curves.curve_list",self.parent.curves.curve_list)
            self.load_curves()
        else:
            msg = MyMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText("Curve with this Re and ncrit already exists!")
            msg.setDetailedText("Curve with Re %s and ncrit %s values already exists. Did you mean to update the existing curve?" % (re_chosen,ncrit_chosen))
            msg.exec_()

        self.connect()



    def remove_curve(self):
        if len(self.parent.curves.curve_list) == 0:
            return
        out_chosen_values = self.get_chosen_values_from_dropdowns()
        re_chosen,ncrit_chosen = out_chosen_values
        self.parent.curves.remove_curve(re_chosen,ncrit_chosen)
        self.load_curves()


    def get_chosen_values_from_dropdowns(self):
        re_chosen = self.picker_reynolds.itemText(self.picker_reynolds.currentIndex())
        ncrit_chosen = self.picker_ncrit.itemText(self.picker_ncrit.currentIndex())
        return float(re_chosen),float(ncrit_chosen)

    def load_curves(self):
        try:
            re_last,ncrit_last = self.get_chosen_values_from_dropdowns()
        except ValueError:
            re_last,ncrit_last = None, None

        if len(self.parent.curves.curve_list) == 0:
            self.disconnect()

            self.picker_reynolds.clear()
            self.picker_ncrit.clear()

            self.connect()
            return

        self.disconnect()

        self.picker_reynolds.clear()
        self.picker_ncrit.clear()
        
        ncrit_list = []
        for curve in self.parent.curves.curve_list:      
            ncrit = curve.ncrit
            if not ncrit in ncrit_list:
                ncrit_list.append(ncrit)

        ncrit_list.sort()
        ncrit_list_str = [str(n) for n in ncrit_list]
        self.picker_ncrit.addItems(ncrit_list_str)

        if ncrit_last in ncrit_list:
            ncrit_index = ncrit_list.index(ncrit_last)
        else:
            ncrit_index = 0
            ncrit_last = ncrit_list[ncrit_index]

        self.picker_ncrit.setCurrentIndex(ncrit_index)

        re_list = []
        for curve in self.parent.curves.curve_list:
            re = curve.Re
            ncrit = curve.ncrit
            if not re in re_list:
                if ncrit == ncrit_last:
                    re_list.append(re)

        re_list.sort()
        re_list_str = [str(r) for r in re_list]
        self.picker_reynolds.addItems(re_list_str)

        if re_last in re_list:
            re_index = re_list.index(re_last)
        else:
            re_index = 0
            re_last = re_list[re_index]

        self.picker_reynolds.setCurrentIndex(re_index)

        self.connect()

        self.load_values_into_table()

    def disconnect(self):
        try:
            self.picker_reynolds.currentIndexChanged.disconnect()
            self.picker_ncrit.currentIndexChanged.disconnect()
        except TypeError:
            pass


    def connect(self):
        self.sig3 = self.picker_reynolds.currentIndexChanged.connect(self.load_curves)
        self.sig4 = self.picker_ncrit.currentIndexChanged.connect(self.load_curves)

    def load_values_into_table(self):
        out_chosen_values = self.get_chosen_values_from_dropdowns()
        if out_chosen_values == None:
            return
        re_chosen,ncrit_chosen = out_chosen_values
        chosen_curve = self.parent.curves.get_curve(re_in=re_chosen,ncrit_in=ncrit_chosen)
        if chosen_curve != None:
            self.current_curve = chosen_curve
            array = [self.current_curve.alpha,self.current_curve.cl,self.current_curve.cd]
            array = transpose(array)
            self.table.clear_table()
            self.table.createTable(array)

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

        for curve in self.parent.curves.get_curves_sorted():
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
        #print("A",self.curve.A,"B",self.curve.B,"A-",self.curve.Am,"B-",self.curve.Bm)
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

        self.forms_dict = {}

        self.settings = {"propeller_mode": False, "tip_loss": False, "hub_loss": False, "new_tip_loss": False,
                         "new_hub_loss": False, "cascade_correction": False, "skewed_wake_correction": False,
                         "rotational_augmentation_correction": False, "rotational_augmentation_correction_method": 1,
                         "mach_number_correction": False, "max_iterations": 100, "convergence_limit": 0.001,
                         "rho": 1.225, "method": 10, 
                         #"linspace_interp": False, "num_interp": 25,
                         "variable_selection": 0,
                         "constant_selection": 0,
                         "constant_speed":5,
                         "constant_rpm":1500,
                         "pitch": 0.0,
                         "v_min": 3, "v_max": 20, "v_num": 10,
                         "rpm_min": 100, "rpm_max": 3000, "rpm_num": 10,
                         "tsr_min":1, "tsr_max":10, "tsr_num":10,
                         "J_min":0.1, "J_max":1.5, "J_num":10,
                         "pitch_min":-15, "pitch_max":15, "pitch_num":10,
                         "relaxation_factor": 0.3, "print_all": False, "print_out": False, "reynolds": 50000,
                         "fix_reynolds": False, "yaw_angle": 0}

        self.settings_to_name = {"propeller_mode": "Propeller mode", "print_out": "Print final iteration data",
                                 "tip_loss": "Prandtl tip loss", "hub_loss": "Prandtl hub loss",
                                 "new_tip_loss": "New tip loss", "new_hub_loss": "New hub loss",
                                 "cascade_correction": "Cascade correction", "max_iterations": "Maximum iterations",
                                 "convergence_limit": "Convergence criteria", "rho": "Air density [kg/m^3]",
                                 "method": "Calculation method",
                                 "variable_selection" : "Variable parameter",
                                 "constant_selection" : "Constant variable",
                                 "constant_speed" : "Wind speed",
                                 "constant_rpm" : "RPM",
                                 "pitch": "Pitch",
                                 "v_min": "Min calc. wind speed [m/s]",
                                 "v_max": "Max calc. wind speed [m/s]",
                                 "v_num": "Number of wind speed points",
                                 "rpm_min": "Min calc. RPM [RPM]",
                                 "rpm_max": "Max calc. RPM [RPM]",
                                 "rpm_num": "Number of RPM points",
                                 "tsr_min":"Min TSR",
                                 "tsr_max":"Max TSR",
                                 "tsr_num":"Num TSR",
                                 "J_min":"Min J",
                                 "J_max":"Max J",
                                 "J_num":"Num J",
                                 "pitch_min":"Min pitch",
                                 "pitch_max":"Max pitch",
                                 "pitch_num":"Num pitch",
                                 "relaxation_factor": "Relaxation factor",
                                 "print_all": "Print every iteration [debug]",
                                 #"num_interp": "Number of sections (interp)",
                                 #"linspace_interp": "Custom number of sections",
                                 "rotational_augmentation_correction": "Rot. augmentation cor.",
                                 "rotational_augmentation_correction_method": "Rot. augmentation cor. method",
                                 "fix_reynolds": "Fix Reynolds", "reynolds": "Reynolds",
                                 "mach_number_correction": "Mach number correction",
                                 "yaw_angle": "Yaw angle [°]", "skewed_wake_correction": "Skewed Wake Correction"}

        self.settings_to_tooltip = {"propeller_mode": "Ta vrednost mora biti izbrana le v primeru, če preračunavamo propeler.",
                                 "print_out": "Izpis končnih vrednosti po konvergenci za vsak odsek",
                                 "tip_loss": "Popravek izgub pri vrhu lopatice (po Prandtlu)",
                                 "hub_loss": "Popravek izgub pri pestu (po Prandtlu)",
                                 "new_tip_loss": "Popravek izgub pri vrhu lopatice (po Speri)",
                                 "new_hub_loss": "Popravek izgub pri pestu (po Speri)",
                                 "cascade_correction": "Kaskadni popravki",
                                 "max_iterations": "Maksimalno število iteracij.",
                                 "convergence_limit": "Konvergenčni kriterij.",
                                 "rho": "Gostota zraka [kg/m^3]",
                                 "method": "Metoda za preračun. Privzeta je e) Aerodyn (Buhl).",
                                 "variable_selection": "Izbira spremenljivega parametra",
                                 "constant_selection":"Konstantna spremenljivka",
                                 "constant_speed":"Constant wind speed",
                                 "constant_rpm":"Constant RPM",
                                 "pitch": "Nastavni kot lopatice",
                                 "v_min": "Minimalna hitrost vetra [m/s]",
                                 "v_max": "Maksimalna hitrost vetra [m/s]",
                                 "v_num": "Število računskih točk (linearno razporejenih) od min hitrosti vetra do max hitrosti vetra",
                                 "rpm_min": "Minimalni vrtljaji/min [RPM]",
                                 "rpm_max": "Maksimalni vrtljaji/min [RPM]",
                                 "rpm_num": "Število računskih točk (linearno razporejenih) od min RPM do max RPM",
                                 "tsr_min":"Min TSR",
                                 "tsr_max":"Max TSR",
                                 "tsr_num":"Num TSR",
                                 "J_min":"Min J",
                                 "J_max":"Max J",
                                 "J_num":"Num J",
                                 "pitch_min":"Min pitch",
                                 "pitch_max":"Max pitch",
                                 "pitch_num":"Num pitch",
                                 "relaxation_factor": "Relaksacijski faktor. Privzeta vrednost: 0.3",
                                 "print_all": "Podroben izpis vrednosti po vsaki iteraciji (upočasni izračun)",
                                 #"num_interp": "Se uporabi samo v primeru, če izberemo 'Custom number of sections' opcijo",
                                 #"linspace_interp": "V primeru, da želimo preračunati več/manj odsekov lopatice po radiju (interpolacija geometrije)",
                                 "rotational_augmentation_correction": "Popravek rotacijske augmentacije",
                                 "rotational_augmentation_correction_method": "Izbira metode za popravek rotacijske augmentacije",
                                 "fix_reynolds": "Izračunaj vse odseke pri enem samem Re",
                                 "reynolds": "Se uporabi samo v primeru, če izberemo 'Fix Reynolds' opcijo",
                                 "mach_number_correction": "Popravek Mach števila (uporabno pri propelerjih)",
                                 "yaw_angle": "Kot vetra glede na smer osi rotorja [°]. Če je turbina obrnjena proti vetru, je 0°.",
                                 "skewed_wake_correction": "Popravek nagnjenega zračnega toka za turbino (Skewed wake)"}

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

        hideable_parameters = ["constant_selection","constant_speed","constant_rpm","pitch","v_min","v_max","v_num","rpm_min","rpm_max","rpm_num","tsr_min","tsr_max","tsr_num","J_min","J_max","J_num","pitch_min","pitch_max","pitch_num","reynolds"]

        for key, value in self.settings.items():
            if key == "method":
                form = QComboBox()
                form.addItems([self.methods_to_names[k] for k, v in self.methods_to_names.items()])
                form.setCurrentIndex(7)
            elif key == "rotational_augmentation_correction_method":
                form = QComboBox()
                form.addItems(["1", "2", "3", "4", "5"])
            elif key == "variable_selection":
                form = QComboBox()
                form.addItems(["RPM and v","TSR","J","pitch"])
                form.currentIndexChanged.connect(self.set_parameter_visibility)
            elif key == "constant_selection":
                form = QComboBox()
                form.addItems(["speed","rpm"])
                form.currentIndexChanged.connect(self.set_parameter_visibility)
            elif key == "fix_reynolds":
                form = QCheckBox()
                form.setTristate(False)
                form.stateChanged.connect(self.set_parameter_visibility)
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

            form.setToolTip(self.settings_to_tooltip[key])
            label = self.settings_to_name[key]
            
            self.form_list.append([label, form, key])
            
            if key in hideable_parameters:
                row = QWidget()
                row_layout = QFormLayout()
                row.setLayout(row_layout)
                row_layout.setContentsMargins(0, 0, 0, 0) # tight layout
                row_layout.addRow(label,form)
                self.fbox.addRow(row)
                self.forms_dict[key] = [form,label,row]
            else:
                self.fbox.addRow(label, form)
                self.forms_dict[key] = [form,label]

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

        self.tsr_string = QLabel("0")
        self.J_string = QLabel("0")
        self.fbox.addRow("TSR:", self.tsr_string)
        self.fbox.addRow("J:", self.J_string)

        self.set_parameter_visibility()

    def set_parameter_visibility(self):
        #print("Setting parameter visibility")
        #print(self.forms_dict)
        #print("form v_min",self.forms_dict["v_min"])
        #self.forms_dict["v_min"][2].hide()
        current_index = self.forms_dict["variable_selection"][0].currentIndex()
        if current_index == 0:
            # RPM / speed
            # rpm, speed need to be shown
            # constant_selection, TSR, J and pitch need to be hidden
            self.forms_dict["constant_selection"][2].hide()
            self.forms_dict["constant_speed"][2].hide()
            self.forms_dict["constant_rpm"][2].hide()
            self.forms_dict["pitch"][2].show()
            self.forms_dict["v_min"][2].show()
            self.forms_dict["v_max"][2].show()
            self.forms_dict["v_num"][2].show()
            self.forms_dict["rpm_min"][2].show()
            self.forms_dict["rpm_max"][2].show()
            self.forms_dict["rpm_num"][2].show()
            self.forms_dict["tsr_min"][2].hide()
            self.forms_dict["tsr_max"][2].hide()
            self.forms_dict["tsr_num"][2].hide()
            self.forms_dict["J_min"][2].hide()
            self.forms_dict["J_max"][2].hide()
            self.forms_dict["J_num"][2].hide()
            self.forms_dict["pitch_min"][2].hide()
            self.forms_dict["pitch_max"][2].hide()
            self.forms_dict["pitch_num"][2].hide()
            pass
        elif current_index == 1:
            # variable TSR
            # tsr, constant_selection, appropriate constant need to be shown
            self.forms_dict["constant_selection"][2].show()
            self.forms_dict["constant_speed"][2].show()
            self.forms_dict["constant_rpm"][2].show()
            self.forms_dict["pitch"][2].show()
            self.forms_dict["v_min"][2].hide()
            self.forms_dict["v_max"][2].hide()
            self.forms_dict["v_num"][2].hide()
            self.forms_dict["rpm_min"][2].hide()
            self.forms_dict["rpm_max"][2].hide()
            self.forms_dict["rpm_num"][2].hide()
            self.forms_dict["tsr_min"][2].show()
            self.forms_dict["tsr_max"][2].show()
            self.forms_dict["tsr_num"][2].show()
            self.forms_dict["J_min"][2].hide()
            self.forms_dict["J_max"][2].hide()
            self.forms_dict["J_num"][2].hide()
            self.forms_dict["pitch_min"][2].hide()
            self.forms_dict["pitch_max"][2].hide()
            self.forms_dict["pitch_num"][2].hide()
            pass
        elif current_index == 2:
            #variable J
            # J, constant_selection, appropriate constant need to be shown
            self.forms_dict["constant_selection"][2].show()
            self.forms_dict["constant_speed"][2].show()
            self.forms_dict["constant_rpm"][2].show()
            self.forms_dict["pitch"][2].show()
            self.forms_dict["v_min"][2].hide()
            self.forms_dict["v_max"][2].hide()
            self.forms_dict["v_num"][2].hide()
            self.forms_dict["rpm_min"][2].hide()
            self.forms_dict["rpm_max"][2].hide()
            self.forms_dict["rpm_num"][2].hide()
            self.forms_dict["tsr_min"][2].hide()
            self.forms_dict["tsr_max"][2].hide()
            self.forms_dict["tsr_num"][2].hide()
            self.forms_dict["J_min"][2].show()
            self.forms_dict["J_max"][2].show()
            self.forms_dict["J_num"][2].show()
            self.forms_dict["pitch_min"][2].hide()
            self.forms_dict["pitch_max"][2].hide()
            self.forms_dict["pitch_num"][2].hide()
            pass
        elif current_index == 3:
            # variable pitch
            # pitch, both constants need to be shown
            self.forms_dict["constant_selection"][2].hide()
            self.forms_dict["constant_speed"][2].show()
            self.forms_dict["constant_rpm"][2].show()
            self.forms_dict["pitch"][2].hide()
            self.forms_dict["v_min"][2].hide()
            self.forms_dict["v_max"][2].hide()
            self.forms_dict["v_num"][2].hide()
            self.forms_dict["rpm_min"][2].hide()
            self.forms_dict["rpm_max"][2].hide()
            self.forms_dict["rpm_num"][2].hide()
            self.forms_dict["tsr_min"][2].hide()
            self.forms_dict["tsr_max"][2].hide()
            self.forms_dict["tsr_num"][2].hide()
            self.forms_dict["J_min"][2].hide()
            self.forms_dict["J_max"][2].hide()
            self.forms_dict["J_num"][2].hide()
            self.forms_dict["pitch_min"][2].show()
            self.forms_dict["pitch_max"][2].show()
            self.forms_dict["pitch_num"][2].show()
            pass

        current_index_constant_value = self.forms_dict["constant_selection"][0].currentIndex()
        if current_index_constant_value == 0:
            self.forms_dict["constant_speed"][2].show()
            self.forms_dict["constant_rpm"][2].hide()
        elif current_index_constant_value == 1:
            self.forms_dict["constant_speed"][2].hide()
            self.forms_dict["constant_rpm"][2].show()

        current_state_fix_reynolds = self.forms_dict["fix_reynolds"][0].isChecked()
        if current_state_fix_reynolds == True:
            self.forms_dict["reynolds"][2].show()
        else:
            self.forms_dict["reynolds"][2].hide()

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
        #sender.setStyleSheet("QLineEdit { background-color: %s; color: #000000 }" % color)

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
            self.p = Process(target=calculate_power_3d, args=[self.runner_input,True])
            self.p.start()

    def add_text(self, string):
        if self.buttonEOF.checkState() == 2:
            self.textEdit.moveCursor(QtGui.QTextCursor.End)
        self.textEdit.insertPlainText(string)

    def clear(self):
        self.textEdit.clear()

    def terminate(self):
        try:
            self.p.terminate()
        except:
            pass

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
                    from results import ResultsWindow
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
        self.target_speed.setToolTip("Ciljna hitrost vetra. Pri tej hitrosti se bo izvajala BEM analiza.")

        self.target_rpm = QLineEdit()
        self.target_rpm.setValidator(self.validator)
        self.target_rpm.textChanged.connect(self.check_state)
        self.target_rpm.textChanged.emit(self.target_rpm.text())
        self._target_rpm = QLabel("Target rpm [RPM]")
        self.form_list.append([self._target_rpm, self.target_rpm])
        self.target_rpm.setToolTip("Ciljni vrtljaji rotorja turbine. Pri teh vrtljajih se bo izvajala BEM analiza.")

        self.min_bound = QLineEdit()
        self.min_bound.setValidator(self.validator)
        self.min_bound.textChanged.connect(self.check_state)
        self.min_bound.textChanged.emit(self.min_bound.text())
        self._min_bound = QLabel("Minimum boundary")
        self.form_list.append([self._min_bound, self.min_bound])
        self.min_bound.setToolTip("Minimalni kot (theta) lopatice.")

        self.max_bound = QLineEdit()
        self.max_bound.setValidator(self.validator)
        self.max_bound.textChanged.connect(self.check_state)
        self.max_bound.textChanged.emit(self.max_bound.text())
        self._max_bound = QLabel("Maximum boundary")
        self.form_list.append([self._max_bound, self.max_bound])
        self.max_bound.setToolTip("Maksimalni kot (theta) lopatice.")

        self.mut_coeff = QLineEdit()
        self.mut_coeff.setValidator(self.validator)
        self.mut_coeff.textChanged.connect(self.check_state)
        self.mut_coeff.textChanged.emit(self.mut_coeff.text())
        self._mut_coeff = QLabel("Mutation coefficient")
        self.form_list.append([self._mut_coeff, self.mut_coeff])
        self.mut_coeff.setToolTip("Mutacijski koeficient nastavlja jakost naključnih mutacij, ki se zgodijo pri vsaki novi generaciji (iteraciji).")

        self.population = QLineEdit()
        self.population.setValidator(self.validator)
        self.population.textChanged.connect(self.check_state)
        self.population.textChanged.emit(self.population.text())
        self._population = QLabel("Population")
        self.form_list.append([self._population, self.population])
        self.population.setToolTip("Število posameznikov v algoritmu diferencialne evolucije.")

        self.num_iter = QLineEdit()
        self.num_iter.setValidator(self.validator)
        self.num_iter.textChanged.connect(self.check_state)
        self.num_iter.textChanged.emit(self.num_iter.text())
        self._num_iter = QLabel("Number of iterations")
        self.form_list.append([self._num_iter, self.num_iter])
        self.num_iter.setToolTip("Število generacij (iteracij). Konvergenčni kriterij je pri tovrstnih algoritmih težko določljiv, zato izberemo fiksno vrednost.")

        self._opt_variable = QLabel("Optimization variable")
        self.opt_variable = QComboBox()
        self.opt_variable.addItems(["max(dQ) (torque->wind turbine)", "max(dT) (thrust->propeller)", "max(weight_dq*dQ-weight_dt*dT)", "max(weight_dt*dT-weight_dq*dQ)"])
        self.opt_variable.setCurrentIndex(0)
        self.form_list.append([self._opt_variable, self.opt_variable])
        self.opt_variable.setToolTip("Optimizacija naj poteka za to izbrano spremenljivko. V primeru vetrne turbine max(dQ).")

        self.weight_dq = QLineEdit()
        self.weight_dq.setValidator(self.validator)
        self.weight_dq.textChanged.connect(self.check_state)
        self.weight_dq.textChanged.emit(self.weight_dq.text())
        self._weight_dq = QLabel("weight_dq")
        self.form_list.append([self._weight_dq, self.weight_dq])
        self.weight_dq.setToolTip("Utež dq v primeru dvojne optimizacije.")

        self.weight_dt = QLineEdit()
        self.weight_dt.setValidator(self.validator)
        self.weight_dt.textChanged.connect(self.check_state)
        self.weight_dt.textChanged.emit(self.weight_dt.text())
        self._weight_dt = QLabel("weight_dt")
        self.form_list.append([self._weight_dt, self.weight_dt])
        self.weight_dt.setToolTip("Utež dt v primeru dvojne optimizacije.")

        self.pitch_optimization = QCheckBox()
        self.pitch_optimization.setChecked(False)
        self._pitch_optimization = QLabel("Pitch optimization")
        self.form_list.append([self._pitch_optimization, self.pitch_optimization])
        self.pitch_optimization.setToolTip("To možnost izberemo samo, kadar nas zanima optimalni nastavni kot lopatice.")

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
        self.win.setWindowTitle("Live Optimization Visualizer")
        self.manager_pyqtgraph = Manager()
        self.queue_pyqtgraph = self.manager_pyqtgraph.list()
        self.queue_pyqtgraph.append([[0],[0],0,0])

        self.tsr_string = QLabel("0")
        self.J_string = QLabel("0")
        self.fbox.addRow("TSR:", self.tsr_string)
        self.fbox.addRow("J:", self.J_string)

    def update_tsr_and_j(self):
        try:
            s = self.get_settings()
            R = float(self.main.wind_turbine_properties.R.text())
            tsr = 2 * np.pi * float(s["target_rpm"]) * R / 60 / float(s["target_speed"])
            self.tsr_string.setText("%.2f" % (tsr))
            J = float(s["target_speed"]) / (float(s["target_rpm"]) / 60 * 2 * R)
            self.J_string.setText("%.2f" % (J))
        except:
            pass

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
        self.update_tsr_and_j()

        sender = self.sender()
        validator = sender.validator()
        state = validator.validate(sender.text(), 0)[0]
        if state == QtGui.QValidator.Acceptable:
            color = "#edf5e1"  # green
        elif state == QtGui.QValidator.Intermediate:
            color = "#fff79a"  # yellow
        else:
            color = "#f6989d"  # red
        #sender.setStyleSheet("QLineEdit { background-color: %s; color: #000000 }" % color)

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
        try:
            self.p.terminate()
        except:
            pass
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
        out["weight_dt"] = self.weight_dt.text()
        out["weight_dq"] = self.weight_dq.text()

        if int(self.opt_variable.currentIndex()) == 0:
            out["optimization_variable"] = "dQ"
        elif int(self.opt_variable.currentIndex()) == 1:
            out["optimization_variable"] = "dT"
        elif int(self.opt_variable.currentIndex()) == 2:
            out["optimization_variable"] = "dQ-dT"
        elif int(self.opt_variable.currentIndex()) == 3:
            out["optimization_variable"] = "dT-dQ"

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
        self.weight_dt.setText(str(inp_dict["weight_dt"]))
        self.weight_dq.setText(str(inp_dict["weight_dq"]))


class ThreadGetter(QThread):
    def __init__(self, parent):
        super(ThreadGetter, self).__init__(parent)
        self.dataCollectionTimer = QtCore.QTimer()
        self.dataCollectionTimer.moveToThread(self)
        self.dataCollectionTimer.timeout.connect(self.updateInProc)

    def run(self):
        self.dataCollectionTimer.start(2)  # 0 causes freeze
        self.loop = QtCore.QEventLoop()
        self.loop.exec_()

    def updateInProc(self):
        if len(self.parent().return_print) > 0:
            t = self.parent().return_print.pop(0)
            self.parent().emitter_add.emit(str(t))
        if self.parent().end_of_file.value == True and len(self.parent().return_print) == 0:
            self.parent().emitter_done.emit()


class DataCaptureThread(QThread):

    def __init__(self, parent, *args, **kwargs):
        QThread.__init__(self, parent, *args, **kwargs)
        self.dataCollectionTimer = QtCore.QTimer()
        self.dataCollectionTimer.moveToThread(self)
        self.dataCollectionTimer.timeout.connect(self.updateInProc)

    def run(self):
        self.dataCollectionTimer.start(50)  # 0 causes freeze
        self.loop = QtCore.QEventLoop()
        self.loop.exec_()

    def updateInProc(self):
        if len(self.parent().parent.queue_pyqtgraph) > 0:
            item = self.parent().parent.queue_pyqtgraph[0]
            x = item[0]
            y = item[1]
            best_x = [item[2]]
            best_y = [item[3]]
            self.parent().curve.setData(x, y)
            self.parent().curve_red.setData(best_x, best_y)


class XFoilThread(QThread):

    progressSignal = QtCore.Signal(int)
    completeSignal = QtCore.Signal(str)

    def __init__(self, parent, *args, **kwargs):
        QThread.__init__(self, parent)
        self.parent = parent

    def set_params(self,dat_path):
        self.dat_path = dat_path

    def run(self):
        out = generate_polars_data(self.dat_path)
        self.parent.xfoil_generated_data = out
        self.completeSignal.emit("Done")

class ScrapeThread(QThread):

    progressSignal = QtCore.Signal(int)
    completeSignal = QtCore.Signal(str)

    def __init__(self, parent, *args, **kwargs):
        QThread.__init__(self, parent)
        self.parent = parent

    def set_params(self,link):
        self.link = link

    def run(self):
        data = scrape_data(self.link.text())
        x,y=get_x_y_from_link(self.link.text())
        out = [data,x,y]
        self.parent.scraping_generated_data = out
        self.completeSignal.emit("Done")

        


class PyQtGraphWindow(QMainWindow):
    def __init__(self, parent):
        super(PyQtGraphWindow, self).__init__(parent)
        self.obj = pg.PlotWidget()
        self.setCentralWidget(self.obj)
        self.curve = self.obj.plot(pen=None, symbol='o', symbolPen=None, symbolSize=4, symbolBrush=('g'))
        self.curve_red = self.obj.plot(pen=None, symbol='o', symbolPen=None, symbolSize=5, symbolBrush=('r'))
        self.obj.setLabel("left", "Optimization variable")
        self.obj.setLabel("bottom", "Theta [°]")
        self.parent = parent
        self.thread = DataCaptureThread(self)

    def start_update(self):
        self.thread.start()

    def stop_update(self):
        self.thread.quit()


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
        self.close()

    def closeEvent(self, event):
        self.emitter.disconnect()
        event.accept()


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

    def closeEvent(self, event):
        self.figure.clear()
        plt.close(self.figure)
        event.accept()  # let the window close


class PrintoutWindow(QMainWindow):
    def __init__(self, parent):
        super(PrintoutWindow, self).__init__(parent)
        self.setWindowTitle("Progress")
        self.setGeometry(50, 50, 500, 300)
        self.parent = parent
        sys.stdout = Stream(newText=self.onUpdateText)
        sys.stderr = Stream(newText=self.onUpdateText)
        self.process  = QtGui.QTextEdit()
        self.setCentralWidget(self.process)
        self.show()
 
    def onUpdateText(self, text):
        cursor = self.process.textCursor()
        cursor.movePosition(QtGui.QTextCursor.End)
        cursor.insertText(text)
        self.process.setTextCursor(cursor)
        self.process.ensureCursorVisible()

    def closeEvent(self, event):
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        QMainWindow.closeEvent(self, event)


class Stream(QtCore.QObject):
    newText = QtCore.pyqtSignal(str)

    def write(self, text):
        self.newText.emit(str(text))


class TabWidget(QTabWidget):
    def __init__(self, parent=None):
        super(TabWidget, self).__init__(parent)
        self.tabs = []

    def add_tab(self, widget, tab_name, tooltip=None):
        self.tabs.append([widget, tab_name])
        self.addTab(widget, tab_name)
        if tooltip != None:
            self.setTabToolTip(len(self.tabs)-1,tooltip)


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
    QLocale.setDefault(QLocale(QLocale.English)) # da je pika decimalno mesto
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
    main.setWindowIcon(app_icon)
    if quick_results:
        main.analysis.run()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()

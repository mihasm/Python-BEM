__author__ = "Miha Smrekar"
__credits__ = ["Miha Smrekar"]
__license__ = "GPL"
__version__ = "0.3.0"
__maintainer__ = "Miha Smrekar"
__email__ = "miha.smrekar9@gmail.com"
__status__ = "Development"

import sys

import numpy
from PyQt5 import QtWidgets
from PyQt5.QtCore import QThread, QTextStream, pyqtSignal, QProcess, QRect
from PyQt5 import QtCore, QtGui
from PyQt5.QtWidgets import (
    QComboBox,
    QMainWindow,
    QPushButton,
    QTextEdit,
    QWidget,
    QFormLayout,
    QLabel,
    QLineEdit,
    QGridLayout,
    QCheckBox,
    QStyleFactory,
    QMessageBox,
    QAction,
    QFileDialog,
)
from PyQt5.QtGui import QPalette, QColor
from numpy import array
from scipy import interpolate

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt

from cp_curve import calculate_power_3d
from results import ResultsWindow
from table import Table
import time
from turbine_data import SET_INIT
from optimisation import Optimizer
from utils import interpolate_geom, to_float, fltr

from multiprocessing import Process, Manager
import multiprocessing
import json

TITLE_STR = "BEM analiza v%s" % __version__

METHODS_STRINGS = {"0":"Original",
                   "1":"b) Spera",
                   "2":"Wiley: Strip theory, incl. wake rot.",
                   "3":"Grant Ingram (without Ct corr.)",
                   "4":"f) Glauert empirical",
                   "5":"Propx",
                   "6":"e) Aerodyn (Buhl)",
                   "7":"QBlade (Buhl)",
                   "8":"d) Shen",
                   "9":"a) Glauert",
                   "10":"Wilson and Walker",
                   "11":"Classical brake state model",
                   "12":"Advanced brake state model",
                   "13":"c) Modified ABS model"}

class MainWindow(QMainWindow):
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
        getSettings=QAction("Get settings",self)
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

        self.curve_manager = CurveManager(self)
        #self.curve_manager.add_foil("s826")
        self.tab_widget.add_tab(self.curve_manager,"Curve Manager")

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

        self.show()

    def set_title(self):
        s = self.wind_turbine_properties.name.text()
        if s == "":
            self.setWindowTitle(TITLE_STR)
        else:
            self.setWindowTitle(TITLE_STR+" - "+s)

    def file_save(self):
        name = QFileDialog.getSaveFileName(self, 'Save File')[0]
        if name != "":
            file = open(name, 'w')
            d = self.get_all_settings()
            d_to_save = fltr(d,(float, int, list, str, bool,numpy.ndarray))
            json_d = json.dumps(d_to_save)
            file.write(json_d)
            file.close()

    def file_load(self):
        file_path = QFileDialog.getOpenFileName(self, "Load File")[0]
        if file_path != "":
            with open(file_path, "r") as fp:
                data = json.load(fp)
            self.set_all_settings(data)
        self.analysis.clear()
        self.optimization.clear()
        self.set_title()

    def get_all_settings(self):
        valid_foils = list(self.curve_manager.get_settings()["curves"].keys()) + ["transition","Transition"]
        try:
            properties = self.wind_turbine_properties.get_settings()
            #curves = self.curves.get_settings()
            settings = self.analysis.get_settings()
            opt_settings = self.optimization.get_settings()
            curve_manager_settings = self.curve_manager.get_settings()
            out = {**properties, **settings, **opt_settings, **curve_manager_settings}#, **curves}
            _r = out["r"]
            _c = out["c"]
            _theta = out["theta"]
            _foils = out["foils"]
            r, c, theta, foils, dr = interpolate_geom(
                _r, _c, _theta, _foils, out["num_interp"], out["linspace_interp"]
            )
            out["r"], out["c"], out["theta"], out["foils"], out["dr"] = r, c, theta, foils, dr
            out["r_in"], out["c_in"], out["theta_in"], out["foils_in"] = _r, _c, _theta, _foils
            print(out)
            return out
        except Exception as e:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText("Error while getting settings")
            msg.setDetailedText(str(e))
            msg.exec_()
            print(e)
            return None

    def set_all_settings(self, inp_dict):
        try:
            self.analysis.set_settings(inp_dict)
        except:
            print("Error setting analysis settings!")
        try:
            self.optimization.set_settings(inp_dict)
        except:
            print("Error setting optimization settings!")
        try:
            self.wind_turbine_properties.set_settings(inp_dict)
        except:
            print("Error setting wind turbine properties settings!")
        try:
            self.curve_manager.set_settings(inp_dict)
        except:
            #print(e)
            print("Error setting curve manager settings!")
            #raise

    def get_input_params(self):
        settings = self.get_all_settings()
        if settings == None:
            return None
        self.return_print = self.manager.list([])
        self.return_results = self.manager.list([])
        self.end_of_file = False
        inp_params = {
            **settings,
            "return_print": self.return_print,
            "return_results": self.return_results,
            "EOF": self.end_of_file
        }
        return inp_params

    def set_buttons_running(self):
        self.analysis.buttonRun.setEnabled(False)
        self.optimization.buttonAngles.setEnabled(False)
        self.optimization.buttonPitch.setEnabled(False)
        self.analysis.buttonStop.setEnabled(True)
        self.optimization.buttonStop.setEnabled(True)

    def set_buttons_await(self):
        self.analysis.buttonRun.setEnabled(True)
        self.optimization.buttonAngles.setEnabled(True)
        self.optimization.buttonPitch.setEnabled(True)
        self.analysis.buttonStop.setEnabled(False)
        self.optimization.buttonStop.setEnabled(False)


class WindTurbineProperties(QWidget):
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

    def get_settings(self):
        out_properties = {
            "Rhub": to_float(self.Rhub.text()),
            "R": to_float(self.R.text()),
            "B": int(self.B.text()),
            "turbine_name":self.name.text(),
        }
        geom_array = self.table_properties.get_values()
        r, c, theta, foils = [], [], [], []
        for row in geom_array:
            if row[0] != "" and row[1] != "" and row[2] != "":
                r.append(to_float(row[0]))
                c.append(to_float(row[1]))
                theta.append(to_float(row[2]))
                foils.append(row[3])
        out_properties["r"] = numpy.array(r)
        out_properties["c"] = numpy.array(c)
        out_properties["theta"] = numpy.array(theta)
        out_properties["foils"] = foils
        return out_properties

    def set_settings(self, dict_settings):
        if "Rhub" in dict_settings:
            t = str(dict_settings["Rhub"])
            self.Rhub.setText(t)
        if "R" in dict_settings:
            t = str(dict_settings["R"])
            self.R.setText(t)
        if "B" in dict_settings:
            t = str(dict_settings["B"])
            self.B.setText(t)
        if "r" in dict_settings and "c" in dict_settings and "theta" in dict_settings and "foils" in dict_settings:
            _array = []
            for r in range(len(dict_settings["r"])):
                _r = dict_settings["r"][r]
                _c = dict_settings["c"][r]
                _theta = dict_settings["theta"][r]
                _f = dict_settings["foils"][r]
                _array.append(
                    [
                        _r,
                        _c,
                        _theta,
                        _f
                    ]
                )
            self.table_properties.createTable(_array)
        if "turbine_name" in dict_settings:
            t = str(dict_settings["turbine_name"])
            self.name.setText(t)
        else:
            self.name.setText("")


class CurveManager(QWidget):
    #popup_close = pyqtSignal(str)
    emitter = pyqtSignal(str)
    def __init__(self, parent=None):
        super(CurveManager, self).__init__(parent)
        
        self.main = self.parent()

        self.grid = QGridLayout()
        self.setLayout(self.grid)

        self.tab_widget = TabWidget(self)
        self.grid.addWidget(self.tab_widget,2,0)

        self.upper_widget = QWidget()
        self.upper_layout = QGridLayout()
        self.upper_widget.setLayout(self.upper_layout)
        self.grid.addWidget(self.upper_widget,1,0)

        self.button_add_foil = QPushButton("Add airfoil")
        self.button_add_foil.clicked.connect(self.add_foil_popup)
        self.button_remove_foil = QPushButton("Remove foil")
        self.button_remove_foil.clicked.connect(self.tab_widget.remove_current_tab)
        self.button_rename_foil = QPushButton("Rename foil")
        self.button_rename_foil.clicked.connect(self.rename_foil_popup)
        self.button_get_settings = QPushButton("get settings")
        self.button_get_settings.clicked.connect(self.get_settings)

        self.upper_layout.addWidget(self.button_add_foil,0,1)
        self.upper_layout.addWidget(self.button_remove_foil,0,2)
        self.upper_layout.addWidget(self.button_rename_foil,0,3)
        self.upper_layout.addWidget(self.button_get_settings,0,4)

    def add_foil_popup(self):
        self.emitter.connect(self.add_foil)
        self.p = PopupText(self,"foil name","airfoil_name",self.emitter,"Add foil")
        self.p.setGeometry(QRect(100, 100, 400, 200))
        self.p.show()

    def add_foil(self,string):
        c = Curves(self)
        self.tab_widget.add_tab(c,string)

    def rename_foil_popup(self):
        self.emitter.connect(self.rename_foil)
        self.p = PopupText(self,"foil name",self.tab_widget.current_tab_name(),self.emitter,"Rename foil")
        self.p.setGeometry(QRect(100, 100, 400, 200))
        self.p.show()
    
    def rename_foil(self,string):
        self.tab_widget.rename_current_tab(string)
        #self.tab_widget.tabs

    def get_settings(self):
        out = {}
        i=0

        #TODO Dont rely on name being set correctly in n!
        for w,n in self.tab_widget.tabs:
            out[n] = w.get_settings()
            i+=1

        return {"curves":out}

    def set_settings(self,dict_settings):
        self.tab_widget.remove_all_tabs()
        if "curves" in dict_settings:
            if len(dict_settings["curves"]) > 0:
                for c_name,c_dict in dict_settings["curves"].items():
                    curve_widget = Curves()
                    curve_widget.set_settings(c_dict)
                    self.tab_widget.add_tab(curve_widget,c_name)


class PopupText(QWidget):
    def __init__(self, parent=None,message="message",default_str="", emitter=None, title="Text popup"):
        QWidget.__init__(self)

        #self.setTitle(title)
        self.emitter = emitter

        self.layout = QGridLayout()
        self.setLayout(self.layout)

        self.message = QLabel(message)
        #
        self.layout.addWidget(self.message,0,0)

        self.inp = QLineEdit()
        self.inp.setText(default_str)
        self.layout.addWidget(self.inp,1,0)

        self.button = QPushButton("OK")
        self.button.clicked.connect(self.send_signal)
        self.layout.addWidget(self.button,2,0)

    def send_signal(self):
        self.emitter.emit(self.inp.text())
        self.emitter.disconnect()
        self.close()


class Curves(QWidget):
    def __init__(self, parent=None):
        super(Curves, self).__init__(parent)
        grid = QGridLayout()
        self.setLayout(grid)

        self.table_dat = Table()
        self.table_dat.createEmpty(2, 50)
        self.table_dat.set_labels(["x", "y"])
        grid.addWidget(self.table_dat, 1, 1)

        self.plt = plt.figure(figsize=(10,5))
        self.canvas = FigureCanvas(self.plt)
        toolbar = NavigationToolbar(self.canvas,self)
        grid.addWidget(self.canvas,1,2)
        self.ax = self.plt.add_subplot(111)
        grid.addWidget(toolbar,2,2)
        #self.ax.plot([1,2,3],[4,5,6])
        #self.ax = self.plt.add_subplot()

        #self.table_cd = Table()
        #self.table_cd.createEmpty(2, 50)
        #self.table_cd.set_labels(["x", "y"])
        #grid.addWidget(self.table_cd, 1, 2)

        #self.max_thickness_str = QLabel("Maximum thickness [ratio of chord length]:")
        #self.max_thickness = QLineEdit()
        self.buttonRefresh = QPushButton("Refresh curve")
        #grid.addWidget(self.max_thickness_str,2,1)
        grid.addWidget(self.buttonRefresh,2,1)
        self.buttonRefresh.clicked.connect(self.draw_curve)

        self.table_dat.tableWidget.cellChanged.connect(self.draw_curve)

    def draw_curve(self):
        #print("draw_curve")
        self.ax.clear()
        #self.ax = self.plt.add_subplot(111)
        x_values = []
        y_values = []
        array_dat = self.table_dat.get_values()
        #print(y_values)
        for r in array_dat:
            if r[0] != "" and r[1] != "":
                try:
                    _x = to_float(r[0])
                    _y = to_float(r[1])
                except:
                    print("Error drawing airfoil because _x or _y isn't a float.")
                    return
                x_values.append(_x)
                y_values.append(_y)
        #print(y_values)
        self.ax.set_xlim(0,1)
        self.ax.set_ylim(-0.5,0.5)
        self.ax.plot(x_values,y_values)
        self.plt.canvas.draw()

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
        y_max = numpy.max(y)
        y_min = numpy.min(y)
        thickness = (abs(y_max)+abs(y_min))/1
        print("Thickness:",thickness)
        return thickness


    def get_settings(self):
        x = []
        y = []
        #AoA_cD = []
        #cD = []

        array_dat = self.table_dat.get_values()
        for r in array_dat:
            if r[0] != "" and r[1] != "":
                _x = to_float(r[0])
                _y = to_float(r[1])
                #if _AoA > 90:
                #    _AoA = _AoA - 180
                x.append(_x)
                y.append(_y)

        #array_cd = self.table_cd.get_values()
        #for r in array_cd:
        #    if r[0] != "" and r[1] != "":
        #        _AoA = to_float(r[0])
        #        _cD = to_float(r[1])
        #        if _AoA > 90:
        #            _AoA = _AoA - 180
        #        AoA_cD.append(_AoA)
        #        cD.append(_cD)

        #f_c_L = interpolate.interp1d(
        #    AoA_cL, cL, fill_value=(cL[0], cL[-1]), bounds_error=False
        #)
        #f_c_D = interpolate.interp1d(
        #    AoA_cD, cD, fill_value=(cD[0], cD[-1]), bounds_error=False
        #)
        #inverse_f_c_L = interpolate.interp1d(
        #    cL, AoA_cL, fill_value=(AoA_cL[0], AoA_cL[-1]), bounds_error=False
        #)
        #alpha_zero = inverse_f_c_L(0.0)

        #max_thickness = float(self.max_thickness.text())
        #dat = self.generate_dat()
        return {
            #"f_c_L": f_c_L,
            #"f_c_D": f_c_D,
            #"AoA_cL": AoA_cL,
            #"AoA_cD": AoA_cD,
            #"cL": cL,
            #"cD": cD,
            "x":x,
            "y":y,
            "max_thickness": self.get_max_thickness(),
            #"dat":dat
            #"inverse_f_c_L": inverse_f_c_L,
            #"alpha_zero":alpha_zero,
        }


    def set_settings(self, dict_settings):
        #array_cl = []
        #for r in range(len(dict_settings["AoA_cL"])):
        #    array_cl.append(
        #        [str(dict_settings["AoA_cL"][r]), str(dict_settings["cL"][r])]
        #    )
        #array_cd = []
        #for r in range(len(dict_settings["AoA_cD"])):
        #    array_cd.append(
        #        [str(dict_settings["AoA_cD"][r]), str(dict_settings["cD"][r])]
        #    )
        array_dat = []
        for r in range(len(dict_settings["x"])):
            array_dat.append(
                [str(dict_settings["x"][r]), str(dict_settings["y"][r])])

        self.table_dat.createTable(array_dat)
        #self.table_cd.createTable(array_cd)
        #print(dict_settings)
        #self.max_thickness.setText(str(dict_settings["max_thickness"]))
        self.draw_curve()


class Analysis(QWidget):
    def __init__(self, parent=None):
        super(Analysis, self).__init__(parent)

        self.main = self.parent()

        self.settings = {
            "tip_loss": False,
            "hub_loss": False,
            "new_tip_loss": False,
            "new_hub_loss": False,
            "cascade_correction": False,
            "rotational_augmentation_correction": False,
            "rotational_augmentation_correction_method": 1,
            "max_iterations": 100,
            "convergence_limit": 0.001,
            "rho": 1.225,
            "method": 10,
            "linspace_interp": False,
            "num_interp": 25,
            "v_min": 3,
            "v_max": 20,
            "v_num": 10,
            "rpm_min": 100,
            "rpm_max": 3000,
            "rpm_num": 10,
            "relaxation_factor": 0.3,
            "print_all": False,
            "print_out": False,
            "reynolds":50000,
            "fix_reynolds":False
        }

        self.settings_to_name = {
            "print_out": "Print final iteration data",
            "tip_loss": "Prandtl tip loss",
            "hub_loss": "Prandtl hub loss",
            "new_tip_loss": "New tip loss",
            "new_hub_loss": "New hub loss",
            "cascade_correction": "Cascade correction",
            "max_iterations": "Maximum iterations",
            "convergence_limit": "Convergence criteria",
            "rho": "Air density [kg/m^3]",
            "method": "Calculation method",
            "v_min": "Min calc. wind speed [m/s]",
            "v_max": "Max calc. wind speed [m/s]",
            "v_num": "Number of wind speed points",
            "rpm_min": "Min calc. RPM [RPM]",
            "rpm_max": "Max calc. RPM [RPM]",
            "rpm_num": "Number of RPM points",
            "relaxation_factor": "Relaxation factor",
            "print_all": "Print every iteration [debug]",
            "num_interp": "Number of sections (interp)",
            "linspace_interp": "Custom number of sections",
            "rotational_augmentation_correction": "Rot. augmentation cor.",
            "rotational_augmentation_correction_method": "Rot. augmentation cor. method",
            "fix_reynolds":"Fix Reynolds",
            "reynolds":"Reynolds"
        }

        self.methods_to_names = METHODS_STRINGS

        self.name_to_methods = {v: k for k, v in self.methods_to_names.items()}
        self.name_to_settings = {v: k for k, v in self.settings_to_name.items()}

        self.grid = QGridLayout()
        self.setLayout(self.grid)

        self.left = QWidget()
        self.fbox = QFormLayout()
        self.left.setLayout(self.fbox)

        self.scroll_area = QtWidgets.QScrollArea()
        self.scroll_widget = QtWidgets.QWidget()
        self.scroll_widget_layout = QtWidgets.QVBoxLayout()

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
                form.addItems([self.methods_to_names[k] for k,v in self.methods_to_names.items()])
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
            key_orig = key
            key = self.settings_to_name[key]
            self.fbox.addRow(key, form)
            self.form_list.append([key, form, key_orig])

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

        # self.manager = Manager()
        # self.return_print = self.manager.list()
        # self.return_results = self.manager.list()

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
        sender.setStyleSheet("QLineEdit { background-color: %s }" % color)

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

    def run(self):
        self.clear()
        check = self.check_forms()
        if check != True:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText("Input validation error.")
            msg.setDetailedText(check)
            msg.exec_()
            return

        self.main.emitter_add.connect(self.add_text)
        self.main.emitter_done.connect(self.done)

        if not self.main.running:
            self.main.set_buttons_running()
            self.main.running = True
            self.runner_input = self.main.get_input_params()
            if self.runner_input == None:
                print("No settings fetched... Terminating.")
                self.terminate()
                return self.done(True)
            self.main.getter.start()
            self.p = Process(target=calculate_power_3d, args=[self.runner_input])
            self.p.start()
        else:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText("Cannot run while existing operation is running")
            msg.setInformativeText(
                "The program detected that an existing operation is running."
            )
            msg.setWindowTitle("Runtime error")
            msg.setDetailedText(
                "Currently tha value MainWindow.running is %s, \
                it should be False."
                % str(self.main.running)
            )
            msg.exec_()

    def add_text(self, string):
        self.textEdit.insertPlainText(string)
        if self.buttonEOF.checkState() == 2:
            self.textEdit.moveCursor(QtGui.QTextCursor.End)

    def done(self, terminated=False):
        self.main.emitter_add.disconnect()
        self.main.emitter_done.disconnect()
        self.main.set_buttons_await()
        self.main.running = False
        self.main.getter.__del__()
        if not terminated:
            self.p.join()
            if len(self.main.return_results) > 0:
                results = self.main.return_results[0]
                if "v" in results:
                    if len(results["v"]) > 0:
                        inp_params = self.runner_input
                        r = ResultsWindow(
                            self,
                            self.main.screen_width,
                            self.main.screen_width,
                            results,
                            inp_params,
                        )
                    else:
                        print("Not enough points to print results...")
                else:
                    print("No results to print...")

    def clear(self):
        self.textEdit.clear()

    def terminate(self):
        if hasattr(self, "p"):
            if self.p.is_alive():
                self.p.terminate()
                self.main.running = False
                self.main.getter.__del__()
                self.done(True)


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

        self.delta_start = QLineEdit()
        self.delta_start.setValidator(self.validator)
        self.delta_start.textChanged.connect(self.check_state)
        self.delta_start.textChanged.emit(self.delta_start.text())
        self._delta_start = QLabel("delta_start")
        self.form_list.append([self._delta_start, self.delta_start])

        self.decrease_factor = QLineEdit()
        self.decrease_factor.setValidator(self.validator)
        self.decrease_factor.textChanged.connect(self.check_state)
        self.decrease_factor.textChanged.emit(self.decrease_factor.text())
        self._decrease_factor = QLabel("decrease_factor")
        self.form_list.append([self._decrease_factor, self.decrease_factor])

        self.min_delta = QLineEdit()
        self.min_delta.setValidator(self.validator)
        self.min_delta.textChanged.connect(self.check_state)
        self.min_delta.textChanged.emit(self.min_delta.text())
        self._min_delta = QLabel("min_delta")
        self.form_list.append([self._min_delta, self.min_delta])

        self.min_add_angle = QLineEdit()
        self.min_add_angle.setValidator(self.validator)
        self.min_add_angle.textChanged.connect(self.check_state)
        self.min_add_angle.textChanged.emit(self.min_add_angle.text())
        self._min_add_angle = QLabel("min_add_angle")
        self.form_list.append([self._min_add_angle, self.min_add_angle])

        self.max_add_angle = QLineEdit()
        self.max_add_angle.setValidator(self.validator)
        self.max_add_angle.textChanged.connect(self.check_state)
        self.max_add_angle.textChanged.emit(self.max_add_angle.text())
        self._max_add_angle = QLabel("max_add_angle")
        self.form_list.append([self._max_add_angle, self.max_add_angle])

        self.angle_step = QLineEdit()
        self.angle_step.setValidator(self.validator)
        self.angle_step.textChanged.connect(self.check_state)
        self.angle_step.textChanged.emit(self.angle_step.text())
        self._angle_step = QLabel("angle_step")
        self.form_list.append([self._angle_step, self.angle_step])

        self.buttonAngles = QPushButton("Run angle optimization")
        self.buttonAngles.clicked.connect(self.run)

        self.buttonPitch = QPushButton("Run pitch optimization")
        self.buttonPitch.clicked.connect(self.run_pitch)

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
        self.fbox.addRow(self._target_speed, self.target_speed)
        self.fbox.addRow(self._target_rpm, self.target_rpm)

        rpmtext = QTextEdit("Target RPM can be empty for pitch calc.")
        rpmtext.setReadOnly(True)
        self.fbox.addRow(QLabel(" "), rpmtext)

        self.fbox.addRow(QLabel("--------"))
        self.fbox.addRow(self._delta_start, self.delta_start)
        self.fbox.addRow(self._min_delta, self.min_delta)
        self.fbox.addRow(self._decrease_factor, self.decrease_factor)
        self.fbox.addRow(self.buttonAngles)

        self.fbox.addRow(QLabel("--------"))
        self.fbox.addRow(self._min_add_angle, self.min_add_angle)
        self.fbox.addRow(self._max_add_angle, self.max_add_angle)
        self.fbox.addRow(self._angle_step, self.angle_step)
        self.fbox.addRow(self.buttonPitch)
        self.fbox.addRow(self.buttonClear, self.buttonStop)
        self.fbox.addRow(self.buttonEOFdescription, self.buttonEOF)


    def check_forms_angles(self):
        out = ""
        _needed_vars = [
            [self._target_speed, self.target_speed],
            [self._target_rpm, self.target_rpm],
            [self._delta_start, self.delta_start],
            [self._decrease_factor, self.decrease_factor],
            [self._min_delta, self.min_delta],
        ]
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

    def check_forms_pitch(self):
        out = ""
        _needed_vars = [
            [self._target_speed, self.target_speed],
            [self._min_add_angle, self.min_add_angle],
            [self._max_add_angle, self.max_add_angle],
            [self._angle_step, self.angle_step],
        ]
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
        sender.setStyleSheet("QLineEdit { background-color: %s }" % color)

    def run(self, run_pitch=False):
        self.clear()
        if run_pitch:
            check = self.check_forms_pitch()
        else:
            check = self.check_forms_angles()
        check_analysis = self.main.analysis.check_forms()
        if check != True or check_analysis != True:
            if check == True:
                check = ""
            if check_analysis == True:
                check_analysis = ""
            check = check + check_analysis
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText("Input validation error")
            msg.setDetailedText(check)
            msg.exec_()
            return

        self.main.emitter_add.connect(self.add_text)
        self.main.emitter_done.connect(self.done)

        if not self.main.running:
            self.main.set_buttons_running()
            self.main.running = True
            self.runner_input = self.main.get_input_params()
            self.main.getter.start()
            self.o = Optimizer(self.runner_input)
            if run_pitch:
                self.p = Process(target=self.o.optimize_pitch, args=[self.runner_input])
            else:
                self.p = Process(
                    target=self.o.optimize_angles, args=[self.runner_input]
                )
            self.p.start()
        else:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText("Cannot run while existing operation is running")
            msg.setInformativeText(
                "The program detected that an existing operation is running."
            )
            msg.setWindowTitle("Runtime error")
            msg.setDetailedText(
                "Currently tha value MainWindow.running is %s, \
                it should be False."
                % str(self.main.running)
            )

    def clear(self):
        self.textEdit.clear()

    def add_text(self, string):
        self.textEdit.insertPlainText(string)
        if self.buttonEOF.checkState() == 2:
            self.textEdit.moveCursor(QtGui.QTextCursor.End)

    def run_pitch(self):
        self.run(True)

    def terminate(self):
        if hasattr(self, "p"):
            if self.p.is_alive():
                self.p.terminate()
                self.main.running = False
                self.main.getter.__del__()
                self.done(True)

    def done(self, terminated=False):
        self.main.emitter_add.disconnect()
        self.main.emitter_done.disconnect()
        self.main.set_buttons_await()
        self.main.running = False
        self.main.getter.__del__()
        if not terminated:
            self.p.join()

    def get_settings(self):
        out = {}
        if self.target_rpm.text() == "":
            out["target_rpm"] = None
        else:
            out["target_rpm"] = self.target_rpm.text()
        out["target_speed"] = self.target_speed.text()
        out["delta_start"] = self.delta_start.text()
        out["decrease_factor"] = self.decrease_factor.text()
        out["min_delta"] = self.min_delta.text()
        out["min_add_angle"] = self.min_add_angle.text()
        out["max_add_angle"] = self.max_add_angle.text()
        out["angle_step"] = self.angle_step.text()
        for k, v in out.items():
            if v == "":
                v = None
            elif v == None:
                pass
            else:
                v = to_float(v)
            out[k] = v
        return out

    def set_settings(self, inp_dict):
        #print(inp_dict)
        self.target_rpm.setText(str(inp_dict["target_rpm"]))
        self.target_speed.setText(str(inp_dict["target_speed"]))
        self.delta_start.setText(str(inp_dict["delta_start"]))
        self.decrease_factor.setText(str(inp_dict["decrease_factor"]))
        self.min_delta.setText(str(inp_dict["min_delta"]))
        self.min_add_angle.setText(str(inp_dict["min_add_angle"]))
        self.max_add_angle.setText(str(inp_dict["max_add_angle"]))
        self.angle_step.setText(str(inp_dict["angle_step"]))


class ThreadGetter(QThread):
    def __init__(self, parent):
        super(ThreadGetter, self).__init__(parent)

    def __del__(self):
        self.wait()

    def run(self):
        print("Running Getter.")
        while True:
            if len(self.parent().return_print) > 0:
                t = self.parent().return_print.pop(0)
                self.parent().emitter_add.emit(str(t))
                if "!!!!EOF!!!!" in t:
                    self.parent().emitter_done.emit()
                    break
            if self.parent().running == False:
                break
        print("Getter finished.")
        return


class TabWidget(QtWidgets.QTabWidget):
    def __init__(self, parent=None):
        super(TabWidget, self).__init__(parent)
        self.tabs = []

    def add_tab(self, widget, tab_name):
        for t,n in self.tabs:
            if n == tab_name:
                print("n",n,"tab_name",tab_name)
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Warning)
                msg.setText("Tab with same name already exists!")
                msg.exec_()
                return
        self.tabs.append([widget,tab_name])
        self.addTab(widget, tab_name)
        return

    def remove_tab(self,index):
        self.removeTab(index)
        del self.tabs[index]

    def remove_all_tabs(self):
        while len(self.tabs) > 0:
            self.remove_tab(0)

    def remove_current_tab(self):
        self.remove_tab(self.currentIndex())

    def rename_current_tab(self,string):
        self.setTabText(self.currentIndex(),string)
        self.tabs[self.currentIndex()][1] = string

    def current_tab_name(self):
        return self.tabText(self.currentIndex())


WHITE =     QColor(255, 255, 255)
BLACK =     QColor(0, 0, 0)
RED =       QColor(255, 0, 0)
PRIMARY =   QColor(53, 53, 53)
SECONDARY = QColor(35, 35, 35)
TERTIARY =  QColor(42, 130, 218)
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
        self.setColor(QPalette.Window,          PRIMARY)
        self.setColor(QPalette.WindowText,      WHITE)
        self.setColor(QPalette.Base,            SECONDARY)
        self.setColor(QPalette.AlternateBase,   PRIMARY)
        self.setColor(QPalette.ToolTipBase,     WHITE)
        self.setColor(QPalette.ToolTipText,     WHITE)
        self.setColor(QPalette.Text,            LIGHT_GRAY)
        self.setColor(QPalette.Button,          PRIMARY)
        self.setColor(QPalette.ButtonText,      WHITE)
        self.setColor(QPalette.BrightText,      RED)
        self.setColor(QPalette.Link,            TERTIARY)
        self.setColor(QPalette.Highlight,       TERTIARY)
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

if __name__ == "__main__":
    if sys.platform.startswith("win"):
        # On Windows calling this function is necessary.
        multiprocessing.freeze_support()
    app = QtWidgets.QApplication([])
    app.setStyle("Fusion")
    if sys.platform.startswith("darwin"):
        #dark theme fix on OSX
        palette = QDarkPalette()
        palette.set_app(app)
        palette.set_stylesheet(app)
    screen = app.primaryScreen()
    size = screen.size()
    main = MainWindow(size.width(), size.height())
    sys.exit(app.exec_())

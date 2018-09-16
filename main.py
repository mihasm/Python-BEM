__author__ = "Miha Smrekar"
__credits__ = ["Miha Smrekar"]
__license__ = "GPL"
__version__ = "0.2.3"
__maintainer__ = "Miha Smrekar"
__email__ = "miha.smrekar9@gmail.com"
__status__ = "Development"

import sys

import numpy
from PyQt5 import QtWidgets
from PyQt5.QtCore import QThread, QTextStream, pyqtSignal, QProcess
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
    QMessageBox
)
from numpy import array
from scipy import interpolate

from cp_curve import calculate_power_3d
from results import ResultsWindow
from table import Table
import time
from turbine_data import SET_INIT
from optimisation import Optimizer

from multiprocessing import Process, Manager
import multiprocessing


class MainWindow(QMainWindow):
    emitter_add = pyqtSignal(str)
    emitter_done = pyqtSignal()

    def __init__(self, width, height):
        super().__init__()
        self.screen_width = width
        self.screen_height = height
        self.setGeometry(width * 0.125, height * 0.125, width * 0.75, height * 0.75)
        self.setWindowTitle("BEM analiza v0.2.3")
        self.tab_widget = TabWidget(self)
        self.setCentralWidget(self.tab_widget)

        self.wind_turbine_properties = WindTurbineProperties(self)
        self.tab_widget.add_tab(self.wind_turbine_properties, "Turbine info")

        self.curves = Curves(self)
        self.tab_widget.add_tab(self.curves, "Cl and Cd")

        self.analysis = Analysis(self)
        self.tab_widget.add_tab(self.analysis, "Analysis")

        self.getter = ThreadGetter(self)

        self.optimization = Optimization(self)
        self.tab_widget.add_tab(self.optimization,"Optimization")

        self.running = False
        self.manager = Manager()

        self.show()

    def get_all_settings(self):
        properties = self.wind_turbine_properties.get_properties()
        curves = self.curves.get_curves()
        settings = self.analysis.get_settings()
        opt_settings = self.optimization.get_settings()
        out = {**properties, **curves, **settings, **opt_settings}
        return out

    def get_input_params(self):
        settings = self.get_all_settings()
        inp_params = {
            **settings,
            "return_print": self.return_print,
            "return_results": self.return_results,
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

        grid = QGridLayout()
        self.setLayout(grid)

        left = QWidget()
        fbox = QFormLayout()
        left.setLayout(fbox)

        self.table_properties = Table()
        self.table_properties.createEmpty(4, 30)
        self.table_properties.set_labels(["r [m]", "c [m]", "theta [deg]", "dr [m]"])

        grid.addWidget(left, 1, 1)
        grid.addWidget(self.table_properties, 1, 2)

        _name = QLabel("Turbine Name")
        self.name = QLineEdit()
        self.name.setText("test")
        fbox.addRow(_name, self.name)

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

        self.set_properties(SET_INIT)

    def get_properties(self):
        out_properties = {
            "Rhub": float(self.Rhub.text()),
            "R": float(self.R.text()),
            "B": int(self.B.text()),
        }
        geom_array = self.table_properties.get_values()
        r, c, theta, dr = [], [], [], []
        for row in geom_array:
            if row[0] != "" and row[1] != "" and row[2] != "" and row[3] != "":
                r.append(float(row[0]))
                c.append(float(row[1]))
                theta.append(float(row[2]))
                dr.append(float(row[3]))
        out_properties["r"] = numpy.array(r)
        out_properties["c"] = numpy.array(c)
        out_properties["theta"] = numpy.array(theta)
        out_properties["dr"] = numpy.array(dr)
        return out_properties

    def set_properties(self, dict_settings):
        if "Rhub" in dict_settings:
            t = str(dict_settings["Rhub"])
            self.Rhub.setText(t)
        if "R" in dict_settings:
            t = str(dict_settings["R"])
            self.R.setText(t)
        if "B" in dict_settings:
            t = str(dict_settings["B"])
            self.B.setText(t)
        if (
            "r" in dict_settings
            and "c" in dict_settings
            and "theta" in dict_settings
            and "dr" in dict_settings
        ):
            _array = []
            for r in range(len(dict_settings["r"])):
                _r = dict_settings["r"][r]
                _c = dict_settings["c"][r]
                _theta = dict_settings["theta"][r]
                _dr = dict_settings["dr"][r]
                _array.append(
                    [
                        dict_settings["r"][r],
                        dict_settings["c"][r],
                        dict_settings["theta"][r],
                        dict_settings["dr"][r],
                    ]
                )
            self.table_properties.createTable(_array)


class Curves(QWidget):
    def __init__(self, parent=None):
        super(Curves, self).__init__(parent)
        grid = QGridLayout()
        self.setLayout(grid)

        self.table_cl = Table()
        self.table_cl.createEmpty(2, 50)
        self.table_cl.set_labels(["AoA", "Cl"])
        grid.addWidget(self.table_cl, 1, 1)

        self.table_cd = Table()
        self.table_cd.createEmpty(2, 50)
        self.table_cd.set_labels(["AoA", "Cd"])
        grid.addWidget(self.table_cd, 1, 2)

        self.set_curves(SET_INIT)

    def get_curves(self):
        AoA_cL = []
        cL = []
        AoA_cD = []
        cD = []

        array_cl = self.table_cl.get_values()
        for r in array_cl:
            if r[0] != "" and r[1] != "":
                AoA_cL.append(float(r[0]))
                cL.append(float(r[1]))

        array_cd = self.table_cd.get_values()
        for r in array_cd:
            if r[0] != "" and r[1] != "":
                AoA_cD.append(float(r[0]))
                cD.append(float(r[1]))

        f_c_L = interpolate.interp1d(
            AoA_cL, cL, fill_value=(cL[0], cL[-1]), bounds_error=False
        )
        f_c_D = interpolate.interp1d(
            AoA_cD, cD, fill_value=(cD[0], cD[-1]), bounds_error=False
        )
        return {
            "f_c_L": f_c_L,
            "f_c_D": f_c_D,
            "AoA_cL": AoA_cL,
            "AoA_cD": AoA_cD,
            "cL": cL,
            "cD": cD,
        }

    def set_curves(self, dict_settings):
        array_cl = []
        for r in range(len(dict_settings["AoA_cL"])):
            array_cl.append(
                [str(dict_settings["AoA_cL"][r]), str(dict_settings["cL"][r])]
            )
        array_cd = []
        for r in range(len(dict_settings["AoA_cD"])):
            array_cd.append(
                [str(dict_settings["AoA_cD"][r]), str(dict_settings["cD"][r])]
            )
        self.table_cl.createTable(array_cl)
        self.table_cd.createTable(array_cd)


class Analysis(QWidget):

    def __init__(self, parent=None):
        super(Analysis, self).__init__(parent)

        self.settings = {
            "tip_loss": False,
            "hub_loss": False,
            "new_tip_loss": False,
            "new_hub_loss": False,
            "cascade_correction": False,
            "max_iterations": 100,
            "convergence_limit": 0.001,
            "rho": 1.225,
            "method": 10,
            "v_min": 3,
            "v_max": 20,
            "v_num": 10,
            "rpm_min": 100,
            "rpm_max": 3000,
            "rpm_num": 10,
            "relaxation_factor": 0.3,
            "print_all": False,
            "print_out": False,
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
        }

        self.name_to_settings = {}
        for k, v in self.settings_to_name.items():
            self.name_to_settings[v] = k

        self.grid = QGridLayout()
        self.setLayout(self.grid)

        self.left = QWidget()
        self.fbox = QFormLayout()
        self.left.setLayout(self.fbox)

        self.grid.addWidget(self.left, 1, 1)

        self.textEdit = QTextEdit()
        self.textEdit.setReadOnly(True)
        self.grid.addWidget(self.textEdit, 1, 2)

        self.buttonRun = QPushButton("Run")

        self.form_list = []
        self.validator = QtGui.QDoubleValidator()

        for key, value in self.settings.items():
            if key == "method":
                form = QComboBox()
                form.addItems(["1", "4", "5", "6", "7", "8", "9", "10"])
            elif isinstance(value, bool):
                form = QCheckBox()
                form.setTristate(value)
            else:
                form = QLineEdit()
                form.setValidator(self.validator)
                form.textChanged.connect(self.check_state)
                form.textChanged.emit(form.text())
                form.insert(str(value))

            key = self.settings_to_name[key]
            self.fbox.addRow(key, form)
            self.form_list.append([key, form])

        self.emptyLabel = QLabel(" ")
        self.buttonRun.clicked.connect(self.run)
        self.buttonClear = QPushButton("Clear screen")
        self.buttonClear.clicked.connect(self.clear)
        self.buttonEOF = QCheckBox()
        self.buttonEOFdescription = QLabel("Scroll to end of screen")
        self.buttonStop = QPushButton("Stop")
        self.buttonStop.clicked.connect(self.terminate)

        self.fbox.addRow(self.emptyLabel, self.buttonRun)
        self.fbox.addRow(self.buttonClear, self.buttonStop)
        self.fbox.addRow(self.buttonEOFdescription,self.buttonEOF)

        #self.manager = Manager()
        #self.return_print = self.manager.list()
        #self.return_results = self.manager.list()

    def check_forms(self):
        out = ""
        for n,f in self.form_list:
            if isinstance(f,QLineEdit):
                state = self.validator.validate(f.text(), 0)[0]
                if state == QtGui.QValidator.Acceptable:
                    pass
                elif state == QtGui.QValidator.Intermediate:
                    out+=("Form %s appears not to be valid.\n" % n)
                else:
                    out+=("Form %s is not of the valid type.\n" % n)
        if out == "":
            return True
        return out

    def check_state(self, *args, **kwargs):
        sender = self.sender()
        validator = sender.validator()
        state = validator.validate(sender.text(), 0)[0]
        if state == QtGui.QValidator.Acceptable:
            color = '#edf5e1' # green
        elif state == QtGui.QValidator.Intermediate:
            color = '#fff79a' # yellow
        else:
            color = '#f6989d' # red
        sender.setStyleSheet('QLineEdit { background-color: %s }' % color)

    def get_settings(self):
        out_settings = {}
        for name, value in self.form_list:
            name = self.name_to_settings[name]
            if isinstance(value, QCheckBox):
                value = bool(value.checkState())
            if isinstance(value, QLineEdit):
                value = float(value.text())
            if isinstance(value, QComboBox):
                value = int(value.currentText())
            out_settings[name] = value
        return out_settings


    def run(self):
        check = self.check_forms()
        if check != True:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText("Input validation error.")
            msg.setDetailedText(check)
            msg.exec_()
            return

        self.parent().parent().parent().emitter_add.connect(self.add_text)
        self.parent().parent().parent().emitter_done.connect(self.done)

        if not self.parent().parent().parent().running:
            self.parent().parent().parent().return_print = self.parent().parent().parent().manager.list([])
            self.parent().parent().parent().return_results = self.parent().parent().parent().manager.list([])
            self.parent().parent().parent().set_buttons_running()
            self.parent().parent().parent().running = True
            self.runner_input  = self.parent().parent().parent().get_input_params()
            self.parent().parent().parent().getter.start()
            self.p = Process(target=calculate_power_3d, kwargs=self.runner_input)
            self.p.start()
        else:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText("Cannot run while existing operation is running")
            msg.setInformativeText("The program detected that an existing operation is running.")
            msg.setWindowTitle("Runtime error")
            msg.setDetailedText("Currently tha value MainWindow.running is %s, \
                it should be False." % str(self.parent().parent().parent().running))
            msg.exec_()


    def add_text(self, string):
        self.textEdit.insertPlainText(string)
        if self.buttonEOF.checkState() == 2:
            self.textEdit.moveCursor(QtGui.QTextCursor.End)

    def done(self, terminated=False):
        self.parent().parent().parent().emitter_add.disconnect()
        self.parent().parent().parent().emitter_done.disconnect()
        self.parent().parent().parent().set_buttons_await()
        self.parent().parent().parent().running = False
        self.parent().parent().parent().getter.__del__()
        if not terminated:
            self.p.join()
            results = self.parent().parent().parent().return_results[0]
            inp_params = self.runner_input
            r = ResultsWindow(
                self,
                self.parent().parent().parent().screen_width,
                self.parent().parent().parent().screen_width,
                results,
                inp_params,
            )

    def clear(self):
        self.textEdit.clear()

    def terminate(self):
        if hasattr(self,"p"):
            if self.p.is_alive():
                self.p.terminate()
                self.parent().parent().parent().running = False
                self.parent().parent().parent().getter.__del__()
                self.done(True)


class Optimization(QWidget):
    def __init__(self, parent=None):
        super(Optimization, self).__init__(parent)
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
        self.form_list.append([self._target_speed,self.target_speed])

        self.target_rpm = QLineEdit()
        self.target_rpm.setValidator(self.validator)
        self.target_rpm.textChanged.connect(self.check_state)
        self.target_rpm.textChanged.emit(self.target_rpm.text())
        self._target_rpm = QLabel("Target rpm [RPM]")
        self.form_list.append([self._target_rpm,self.target_rpm])

        
        self.delta_start     = QLineEdit()
        self.delta_start.setValidator(self.validator)
        self.delta_start.textChanged.connect(self.check_state)
        self.delta_start.textChanged.emit(self.delta_start.text())
        self._delta_start    = QLabel("delta_start")
        self.form_list.append([self._delta_start,self.delta_start])

        self.decrease_factor = QLineEdit()
        self.decrease_factor.setValidator(self.validator)
        self.decrease_factor.textChanged.connect(self.check_state)
        self.decrease_factor.textChanged.emit(self.decrease_factor.text())
        self._decrease_factor= QLabel("decrease_factor")
        self.form_list.append([self._decrease_factor,self.decrease_factor])

        self.min_delta       = QLineEdit()
        self.min_delta.setValidator(self.validator)
        self.min_delta.textChanged.connect(self.check_state)
        self.min_delta.textChanged.emit(self.min_delta.text())
        self._min_delta      = QLabel("min_delta")
        self.form_list.append([self._min_delta,self.min_delta])

        self.min_add_angle   = QLineEdit()
        self.min_add_angle.setValidator(self.validator)
        self.min_add_angle.textChanged.connect(self.check_state)
        self.min_add_angle.textChanged.emit(self.min_add_angle.text())
        self._min_add_angle  = QLabel("min_add_angle")
        self.form_list.append([self._min_add_angle,self.min_add_angle])

        self.max_add_angle   = QLineEdit()
        self.max_add_angle.setValidator(self.validator)
        self.max_add_angle.textChanged.connect(self.check_state)
        self.max_add_angle.textChanged.emit(self.max_add_angle.text())
        self._max_add_angle  = QLabel("max_add_angle")
        self.form_list.append([self._max_add_angle,self.max_add_angle])

        self.angle_step      = QLineEdit()
        self.angle_step.setValidator(self.validator)
        self.angle_step.textChanged.connect(self.check_state)
        self.angle_step.textChanged.emit(self.angle_step.text())
        self._angle_step     = QLabel("angle_step")
        self.form_list.append([self._angle_step,self.angle_step])

        self.buttonAngles = QPushButton("Run angle optimization")
        self.buttonAngles.clicked.connect(self.run)

        self.buttonPitch = QPushButton("Run pitch optimization")
        self.buttonPitch.clicked.connect(self.run_pitch)

        self.buttonStop = QPushButton("Stop")
        self.buttonStop.clicked.connect(self.terminate)

        self.buttonClear = QPushButton("Clear screen")
        self.buttonClear.clicked.connect(self.clear)

        self.buttonEOF = QCheckBox()
        self.buttonEOFdescription = QLabel("Scroll to end of screen")

        self.grid = QGridLayout()
        self.setLayout(self.grid)
        self.grid.addWidget(self.left,1,1)
        self.grid.addWidget(self.textEdit,1,2)

        self.fbox = QFormLayout()
        self.left.setLayout(self.fbox)
        self.fbox.addRow(self._target_speed,self.target_speed)
        self.fbox.addRow(self._target_rpm,self.target_rpm)

        rpmtext = QTextEdit("Target RPM can be empty for pitch calc.")
        rpmtext.setReadOnly(True)
        self.fbox.addRow(QLabel(" "),rpmtext)


        self.fbox.addRow(QLabel("--------"))
        self.fbox.addRow(self._delta_start,self.delta_start)
        self.fbox.addRow(self._min_delta,self.min_delta)
        self.fbox.addRow(self._decrease_factor,self.decrease_factor)
        self.fbox.addRow(self.buttonAngles)

        self.fbox.addRow(QLabel("--------"))
        self.fbox.addRow(self._min_add_angle,self.min_add_angle)
        self.fbox.addRow(self._max_add_angle,self.max_add_angle)
        self.fbox.addRow(self._angle_step,self.angle_step)
        self.fbox.addRow(self.buttonPitch)
        self.fbox.addRow(self.buttonClear,self.buttonStop)
        self.fbox.addRow(self.buttonEOFdescription,self.buttonEOF)

        self.set_settings(SET_INIT)

    def check_forms_angles(self):
        out = ""
        _needed_vars=[[self._target_speed,self.target_speed],[self._target_rpm,self.target_rpm],[self._delta_start,self.delta_start],[self._decrease_factor,self.decrease_factor],[self._min_delta,self.min_delta]]
        for n,f in _needed_vars:
            if isinstance(f,QLineEdit):
                state = self.validator.validate(f.text(), 0)[0]
                if state == QtGui.QValidator.Acceptable:
                    pass
                elif state == QtGui.QValidator.Intermediate:
                    out+=("Form %s appears not to be valid.\n" % n.text())
                else:
                    out+=("Form %s is not of the valid type.\n" % n.text())
        if out == "":
            return True
        return out

    def check_forms_pitch(self):
        out = ""
        _needed_vars=[[self._target_speed,self.target_speed],[self._min_add_angle,self.min_add_angle],[self._max_add_angle,self.max_add_angle],[self._angle_step,self.angle_step]]
        for n,f in _needed_vars:
            if isinstance(f,QLineEdit):
                state = self.validator.validate(f.text(), 0)[0]
                if state == QtGui.QValidator.Acceptable:
                    pass
                elif state == QtGui.QValidator.Intermediate:
                    out+=("Form %s appears not to be valid.\n" % n.text())
                else:
                    out+=("Form %s is not of the valid type.\n" % n.text())
        if out == "":
            return True
        return out

    def check_state(self, *args, **kwargs):
        sender = self.sender()
        validator = sender.validator()
        state = validator.validate(sender.text(), 0)[0]
        if state == QtGui.QValidator.Acceptable:
            color = '#edf5e1' # green
        elif state == QtGui.QValidator.Intermediate:
            color = '#fff79a' # yellow
        else:
            color = '#f6989d' # red
        sender.setStyleSheet('QLineEdit { background-color: %s }' % color)

    def run(self,run_pitch=False):
        if run_pitch:
            check = self.check_forms_pitch()
        else:
            check = self.check_forms_angles()
        check_analysis = self.parent().parent().parent().analysis.check_forms()
        if check != True or check_analysis != True:
            if check == True:
                check = ""
            if check_analysis == True:
                check_analysis = ""
            check = check+check_analysis
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText("Input validation error")
            msg.setDetailedText(check)
            msg.exec_()
            return

        self.parent().parent().parent().emitter_add.connect(self.add_text)
        self.parent().parent().parent().emitter_done.connect(self.done)

        if not self.parent().parent().parent().running:
            self.parent().parent().parent().return_print = self.parent().parent().parent().manager.list([])
            self.parent().parent().parent().return_results = self.parent().parent().parent().manager.list([])
            self.parent().parent().parent().set_buttons_running()
            self.parent().parent().parent().running = True
            self.runner_input  = self.parent().parent().parent().get_input_params()
            self.parent().parent().parent().getter.start()
            self.o = Optimizer(**self.runner_input)
            if run_pitch:
                self.p = Process(target=self.o.optimize_pitch, kwargs=self.runner_input)
            else:
                self.p = Process(target=self.o.optimize_angles, kwargs=self.runner_input)
            self.p.start()
        else:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText("Cannot run while existing operation is running")
            msg.setInformativeText("The program detected that an existing operation is running.")
            msg.setWindowTitle("Runtime error")
            msg.setDetailedText("Currently tha value MainWindow.running is %s, \
                it should be False." % str(self.parent().parent().parent().running))

    def clear(self):
        self.textEdit.clear()

    def add_text(self, string):
        self.textEdit.insertPlainText(string)
        if self.buttonEOF.checkState() == 2:
            self.textEdit.moveCursor(QtGui.QTextCursor.End)

    def run_pitch(self):
        self.run(True)

    def terminate(self):
        if hasattr(self,"p"):
            if self.p.is_alive():
                self.p.terminate()
                self.parent().parent().parent().running = False
                self.parent().parent().parent().getter.__del__()
                self.done(True)

    def done(self, terminated=False):
        self.parent().parent().parent().emitter_add.disconnect()
        self.parent().parent().parent().emitter_done.disconnect()
        self.parent().parent().parent().set_buttons_await()
        self.parent().parent().parent().running = False
        self.parent().parent().parent().getter.__del__()
        if not terminated:
            self.p.join()

    def get_settings(self):
        out = {}
        if self.target_rpm.text() == "":
            out["target_rpm"] = None
        else:
            out["target_rpm"] = self.target_rpm.text()
        out["target_speed"] = self.target_speed.text()
        out["delta_start"]=self.delta_start.text()
        out["decrease_factor"]=self.decrease_factor.text()
        out["min_delta"]=self.min_delta.text()
        out["min_add_angle"]=self.min_add_angle.text()
        out["max_add_angle"]=self.max_add_angle.text()
        out["angle_step"]=self.angle_step.text()
        for k,v in out.items():
            if v == "":
                v = None
            elif v == None:
                pass
            else:
                v = float(v)
            out[k] = v
        #print(out)
        return out

    def set_settings(self,inp_dict):
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
        #print("This was already in:",self.parent().return_print)
        while True:
            if len(self.parent().return_print) > 0:
                t = self.parent().return_print.pop(0)
                self.parent().emitter_add.emit(str(t))
                if t == "!!!!EOF!!!!":
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
        self.tabs.append(widget)
        self.addTab(self.tabs[-1], tab_name)


if __name__ == "__main__":
    if sys.platform.startswith('win'):
        # On Windows calling this function is necessary.
        multiprocessing.freeze_support()
    app = QtWidgets.QApplication([])
    app.setStyle(QStyleFactory.create('Fusion'))
    screen = app.primaryScreen()
    size = screen.size()
    main = MainWindow(size.width(), size.height())
    sys.exit(app.exec_())

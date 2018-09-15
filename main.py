__author__ = "Miha Smrekar"
__credits__ = ["Miha Smrekar"]
__license__ = "GPL"
__version__ = "0.2"
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
)
from numpy import array
from scipy import interpolate

from cp_curve import calculate_power_3d
from results import ResultsWindow
from table import Table
import time

from multiprocessing import Process, Manager

SET_INIT = {
    "Rhub": 0.1,
    "R": 0.776,
    "B": 5,
    "r": array(
        [
            0.15,
            0.18,
            0.21,
            0.24,
            0.27,
            0.3,
            0.33,
            0.36,
            0.39,
            0.42,
            0.45,
            0.48,
            0.51,
            0.54,
            0.57,
            0.6,
            0.63,
            0.66,
            0.69,
            0.72,
            0.75,
            0.76,
        ]
    ),
    "c": array(
        [
            0.12429,
            0.14448,
            0.15224,
            0.14968,
            0.14213,
            0.13129,
            0.12198,
            0.1136,
            0.10666,
            0.10015,
            0.09487,
            0.08959,
            0.08441,
            0.08029,
            0.07728,
            0.07371,
            0.07097,
            0.06797,
            0.06637,
            0.06337,
            0.06072,
            0.05984,
        ]
    ),
    "theta": array(
        [
            29.24,
            33.41,
            34.9,
            34.6,
            32.39,
            20.96,
            21.71,
            13.33,
            11.04,
            8.18,
            7.05,
            5.88,
            3.41,
            1.09,
            0.57,
            -0.2,
            -0.66,
            -2.05,
            -2.29,
            -2.6,
            -2.92,
            -2.37,
        ]
    ),
    "dr": array(
        [
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.03,
            0.01,
        ]
    ),
    "AoA_cL": [
        "-22.5",
        "-22",
        "-21.5",
        "-21",
        "-20.5",
        "-20",
        "-19.5",
        "-19",
        "-18.5",
        "-18",
        "-17.5",
        "-17",
        "-16.5",
        "-16",
        "-15.5",
        "-15",
        "-14.5",
        "-14",
        "-13.5",
        "-13",
        "-12.5",
        "-12",
        "-11.5",
        "-10.5",
        "-10",
        "-9.5",
        "-9",
        "-8.5",
        "-8",
        "-7.5",
        "-7",
        "-6.5",
        "-6",
        "-5.5",
        "-5",
        "-4.5",
        "-4",
        "-3.5",
        "-3",
        "-2.5",
        "-2",
        "-1.5",
        "-1",
        "0",
        "0.5",
        "1",
        "1.5",
        "2",
        "2.5",
        "3.5",
        "4",
        "4.5",
        "5",
        "5.5",
        "6",
        "6.5",
        "7",
        "7.5",
        "8",
        "8.5",
        "9",
        "9.5",
        "10",
        "10.5",
        "11",
        "11.5",
        "12",
        "12.5",
        "13",
        "13.5",
        "14",
        "14.5",
        "15",
        "15.5",
        "16",
        "16.5",
        "17",
        "17.5",
        "18",
        "18.5",
        "19",
        "19.5",
        "20",
        "20.5",
        "21",
        "21.5",
        "22",
        "22.5",
        "23",
    ],
    "AoA_cD": [
        "-22.5",
        "-22",
        "-21.5",
        "-21",
        "-20.5",
        "-20",
        "-19.5",
        "-19",
        "-18.5",
        "-18",
        "-17.5",
        "-17",
        "-16.5",
        "-16",
        "-15.5",
        "-15",
        "-14.5",
        "-14",
        "-13.5",
        "-13",
        "-12.5",
        "-12",
        "-11.5",
        "-10.5",
        "-10",
        "-9.5",
        "-9",
        "-8.5",
        "-8",
        "-7.5",
        "-7",
        "-6.5",
        "-6",
        "-5.5",
        "-5",
        "-4.5",
        "-4",
        "-3.5",
        "-3",
        "-2.5",
        "-2",
        "-1.5",
        "-1",
        "0",
        "0.5",
        "1",
        "1.5",
        "2",
        "2.5",
        "3.5",
        "4",
        "4.5",
        "5",
        "5.5",
        "6",
        "6.5",
        "7",
        "7.5",
        "8",
        "8.5",
        "9",
        "9.5",
        "10",
        "10.5",
        "11",
        "11.5",
        "12",
        "12.5",
        "13",
        "13.5",
        "14",
        "14.5",
        "15",
        "15.5",
        "16",
        "16.5",
        "17",
        "17.5",
        "18",
        "18.5",
        "19",
        "19.5",
        "20",
        "20.5",
        "21",
        "21.5",
        "22",
        "22.5",
        "23",
    ],
    "cL": [
        "-0.7256",
        "-0.709",
        "-0.6931",
        "-0.6773",
        "-0.6615",
        "-0.6454",
        "-0.6293",
        "-0.613",
        "-0.5967",
        "-0.5801",
        "-0.5632",
        "-0.5496",
        "-0.5375",
        "-0.525",
        "-0.5122",
        "-0.4998",
        "-0.4879",
        "-0.4764",
        "-0.4657",
        "-0.4555",
        "-0.4461",
        "-0.4381",
        "-0.44",
        "-0.7456",
        "-0.7976",
        "-0.7721",
        "-0.7264",
        "-0.6877",
        "-0.6486",
        "-0.6122",
        "-0.5484",
        "-0.4816",
        "-0.4063",
        "-0.3301",
        "-0.2503",
        "-0.1748",
        "-0.095",
        "-0.025",
        "0.0543",
        "0.1214",
        "0.1899",
        "0.2577",
        "0.3317",
        "0.3828",
        "0.4392",
        "0.4958",
        "0.5526",
        "0.6092",
        "0.6658",
        "0.7753",
        "0.8283",
        "0.8781",
        "0.9256",
        "0.9635",
        "0.9946",
        "1.0353",
        "1.0764",
        "1.1103",
        "1.1517",
        "1.1913",
        "1.2229",
        "1.2573",
        "1.2778",
        "1.2813",
        "1.292",
        "1.3016",
        "1.3101",
        "1.3165",
        "1.3134",
        "1.295",
        "1.2607",
        "1.2106",
        "1.1565",
        "1.1301",
        "1.0992",
        "1.0525",
        "0.9139",
        "0.9223",
        "0.9316",
        "0.9415",
        "0.9517",
        "0.9622",
        "0.9727",
        "0.9833",
        "0.9939",
        "1.0044",
        "1.0137",
        "1.0225",
        "1.033",
    ],
    "cD": [
        "0.27666",
        "0.26956",
        "0.2632",
        "0.25678",
        "0.25043",
        "0.24415",
        "0.23798",
        "0.2319",
        "0.22592",
        "0.22012",
        "0.21539",
        "0.21116",
        "0.20595",
        "0.19986",
        "0.19321",
        "0.18622",
        "0.17896",
        "0.17146",
        "0.16379",
        "0.15574",
        "0.14752",
        "0.13893",
        "0.12623",
        "0.0473",
        "0.03114",
        "0.02624",
        "0.02735",
        "0.02669",
        "0.02768",
        "0.0283",
        "0.03041",
        "0.02477",
        "0.02581",
        "0.02416",
        "0.02044",
        "0.01818",
        "0.01565",
        "0.0147",
        "0.01362",
        "0.01264",
        "0.0117",
        "0.01043",
        "0.00939",
        "0.00853",
        "0.00865",
        "0.00883",
        "0.00906",
        "0.00943",
        "0.00979",
        "0.01066",
        "0.01114",
        "0.01174",
        "0.01254",
        "0.014",
        "0.01613",
        "0.01735",
        "0.01837",
        "0.01978",
        "0.0205",
        "0.02128",
        "0.02248",
        "0.02344",
        "0.02549",
        "0.02892",
        "0.03215",
        "0.03591",
        "0.04016",
        "0.04507",
        "0.0513",
        "0.05917",
        "0.0693",
        "0.08243",
        "0.09766",
        "0.11011",
        "0.12521",
        "0.14697",
        "0.21055",
        "0.22036",
        "0.22996",
        "0.23939",
        "0.24868",
        "0.25783",
        "0.26684",
        "0.27568",
        "0.28434",
        "0.29282",
        "0.30075",
        "0.30803",
        "0.31537",
    ],
    "print_out": False,
    "tip_loss": False,
    "hub_loss": False,
    "new_tip_loss": False,
    "new_hub_loss": False,
    "cascade_correction": False,
    "max_iterations": 100.0,
    "convergence_limit": 0.001,
    "rho": 1.225,
    "method": 10.0,
    "v_min": 3.0,
    "v_max": 20.0,
    "v_num": 10.0,
    "rpm_min": 100.0,
    "rpm_max": 3000.0,
    "rpm_num": 10.0,
}

SET_INIT["r"] = numpy.array(SET_INIT["r"])
SET_INIT["c"] = numpy.array(SET_INIT["c"])
SET_INIT["dr"] = numpy.array(SET_INIT["dr"])
SET_INIT["theta"] = numpy.array(SET_INIT["theta"])


class MainWindow(QMainWindow):
    def __init__(self, width, height):
        super().__init__()
        self.screen_width = width
        self.screen_height = height
        self.setGeometry(width * 0.125, height * 0.125, width * 0.75, height * 0.75)
        self.setWindowTitle("BEM analiza v0.2")
        self.tab_widget = TabWidget(self)
        self.setCentralWidget(self.tab_widget)

        self.wind_turbine_properties = WindTurbineProperties(self)
        self.tab_widget.add_tab(self.wind_turbine_properties, "Turbine info")

        self.curves = Curves(self)
        self.tab_widget.add_tab(self.curves, "Cl and Cd")

        self.analysis = Analysis(self)
        self.tab_widget.add_tab(self.analysis, "Analysis")

        self.analysis.getter = ThreadGetter(self)

        self.show()

    def get_all_settings(self):
        properties = self.wind_turbine_properties.get_properties()
        curves = self.curves.get_curves()
        settings = self.analysis.get_settings()
        out = {**properties, **curves, **settings}
        return out


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
    emitter_add = pyqtSignal(str)
    emitter_done = pyqtSignal()

    def __init__(self, parent=None):
        super(Analysis, self).__init__(parent)

        self.running = False
        self.running_getter = False

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
        self.grid.addWidget(self.textEdit, 1, 2)

        self.buttonSettings = QPushButton("Get Settings")
        self.buttonRun = QPushButton("Run")

        self.form_list = []
        for key, value in self.settings.items():
            if key == "method":
                form = QComboBox()
                form.addItems(["1", "4", "5", "6", "7", "8", "9", "10"])
            elif isinstance(value, bool):
                form = QCheckBox()
                form.setTristate(value)
            else:
                form = QLineEdit()
                form.insert(str(value))
            key = self.settings_to_name[key]
            self.fbox.addRow(key, form)
            self.form_list.append([key, form])

        self.buttonSettings.clicked.connect(self.parent().get_all_settings)
        self.buttonRun.clicked.connect(self.run)
        self.buttonClear = QPushButton("Clear screen")
        self.buttonClear.clicked.connect(self.clear)
        self.buttonEOF = QCheckBox("Jump to end of screen")
        self.buttonStop = QPushButton("Stop")
        self.buttonStop.clicked.connect(self.terminate)
        self.buttonStop.setEnabled(False)

        self.fbox.addRow(self.buttonEOF)
        self.fbox.addRow(self.buttonSettings, self.buttonClear)
        self.fbox.addRow(self.buttonRun, self.buttonStop)

        self.manager = Manager()
        self.return_print = self.manager.list()
        self.return_results = self.manager.list()

        self.emitter_add.connect(self.add_text)
        self.emitter_done.connect(self.done)

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
        self.buttonRun.setEnabled(False)
        self.buttonStop.setEnabled(True)
        self.running = True
        self.getter.start()
        inp_params = self.parent().parent().parent().get_all_settings()
        self.runner_input = {
            **inp_params,
            "return_print": self.return_print,
            "return_results": self.return_results,
        }
        self.p = Process(target=calculate_power_3d, kwargs=self.runner_input)
        self.p.start()

    def add_text(self, string):
        self.textEdit.insertPlainText(string)
        if self.buttonEOF.checkState() == 2:
            self.textEdit.moveCursor(QtGui.QTextCursor.End)

    def done(self, terminated=False):
        self.buttonRun.setEnabled(True)
        self.buttonStop.setEnabled(False)
        self.running = False
        if not terminated:
            self.p.join()
            results = self.return_results[0]
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
        self.p.terminate()
        self.done(True)


class ThreadGetter(QThread):
    def __init__(self, parent):
        super(ThreadGetter, self).__init__(parent)

    def __del__(self):
        self.wait()

    def run(self):
        while True:
            if len(self.parent().analysis.return_print) > 0:
                t = self.parent().analysis.return_print.pop(0)
                self.parent().analysis.emitter_add.emit(str(t))
                if t == "!!!!EOF!!!!":
                    self.parent().analysis.emitter_done.emit()
                    break
            if self.parent().analysis.running == False:
                break
        return


class TabWidget(QtWidgets.QTabWidget):
    def __init__(self, parent=None):
        super(TabWidget, self).__init__(parent)
        self.tabs = []

    def add_tab(self, widget, tab_name):
        self.tabs.append(widget)
        self.addTab(self.tabs[-1], tab_name)


if __name__ == "__main__":
    app = QtWidgets.QApplication([])
    screen = app.primaryScreen()
    size = screen.size()
    main = MainWindow(size.width(), size.height())
    sys.exit(app.exec_())

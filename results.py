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

import json
import warnings

import matplotlib as mpl
import matplotlib.cbook
import matplotlib.pyplot as plt
import numpy as np
from PyQt5 import QtWidgets
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import (QComboBox, QMainWindow, QWidget, QAction, QFileDialog)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from UI.Table import Table
from calculation import OUTPUT_VARIABLES_LIST
from utils import dict_to_ar, greek_letters_to_string, transpose, filter_3d_results

warnings.filterwarnings("ignore", category=matplotlib.cbook.MatplotlibDeprecationWarning)


class ResultsWindow(QMainWindow):
    """

    """

    def __init__(self, parent, width, height, results_3d, input_data):
        super(ResultsWindow, self).__init__(parent)
        self.input_data = input_data
        self.results_3d = results_3d
        self.title = "Results"

        self.statusBar()
        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu('&File')

        saveFile = QAction("&Export to JSON", self)
        saveFile.setShortcut("Ctrl+S")
        saveFile.setStatusTip('Export to JSON')
        saveFile.triggered.connect(self.export)
        fileMenu.addAction(saveFile)

        self.screen_width = width
        self.screen_height = height
        self.setWindowTitle(self.title)
        self.setGeometry(int(self.screen_width * 0.125), int(self.screen_height * 0.125), int(self.screen_width * 0.75),
                         int(self.screen_width * 0.75 * 0.4))
        self.tab_widget = TabWidget(self)
        self.setCentralWidget(self.tab_widget)

        # list_of_variable_parameter = []
        if input_data["variable_selection"] == 0:
            list_of_variable_parameter = np.array(results_3d["TSR"])
            variable_parameter_title = r"$\lambda$"
        elif input_data["variable_selection"] == 1:
            list_of_variable_parameter = np.array(results_3d["TSR"])
            variable_parameter_title = r"$\lambda$"
        elif input_data["variable_selection"] == 2:
            list_of_variable_parameter = np.array(results_3d["J"])
            variable_parameter_title = "J"
        elif input_data["variable_selection"] == 3:
            list_of_variable_parameter = np.array(results_3d["pitch"])
            variable_parameter_title = "pitch"
        elif input_data["variable_selection"] == 4:
            list_of_variable_parameter = np.array(results_3d["pitch"])
            variable_parameter_title = "TSR"

        ########### CP(lambda) CURVE ###############
        if input_data["turbine_type"] == 0:  # wind turbine
            f2 = self.tab_widget.add_tab_figure("Cp curve")

            mat = np.array(
                [results_3d["TSR"], results_3d["cp"], results_3d["pitch"], results_3d["blade_stall_percentage"]])
            # transpose
            mat = mat.transpose()
            # sort by multiple columns
            r = np.core.records.fromarrays([mat[:, 2], mat[:, 0]], names='a,b')
            mat = mat[r.argsort()]
            # split at different column values
            dif = np.diff(mat[:, 2])
            ind_split = np.where(dif != 0)[0]
            splitted = np.split(mat, ind_split + 1)

            ax_cp = f2.add_subplot(111)
            ax_cp.set_title(r"$C_{P}-\lambda$")
            ax_cp.set_xlabel(r"$\lambda$")
            ax_cp.set_ylabel(r"$C_{P}$")
            ax_cp.grid()

            for a in splitted:
                x = a[:, 0]
                y = a[:, 1]
                label = str(round(a[0, 2], 2))
                stall = (1 - a[:, 3])
                ax_cp.plot(x, y, label=label)
                sc = ax_cp.scatter(x, y, c=stall, cmap="RdYlGn")

            norm = mpl.colors.Normalize(vmin=0, vmax=1)
            sm = plt.cm.ScalarMappable(cmap="RdYlGn", norm=norm)
            sm.set_array([])

            cbar = plt.colorbar(sm, ticks=[0, 0.5, 1])
            cbar.set_label("Stalled area")
            cbar.ax.set_yticklabels(["100 %", "50 %", "0 %"])

            if len(splitted) > 1:
                ax_cp.legend(title="Pitch [°]")
        ############################################

        ############ Propeller curves ##############
        if input_data["turbine_type"] == 1:  # propeller
            f3 = self.tab_widget.add_tab_figure("Prop curves")
            J, eff, ct, cp, cq = results_3d["J"], results_3d["eff"], results_3d["ct"], results_3d["cp"], results_3d[
                "cq"]

            ax3 = f3.add_subplot(111)
            ax3.set_title("Propeller curves")
            ax3.set_xlabel("J")
            ax3.set_ylabel(r"$\eta$")
            ax3.grid()
            ax3.plot(J, eff, "-g", label=r"$\eta$")

            ax3_2 = ax3.twinx()
            ax3_2.plot(J, ct, "-r", label=r"$C_T$")
            ax3_2.plot(J, cp, "-b", label=r"$C_P$")
            ax3_2.plot(J, cq, "-", color="orange", label=r"$C_Q$")
            ax3_2.set_ylabel(r"$C_T, C_P, C_Q$")

            f3.legend(loc="center right")

        if input_data["turbine_type"] == 1:  # propeller
            f4 = self.tab_widget.add_tab_figure("Propeller efficency")
            mat = np.array([results_3d["J"], results_3d["eff"], results_3d["pitch"]])
            # transpose
            mat = mat.transpose()
            # sort by multiple columns
            r = np.core.records.fromarrays([mat[:, 2], mat[:, 0]], names='a,b')
            mat = mat[r.argsort()]
            # split at different column values
            dif = np.diff(mat[:, 2])
            ind_split = np.where(dif != 0)[0]
            splitted = np.split(mat, ind_split + 1)

            ax_cp = f4.add_subplot(111)
            ax_cp.set_title(r"Propeller efficiency")
            ax_cp.set_xlabel(r"J")
            ax_cp.set_ylabel(r"$\eta$")
            ax_cp.grid()

            for a in splitted:
                x = a[:, 0]
                y = a[:, 1]
                label = str(round(a[0, 2], 2))  # label by pitch
                ax_cp.plot(x, y, label=label)

            if len(splitted) > 1:
                ax_cp.legend(title="Pitch [°]")
        ############################################

        ############ geometry check ################
        array_geom = []
        for r in range(len(input_data["r"])):
            _r = input_data["r"][r]
            _c = input_data["c"][r]
            _theta = input_data["theta"][r]
            _dr = input_data["dr"][r]
            array_geom.append([input_data["r"][r], input_data["c"]
            [r], input_data["theta"][r], input_data["dr"][r], ])
        t_geom = Table()
        t_geom.createTable(array_geom)
        t_geom.set_labels(["r", "c", "theta", "dr"])
        self.tab_widget.add_tab_widget(t_geom, "Geometry check")
        ############################################

        ########### data ###########################
        data = dict_to_ar(results_3d)
        names = []
        symbols_and_units = []
        for i in data[0]:
            names.append(OUTPUT_VARIABLES_LIST[i]["name"])
            symbol = OUTPUT_VARIABLES_LIST[i]["symbol"]
            symbol = greek_letters_to_string(symbol)
            symbol = symbol.replace("$", "")
            unit = OUTPUT_VARIABLES_LIST[i]["unit"]
            unit = unit.replace("$", "")
            if unit == "":
                symbols_and_units.append(symbol)
            else:
                symbols_and_units.append(symbol + " [" + unit + "]")
        del data[0]
        data.insert(0, names)
        data.insert(1, symbols_and_units)
        t = Table()
        t.createTable(data)
        self.tab_widget.add_tab_widget(t, "Data")
        ############################################

        ########## CUSTOM GRAPH ####################
        self.custom_graph = CustomGraphWidget(self)
        self.custom_graph.set_data(results_3d)
        self.tab_widget.add_tab_widget(self.custom_graph, "Custom graph")
        ############################################

        self.show()

    def closeEvent(self, event):
        """

        :param event:
        """
        for f in self.tab_widget.figures:
            f.clear()
            plt.close(f)
        event.accept()  # let the window close

    def set_menubar(self):
        """

        """
        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu("File")
        exitButton = QAction(QIcon("exit24.png"), "Exit", self)
        exitButton.setShortcut("Ctrl+Q")
        exitButton.setStatusTip("Exit application")
        exitButton.triggered.connect(exit)
        fileMenu.addAction(exitButton)

    def export(self):
        """

        """
        name = QFileDialog.getSaveFileName(self, 'Save File', "", "JSON (*.json)")[0]
        if name != "":
            d = self.results_3d
            d_out = filter_3d_results(d)
            print(d_out)
            json_d = json.dumps(d_out)
            file = open(name, 'w')
            file.write(json_d)
            file.close()


class TabWidget(QtWidgets.QTabWidget):
    """

    """

    def __init__(self, parent=None):
        super(TabWidget, self).__init__(parent)
        self.tabs = []
        self.figures = []
        self.canvas = []
        self.toolbars = []

    def add_tab_figure(self, tab_name):
        """

        :param tab_name:
        :return:
        """
        self.tabs.append(QtWidgets.QWidget())
        self.addTab(self.tabs[-1], tab_name)
        self.figures.append(plt.figure(figsize=(10, 5)))
        self.canvas.append(FigureCanvas(self.figures[-1]))
        toolbar = NavigationToolbar(self.canvas[-1], self)
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.canvas[-1])
        layout.addWidget(toolbar)
        self.tabs[-1].setLayout(layout)
        return self.figures[-1]

    def add_tab_widget(self, widget, tab_name):
        """

        :param widget:
        :param tab_name:
        """
        self.tabs.append(widget)
        self.addTab(self.tabs[-1], tab_name)


class CustomGraphWidget(QWidget):
    """

    """

    def __init__(self, parent=None):
        super(CustomGraphWidget, self).__init__(parent)
        self.parent = parent

        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.ax = self.figure.add_subplot(111)

        self.layout = QtWidgets.QGridLayout()
        self.layout.addWidget(self.canvas, 1, 1)
        self.layout.addWidget(self.toolbar, 2, 1)

        self.table = Table()
        self.layout.addWidget(self.table, 1, 2)

        self.comboBox_x = QComboBox(self)
        self.comboBox_y = QComboBox(self)

        self.layout.addWidget(self.comboBox_x, 3, 1)
        self.layout.addWidget(self.comboBox_y, 4, 1)

        self.inverse_list = {}
        list_of_options = []

        for k, v in OUTPUT_VARIABLES_LIST.items():
            if v["unit"] != "":
                var_name = v["name"] + " " + v["symbol"] + " [" + v["unit"] + "]"
            else:
                var_name = v["name"] + " " + v["symbol"]
            var_name = greek_letters_to_string(var_name)
            var_name = var_name.replace("$", "")
            list_of_options.append(var_name)
            self.inverse_list[var_name] = k

        list_of_options.sort()
        for l in list_of_options:
            self.comboBox_x.addItem(l)
            self.comboBox_y.addItem(l)

        # self.button_draw = QPushButton("Draw")
        # self.layout.addWidget(self.button_draw)
        # self.button_draw.clicked.connect(self.draw_graph)

        self.comboBox_x.currentTextChanged.connect(self.draw_graph)
        self.comboBox_y.currentTextChanged.connect(self.draw_graph)

        self.setLayout(self.layout)

    def draw_graph(self):
        """

        """
        self.ax.clear()

        # set variable

        print("Current variable parameter:", self.parent.input_data["variable_selection"])

        if self.parent.input_data["variable_selection"] == 0:
            list_of_variable_parameter = np.array(self.results["TSR"])
            variable_parameter_title = r"$\lambda$"
        elif self.parent.input_data["variable_selection"] == 1:
            list_of_variable_parameter = np.array(self.results["TSR"])
            variable_parameter_title = r"$\lambda$"
        elif self.parent.input_data["variable_selection"] == 2:
            list_of_variable_parameter = np.array(self.results["J"])
            variable_parameter_title = "J"
        elif self.parent.input_data["variable_selection"] == 3:
            list_of_variable_parameter = np.array(self.results["pitch"])
            variable_parameter_title = "pitch"
        elif self.parent.input_data["variable_selection"] == 4:
            list_of_variable_parameter_1 = np.array(self.results["pitch"])
            variable_parameter_title_1 = "pitch"
            list_of_variable_parameter_2 = np.array(self.results["TSR"])
            variable_parameter_title_2 = "TSR"

        x_data = self.inverse_list[str(self.comboBox_x.currentText())]
        y_data = self.inverse_list[str(self.comboBox_y.currentText())]
        x_ar = np.array(self.results[x_data])
        y_ar = np.array(self.results[y_data])

        # set labels
        xlabel = OUTPUT_VARIABLES_LIST[x_data]["symbol"]
        if OUTPUT_VARIABLES_LIST[x_data]["unit"] != "":
            xlabel = xlabel + " [" + OUTPUT_VARIABLES_LIST[x_data]["unit"] + "]"

        ylabel = OUTPUT_VARIABLES_LIST[y_data]["symbol"]
        if OUTPUT_VARIABLES_LIST[y_data]["unit"] != "":
            ylabel = ylabel + " [" + OUTPUT_VARIABLES_LIST[y_data]["unit"] + "]"

        is_3d = False

        type_x = OUTPUT_VARIABLES_LIST[x_data]["type"]
        type_y = OUTPUT_VARIABLES_LIST[y_data]["type"]

        try:
            self.figure.delaxes(self.ax)
        except Exception as e:
            print(e)
            print("Did not find axes, not deleting...")
            pass

        if type_x == "array" and type_y == "array":
            print("case 1")
            self.ax = self.figure.add_subplot(111)
            x_ar = np.transpose(x_ar)
            y_ar = np.transpose(y_ar)
            self.ax.plot(x_ar, y_ar, label="2d plot")
            self.ax.legend(np.round(list_of_variable_parameter, 2), title=variable_parameter_title)
            # data_table = np.column_stack((x_ar,y_ar))
            # data_table=transpose(data_table)
            self.table.createTable(y_ar)
            # self.table.clear_table()


        elif type_x == "array" and type_y == "float":
            print("case 2")
            self.ax = self.figure.add_subplot(111)
            # x_ar=np.transpose(x_ar)
            # y_ar=np.transpose(y_ar)

            self.ax.plot(x_ar, y_ar, label="2d plot")
            self.ax.legend(self.parent.input_data["r"], title="r [m]")
            data_table = np.column_stack((y_ar, x_ar))
            self.table.createTable(data_table)


        elif type_x == "float" and type_y == "array":
            print("case 3")
            self.ax = self.figure.add_subplot(111)
            x_ar = np.transpose(x_ar)
            self.ax.plot(x_ar, y_ar, label="2d plot")
            self.ax.legend(self.parent.input_data["r"], title="r [m]")
            data_table = list(np.column_stack((x_ar, y_ar)))
            self.table.createTable(data_table)

        elif type_x == "float" and type_y == "float":
            print("case 4")
            self.ax = self.figure.add_subplot(111)
            self.ax.plot(x_ar, y_ar, label="2d plot")
            data_table = [x_ar, y_ar]
            data_table = transpose(data_table)
            self.table.createTable(data_table)

        self.ax.set_title(OUTPUT_VARIABLES_LIST[y_data]["name"] + " vs. " + OUTPUT_VARIABLES_LIST[x_data]["name"])
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)

        self.canvas.draw()

    def set_data(self, results):
        """

        :param results:
        """
        self.results = results
        self.draw_graph()

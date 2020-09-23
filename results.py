import warnings

from PyQt5 import QtWidgets
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import (QComboBox, QMainWindow, QPushButton, QTextEdit, QWidget, QFormLayout,
                             QLabel, QLineEdit, QGridLayout, QCheckBox, QStyleFactory, QMessageBox, QAction, QFileDialog, QSlider)
from calculation import OUTPUT_VARIABLES_LIST
from calculation_runner import max_calculate
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.cbook
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
from utils import sort_xy, dict_to_ar, ErrorMessageBox, MyMessageBox, greek_letters_to_string, transpose, Table

warnings.filterwarnings("ignore",category=matplotlib.cbook.MatplotlibDeprecationWarning)



class ResultsWindow(QMainWindow):
    def __init__(self, parent, width, height, results_3d, input_data):
        super(ResultsWindow, self).__init__(parent)
        self.input_data = input_data

        self.title = "Results"
        self.screen_width = width
        self.screen_height = height
        self.setWindowTitle(self.title)
        self.setGeometry(self.screen_width * 0.125, self.screen_height * 0.125, self.screen_width * 0.75,
                         self.screen_width * 0.75 * 0.4, )
        self.tab_widget = TabWidget(self)
        self.setCentralWidget(self.tab_widget)

        #list_of_variable_parameter = []
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

        # get data from results
        X, Y, Z = results_3d["v"], results_3d["rpm"], results_3d["power"]
        max_x, max_y, max_z = max_calculate(X, Y, Z)

        ########### CP(lambda) CURVE ###############
        f2 = self.tab_widget.add_tab_figure("Cp curve")
        TSR, CP = sort_xy(results_3d["TSR"], results_3d["cp"])
        ax_cp = self.tab_widget.add_2d_plot_to_figure(
            f2, TSR, CP, 121, "", r"$\lambda$", r"$C_P$", look="-b", label=r"$C_P$ curve")
        ax_cp.legend(fontsize=18)
        ############################################

        ########### CT(lambda) CURVE ###############
        ax_ct = self.tab_widget.add_2d_plot_to_figure(f2, results_3d["TSR"], results_3d["ct"], 122, "", r"$\lambda$", r"$C_T$",
                                              look="b-", label=r"$C_T$ curve")
        ax_ct.legend(fontsize=18)
        ############################################

        ########### Ct(a) curve ####################
        f3 = self.tab_widget.add_tab_figure("Ct_r(a) curve")
        ax_ct_r = self.tab_widget.add_2d_plot_to_figure(f3, results_3d["a"], results_3d["Ct"], 111, "", "a", r"$C_{T_r}$",
                                                        look="o", x_min=0, x_max=1, label="Ct curve")
        leg = ax_ct_r.legend(
            np.round(np.array(results_3d["TSR"]), 2), fontsize=20)
        leg.set_title(r"$\lambda$", prop={'size': 20})
        ############################################

        ########## 3D Power surface plot ###########
        # noinspection PyBroadException
        if len(set(X)) >= 3 and len(set(Y)) >= 3:
            f4 = self.tab_widget.add_tab_figure("3D power")
            self.tab_widget.add_surface_plot(
                f4, X, Y, Z, 111, "Power (windspeed,RPM)", "windspeed[m/s]", "RPM", "Power [W]", label="Power")
        ############################################

        ########## PROPELLER PLOTS #################
        if input_data["propeller_mode"] == True:
            # ct_p(J) curve
            f5 = self.tab_widget.add_tab_figure("ct(J) curve (propeller)")
            self.tab_widget.add_2d_plot_to_figure(f5, results_3d["J"], results_3d["ct"], 111, "Propeler curves",
                                                  "J = 1/lambda", "ct", look="b-", label="Ct (Thrust)")

            # cp_p(J) curve
            self.tab_widget.add_2d_plot_to_figure(f5, results_3d['J'], results_3d['cp'], 111, None, None,
                                                  None, look='r-', label="Cp (Power)")

            # eff_p(J) curve
            ax1 = self.tab_widget.add_2d_plot_to_figure(f5, results_3d['J'], results_3d['eff'], 111, None, None,
                                                        None, look='g-', label="Eff (Efficiency)")
            ax1.legend()

        ############################################

        ########## CL plot #########################
        f7 = self.tab_widget.add_tab_figure("check CL")
        self.tab_widget.add_3d_scatter_plot(f7, np.array(results_3d["Re"]).flatten(),
                                            np.degrees(
                                                np.array(results_3d["alpha"]).flatten()),
                                            np.array(results_3d["cL"]).flatten(), 111, "Title", "Re", "alpha[deg]",
            "cL", label="Used lift coefficient")
        ############################################

        ########## U4 (r) ##########################
        f8 = self.tab_widget.add_tab_figure("Windspeed U3 (radius)")
        for _u in results_3d['U3']:
            ax2 = self.tab_widget.add_2d_plot_to_figure(
                f8, _u, input_data['r'], 111, '', r'$v_4$ [m/s]', 'r [m]', look="-", c=np.random.rand(3,), label="Windspeed close behind"+str(_u))

        leg_hitrosti = ax2.legend(np.round(list_of_variable_parameter, 2))
        leg_hitrosti.set_title(variable_parameter_title, prop={'size': 20})
        ############################################

        ########## U4 (r) ##########################
        f8 = self.tab_widget.add_tab_figure("Windspeed U4 (radius)")
        for _u in results_3d['U4']:
            ax2 = self.tab_widget.add_2d_plot_to_figure(
                f8, _u, input_data['r'], 111, '', r'$v_4$ [m/s]', 'r [m]', look="-", c=np.random.rand(3,), label="Windspeed far away behind"+str(_u))
        leg_hitrosti = ax2.legend(np.round(list_of_variable_parameter, 2))
        leg_hitrosti.set_title(variable_parameter_title, prop={'size': 20})
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
            symbol = symbol.replace("$","")
            unit = OUTPUT_VARIABLES_LIST[i]["unit"]
            unit = unit.replace("$","")
            if unit == "":
                symbols_and_units.append(symbol)
            else:
                symbols_and_units.append(symbol+" ["+unit+"]")
        del data[0]
        data.insert(0,names)
        data.insert(1,symbols_and_units)
        #data.insert(2,units)
        t = Table()
        t.createTable(data)

        self.tab_widget.add_tab_widget(t, "Data")
        ############################################

        ########## CUSTOM GRAPH ####################
        self.custom_graph = CustomGraphWidget(self)
        self.custom_graph.set_data(results_3d)
        self.tab_widget.add_tab_widget(self.custom_graph,"Custom graph")
        ############################################

        self.show()

    def closeEvent(self, event):
        for f in self.tab_widget.figures:
            f.clear()
            plt.close(f)
        event.accept()  # let the window close

    def set_menubar(self):
        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu("File")
        exitButton = QAction(QIcon("exit24.png"), "Exit", self)
        exitButton.setShortcut("Ctrl+Q")
        exitButton.setStatusTip("Exit application")
        exitButton.triggered.connect(exit)
        fileMenu.addAction(exitButton)


class TabWidget(QtWidgets.QTabWidget):
    def __init__(self, parent=None):
        super(TabWidget, self).__init__(parent)
        self.tabs = []
        self.figures = []
        self.canvas = []
        self.toolbars = []

    def add_tab_figure(self, tab_name):
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

    def add_2d_plot_to_figure(self, f, x, y, whi, title=None, x_name=None, y_name=None, x_min=None, x_max=None,
                              y_min=None, y_max=None, look=None, legend=False, **kwargs):
        ax = f.add_subplot(whi)
        if look:
            ax.plot(x, y, look, **kwargs)
        else:
            ax.plot(x, y, **kwargs)
        if title:
            plt.title(title, fontsize=25)
        if x_name:
            ax.set_xlabel(x_name, fontsize=20)
        if y_name:
            ax.set_ylabel(y_name, fontsize=20)
        if x_min != None and x_max != None:
            ax.set_xlim(x_min, x_max)
        if y_min != None and y_max != None:
            ax.set_ylim(y_min, y_max)
        ax.xaxis.set_tick_params(labelsize=15)
        ax.yaxis.set_tick_params(labelsize=15)
        self.canvas[-1].draw()
        return ax

    def add_surface_plot(self, f, x, y, z, whi, title=None, x_name=None, y_name=None, z_name=None, x_min=None,
                         x_max=None, y_min=None, y_max=None, z_min=None, z_max=None, legend=False, **kwargs):
        ax = f.add_subplot(whi, projection="3d")
        p0 = ax.plot_trisurf(x, y, z, cmap=plt.cm.CMRmap, **kwargs)
        cbar = plt.colorbar(p0)

        if title:
            plt.title(title)
        if x_name:
            ax.set_xlabel(x_name)
        if y_name:
            ax.set_ylabel(y_name)
        if z_name:
            ax.set_zlabel(z_name)
            cbar.set_label(z_name)
        if x_min != None and x_max != None:
            ax.set_xlim(x_min, x_max)
        if y_min != None and y_max != None:
            ax.set_ylim(y_min, y_max)
        if z_min != None and z_max != None:
            ax.set_zlim(z_min, z_max)
        self.canvas[-1].draw()

    def add_3d_scatter_plot(self, f, x, y, z, whi, title=None, x_name=None, y_name=None, z_name=None, x_min=None,
                            x_max=None, y_min=None, y_max=None, z_min=None, z_max=None, legend=False, label=None):
        ax = f.add_subplot(whi, projection="3d")
        p0 = ax.scatter(x, y, z, cmap=plt.cm.CMRmap, label=label)
        # cbar = plt.colorbar(p0)

        if title:
            plt.title(title)
        if x_name:
            ax.set_xlabel(x_name)
        if y_name:
            ax.set_ylabel(y_name)
        if z_name:
            ax.set_zlabel(z_name)  # cbar.set_label(z_name)
        if x_min != None and x_max != None:
            ax.set_xlim(x_min, x_max)
        if y_min != None and y_max != None:
            ax.set_ylim(y_min, y_max)
        if z_min != None and z_max != None:
            ax.set_zlim(z_min, z_max)
        self.canvas[-1].draw()

    def add_tab_widget(self, widget, tab_name):
        self.tabs.append(widget)
        self.addTab(self.tabs[-1], tab_name)

    def add_runner_widget(self, input_data, tab_name="Runner"):
        r = RunnerWidget(input_data=input_data)
        self.tabs.append(r)
        self.addTab(self.tabs[-1], tab_name)

class CustomGraphWidget(QWidget):
    def __init__(self, parent=None):
        super(CustomGraphWidget, self).__init__(parent)
        self.parent = parent

        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.ax = self.figure.add_subplot(111)

        self.layout = QtWidgets.QGridLayout()
        self.layout.addWidget(self.canvas,1,1)
        self.layout.addWidget(self.toolbar,2,1)
        
        self.table=Table()
        self.layout.addWidget(self.table,1,2)

        self.comboBox_x = QComboBox(self)
        self.comboBox_y = QComboBox(self)

        self.layout.addWidget(self.comboBox_x,3,1)
        self.layout.addWidget(self.comboBox_y,4,1)

        self.inverse_list = {}
        list_of_options = []

        for k,v in OUTPUT_VARIABLES_LIST.items():
            if v["unit"] != "":
                var_name = v["name"]+" "+v["symbol"]+" ["+v["unit"]+"]"
            else:
                var_name = v["name"]+" "+v["symbol"]
            var_name = greek_letters_to_string(var_name)
            var_name = var_name.replace("$","")
            list_of_options.append(var_name)
            self.inverse_list[var_name]=k

        list_of_options.sort()
        for l in list_of_options:
            self.comboBox_x.addItem(l)
            self.comboBox_y.addItem(l)

        #self.button_draw = QPushButton("Draw")
        #self.layout.addWidget(self.button_draw)
        #self.button_draw.clicked.connect(self.draw_graph)

        self.comboBox_x.currentTextChanged.connect(self.draw_graph)
        self.comboBox_y.currentTextChanged.connect(self.draw_graph)

        self.setLayout(self.layout)

    def draw_graph(self):
        self.ax.clear()

        # set variable

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

        # get x,y data

        x_data = self.inverse_list[str(self.comboBox_x.currentText())]
        y_data = self.inverse_list[str(self.comboBox_y.currentText())]
        x_ar = np.array(self.results[x_data])
        y_ar = np.array(self.results[y_data])

        if not x_data in OUTPUT_VARIABLES_LIST:
            raise Exception("Variable %s not defined in OUTPUT_VARIABLES_LIST" % x_data)

        if not y_data in OUTPUT_VARIABLES_LIST:
            raise Exception("Variable %s not defined in OUTPUT_VARIABLES_LIST" % x_data)

        # set labels

        xlabel = OUTPUT_VARIABLES_LIST[x_data]["symbol"]
        if OUTPUT_VARIABLES_LIST[x_data]["unit"] != "":
            xlabel = xlabel+" ["+OUTPUT_VARIABLES_LIST[x_data]["unit"]+"]"

        ylabel = OUTPUT_VARIABLES_LIST[y_data]["symbol"]
        if OUTPUT_VARIABLES_LIST[y_data]["unit"] != "":
            ylabel = ylabel+" ["+OUTPUT_VARIABLES_LIST[y_data]["unit"]+"]"

        is_3d = False

        type_x = OUTPUT_VARIABLES_LIST[x_data]["type"]
        type_y = OUTPUT_VARIABLES_LIST[y_data]["type"]

        try:
            self.figure.delaxes(self.ax)
        except KeyError:
            print("Did not find axes, not deleting...")
            pass

        if type_x == "array" and type_y == "array":
            print("case 1")
            self.ax = self.figure.add_subplot(111)
            x_ar=np.transpose(x_ar)
            y_ar=np.transpose(y_ar)
            self.ax.plot(x_ar,y_ar,label="2d plot")
            self.ax.legend(np.round(list_of_variable_parameter,2),title=variable_parameter_title)
            #data_table = np.column_stack((x_ar,y_ar))
            #data_table=transpose(data_table)
            #self.table.createTable(data_table)
            self.table.clear_table()

        elif type_x == "array" and type_y == "float":
            print("case 2")
            self.ax = self.figure.add_subplot(111)
            #x_ar=np.transpose(x_ar)
            #y_ar=np.transpose(y_ar)

            self.ax.plot(x_ar,y_ar,label="2d plot")
            self.ax.legend(self.parent.input_data["r"],title="r [m]")
            data_table = np.column_stack((y_ar,x_ar))
            self.table.createTable(data_table)
            

        elif type_x == "float" and type_y == "array":
            print("case 3")
            self.ax = self.figure.add_subplot(111)
            x_ar=np.transpose(x_ar)
            self.ax.plot(x_ar,y_ar,label="2d plot")
            self.ax.legend(self.parent.input_data["r"],title="r [m]")
            data_table = list(np.column_stack((x_ar,y_ar)))
            self.table.createTable(data_table)

        elif type_x == "float" and type_y == "float":
            print("case 4")
            self.ax = self.figure.add_subplot(111)
            self.ax.plot(x_ar,y_ar,label="2d plot")
            self.ax.legend(self.parent.input_data["r"],title="r [m]")
            data_table = [x_ar,y_ar]
            data_table=transpose(data_table)
            self.table.createTable(data_table)
            

        self.ax.set_title(OUTPUT_VARIABLES_LIST[y_data]["name"]+" vs. "+OUTPUT_VARIABLES_LIST[x_data]["name"])
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)
        

        self.canvas.draw()

    def set_data(self,results):
        self.results = results
        self.draw_graph()

import numpy
import warnings
from PyQt5 import QtWidgets
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import (QComboBox, QMainWindow, QPushButton, QTextEdit, QWidget, QFormLayout,
                             QLabel, QLineEdit, QGridLayout, QCheckBox, QStyleFactory, QMessageBox, QAction, QFileDialog, QSlider)
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure


from calculation_runner import max_calculate
from table import Table
from utils import sort_xy, dict_to_ar, ErrorMessageBox, MyMessageBox

from calculation import OUTPUT_VARIABLES_LIST

import matplotlib.cbook
warnings.filterwarnings("ignore",category=matplotlib.cbook.MatplotlibDeprecationWarning)



class ResultsWindow(QMainWindow):
    def __init__(self, parent, width, height, results_3d, input_data):
        super(ResultsWindow, self).__init__(parent)
        self.title = "Results"
        self.screen_width = width
        self.screen_height = height
        self.setWindowTitle(self.title)
        self.setGeometry(self.screen_width * 0.125, self.screen_height * 0.125, self.screen_width * 0.75,
                         self.screen_width * 0.75 * 0.4, )
        self.tab_widget = TabWidget(self)
        self.setCentralWidget(self.tab_widget)

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
        f8 = self.tab_widget.add_tab_figure("Windspeed (radius)")
        for _u in results_3d['U4']:
            ax2 = self.tab_widget.add_2d_plot_to_figure(
                f8, _u, input_data['r'], 111, '', r'$v_4$ [m/s]', 'r [m]', look="-", c=numpy.random.rand(3,), label="Windspeed far away behind"+str(_u))
        leg_hitrosti = ax2.legend(np.round(np.array(results_3d["TSR"]), 2))
        leg_hitrosti.set_title(r"$\lambda$", prop={'size': 20})
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
        t = Table()
        t.createTable(data)

        self.tab_widget.add_tab_widget(t, "Data")
        ############################################

        ########## CUSTOM GRAPH ####################
        self.custom_graph = CustomGraphWidget(self)
        self.tab_widget.add_tab_widget(self.custom_graph,"Custom graph")
        self.custom_graph.set_data(results_3d)
        ############################################

        self.show()

    def closeEvent(self, event):
        for f in self.tab_widget.figures:
            f.clear()
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

        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.ax = self.figure.add_subplot(111)

        self.layout = QtWidgets.QVBoxLayout()
        self.layout.addWidget(self.toolbar)
        self.layout.addWidget(self.canvas)

        self.comboBox_x = QComboBox(self)
        self.comboBox_y = QComboBox(self)

        self.layout.addWidget(self.comboBox_x)
        self.layout.addWidget(self.comboBox_y)

        self.inverse_list = {}
        list_of_options = []
        for k,v in OUTPUT_VARIABLES_LIST.items():
            list_of_options.append(v["name"])
            self.inverse_list[v["name"]]=k
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

        x_data = self.inverse_list[str(self.comboBox_x.currentText())]
        y_data = self.inverse_list[str(self.comboBox_y.currentText())]
        x_ar = np.array(self.results[x_data])
        y_ar = np.array(self.results[y_data])

        legend_r = False

        if OUTPUT_VARIABLES_LIST[x_data]["type"] == "array" and OUTPUT_VARIABLES_LIST[y_data]["type"] == "array":
            x_ar = numpy.array(x_ar)
            x_ar = numpy.transpose(x_ar)
            y_ar = numpy.array(y_ar)
            y_ar = numpy.transpose(y_ar)
            legend_r = True

        elif OUTPUT_VARIABLES_LIST[x_data]["type"] == "array" or OUTPUT_VARIABLES_LIST[y_data]["type"] == "array":
            legend_r = True

        draw_3d = True

        self.figure.delaxes(self.ax)
        
        if legend_r:
            if draw_3d:
                self.ax = self.figure.add_subplot(111, projection="3d", label="3d plot")
                #ax = f.add_subplot(whi, projection="3d")
                num_columns = None
                if len(x_ar.shape) == 1:
                    num_columns = y_ar.shape[1]
                    x_ar = np.array(x_ar)
                    x_ar = np.tile(x_ar,(num_columns,1))
                    x_ar = np.transpose(x_ar)
                if len(y_ar.shape) == 1:
                    num_columns = x_ar.shape[1]
                    y_ar = np.array(y_ar)
                    y_ar = np.tile(y_ar,(num_columns,1))
                    y_ar = np.transpose(y_ar)

                if num_columns == None:
                    num_columns = x_ar.shape[1]

                z_ar = np.tile(np.array(self.results["TSR"]),(num_columns,1))
                z_ar = np.transpose(z_ar)
                try:
                    #p0 = self.ax.plot_trisurf(z_ar.flatten(),x_ar.flatten(),y_ar.flatten(), cmap=plt.cm.CMRmap)
                    p0 = self.ax.scatter(z_ar.flatten(),x_ar.flatten(),y_ar.flatten(),c=y_ar.flatten(),cmap=plt.cm.CMRmap, label="3d plot")
                    #cbar = plt.colorbar(p0,cax=self.ax)
                except:
                    msg = MyMessageBox()
                
                xlabel = r"$\lambda$"
                ylabel = OUTPUT_VARIABLES_LIST[x_data]["symbol"]
                if OUTPUT_VARIABLES_LIST[x_data]["unit"] != "":
                    ylabel = ylabel+" ["+OUTPUT_VARIABLES_LIST[x_data]["unit"]+"]"

                zlabel = OUTPUT_VARIABLES_LIST[y_data]["symbol"]
                if OUTPUT_VARIABLES_LIST[x_data]["unit"] != "":
                    zlabel = zlabel+" ["+OUTPUT_VARIABLES_LIST[y_data]["unit"]+"]"
                
                self.ax.set_xlabel(xlabel)
                self.ax.set_ylabel(ylabel)
                self.ax.set_zlabel(zlabel)
                
            else:
                self.ax = self.figure.add_subplot(111)
                self.ax.plot(x_ar,y_ar,label="2d plot")
                self.ax.legend(numpy.round(self.results["TSR"],2),title=r"$\lambda$")
        else:
            self.ax = self.figure.add_subplot(111)
            self.ax.plot(x_ar,y_ar,label="2d plot")

            self.ax.set_title(OUTPUT_VARIABLES_LIST[y_data]["name"]+" vs. "+OUTPUT_VARIABLES_LIST[x_data]["name"])

            xlabel = OUTPUT_VARIABLES_LIST[x_data]["symbol"]
            if OUTPUT_VARIABLES_LIST[x_data]["unit"] != "":
                xlabel = xlabel+" ["+OUTPUT_VARIABLES_LIST[x_data]["unit"]+"]"

            ylabel = OUTPUT_VARIABLES_LIST[y_data]["symbol"]
            if OUTPUT_VARIABLES_LIST[x_data]["unit"] != "":
                ylabel = ylabel+" ["+OUTPUT_VARIABLES_LIST[y_data]["unit"]+"]"
            
            self.ax.set_xlabel(xlabel)
            self.ax.set_ylabel(ylabel)

        self.canvas.draw()

    def set_data(self,results):
        self.results = results
        self.draw_graph()

__author__ = "Miha Smrekar"
__credits__ = ["Miha Smrekar"]
__license__ = "GPL"
__version__ = "0.3.0"
__maintainer__ = "Miha Smrekar"
__email__ = "miha.smrekar9@gmail.com"
__status__ = "Development"


import numpy
from PyQt5 import QtWidgets
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import QMainWindow, QAction
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt

from cp_curve import max_calculate
from table import Table
from utils import sort_xy, dict_to_ar

#### TO REMOVE
from comparison import TSR_exp, CP_exp


class ResultsWindow(QMainWindow):
    def __init__(self, parent, width, height, results_3d, input_data):
        super(ResultsWindow, self).__init__(parent)
        self.title = "Rezultati"
        # self.left = 0
        # self.top = 0
        self.screen_width = width
        self.screen_height = height
        self.setWindowTitle(self.title)
        self.setGeometry(
            self.screen_width * 0.125,
            self.screen_height * 0.125,
            self.screen_width * 0.75,
            self.screen_height * 0.75,
        )

        self.tab_widget = TabWidget(self)
        self.setCentralWidget(self.tab_widget)

        X, Y, Z = results_3d["v"], results_3d["rpm"], results_3d["power"]
        max_x, max_y, max_z = max_calculate(X, Y, Z)

        array_geom = []
        for r in range(len(input_data["r"])):
            _r = input_data["r"][r]
            _c = input_data["c"][r]
            _theta = input_data["theta"][r]
            _dr = input_data["dr"][r]
            array_geom.append(
                [
                    input_data["r"][r],
                    input_data["c"][r],
                    input_data["theta"][r],
                    input_data["dr"][r],
                ]
            )

        t_geom = Table()
        t_geom.createTable(array_geom)
        t_geom.set_labels(["r", "c", "theta", "dr"])
        self.tab_widget.add_tab_widget(t_geom, "Geometry check")

        # alpha = numpy.linspace(-90, 90, 100)
        # f_c_L = input_data["f_c_L"]
        # f_c_D = input_data["f_c_D"]
        # cL = f_c_L(alpha)
        # cD = f_c_D(alpha)
        #AoA_cL = input_data["AoA_cL"]
        #AoA_cD = input_data["AoA_cD"]
        #cL = input_data["cL"]
        #cD = input_data["cD"]

        #f = self.tab_widget.add_tab_figure("Cl/Cd check")
        #self.tab_widget.add_2d_plot_to_figure(
        #    f, AoA_cL, cL, 121, "cL", "alpha", "cL", look="og"
        #)
        #self.tab_widget.add_2d_plot_to_figure(
        ##    f, AoA_cD, cD, 122, "cD", "alpha", "cD", look="or"
        #)

        f2 = self.tab_widget.add_tab_figure("Moč in Cp")
        self.tab_widget.add_2d_plot_to_figure(
            f2, max_x, max_z, 121, "Moč vs veter", "veter [m/s]", "moč [W]"
        )

        TSR, CP = sort_xy(results_3d["TSR"], results_3d["cp_w"])
        self.tab_widget.add_2d_plot_to_figure(
            f2, TSR, CP, 122, "Cp krivulja", "lambda", "Cp", look="o"
        )
        
        #for comparison Karlsen et al., S826, 0.45m, 3B
        self.tab_widget.add_2d_plot_to_figure(
            f2,TSR_exp,CP_exp, 122, "Cp_exp", "lambda", "Cp",look="-r")
        

        # noinspection PyBroadException
        try:
            if len(X) >= 3 and len(Y) >= 3:
                f3 = self.tab_widget.add_tab_figure("3D moč")
                self.tab_widget.add_surface_plot(
                    f3, X, Y, Z, 111, "Moč (veter,rpm)", "veter[m/s]", "rpm", "moč [W]"
                )
        except:
            print("Could not create 3D surface plot...")

        f4 = self.tab_widget.add_tab_figure("Ct krivulja")
        self.tab_widget.add_2d_plot_to_figure(
            f4,
            results_3d["a"],
            results_3d["Ct"],
            111,
            "Ct(a) krivulja",
            "a",
            "Ct",
            look="o",
            x_min=0,
            x_max=1,
            y_min=0,
            y_max=1
        )

        f5 = self.tab_widget.add_tab_figure("ct/J krivulja (propeler)")
        self.tab_widget.add_2d_plot_to_figure(
            f5,
            results_3d["J"],
            results_3d["ct_p"],
            111,
            "ct_p(J) krivulja",
            "J = 1/lambda",
            "ct",
            look="-",
            x_min=0,
            x_max=1,
            y_min=0,
            y_max=0.1)

        f6 = self.tab_widget.add_tab_figure("cp/J krivulja (propeler)")
        self.tab_widget.add_2d_plot_to_figure(
            f6,
            results_3d['J'],
            results_3d['cp_p'],
            111,
            'cp(J) krivulja',
            'J',
            'cp',
            look='-',
            x_min=0,
            x_max=1,
            y_min=0,
            y_max=1)

        data = dict_to_ar(results_3d)
        t = Table()
        t.createTable(data)

        self.tab_widget.add_tab_widget(t, "Data")

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

    def add_2d_plot_to_figure(
            self,
            f,
            x,
            y,
            whi,
            title=None,
            x_name=None,
            y_name=None,
            x_min=None,
            x_max=None,
            y_min=None,
            y_max=None,
            look=None,
            legend=False,
    ):
        ax = f.add_subplot(whi)
        if look:
            ax.plot(x, y, look)
        else:
            ax.plot(x, y)
        if title:
            plt.title(title)
        if x_name:
            ax.set_xlabel(x_name)
        if y_name:
            ax.set_ylabel(y_name)
        if x_min != None and x_max != None:
            ax.set_xlim(x_min, x_max)
        if y_min != None and y_max != None:
            ax.set_ylim(y_min, y_max)
        self.canvas[-1].draw()

    def add_surface_plot(
            self,
            f,
            x,
            y,
            z,
            whi,
            title=None,
            x_name=None,
            y_name=None,
            z_name=None,
            x_min=None,
            x_max=None,
            y_min=None,
            y_max=None,
            z_min=None,
            z_max=None,
            legend=False,
    ):
        ax = f.add_subplot(whi, projection="3d")
        p0 = ax.plot_trisurf(x, y, z, cmap=plt.cm.CMRmap)
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

    def add_tab_widget(self, widget, tab_name):
        self.tabs.append(widget)
        self.addTab(self.tabs[-1], tab_name)

    def add_runner_widget(self, input_data, tab_name="Runner"):
        # print(input_data)
        r = RunnerWidget(input_data=input_data)
        self.tabs.append(r)
        self.addTab(self.tabs[-1], tab_name)

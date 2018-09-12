import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from PyQt5 import QtGui, QtCore, QtWidgets
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import QMainWindow, QApplication, QPushButton, QWidget, QAction, QTabWidget, QVBoxLayout
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
import random
from table import Table


class MainWindow(QMainWindow):

    def __init__(self, width, height):
        super().__init__()
        self.title = 'BEM analiza v0.1'
        # self.left = 0
        # self.top = 0
        self.screen_width = width
        self.screen_height = height
        self.setWindowTitle(self.title)
        self.setGeometry(self.screen_width * 0.125, self.screen_height * 0.125, self.screen_width * 0.75,
                         self.screen_height * 0.75)

        self.tab_widget = TabWidget(self)
        self.setCentralWidget(self.tab_widget)

        self.set_menubar()

        self.show()

    def set_menubar(self):
        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu('File')
        exitButton = QAction(QIcon('exit24.png'), 'Exit', self)
        exitButton.setShortcut('Ctrl+Q')
        exitButton.setStatusTip('Exit application')
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
                              y_min=None, y_max=None, look=None, legend=False):
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

    def add_surface_plot(self, f, x, y, z, whi, title=None, x_name=None, y_name=None, z_name=None, x_min=None,
                         x_max=None, y_min=None, y_max=None, z_min=None, z_max=None, legend=False):
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

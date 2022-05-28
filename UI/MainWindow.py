import json
import os
from multiprocessing import Manager

import numpy as np
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import QMainWindow, QAction, QFileDialog, QApplication

from UI.AirfoilManager import AirfoilManager
from UI.Analysis import Analysis
from UI.Optimization import Optimization
from UI.WindTurbineProperties import WindTurbineProperties
from UI.helpers import ThreadGetter, TabWidget
from main import TITLE_STR, application_path
from turbine_data import SET_INIT
from utils import create_folder, fltr, ErrorMessageBox


class MainWindow(QMainWindow):
    """
    Main class that sets up the UI layout.
    All QWidgets within the MainWindow class are displayed using the custom TabWidget.
    The MainWindow class also holds the functions for saving and loading data from files.

    The MainWindow class stores references to all other subclasses used in the program.
    There are four subclasses in MainWindow called AirfoilManager, WindTurbineProperties, Analysis, and ThreadGetter.
    """
    emitter_add = pyqtSignal(str)
    emitter_done = pyqtSignal()

    def __init__(self, width, height):
        super().__init__()

        self.screen_width = width
        self.screen_height = height

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

        # self.setGeometry(width * 0.125, height * 0.125, width * 0.75, height * 0.75)
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
        self.tab_widget.add_tab(self.optimization, "Optimization",
                                "Optimizacija lopatice s pomočjo algoritma diferencialne evolucije")

        self.running = False
        self.manager = Manager()
        self.set_all_settings(SET_INIT)

        create_folder(os.path.join(application_path, "foils"))  # Used by XFoil

        self.set_process_stopped()

        self.show()

    def set_title(self):
        """
        Sets the title of the UI window, based on the name of the wind turbine.
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
        name = QFileDialog.getSaveFileName(self, 'Save File', "", "BEM (*.bem)")[0]
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
        file_path = QFileDialog.getOpenFileName(self, "Load File", "", "BEM (*.bem)")[0]
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
            # pprint(out)
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

    def wheelEvent(self, event):
        """

        :param event:
        """
        modifiers = QtWidgets.QApplication.keyboardModifiers()
        if modifiers == QtCore.Qt.ControlModifier:
            font = QApplication.instance().font()
            size = font.pointSizeF()
            if event.angleDelta().y() > 0:
                size = size + 1
            else:
                size = size - 1
            font.setPointSize(size)
            QApplication.instance().setFont(font)
            for w in QApplication.allWidgets():
                w.setFont(font)

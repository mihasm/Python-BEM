import numpy as np
import pyqtgraph as pg
from PyQt5.QtWidgets import QWidget, QGridLayout, QFormLayout, QScrollArea, QVBoxLayout, QPushButton, QLineEdit, \
    QCheckBox, QLabel, QComboBox, QFileDialog
from matplotlib import cm
from matplotlib import pyplot as plt

from UI.Curve import Curve
from UI.CurveEditor import CurveEditor
from UI.CurveViewer import CurveViewer
from UI.Curves import Curves
from UI.Table import Table
from UI.helpers import XFoilThread, ScrapeThread, MatplotlibWindow, PrintoutWindow, XfoilOptionsWindow, MyMessageBox, \
    ErrorMessageBox
from utils import import_dat, import_nrel_dat, interp_at, generate_dat, to_float, get_centroid_coordinates


class Airfoils(QWidget):
    """

    """
    def __init__(self, airfoil_name, parent=None):
        super(Airfoils, self).__init__(parent)

        self.thread = None
        self.xfoil_window = None

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

        self.table_dat = Table(True)
        self.table_dat.createEmpty(2, 50)
        self.table_dat.set_labels(["x", "y"])
        self.grid.addWidget(self.table_dat, 1, 2)

        self.plt = plt.figure(figsize=(3, 3))
        self.ax = self.plt.add_subplot(111)

        self.airfoil_graph = pg.PlotWidget()
        self.grid.addWidget(self.airfoil_graph, 1, 3)
        self.airfoil_graph.getViewBox().setAspectLocked(lock=True, ratio=1)

        self.buttonRefresh = QPushButton("Refresh curve")
        self.grid.addWidget(self.buttonRefresh, 3, 3)
        self.buttonRefresh.setToolTip("Refresh graph curve based on table data")
        self.buttonRefresh.clicked.connect(self.refresh)

        self.fbox.addRow(QLabel("———— 1. Import or generate ————"))

        self.link = QLineEdit("link (airfoiltools.com)")
        self.fbox.addRow(self.link)
        self.link.setToolTip(
            "Link to the airfoiltools.com profile web page")

        self.button_generate_curves_link = QPushButton("Scrape curves from link")
        self.button_generate_curves_link.clicked.connect(self.generate_curves_link)
        self.fbox.addRow(self.button_generate_curves_link)
        self.button_generate_curves_link.setToolTip(
            "Download cL/cD/profile from airfoiltools.com")

        self.button_generate_curves_xfoil = QPushButton("Generate XFOIL curves")
        self.button_generate_curves_xfoil.clicked.connect(self.open_xfoil_options)
        self.fbox.addRow(self.button_generate_curves_xfoil)
        self.button_generate_curves_xfoil.setToolTip(
            "Generate cL/cD using profile xy points with XFOIL")

        self.button_import_dat_from_file = QPushButton("Import .dat")
        self.button_import_dat_from_file.clicked.connect(self.dat_importer)
        self.fbox.addRow(self.button_import_dat_from_file)
        self.button_import_dat_from_file.setToolTip("Import .dat file")

        self.button_import_nrel_dat_from_file = QPushButton("Import NREL .dat")
        self.button_import_nrel_dat_from_file.clicked.connect(self.nrel_dat_importer)
        self.fbox.addRow(self.button_import_nrel_dat_from_file)
        self.button_import_nrel_dat_from_file.setToolTip("Import NREL .dat file (cl and cd curves)")

        self.fbox.addRow(QLabel("———— 2. Prepare ————"))

        self.button_open_viewer = QPushButton("Curve Extrapolation Editor")
        self.button_open_viewer.clicked.connect(self.open_viewer)
        self.fbox.addRow(self.button_open_viewer)
        self.button_open_viewer.setToolTip(
            "Modify extrapolation coefficients for Montgomerie method")

        self.button_curve_editor = QPushButton("Curve Data Editor")
        self.button_curve_editor.clicked.connect(self.open_curve_editor)
        self.fbox.addRow(self.button_curve_editor)
        self.button_curve_editor.setToolTip("Manually edit cL/cD curve data")

        self.fbox.addRow(QLabel("———— 3. Verify ————"))

        self.button_visualize = QPushButton("Show cL/cD(Re,alpha) curves")
        self.button_visualize.clicked.connect(self.visualize)
        self.fbox.addRow(self.button_visualize)
        self.button_visualize.setToolTip(
            "Show 3D graph of extrapolated curves")

        self.fbox.addRow(QLabel("———— 4. Analysis settings ————"))

        self.extrapolation_bool = QCheckBox("Use extrapolation")
        self.fbox.addRow(self.extrapolation_bool)
        self.extrapolation_bool.setToolTip(
            "Toggle usage of extrapolated curves in the analysis")

        self.ncrit_selection = QComboBox()
        self._ncrit_selection = QLabel("Ncrit")
        self.fbox.addRow(self._ncrit_selection, self.ncrit_selection)
        self.ncrit_selection.setToolTip(
            "Use curves with this N value for calculation (e^N transition theory)")

        self.centroid_widget = QWidget()
        self.centroid_grid = QGridLayout()
        self.centroid_widget.setLayout(self.centroid_grid)
        self.centroid_label = QLabel("Center:")
        self.grid.addWidget(self.centroid_widget, 2, 2)
        self.centroid_label.setToolTip(
            "Center of rotation for 3D export to Solidworks. Doesn't affect analysis.")

        self.centroid_x_edit = QLineEdit()
        self.centroid_y_edit = QLineEdit()

        self.centroid_grid.addWidget(self.centroid_label, 1, 1)
        self.centroid_grid.addWidget(self.centroid_x_edit, 1, 2)
        self.centroid_grid.addWidget(self.centroid_y_edit, 1, 3)

        self.get_centroid_button = QPushButton("Calculate centroid")
        self.get_centroid_button.clicked.connect(self.calculate_centroid)
        self.grid.addWidget(self.get_centroid_button, 3, 2)
        self.get_centroid_button.setToolTip(
            "Recalculate the centroid for the profile")

        self.grid.setColumnStretch(1, 1)
        self.grid.setColumnStretch(2, 1)
        self.grid.setColumnStretch(3, 2)

        self.window = None
        self.curve_editor = CurveEditor(self)

    def dat_importer(self):
        """
        Loads the wind turbine data from a file. Also clears the calculation text areas and sets the appropriate title.
        """
        file_path = QFileDialog.getOpenFileName(self, "Import .dat file", "", "dat (*.dat)")[0]
        if file_path != "":
            x, y = import_dat(file_path)
            self.set_x_y(x, y)
            self.refresh()

    def nrel_dat_importer(self):
        """
        Loads the wind turbine data from a file. Also clears the calculation text areas and sets the appropriate title.
        """
        file_path = QFileDialog.getOpenFileName(self, "Import .dat file", "", "dat (*.dat)")[0]
        if file_path != "":
            data = import_nrel_dat(file_path)
            self.populate_curve_list(data)

    def visualize(self):
        """

        """
        print("Visualizing")
        data = self.gather_curves()
        print(data)
        data = data[np.in1d(data[:, 1], float(self.ncrit_selection.currentText()))]  # current Ncrit

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
        """

        """
        print("opening viewer")
        self.viewer.show()
        self.viewer.generate_views()

    def open_curve_editor(self):
        """

        """
        print("opening curve editor")
        self.curve_editor.show()
        self.curve_editor.load_curves()

    def open_xfoil_options(self):
        """

        """
        if self.xfoil_window == None:
            self.xfoil_window = XfoilOptionsWindow(self)
            self.xfoil_window.setWindowTitle("XFOIL runner: " + self.airfoil_name)
        else:
            self.xfoil_window.show()
            self.xfoil_window.activateWindow()

    def generate_curves_xfoil(self):
        """

        """
        print("Generating xfoil curves")
        x, y = self.get_x_y()
        generate_dat(self.airfoil_name, x, y)

        dat_path = self.airfoil_name + ".dat"
        alpha_from = to_float(self.xfoil_window.alpha_from.text())
        alpha_to = to_float(self.xfoil_window.alpha_to.text())
        alpha_num = int(self.xfoil_window.alpha_num.text())
        reynolds_from = to_float(self.xfoil_window.reynolds_from.text())
        reynolds_to = to_float(self.xfoil_window.reynolds_to.text())
        reynolds_num = int(self.xfoil_window.reynolds_num.text())
        ncrit = to_float(self.xfoil_window.ncrit.text())

        if self.thread != None:
            self.xfoil_window.button_run_xfoil.setDisabled(True)
            self.xfoil_window.button_stop_xfoil.setEnabled(True)
            msg = MyMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText("XFoil thread still running! Press Stop again.")
            msg.setDetailedText("XFoil thread was started and was not set back to None. Try pressing Stop or saving and restarting.")
            msg.exec_()
            return

        if self.window != None:
            self.window.close()

        self.window = PrintoutWindow(self)
        self.thread = XFoilThread(self)
        self.thread.set_params(dat_path,
                               alpha_from, alpha_to, alpha_num,
                               reynolds_from, reynolds_to, reynolds_num,
                               ncrit)

        self.xfoil_window.button_run_xfoil.setDisabled(True)
        self.xfoil_window.button_stop_xfoil.setEnabled(True)
        self.thread.completeSignal.connect(self.xfoil_completion)
        self.thread.start()

    def stop_xfoil(self):
        """

        """
        self.thread.terminate()
        self.xfoil_window.button_run_xfoil.setEnabled(True)
        self.xfoil_window.button_stop_xfoil.setDisabled(True)
        self.thread = None

    def xfoil_completion(self, nothing_important):
        """

        :param nothing_important:
        """
        self.populate_curve_list(self.xfoil_generated_data)
        self.refresh()
        self.xfoil_window.button_run_xfoil.setEnabled(True)
        self.xfoil_window.button_stop_xfoil.setDisabled(True)
        self.thread = None

    def generate_curves_link(self):
        """

        """
        print("Scraping from link...")
        if self.window != None:
            self.window.close()
        print("Open window")
        self.window = PrintoutWindow(self)
        print("Thread")
        self.thread = ScrapeThread(self)
        print("Params")
        self.thread.set_params(self.link)
        print("Connecting")
        self.thread.completeSignal.connect(self.generate_curves_link_finished)
        print("Starting")
        self.thread.start()

    def generate_curves_link_finished(self, nothing_important):
        """

        :param nothing_important:
        """
        print("Finished, populating...")
        self.table_dat.clear_table()
        self.populate_curve_list(self.scraping_generated_data[0])
        self.set_x_y(self.scraping_generated_data[1], self.scraping_generated_data[2])
        self.refresh()

    def populate_curve_list(self, data):
        """

        :param data:
        """
        print(data)
        self.curves.curve_list = []
        x, y = self.get_x_y()
        ncrit_list = np.unique(data[:, 1])
        for ncrit_selected in ncrit_list:
            rows_with_ncrit = data[np.in1d(data[:, 1], ncrit_selected)]
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
        """

        :return:
        """
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

        self.airfoil_graph.plot(x_values, y_values)

        try:
            centroid_x = float(self.centroid_x_edit.text())
            centroid_y = float(self.centroid_y_edit.text())
            self.ax.plot(centroid_x, centroid_y, "r+")
        except:
            msg = ErrorMessageBox()

        self.refresh_ncrits_combobox()

    def get_max_thickness(self):
        """

        :return:
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
        if len(y) > 0:
            y_max = np.max(y)
            y_min = np.min(y)
            thickness = (abs(y_max) + abs(y_min)) / 1
            return thickness
        return 0.1  # default value

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

    def set_x_y(self, x, y):
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
        """

        :return:
        """
        foil_x, foil_y = self.get_x_y()
        x, y = get_centroid_coordinates(foil_x, foil_y)
        # print("x:",x,"y:",y)
        self.centroid_x_edit.setText(str(round(x, 6)))
        self.centroid_y_edit.setText(str(round(y, 6)))
        return x, y

    def get_ncrits(self):
        """

        :return:
        """
        curves = self.gather_curves()
        if len(curves) > 0:
            ncrit_list = np.unique(curves[:, 1])
            return ncrit_list

    def refresh_ncrits_combobox(self):
        """

        """
        self.ncrit_selection.clear()
        ncrits = self.get_ncrits()
        if ncrits is not None:
            self.ncrit_selection.addItems([str(n) for n in list(ncrits)])

    def gather_curves(self):
        """

        :return:
        """
        return self.curves.gather_curves(self.extrapolation_bool.checkState())

    def get_settings(self):
        """

        :return:
        """
        out = {}

        x, y = self.get_x_y()
        try:
            centroid_x = float(self.centroid_x_edit.text())
            centroid_y = float(self.centroid_y_edit.text())
        except:
            centroid_x, centroid_y = 0.0, 0.0

        try:
            ncrit_selected = float(self.ncrit_selection.currentText())
        except:
            ncrit_selected = 0.0

        try:
            extrapolation_bool = self.extrapolation_bool.checkState()
        except:
            extrapolation_bool = True
        out = {"x": x,
               "y": y,
               "max_thickness": self.get_max_thickness(),
               "link": self.link.text(),
               "interp_function_cl": self.interp_function_cl,
               "interp_function_cd": self.interp_function_cd,
               "curves": self.curves.save_curves(),
               "gathered_curves": self.gather_curves(),
               "centroid_x": centroid_x,
               "centroid_y": centroid_y,
               "ncrit_selected": ncrit_selected,
               "extrapolation_bool": extrapolation_bool,
               "stall_angles": self.curves.get_stall_angles()}
        return out

    def set_settings(self, dict_settings):
        """

        :param dict_settings:
        """

        if "x" in dict_settings and "y" in dict_settings:
            x, y = dict_settings["x"], dict_settings["y"]
            self.set_x_y(x, y)
        
        if "link" in dict_settings.keys():
            self.link.setText(dict_settings["link"])
        
        if "curves" in dict_settings.keys():
            self.curves.load_curves(dict_settings["curves"])
        
        if "extrapolation_bool" in dict_settings.keys():
            self.extrapolation_bool.setChecked(dict_settings["extrapolation_bool"])
        
        self.refresh()

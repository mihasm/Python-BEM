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

from multiprocessing import Manager
from multiprocessing.context import Process

import numpy as np
from PyQt5 import QtCore, QtGui
from PyQt5.QtWidgets import QWidget, QGridLayout, QFormLayout, QScrollArea, QVBoxLayout, QTextEdit, QPushButton, \
    QComboBox, QCheckBox, QLineEdit, QLabel, QMessageBox

from UI.helpers import PyQtGraphWindow, OptimizationPreviewWindow, MyMessageBox
from UI.Table import Table
from optimization import optimize
from utils import to_float


class OptimizationResultWindow(QWidget):
    def __init__(self, optimization_widget, geometry, summary_lines):
        super(OptimizationResultWindow, self).__init__(None)
        self.optimization_widget = optimization_widget
        self.main = optimization_widget.main
        self.geometry = geometry

        self.setWindowFlags(QtCore.Qt.Window)
        self.setWindowTitle("Optimized geometry")
        self.setMinimumSize(720, 520)

        layout = QGridLayout()
        self.setLayout(layout)

        self.summary = QTextEdit()
        self.summary.setReadOnly(True)
        self.summary.setPlainText("\n".join(summary_lines))
        layout.addWidget(self.summary, 0, 0)

        self.table = Table(True)
        self.table.set_labels(["r [m]", "c [m]", "theta [deg]", "airfoil"])
        rows = []
        for r, c, theta, foil in zip(geometry["r_in"], geometry["c_in"], geometry["theta_in"], geometry["foils_in"]):
            rows.append([round(float(r), 8), round(float(c), 8), round(float(theta), 8), str(foil)])
        self.table.createTable(rows)
        layout.addWidget(self.table, 1, 0)

        self.button_import = QPushButton("Import this geometry back into the software")
        self.button_import.clicked.connect(self.import_geometry)
        layout.addWidget(self.button_import, 2, 0)

        self.show()
        self.raise_()
        self.activateWindow()

    def import_geometry(self):
        settings = self.main.get_all_settings()
        if settings is None:
            return
        settings.update(self.geometry)
        if self.geometry.get("pitch") is not None:
            settings["pitch"] = self.geometry.get("pitch")
        self.main.set_all_settings(settings)
        self.main.analysis.clear()
        self.optimization_widget.clear()
        self.main.tab_widget.setCurrentIndex(1)
        QMessageBox.information(self, "Geometry imported", "The optimized geometry was imported into Turbine info.")


class Optimization(QWidget):
    """

    """
    def __init__(self, parent=None):
        super(Optimization, self).__init__(parent)
        self.main = self.parent()

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

        self.right_output_text_area = QTextEdit()
        self.right_output_text_area.setReadOnly(True)
        self.grid.addWidget(self.right_output_text_area, 1, 2)

        self.validator = QtGui.QDoubleValidator()

        self.form_list = []

        self.target_speed = QLineEdit()
        self.target_speed.setValidator(self.validator)
        self.target_speed.textChanged.connect(self.check_state)
        self.target_speed.textChanged.emit(self.target_speed.text())
        self._target_speed = QLabel("Target speed [m/s]")
        self.form_list.append([self._target_speed, self.target_speed])
        self.target_speed.setToolTip("Ciljna hitrost vetra. Pri tej hitrosti se bo izvajala BEM analiza.")

        self.target_rpm = QLineEdit()
        self.target_rpm.setValidator(self.validator)
        self.target_rpm.textChanged.connect(self.check_state)
        self.target_rpm.textChanged.emit(self.target_rpm.text())
        self._target_rpm = QLabel("Target rpm [RPM]")
        self.form_list.append([self._target_rpm, self.target_rpm])
        self.target_rpm.setToolTip("Ciljni vrtljaji rotorja turbine. Pri teh vrtljajih se bo izvajala BEM analiza.")

        self.mut_coeff = QLineEdit()
        self.mut_coeff.setValidator(self.validator)
        self.mut_coeff.textChanged.connect(self.check_state)
        self.mut_coeff.textChanged.emit(self.mut_coeff.text())
        self._mut_coeff = QLabel("Mutation coefficient")
        self.form_list.append([self._mut_coeff, self.mut_coeff])
        self.mut_coeff.setToolTip(
            "Mutacijski koeficient nastavlja jakost naključnih mutacij, ki se zgodijo pri vsaki novi generaciji (iteraciji).")

        self.population = QLineEdit()
        self.population.setValidator(self.validator)
        self.population.textChanged.connect(self.check_state)
        self.population.textChanged.emit(self.population.text())
        self._population = QLabel("Population")
        self.form_list.append([self._population, self.population])
        self.population.setToolTip("Število posameznikov v algoritmu diferencialne evolucije.")

        self.num_iter = QLineEdit()
        self.num_iter.setValidator(self.validator)
        self.num_iter.textChanged.connect(self.check_state)
        self.num_iter.textChanged.emit(self.num_iter.text())
        self._num_iter = QLabel("Number of iterations")
        self.form_list.append([self._num_iter, self.num_iter])
        self.num_iter.setToolTip(
            "Število generacij (iteracij). Konvergenčni kriterij je pri tovrstnih algoritmih težko določljiv, zato izberemo fiksno vrednost.")

        self.optimization_mode_selection = QComboBox()
        self.optimization_mode_selection.addItems(["Section-by-section", "Whole-blade"])
        self.optimization_mode_selection.currentIndexChanged.connect(self.update_mode_visibility)
        self._optimization_mode = QLabel("Optimization mode")
        self.form_list.append([self._optimization_mode, self.optimization_mode_selection])

        self.whole_blade_objective_selection = QComboBox()
        self.whole_blade_objective_selection.addItems([
            "Min torque at target thrust",
            "Max thrust at power limit",
            "Max Cp",
            "Max efficiency at target thrust",
            "Weighted custom objective",
        ])
        self.whole_blade_objective_selection.currentIndexChanged.connect(self.update_mode_visibility)
        self._whole_blade_objective = QLabel("Whole-blade objective")
        self.form_list.append([self._whole_blade_objective, self.whole_blade_objective_selection])

        self.whole_blade_solver_selection = QComboBox()
        self.whole_blade_solver_selection.addItems(["SLSQP", "Differential evolution"])
        self._whole_blade_solver = QLabel("Whole-blade solver")
        self.form_list.append([self._whole_blade_solver, self.whole_blade_solver_selection])

        self.whole_blade_twist_cp = QLineEdit("5")
        self.whole_blade_twist_cp.setValidator(self.validator)
        self._whole_blade_twist_cp = QLabel("Twist control points")
        self.form_list.append([self._whole_blade_twist_cp, self.whole_blade_twist_cp])

        self.whole_blade_chord_cp = QLineEdit("4")
        self.whole_blade_chord_cp.setValidator(self.validator)
        self._whole_blade_chord_cp = QLabel("Chord control points")
        self.form_list.append([self._whole_blade_chord_cp, self.whole_blade_chord_cp])

        self.whole_blade_root_cp_lock = QLineEdit("1")
        self.whole_blade_root_cp_lock.setValidator(self.validator)
        self._whole_blade_root_cp_lock = QLabel("Locked root CPs")
        self.form_list.append([self._whole_blade_root_cp_lock, self.whole_blade_root_cp_lock])

        self.whole_blade_tip_cp_lock = QLineEdit("1")
        self.whole_blade_tip_cp_lock.setValidator(self.validator)
        self._whole_blade_tip_cp_lock = QLabel("Locked tip CPs")
        self.form_list.append([self._whole_blade_tip_cp_lock, self.whole_blade_tip_cp_lock])

        self.whole_blade_target_thrust = QLineEdit()
        self.whole_blade_target_thrust.setValidator(self.validator)
        self._whole_blade_target_thrust = QLabel("Target thrust [N]")
        self.form_list.append([self._whole_blade_target_thrust, self.whole_blade_target_thrust])

        self.whole_blade_max_power = QLineEdit()
        self.whole_blade_max_power.setValidator(self.validator)
        self._whole_blade_max_power = QLabel("Max power [W]")
        self.form_list.append([self._whole_blade_max_power, self.whole_blade_max_power])

        self.whole_blade_max_stall_percentage = QLineEdit()
        self.whole_blade_max_stall_percentage.setValidator(self.validator)
        self._whole_blade_max_stall_percentage = QLabel("Max stall fraction")
        self.form_list.append([self._whole_blade_max_stall_percentage, self.whole_blade_max_stall_percentage])

        self.whole_blade_monotonic_chord = QCheckBox()
        self.whole_blade_monotonic_chord.setChecked(True)
        self._whole_blade_monotonic_chord = QLabel("Monotonic chord after max")
        self.form_list.append([self._whole_blade_monotonic_chord, self.whole_blade_monotonic_chord])

        self.input_variables_widget = QWidget()
        self.input_variables_layout = QGridLayout()
        self.input_variables_widget.setLayout(self.input_variables_layout)
        self.button_add_input_variable = QPushButton("Add input variable")
        self.button_add_input_variable.clicked.connect(self.add_input_variable)
        self.input_variable_selection = QComboBox()
        self.input_variable_selection.addItems(["_theta", "_c", "pitch"])
        self.list_input_variables = []

        self.output_variables_widget = QWidget()
        self.output_variables_layout = QGridLayout()
        self.output_variables_widget.setLayout(self.output_variables_layout)
        self.button_add_output_variable = QPushButton("Add output variable")
        self.button_add_output_variable.clicked.connect(self.add_output_variable)
        self.output_variable_selection = QComboBox()
        self.output_variable_selection.addItems([
            "dQ", "dT", "a", "a'", "Cl", "Cd", "dFn", "dFt", "U4", "alpha", "phi",
            "thrust", "power", "cp", "ct", "cq", "eff", "M", "Q", "blade_stall_percentage",
        ])
        self.list_output_variables = []

        self.target_variables_widget = QWidget()
        self.target_variables_layout = QGridLayout()
        self.target_variables_widget.setLayout(self.target_variables_layout)
        self.button_add_target_variable = QPushButton("Add target variable")
        self.button_add_target_variable.clicked.connect(self.add_target_variable)
        self.target_variable_selection = QComboBox()
        self.target_variable_selection.addItems([
            "dQ", "dT", "a", "a'", "Cl", "Cd", "dFn", "dFt", "U4", "alpha", "phi",
            "thrust", "power", "cp", "ct", "cq", "eff", "M", "Q", "blade_stall_percentage",
        ])
        self.list_target_variables = []

        self.buttonOptimization = QPushButton("Run optimization")
        self.buttonOptimization.clicked.connect(self.run)

        self.buttonStop = QPushButton("Stop")
        self.buttonStop.clicked.connect(self.terminate)

        self.buttonClear = QPushButton("Clear screen")
        self.buttonClear.clicked.connect(self.clear)

        self.buttonEOF = QCheckBox()
        self.buttonEOF.setChecked(True)
        self.buttonEOFdescription = QLabel("Scroll to end of screen")

        for a, b in self.form_list:
            self.fbox.addRow(a, b)

        self.fbox.addRow(self.button_add_input_variable, self.input_variable_selection)
        self.fbox.addRow(self.input_variables_widget)

        self.fbox.addRow(self.button_add_output_variable, self.output_variable_selection)
        self.fbox.addRow(self.output_variables_widget)

        self.fbox.addRow(self.button_add_target_variable, self.target_variable_selection)
        self.fbox.addRow(self.target_variables_widget)

        self.fbox.addRow(self.buttonOptimization)
        self.fbox.addRow("", QLabel())

        self.fbox.addRow(self.buttonClear, self.buttonStop)
        self.fbox.addRow(self.buttonEOFdescription, self.buttonEOF)

        self.win = PyQtGraphWindow(self)
        self.win.setWindowTitle("Live Optimization Visualizer")
        self.preview_win = None
        self.manager_pyqtgraph = Manager()
        self.queue_pyqtgraph = self.manager_pyqtgraph.list()
        self.queue_pyqtgraph.append([[0], [0], 0, 0, None])

        self.tsr_string = QLabel("0")
        self.J_string = QLabel("0")
        self.fbox.addRow("TSR:", self.tsr_string)
        self.fbox.addRow("J:", self.J_string)

        self.whole_blade_widgets = [
            self._whole_blade_objective, self.whole_blade_objective_selection,
            self._whole_blade_solver, self.whole_blade_solver_selection,
            self._whole_blade_twist_cp, self.whole_blade_twist_cp,
            self._whole_blade_chord_cp, self.whole_blade_chord_cp,
            self._whole_blade_root_cp_lock, self.whole_blade_root_cp_lock,
            self._whole_blade_tip_cp_lock, self.whole_blade_tip_cp_lock,
            self._whole_blade_target_thrust, self.whole_blade_target_thrust,
            self._whole_blade_max_power, self.whole_blade_max_power,
            self._whole_blade_max_stall_percentage, self.whole_blade_max_stall_percentage,
            self._whole_blade_monotonic_chord, self.whole_blade_monotonic_chord,
        ]
        self.section_objective_widgets = [
            self.button_add_output_variable, self.output_variable_selection, self.output_variables_widget,
            self.button_add_target_variable, self.target_variable_selection, self.target_variables_widget,
        ]

        self.result_window = None
        self.update_mode_visibility()

    def refresh_input_variables(self):
        """

        """
        for i in reversed(range(self.input_variables_layout.count())):
            self.input_variables_layout.itemAt(i).widget().setParent(None)

        i = 0
        for var, min_b, max_b in self.list_input_variables:
            name = QLabel(var)
            min_bound = QLineEdit(str(min_b))
            min_bound.textChanged.connect(self.update_input_variables_list)
            max_bound = QLineEdit(str(max_b))
            max_bound.textChanged.connect(self.update_input_variables_list)
            delete_button = QPushButton("X")
            delete_button.clicked.connect(self.delete_input_variable)
            self.input_variables_layout.addWidget(name, i, 0)
            self.input_variables_layout.addWidget(min_bound, i, 1)
            self.input_variables_layout.addWidget(max_bound, i, 2)
            self.input_variables_layout.addWidget(delete_button, i, 3)
            i += 1

    def refresh_output_variables(self):
        """

        """
        for i in reversed(range(self.output_variables_layout.count())):
            self.output_variables_layout.itemAt(i).widget().setParent(None)

        i = 0
        for var, coefficient in self.list_output_variables:
            name = QLabel(var)
            coefficient = QLineEdit(str(coefficient))
            coefficient.textChanged.connect(self.update_output_variables_list)
            delete_button = QPushButton("X")
            delete_button.clicked.connect(self.delete_output_variable)
            self.output_variables_layout.addWidget(name, i, 0)
            self.output_variables_layout.addWidget(coefficient, i, 1)
            self.output_variables_layout.addWidget(delete_button, i, 2)
            i += 1

    def refresh_target_variables(self):
        """

        """
        for i in reversed(range(self.target_variables_layout.count())):
            self.target_variables_layout.itemAt(i).widget().setParent(None)

        i = 0
        for var, target_value, coefficient in self.list_target_variables:
            name = QLabel(var)

            target_value = QLineEdit(str(target_value))
            target_value.textChanged.connect(self.update_target_variables_list)

            coefficient = QLineEdit(str(coefficient))
            coefficient.textChanged.connect(self.update_target_variables_list)

            delete_button = QPushButton("X")
            delete_button.clicked.connect(self.delete_target_variable)

            self.target_variables_layout.addWidget(name, i, 0)
            self.target_variables_layout.addWidget(target_value, i, 1)
            self.target_variables_layout.addWidget(coefficient, i, 2)
            self.target_variables_layout.addWidget(delete_button, i, 3)
            i += 1

    def update_input_variables_list(self):
        """

        """
        for row in range(len(self.list_input_variables)):
            try:
                min_b = float(self.input_variables_layout.itemAtPosition(row, 1).widget().text())
                self.list_input_variables[row][1] = min_b
            except ValueError:
                pass
            try:
                max_b = float(self.input_variables_layout.itemAtPosition(row, 2).widget().text())
                self.list_input_variables[row][2] = max_b
            except ValueError:
                pass

    def update_output_variables_list(self):
        """

        """
        for row in range(len(self.list_output_variables)):
            try:
                coeff = float(self.output_variables_layout.itemAtPosition(row, 1).widget().text())
                self.list_output_variables[row][1] = coeff
            except ValueError:
                pass

    def update_target_variables_list(self):
        """

        """
        for row in range(len(self.list_target_variables)):
            try:
                target_value = float(self.target_variables_layout.itemAtPosition(row, 1).widget().text())
                self.list_target_variables[row][1] = target_value
            except ValueError:
                pass
            try:
                coeff = float(self.target_variables_layout.itemAtPosition(row, 2).widget().text())
                self.list_target_variables[row][2] = coeff
            except ValueError:
                pass

    def delete_input_variable(self):
        """

        """
        button = self.sender()
        index = self.input_variables_layout.indexOf(button)
        row, column, _, _ = self.input_variables_layout.getItemPosition(index)
        del self.list_input_variables[row]
        self.refresh_input_variables()

    def delete_output_variable(self):
        """

        """
        button = self.sender()
        index = self.output_variables_layout.indexOf(button)
        row, column, _, _ = self.output_variables_layout.getItemPosition(index)
        del self.list_output_variables[row]
        self.refresh_output_variables()

    def delete_target_variable(self):
        """

        """
        button = self.sender()
        index = self.target_variables_layout.indexOf(button)
        row, column, _, _ = self.target_variables_layout.getItemPosition(index)
        del self.list_target_variables[row]
        self.refresh_target_variables()

    def add_input_variable(self):
        """

        """
        variable = str(self.input_variable_selection.currentText())
        default_bounds = {
            "_theta": [-45, 45],
            "_c": [0.001, 5],
            "pitch": [-45, 45],
        }
        bounds = default_bounds.get(variable, [0, 0])
        self.list_input_variables.append([variable, bounds[0], bounds[1]])
        self.refresh_input_variables()

    def add_output_variable(self):
        """

        """
        variable = str(self.output_variable_selection.currentText())
        self.list_output_variables.append([variable, 0])
        self.refresh_output_variables()

    def add_target_variable(self):
        """

        """
        variable = str(self.target_variable_selection.currentText())
        self.list_target_variables.append([variable, 0, 1])
        self.refresh_target_variables()

    def set_widgets_visible(self, widgets, visible):
        for widget in widgets:
            widget.setVisible(visible)

    def update_mode_visibility(self):
        is_whole_blade = self.optimization_mode_selection.currentIndex() == 1
        objective_index = self.whole_blade_objective_selection.currentIndex()
        uses_custom_objective = objective_index == 4

        self.set_widgets_visible(self.whole_blade_widgets, is_whole_blade)
        self.set_widgets_visible(self.section_objective_widgets, (not is_whole_blade) or uses_custom_objective)

        show_target_thrust = is_whole_blade and objective_index in [0, 3]
        self.set_widgets_visible([self._whole_blade_target_thrust, self.whole_blade_target_thrust], show_target_thrust)

        show_max_power = is_whole_blade and objective_index == 1
        self.set_widgets_visible([self._whole_blade_max_power, self.whole_blade_max_power], show_max_power)

        if is_whole_blade:
            self.button_add_input_variable.setText("Add design variable")
        else:
            self.button_add_input_variable.setText("Add input variable")

    def update_tsr_and_j(self):
        """

        """
        try:
            s = self.get_settings()
            R = float(self.main.wind_turbine_properties.R.text()) * float(
                self.main.wind_turbine_properties.geometry_scale.text())
            tsr = 2 * np.pi * float(s["target_rpm"]) * R / 60 / float(s["target_speed"])
            self.tsr_string.setText("%.2f" % tsr)
            J = float(s["target_speed"]) / (float(s["target_rpm"]) / 60 * 2 * R)
            self.J_string.setText("%.2f" % J)
        except:
            pass

    def check_forms_angles(self):
        """

        :return:
        """
        out = ""
        _needed_vars = [[self._target_speed, self.target_speed], [self._target_rpm, self.target_rpm], ]
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
        """

        :param args:
        :param kwargs:
        """
        self.update_tsr_and_j()

        sender = self.sender()
        validator = sender.validator()
        state = validator.validate(sender.text(), 0)[0]
        if state == QtGui.QValidator.Acceptable:
            color = "#edf5e1"  # green
        elif state == QtGui.QValidator.Intermediate:
            color = "#fff79a"  # yellow
        else:
            color = "#f6989d"  # red
        # sender.setStyleSheet("QLineEdit { background-color: %s; color: #000000 }" % color)

    def validate_inputs(self):
        """

        :return:
        """
        check = self.check_forms_angles()
        check_analysis = self.main.analysis.check_forms()
        if check != True or check_analysis != True:
            if check == True:
                check = ""
            if check_analysis == True:
                check_analysis = ""
            check = check + check_analysis
            msg = MyMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText("Input validation error")
            msg.setDetailedText(check)
            msg.exec_()
            return False
        return True

    def clear(self):
        """

        """
        self.right_output_text_area.clear()

    def add_text(self, string):
        """

        :param string:
        """
        self.right_output_text_area.insertPlainText(string)
        if self.buttonEOF.checkState() == 2:
            self.right_output_text_area.moveCursor(QtGui.QTextCursor.End)

    def run(self):
        """

        :return:
        """
        self.clear()

        if not self.validate_inputs():
            return

        self.main.emitter_add.connect(self.add_text)
        self.main.emitter_done.connect(self.terminate)

        if not self.main.running:
            self.main.set_process_running()
            self.runner_input = self.main.get_input_params()
            is_whole_blade = self.optimization_mode_selection.currentIndex() == 1
            self.win.set_axis_mode(is_whole_blade)
            if self.preview_win is None:
                self.preview_win = OptimizationPreviewWindow(self)
            self.preview_win.reset()
            self.queue_pyqtgraph[0] = [[0], [0], 0, 0, None]
            self.main.getter.start()
            self.p = Process(target=optimize, args=[self.runner_input, self.queue_pyqtgraph])
            self.p.start()
            self.win.show()
            self.preview_win.show()
            self.win.start_update()

    def terminate(self):
        """

        """
        self.main.set_process_stopped()
        try:
            self.p.terminate()
        except:
            pass
        self.main.getter.quit()

        self.main.emitter_add.disconnect()
        self.main.emitter_done.disconnect()

        self.win.stop_update()
        if self.preview_win is not None:
            self.preview_win.hide()
        self.show_results()

    def show_results(self):
        if len(self.main.return_results) == 0:
            return
        result = self.main.return_results[-1]
        if "optimized_geometry" not in result:
            return

        geometry = result["optimized_geometry"]
        summary_lines = []
        if "thrust" in result:
            summary_lines.append("Thrust: %s N" % round(float(result["thrust"]), 6))
        if "Q" in result:
            summary_lines.append("Torque: %s Nm" % round(float(result["Q"]), 6))
        if "power" in result:
            summary_lines.append("Power: %s W" % round(float(result["power"]), 6))
        if "cp" in result:
            summary_lines.append("Cp: %s" % round(float(result["cp"]), 6))
        if "eff" in result:
            summary_lines.append("Efficiency: %s" % round(float(result["eff"]), 6))
        if geometry.get("pitch") is not None:
            summary_lines.append("Pitch: %s deg" % round(float(geometry["pitch"]), 6))
        summary_lines.append("Rows: %d" % len(geometry.get("r_in", [])))

        self.result_window = OptimizationResultWindow(self, geometry, summary_lines)

    def get_settings(self):
        """

        :return:
        """
        out = {}
        out["target_rpm"] = self.target_rpm.text()
        out["target_speed"] = self.target_speed.text()
        out["mut_coeff"] = self.mut_coeff.text()
        out["population"] = self.population.text()
        out["num_iter"] = self.num_iter.text()
        out["optimization_mode"] = self.optimization_mode_selection.currentIndex()
        out["whole_blade_objective"] = self.whole_blade_objective_selection.currentIndex()
        out["whole_blade_solver"] = self.whole_blade_solver_selection.currentIndex()
        out["whole_blade_twist_cp"] = self.whole_blade_twist_cp.text()
        out["whole_blade_chord_cp"] = self.whole_blade_chord_cp.text()
        out["whole_blade_root_cp_lock"] = self.whole_blade_root_cp_lock.text()
        out["whole_blade_tip_cp_lock"] = self.whole_blade_tip_cp_lock.text()
        out["whole_blade_target_thrust"] = self.whole_blade_target_thrust.text()
        out["whole_blade_max_power"] = self.whole_blade_max_power.text()
        out["whole_blade_max_stall_percentage"] = self.whole_blade_max_stall_percentage.text()
        out["whole_blade_monotonic_chord"] = self.whole_blade_monotonic_chord.isChecked()

        for k, v in out.items():
            if isinstance(v, bool):
                continue
            if k in ["optimization_mode", "whole_blade_objective", "whole_blade_solver"]:
                continue
            if v == "":
                v = None
            elif v == None:
                pass
            else:
                v = to_float(v)
            out[k] = v

        self.update_input_variables_list()
        self.update_output_variables_list()
        self.update_target_variables_list()
        out["optimization_inputs"] = self.list_input_variables
        out["optimization_outputs"] = self.list_output_variables
        out["optimization_targets"] = self.list_target_variables
        return out

    @staticmethod
    def set_line_edit_value(widget, value):
        if value is None:
            widget.setText("")
        else:
            widget.setText(str(value))

    def set_settings(self, inp_dict):
        """

        :param inp_dict:
        """
        if "target_rpm" in inp_dict.keys():
            self.set_line_edit_value(self.target_rpm, inp_dict["target_rpm"])
        if "target_speed" in inp_dict.keys():
            self.set_line_edit_value(self.target_speed, inp_dict["target_speed"])
        if "mut_coeff" in inp_dict.keys():
            self.set_line_edit_value(self.mut_coeff, inp_dict["mut_coeff"])
        if "population" in inp_dict.keys():
            self.set_line_edit_value(self.population, inp_dict["population"])
        if "num_iter" in inp_dict.keys():
            self.set_line_edit_value(self.num_iter, inp_dict["num_iter"])
        if "optimization_mode" in inp_dict.keys():
            self.optimization_mode_selection.setCurrentIndex(int(inp_dict["optimization_mode"]))
        if "whole_blade_objective" in inp_dict.keys():
            self.whole_blade_objective_selection.setCurrentIndex(int(inp_dict["whole_blade_objective"]))
        if "whole_blade_solver" in inp_dict.keys():
            self.whole_blade_solver_selection.setCurrentIndex(int(inp_dict["whole_blade_solver"]))
        if "whole_blade_twist_cp" in inp_dict.keys():
            self.set_line_edit_value(self.whole_blade_twist_cp, inp_dict["whole_blade_twist_cp"])
        if "whole_blade_chord_cp" in inp_dict.keys():
            self.set_line_edit_value(self.whole_blade_chord_cp, inp_dict["whole_blade_chord_cp"])
        if "whole_blade_root_cp_lock" in inp_dict.keys():
            self.set_line_edit_value(self.whole_blade_root_cp_lock, inp_dict["whole_blade_root_cp_lock"])
        if "whole_blade_tip_cp_lock" in inp_dict.keys():
            self.set_line_edit_value(self.whole_blade_tip_cp_lock, inp_dict["whole_blade_tip_cp_lock"])
        if "whole_blade_target_thrust" in inp_dict.keys():
            self.set_line_edit_value(self.whole_blade_target_thrust, inp_dict["whole_blade_target_thrust"])
        if "whole_blade_max_power" in inp_dict.keys():
            self.set_line_edit_value(self.whole_blade_max_power, inp_dict["whole_blade_max_power"])
        if "whole_blade_max_stall_percentage" in inp_dict.keys():
            self.set_line_edit_value(self.whole_blade_max_stall_percentage, inp_dict["whole_blade_max_stall_percentage"])
        if "whole_blade_monotonic_chord" in inp_dict.keys():
            self.whole_blade_monotonic_chord.setChecked(bool(inp_dict["whole_blade_monotonic_chord"]))
        if "optimization_inputs" in inp_dict.keys():
            self.list_input_variables = inp_dict["optimization_inputs"]
        if "optimization_outputs" in inp_dict.keys():
            self.list_output_variables = inp_dict["optimization_outputs"]
        if "optimization_targets" in inp_dict.keys():
            self.list_target_variables = inp_dict["optimization_targets"]
        
        self.refresh_output_variables()
        self.refresh_input_variables()
        self.refresh_target_variables()
        self.update_mode_visibility()

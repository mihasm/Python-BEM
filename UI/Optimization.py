from multiprocessing import Manager
from multiprocessing.context import Process

import numpy as np
from PyQt5 import QtGui
from PyQt5.QtWidgets import QWidget, QTextEdit, QLineEdit, QLabel, QComboBox, QCheckBox, QPushButton, QGridLayout, \
    QFormLayout, QMessageBox

from UI.helpers import PyQtGraphWindow
from optimization import optimize_angles_genetic
from utils import MyMessageBox, to_float


class Optimization(QWidget):
    def __init__(self, parent=None):
        super(Optimization, self).__init__(parent)
        self.main = self.parent()

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
        self.form_list.append([self._target_speed, self.target_speed])
        self.target_speed.setToolTip("Ciljna hitrost vetra. Pri tej hitrosti se bo izvajala BEM analiza.")

        self.target_rpm = QLineEdit()
        self.target_rpm.setValidator(self.validator)
        self.target_rpm.textChanged.connect(self.check_state)
        self.target_rpm.textChanged.emit(self.target_rpm.text())
        self._target_rpm = QLabel("Target rpm [RPM]")
        self.form_list.append([self._target_rpm, self.target_rpm])
        self.target_rpm.setToolTip("Ciljni vrtljaji rotorja turbine. Pri teh vrtljajih se bo izvajala BEM analiza.")

        self.min_bound = QLineEdit()
        self.min_bound.setValidator(self.validator)
        self.min_bound.textChanged.connect(self.check_state)
        self.min_bound.textChanged.emit(self.min_bound.text())
        self._min_bound = QLabel("Minimum boundary")
        self.form_list.append([self._min_bound, self.min_bound])
        self.min_bound.setToolTip("Minimalni kot (theta) lopatice.")

        self.max_bound = QLineEdit()
        self.max_bound.setValidator(self.validator)
        self.max_bound.textChanged.connect(self.check_state)
        self.max_bound.textChanged.emit(self.max_bound.text())
        self._max_bound = QLabel("Maximum boundary")
        self.form_list.append([self._max_bound, self.max_bound])
        self.max_bound.setToolTip("Maksimalni kot (theta) lopatice.")

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

        self._opt_variable = QLabel("Optimization variable")
        self.opt_variable = QComboBox()
        self.opt_variable.addItems(
            ["max(dQ) (torque->wind turbine)", "max(dT) (thrust->propeller)", "max(weight_dq*dQ-weight_dt*dT)",
             "max(weight_dt*dT-weight_dq*dQ)"])
        self.opt_variable.setCurrentIndex(0)
        self.form_list.append([self._opt_variable, self.opt_variable])
        self.opt_variable.setToolTip(
            "Optimizacija naj poteka za to izbrano spremenljivko. V primeru vetrne turbine max(dQ).")

        self.weight_dq = QLineEdit()
        self.weight_dq.setValidator(self.validator)
        self.weight_dq.textChanged.connect(self.check_state)
        self.weight_dq.textChanged.emit(self.weight_dq.text())
        self._weight_dq = QLabel("weight_dq")
        self.form_list.append([self._weight_dq, self.weight_dq])
        self.weight_dq.setToolTip("Utež dq v primeru dvojne optimizacije.")

        self.weight_dt = QLineEdit()
        self.weight_dt.setValidator(self.validator)
        self.weight_dt.textChanged.connect(self.check_state)
        self.weight_dt.textChanged.emit(self.weight_dt.text())
        self._weight_dt = QLabel("weight_dt")
        self.form_list.append([self._weight_dt, self.weight_dt])
        self.weight_dt.setToolTip("Utež dt v primeru dvojne optimizacije.")

        self.pitch_optimization = QCheckBox()
        self.pitch_optimization.setChecked(False)
        self._pitch_optimization = QLabel("Pitch optimization")
        self.form_list.append([self._pitch_optimization, self.pitch_optimization])
        self.pitch_optimization.setToolTip(
            "To možnost izberemo samo, kadar nas zanima optimalni nastavni kot lopatice.")

        self.buttonOptimization = QPushButton("Run optimization")
        self.buttonOptimization.clicked.connect(self.run)

        self.buttonStop = QPushButton("Stop")
        self.buttonStop.clicked.connect(self.terminate)

        self.buttonClear = QPushButton("Clear screen")
        self.buttonClear.clicked.connect(self.clear)

        self.buttonEOF = QCheckBox()
        self.buttonEOF.setChecked(True)
        self.buttonEOFdescription = QLabel("Scroll to end of screen")

        self.grid = QGridLayout()
        self.setLayout(self.grid)
        self.grid.addWidget(self.left, 1, 1)
        self.grid.addWidget(self.textEdit, 1, 2)

        self.fbox = QFormLayout()
        self.left.setLayout(self.fbox)

        for a, b in self.form_list:
            self.fbox.addRow(a, b)

        self.fbox.addRow(self.buttonOptimization)
        self.fbox.addRow("", QLabel())

        self.fbox.addRow(self.buttonClear, self.buttonStop)
        self.fbox.addRow(self.buttonEOFdescription, self.buttonEOF)

        self.win = PyQtGraphWindow(self)
        self.win.setWindowTitle("Live Optimization Visualizer")
        self.manager_pyqtgraph = Manager()
        self.queue_pyqtgraph = self.manager_pyqtgraph.list()
        self.queue_pyqtgraph.append([[0], [0], 0, 0])

        self.tsr_string = QLabel("0")
        self.J_string = QLabel("0")
        self.fbox.addRow("TSR:", self.tsr_string)
        self.fbox.addRow("J:", self.J_string)

    def update_tsr_and_j(self):
        try:
            s = self.get_settings()
            R = float(self.main.wind_turbine_properties.R.text())
            tsr = 2 * np.pi * float(s["target_rpm"]) * R / 60 / float(s["target_speed"])
            self.tsr_string.setText("%.2f" % (tsr))
            J = float(s["target_speed"]) / (float(s["target_rpm"]) / 60 * 2 * R)
            self.J_string.setText("%.2f" % (J))
        except:
            pass

    def check_forms_angles(self):
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
        self.textEdit.clear()

    def add_text(self, string):
        self.textEdit.insertPlainText(string)
        if self.buttonEOF.checkState() == 2:
            self.textEdit.moveCursor(QtGui.QTextCursor.End)

    def run(self):
        self.clear()

        if not self.validate_inputs():
            return

        self.main.emitter_add.connect(self.add_text)
        self.main.emitter_done.connect(self.terminate)

        if not self.main.running:
            self.main.set_process_running()
            self.runner_input = self.main.get_input_params()
            self.main.getter.start()
            self.p = Process(target=optimize_angles_genetic, args=[self.runner_input, self.queue_pyqtgraph])
            self.p.start()
            self.win.show()
            self.win.start_update()

    def terminate(self):
        self.main.set_process_stopped()
        try:
            self.p.terminate()
        except:
            pass
        self.main.getter.quit()

        self.main.emitter_add.disconnect()
        self.main.emitter_done.disconnect()

        self.win.stop_update()

    def get_settings(self):
        out = {}
        out["target_rpm"] = self.target_rpm.text()
        out["target_speed"] = self.target_speed.text()
        out["pitch_optimization"] = bool(self.pitch_optimization.checkState())
        out["min_bound"] = self.min_bound.text()
        out["max_bound"] = self.max_bound.text()
        out["mut_coeff"] = self.mut_coeff.text()
        out["population"] = self.population.text()
        out["num_iter"] = self.num_iter.text()
        out["weight_dt"] = self.weight_dt.text()
        out["weight_dq"] = self.weight_dq.text()

        if int(self.opt_variable.currentIndex()) == 0:
            out["optimization_variable"] = "dQ"
        elif int(self.opt_variable.currentIndex()) == 1:
            out["optimization_variable"] = "dT"
        elif int(self.opt_variable.currentIndex()) == 2:
            out["optimization_variable"] = "dQ-dT"
        elif int(self.opt_variable.currentIndex()) == 3:
            out["optimization_variable"] = "dT-dQ"

        for k, v in out.items():
            if v == "":
                v = None
            elif v == None:
                pass
            else:
                if not k == "optimization_variable":
                    v = to_float(v)
            out[k] = v
        return out

    def set_settings(self, inp_dict):
        self.target_rpm.setText(str(inp_dict["target_rpm"]))
        self.target_speed.setText(str(inp_dict["target_speed"]))
        self.pitch_optimization.setChecked(inp_dict["pitch_optimization"])
        self.min_bound.setText(str(inp_dict["min_bound"]))
        self.max_bound.setText(str(inp_dict["max_bound"]))
        self.mut_coeff.setText(str(inp_dict["mut_coeff"]))
        self.population.setText(str(inp_dict["population"]))
        self.num_iter.setText(str(inp_dict["num_iter"]))
        self.weight_dt.setText(str(inp_dict["weight_dt"]))
        self.weight_dq.setText(str(inp_dict["weight_dq"]))
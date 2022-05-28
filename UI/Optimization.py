from multiprocessing import Manager
from multiprocessing.context import Process

import numpy as np
from PyQt5 import QtGui
from PyQt5.QtWidgets import QWidget, QGridLayout, QFormLayout, QScrollArea, QVBoxLayout, QTextEdit, QPushButton, \
    QComboBox, QCheckBox, QLineEdit, QLabel, QMessageBox

from UI.helpers import PyQtGraphWindow
from optimization import optimize
from utils import MyMessageBox, to_float


class Optimization(QWidget):
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

        self.input_variables_widget =QWidget()
        self.input_variables_layout = QGridLayout()
        self.input_variables_widget.setLayout(self.input_variables_layout)
        self.button_add_input_variable = QPushButton("Add input variable")
        self.button_add_input_variable.clicked.connect(self.add_input_variable)
        self.input_variable_selection = QComboBox()
        self.input_variable_selection.addItems(["_theta","_c"])
        self.list_input_variables = []

        self.output_variables_widget =QWidget()
        self.output_variables_layout = QGridLayout()
        self.output_variables_widget.setLayout(self.output_variables_layout)
        self.button_add_output_variable = QPushButton("Add output variable")
        self.button_add_output_variable.clicked.connect(self.add_output_variable)
        self.output_variable_selection = QComboBox()
        self.output_variable_selection.addItems(["dQ","dT","a","a'","Cl","Cd","dFn","dFt","U4","alpha","phi"])
        self.list_output_variables = []

        self.target_variables_widget =QWidget()
        self.target_variables_layout = QGridLayout()
        self.target_variables_widget.setLayout(self.target_variables_layout)
        self.button_add_target_variable = QPushButton("Add target variable")
        self.button_add_target_variable.clicked.connect(self.add_target_variable)
        self.target_variable_selection = QComboBox()
        self.target_variable_selection.addItems(["dQ","dT","a","a'","Cl","Cd","dFn","dFt","U4","alpha","phi"])
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

        self.fbox.addRow(self.button_add_input_variable,self.input_variable_selection)
        self.fbox.addRow(self.input_variables_widget)

        self.fbox.addRow(self.button_add_output_variable,self.output_variable_selection)
        self.fbox.addRow(self.output_variables_widget)

        self.fbox.addRow(self.button_add_target_variable,self.target_variable_selection)
        self.fbox.addRow(self.target_variables_widget)

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

    def refresh_input_variables(self):
        for i in reversed(range(self.input_variables_layout.count())): 
            self.input_variables_layout.itemAt(i).widget().setParent(None)

        i=0
        for var,min_b,max_b in self.list_input_variables:
            name = QLabel(var)
            min_bound = QLineEdit(str(min_b))
            min_bound.textChanged.connect(self.update_input_variables_list)
            max_bound = QLineEdit(str(max_b))
            max_bound.textChanged.connect(self.update_input_variables_list)
            delete_button = QPushButton("X")
            delete_button.clicked.connect(self.delete_input_variable)
            self.input_variables_layout.addWidget(name,i,0)
            self.input_variables_layout.addWidget(min_bound,i,1)
            self.input_variables_layout.addWidget(max_bound,i,2)
            self.input_variables_layout.addWidget(delete_button,i,3)
            i+=1

    def refresh_output_variables(self):
        for i in reversed(range(self.output_variables_layout.count())): 
            self.output_variables_layout.itemAt(i).widget().setParent(None)

        i=0
        for var,coefficient in self.list_output_variables:
            name = QLabel(var)
            coefficient = QLineEdit(str(coefficient))
            coefficient.textChanged.connect(self.update_output_variables_list)
            delete_button = QPushButton("X")
            delete_button.clicked.connect(self.delete_output_variable)
            self.output_variables_layout.addWidget(name,i,0)
            self.output_variables_layout.addWidget(coefficient,i,1)
            self.output_variables_layout.addWidget(delete_button,i,2)
            i+=1

    def refresh_target_variables(self):
        for i in reversed(range(self.target_variables_layout.count())): 
            self.target_variables_layout.itemAt(i).widget().setParent(None)

        i=0
        for var,target_value,coefficient in self.list_target_variables:
            name = QLabel(var)

            target_value = QLineEdit(str(target_value))
            target_value.textChanged.connect(self.update_target_variables_list)

            coefficient = QLineEdit(str(coefficient))
            coefficient.textChanged.connect(self.update_target_variables_list)

            delete_button = QPushButton("X")
            delete_button.clicked.connect(self.delete_target_variable)
            

            self.target_variables_layout.addWidget(name,i,0)
            self.target_variables_layout.addWidget(target_value,i,1)
            self.target_variables_layout.addWidget(coefficient,i,2)
            self.target_variables_layout.addWidget(delete_button,i,3)
            i+=1

    def update_input_variables_list(self):
        for row in range(len(self.list_input_variables)):
            try:
                min_b = float(self.input_variables_layout.itemAtPosition(row,1).widget().text())
                self.list_input_variables[row][1] = min_b
            except ValueError:
                pass
            try:
                max_b = float(self.input_variables_layout.itemAtPosition(row,2).widget().text())
                self.list_input_variables[row][2] = max_b
            except ValueError:
                pass

    def update_output_variables_list(self):
        for row in range(len(self.list_output_variables)):
            try:
                coeff = float(self.output_variables_layout.itemAtPosition(row,1).widget().text())
                self.list_output_variables[row][1] = coeff
            except ValueError:
                pass

    def update_target_variables_list(self):
        for row in range(len(self.list_target_variables)):
            try:
                target_value = float(self.target_variables_layout.itemAtPosition(row,1).widget().text())
                self.list_target_variables[row][1] = target_value
            except ValueError:
                pass
            try:
                coeff = float(self.target_variables_layout.itemAtPosition(row,2).widget().text())
                self.list_target_variables[row][2] = coeff
            except ValueError:
                pass

    def delete_input_variable(self):
        button = self.sender()
        index = self.input_variables_layout.indexOf(button)
        row,column,_,_ = self.input_variables_layout.getItemPosition(index)
        del self.list_input_variables[row]
        self.refresh_input_variables()

    def delete_output_variable(self):
        button = self.sender()
        index = self.output_variables_layout.indexOf(button)
        row,column,_,_ = self.output_variables_layout.getItemPosition(index)
        del self.list_output_variables[row]
        self.refresh_output_variables()

    def delete_target_variable(self):
        button = self.sender()
        index = self.target_variables_layout.indexOf(button)
        row,column,_,_ = self.target_variables_layout.getItemPosition(index)
        del self.list_target_variables[row]
        self.refresh_target_variables()

    def add_input_variable(self):
        variable = str(self.input_variable_selection.currentText())
        self.list_input_variables.append([variable,0,0])
        self.refresh_input_variables()

    def add_output_variable(self):
        variable = str(self.output_variable_selection.currentText())
        self.list_output_variables.append([variable,0])
        self.refresh_output_variables()

    def add_target_variable(self):
        variable = str(self.target_variable_selection.currentText())
        self.list_target_variables.append([variable,0,1])
        self.refresh_target_variables()

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
        self.right_output_text_area.clear()

    def add_text(self, string):
        self.right_output_text_area.insertPlainText(string)
        if self.buttonEOF.checkState() == 2:
            self.right_output_text_area.moveCursor(QtGui.QTextCursor.End)

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
            self.p = Process(target=optimize, args=[self.runner_input, self.queue_pyqtgraph])
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
        out["mut_coeff"] = self.mut_coeff.text()
        out["population"] = self.population.text()
        out["num_iter"] = self.num_iter.text()
        
        for k, v in out.items():
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
        out["optimization_inputs"]=self.list_input_variables
        out["optimization_outputs"]=self.list_output_variables
        out["optimization_targets"]=self.list_target_variables
        return out

    def set_settings(self, inp_dict):
        self.target_rpm.setText(str(inp_dict["target_rpm"]))
        self.target_speed.setText(str(inp_dict["target_speed"]))
        #self.pitch_optimization.setChecked(inp_dict["pitch_optimization"])
        self.mut_coeff.setText(str(inp_dict["mut_coeff"]))
        self.population.setText(str(inp_dict["population"]))
        self.num_iter.setText(str(inp_dict["num_iter"]))
        self.list_input_variables = inp_dict["optimization_inputs"]
        self.list_output_variables = inp_dict["optimization_outputs"]
        self.refresh_output_variables()
        self.refresh_input_variables()
        self.refresh_target_variables()
        #self.weight_dt.setText(str(inp_dict["weight_dt"]))
        #self.weight_dq.setText(str(inp_dict["weight_dq"]))
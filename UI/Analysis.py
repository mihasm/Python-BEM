from multiprocessing.context import Process

import numpy as np
from PyQt5 import QtGui
from PyQt5.QtWidgets import QWidget, QGridLayout, QFormLayout, QScrollArea, QVBoxLayout, QTextEdit, QPushButton, \
    QComboBox, QCheckBox, QLineEdit, QLabel, QMessageBox

from calculation_runner import calculate_power_3d
from popravki import METHODS_STRINGS
from utils import to_float
from UI.helpers import MyMessageBox


class Analysis(QWidget):
    """

    """
    def __init__(self, parent=None):
        super(Analysis, self).__init__(parent)

        self.main = self.parent()

        self.forms_dict = {}

        self.settings = {"invert_alpha": False,
                         "use_minimization_solver": False,
                         "tip_loss_mode": ["None", "Prandtl", "Shen", "Adkins"],
                         "hub_loss_mode": ["None", "Prandtl", "Shen"],
                         "cascade_correction": False,
                         "skewed_wake_correction": False,
                         "rotational_augmentation_correction": False,
                         "rotational_augmentation_correction_method": ["Snel et al.", "Du & Selig", "Chaviaropoulos and Hansen", "Lindenburg", "Dumitrescu and Cardos", "Snel et al. (prop.)", "Gur & Rosen (prop.)"],
                         "mach_number_correction": False,
                         "rho": 1.225,
                         "kin_viscosity": 1.4207E-5,
                         "method": 6,
                         "variable_selection": ["RPM and v", "TSR", "J", "pitch", "pitch+TSR", "pitch+J"],
                         "constant_selection": ["speed", "rpm"],
                         "constant_speed": 5,
                         "constant_rpm": 1500,
                         "pitch": 0.0,
                         "v_min": 3,
                         "v_max": 20,
                         "v_num": 10,
                         "rpm_min": 100,
                         "rpm_max": 3000,
                         "rpm_num": 10,
                         "tsr_min": 1,
                         "tsr_max": 10,
                         "tsr_num": 10,
                         "J_min": 0.1,
                         "J_max": 1.5,
                         "J_num": 10,
                         "pitch_min": -15,
                         "pitch_max": 15,
                         "pitch_num": 10,
                         "max_iterations": 100,
                         "convergence_limit": 0.001,
                         "relaxation_factor": 0.3,
                         "a_initial": 0.3,
                         "aprime_initial": 0.01,
                         "print_all": False,
                         "print_out": False,
                         "reynolds": 50000,
                         "fix_reynolds": False,
                         "yaw_angle": 0}

        self.settings_to_name = {"print_out": "Print final iteration data",
                                 "invert_alpha": "Invert alpha",
                                 "tip_loss_mode": "Tip loss",
                                 "hub_loss_mode": "Hub loss",
                                 "cascade_correction": "Cascade correction",
                                 "rho": "Air density [kg/m^3]",
                                 "kin_viscosity": "Kinematic viscosity [m^2/s]",
                                 "method": "Calculation method",
                                 "variable_selection": "Variable parameter",
                                 "constant_selection": "Constant variable",
                                 "constant_speed": "Wind speed",
                                 "constant_rpm": "RPM",
                                 "pitch": "Pitch",
                                 "v_min": "Min calc. wind speed [m/s]",
                                 "v_max": "Max calc. wind speed [m/s]",
                                 "v_num": "Number of wind speed points",
                                 "rpm_min": "Min calc. RPM [RPM]",
                                 "rpm_max": "Max calc. RPM [RPM]",
                                 "rpm_num": "Number of RPM points",
                                 "tsr_min": "Min TSR",
                                 "tsr_max": "Max TSR",
                                 "tsr_num": "Num TSR",
                                 "J_min": "Min J",
                                 "J_max": "Max J",
                                 "J_num": "Num J",
                                 "pitch_min": "Min pitch",
                                 "pitch_max": "Max pitch",
                                 "pitch_num": "Num pitch",
                                 "max_iterations": "Maximum iterations",
                                 "convergence_limit": "Convergence criteria",
                                 "relaxation_factor": "Relaxation factor",
                                 "a_initial": "Initial a",
                                 "aprime_initial": "Initial a'",
                                 "print_all": "Print every iteration [debug]",
                                 "rotational_augmentation_correction": "Rot. augmentation cor.",
                                 "rotational_augmentation_correction_method": "Rot. augmentation cor. method",
                                 "fix_reynolds": "Fix Reynolds", "reynolds": "Reynolds",
                                 "mach_number_correction": "Mach number correction",
                                 "yaw_angle": "Yaw angle [°]", "skewed_wake_correction": "Skewed Wake Correction",
                                 "use_minimization_solver": "Use minimization solver"}

        self.settings_to_tooltip = {
            "invert_alpha": "Uporabimo, kadar računamo propeler v obratovalnem stanju vetrnice.",
            "print_out": "Izpis končnih vrednosti po konvergenci za vsak odsek",
            "tip_loss_mode": "Popravek izgub pri vrhu lopatice",
            "hub_loss_mode": "Popravek izgub pri pestu",
            "cascade_correction": "Kaskadni popravki",
            "rho": "Gostota fluida [kg/m^3]",
            "kin_viscosity": "Kinematična viskoznost fluida [m^2/s]",
            "method": "Metoda za preračun. Privzeta je e) Aerodyn (Buhl).",
            "variable_selection": "Izbira spremenljivega parametra",
            "constant_selection": "Konstantna spremenljivka",
            "constant_speed": "Constant wind speed",
            "constant_rpm": "Constant RPM",
            "pitch": "Nastavni kot lopatice",
            "v_min": "Minimalna hitrost fluida [m/s]",
            "v_max": "Maksimalna hitrost fluida [m/s]",
            "v_num": "Število računskih točk (linearno razporejenih) od min hitrosti fluida do max hitrosti fluida",
            "rpm_min": "Minimalni vrtljaji/min [RPM]",
            "rpm_max": "Maksimalni vrtljaji/min [RPM]",
            "rpm_num": "Število računskih točk (linearno razporejenih) od min RPM do max RPM",
            "tsr_min": "Min TSR",
            "tsr_max": "Max TSR",
            "tsr_num": "Num TSR",
            "J_min": "Min J",
            "J_max": "Max J",
            "J_num": "Num J",
            "pitch_min": "Min pitch",
            "pitch_max": "Max pitch",
            "pitch_num": "Num pitch",
            "max_iterations": "Maksimalno število iteracij.",
            "convergence_limit": "Konvergenčni kriterij.",
            "relaxation_factor": "Relaksacijski faktor. Privzeta vrednost: 0.3",
            "print_all": "Podroben izpis vrednosti po vsaki iteraciji (upočasni izračun)",
            "rotational_augmentation_correction": "Popravek rotacijske augmentacije",
            "rotational_augmentation_correction_method": "Izbira metode za popravek rotacijske augmentacije",
            "fix_reynolds": "Izračunaj vse odseke pri enem samem Re",
            "reynolds": "Se uporabi samo v primeru, če izberemo 'Fix Reynolds' opcijo",
            "mach_number_correction": "Popravek Mach števila (uporabno pri propelerjih)",
            "yaw_angle": "Kot fluida glede na smer osi rotorja [°]. Če je turbina obrnjena proti fluidu, je 0°.",
            "skewed_wake_correction": "Popravek nagnjenega zračnega toka za turbino (Skewed wake)",
            "use_minimization_solver": "Uporaba minimizacijskega algoritma za iskanje indukcijskih faktorjev",
            "a_initial": "Začetna vrednost aksialnega indukcijskega faktorja. Vpliva na konvergenco (propeler)!",
            "aprime_initial": "Začetna vrednost tangencialnega indukcijskega faktorja. Vpliva na konvergenco."}

        self.list_settings_for_updating_tsr = ["v_min", "v_max", "v_num", "rpm_min", "rpm_max", "rpm_num"]

        self.methods_to_names = METHODS_STRINGS

        self.name_to_methods = {v: k for k, v in self.methods_to_names.items()}
        self.name_to_settings = {v: k for k, v in self.settings_to_name.items()}

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

        self.buttonRun = QPushButton("Run")

        self.validator = QtGui.QDoubleValidator()

        self.form_list = []

        hideable_parameters = ["constant_selection", "constant_speed", "constant_rpm", "pitch", "v_min", "v_max",
                               "v_num", "rpm_min", "rpm_max", "rpm_num", "tsr_min", "tsr_max", "tsr_num", "J_min",
                               "J_max", "J_num", "pitch_min", "pitch_max", "pitch_num", "reynolds"]

        for key, value in self.settings.items():
            if key == "method":
                form = QComboBox()
                form.addItems([self.methods_to_names[k] for k, v in self.methods_to_names.items()])
                form.setCurrentIndex(0)
            elif isinstance(value, list):
                form = QComboBox()
                form.addItems([l for l in value])
                form.currentIndexChanged.connect(self.set_parameter_visibility)
            elif isinstance(value, bool):
                form = QCheckBox()
                form.setTristate(value)
                form.stateChanged.connect(self.set_parameter_visibility)
            else:
                form = QLineEdit()
                form.setValidator(self.validator)
                form.textChanged.connect(self.check_state)
                form.textChanged.emit(form.text())
                form.insert(str(value))
                if key in self.list_settings_for_updating_tsr:
                    form.textChanged.connect(self.update_tsr_and_j)

            form.setToolTip(self.settings_to_tooltip[key])
            label = self.settings_to_name[key]

            self.form_list.append([label, form, key])

            if key in hideable_parameters:
                row = QWidget()
                row_layout = QFormLayout()
                row.setLayout(row_layout)
                row_layout.setContentsMargins(0, 0, 0, 0)  # tight layout
                row_layout.addRow(label, form)
                self.fbox.addRow(row)
                self.forms_dict[key] = [form, label, row]
            else:
                self.fbox.addRow(label, form)
                self.forms_dict[key] = [form, label]

        self.emptyLabel = QLabel(" ")
        self.buttonRun.clicked.connect(self.run)
        self.buttonClear = QPushButton("Clear screen")
        self.buttonClear.clicked.connect(self.clear)
        self.buttonEOF = QCheckBox()
        self.buttonEOF.setChecked(True)
        self.buttonEOFdescription = QLabel("Scroll to end of screen")
        self.buttonStop = QPushButton("Stop")
        self.buttonStop.clicked.connect(self.terminate)

        self.fbox.addRow(self.emptyLabel, self.buttonRun)
        self.fbox.addRow(self.buttonClear, self.buttonStop)
        self.fbox.addRow(self.buttonEOFdescription, self.buttonEOF)

        self.tsr_string = QLabel("0")
        self.J_string = QLabel("0")
        self.fbox.addRow("TSR:", self.tsr_string)
        self.fbox.addRow("J:", self.J_string)

        self.turbine_type = self.parent().wind_turbine_properties.turbine_type

        self.turbine_type.currentIndexChanged.connect(self.refresh_methods)

        self.set_parameter_visibility()
        self.refresh_methods()

    def set_parameter_visibility(self):
        """

        """
        current_index = self.forms_dict["variable_selection"][0].currentIndex()
        if current_index == 0:
            # RPM / speed
            # rpm, speed need to be shown
            # constant_selection, TSR, J and pitch need to be hidden
            self.forms_dict["constant_selection"][2].hide()
            self.forms_dict["constant_speed"][2].hide()
            self.forms_dict["constant_rpm"][2].hide()
            self.forms_dict["pitch"][2].show()
            self.forms_dict["v_min"][2].show()
            self.forms_dict["v_max"][2].show()
            self.forms_dict["v_num"][2].show()
            self.forms_dict["rpm_min"][2].show()
            self.forms_dict["rpm_max"][2].show()
            self.forms_dict["rpm_num"][2].show()
            self.forms_dict["tsr_min"][2].hide()
            self.forms_dict["tsr_max"][2].hide()
            self.forms_dict["tsr_num"][2].hide()
            self.forms_dict["J_min"][2].hide()
            self.forms_dict["J_max"][2].hide()
            self.forms_dict["J_num"][2].hide()
            self.forms_dict["pitch_min"][2].hide()
            self.forms_dict["pitch_max"][2].hide()
            self.forms_dict["pitch_num"][2].hide()
            pass
        elif current_index == 1:
            # variable TSR
            # tsr, constant_selection, appropriate constant need to be shown
            self.forms_dict["constant_selection"][2].show()
            self.forms_dict["constant_speed"][2].show()
            self.forms_dict["constant_rpm"][2].show()
            self.forms_dict["pitch"][2].show()
            self.forms_dict["v_min"][2].hide()
            self.forms_dict["v_max"][2].hide()
            self.forms_dict["v_num"][2].hide()
            self.forms_dict["rpm_min"][2].hide()
            self.forms_dict["rpm_max"][2].hide()
            self.forms_dict["rpm_num"][2].hide()
            self.forms_dict["tsr_min"][2].show()
            self.forms_dict["tsr_max"][2].show()
            self.forms_dict["tsr_num"][2].show()
            self.forms_dict["J_min"][2].hide()
            self.forms_dict["J_max"][2].hide()
            self.forms_dict["J_num"][2].hide()
            self.forms_dict["pitch_min"][2].hide()
            self.forms_dict["pitch_max"][2].hide()
            self.forms_dict["pitch_num"][2].hide()
            pass
        elif current_index == 2:
            # variable J
            # J, constant_selection, appropriate constant need to be shown
            self.forms_dict["constant_selection"][2].show()
            self.forms_dict["constant_speed"][2].show()
            self.forms_dict["constant_rpm"][2].show()
            self.forms_dict["pitch"][2].show()
            self.forms_dict["v_min"][2].hide()
            self.forms_dict["v_max"][2].hide()
            self.forms_dict["v_num"][2].hide()
            self.forms_dict["rpm_min"][2].hide()
            self.forms_dict["rpm_max"][2].hide()
            self.forms_dict["rpm_num"][2].hide()
            self.forms_dict["tsr_min"][2].hide()
            self.forms_dict["tsr_max"][2].hide()
            self.forms_dict["tsr_num"][2].hide()
            self.forms_dict["J_min"][2].show()
            self.forms_dict["J_max"][2].show()
            self.forms_dict["J_num"][2].show()
            self.forms_dict["pitch_min"][2].hide()
            self.forms_dict["pitch_max"][2].hide()
            self.forms_dict["pitch_num"][2].hide()
            pass
        elif current_index == 3:
            # variable pitch
            # pitch, both constants need to be shown
            self.forms_dict["constant_selection"][2].hide()
            self.forms_dict["constant_speed"][2].show()
            self.forms_dict["constant_rpm"][2].show()
            self.forms_dict["pitch"][2].hide()
            self.forms_dict["v_min"][2].hide()
            self.forms_dict["v_max"][2].hide()
            self.forms_dict["v_num"][2].hide()
            self.forms_dict["rpm_min"][2].hide()
            self.forms_dict["rpm_max"][2].hide()
            self.forms_dict["rpm_num"][2].hide()
            self.forms_dict["tsr_min"][2].hide()
            self.forms_dict["tsr_max"][2].hide()
            self.forms_dict["tsr_num"][2].hide()
            self.forms_dict["J_min"][2].hide()
            self.forms_dict["J_max"][2].hide()
            self.forms_dict["J_num"][2].hide()
            self.forms_dict["pitch_min"][2].show()
            self.forms_dict["pitch_max"][2].show()
            self.forms_dict["pitch_num"][2].show()
            pass
        elif current_index == 4:
            # variable pitch + TSR
            self.forms_dict["constant_selection"][2].show()
            self.forms_dict["constant_speed"][2].show()
            self.forms_dict["constant_rpm"][2].show()
            self.forms_dict["pitch"][2].hide()
            self.forms_dict["v_min"][2].hide()
            self.forms_dict["v_max"][2].hide()
            self.forms_dict["v_num"][2].hide()
            self.forms_dict["rpm_min"][2].hide()
            self.forms_dict["rpm_max"][2].hide()
            self.forms_dict["rpm_num"][2].hide()
            self.forms_dict["tsr_min"][2].show()
            self.forms_dict["tsr_max"][2].show()
            self.forms_dict["tsr_num"][2].show()
            self.forms_dict["J_min"][2].hide()
            self.forms_dict["J_max"][2].hide()
            self.forms_dict["J_num"][2].hide()
            self.forms_dict["pitch_min"][2].show()
            self.forms_dict["pitch_max"][2].show()
            self.forms_dict["pitch_num"][2].show()
            pass
        elif current_index == 5:
            # variable pitch + J
            self.forms_dict["constant_selection"][2].show()
            self.forms_dict["constant_speed"][2].show()
            self.forms_dict["constant_rpm"][2].show()
            self.forms_dict["pitch"][2].hide()
            self.forms_dict["v_min"][2].hide()
            self.forms_dict["v_max"][2].hide()
            self.forms_dict["v_num"][2].hide()
            self.forms_dict["rpm_min"][2].hide()
            self.forms_dict["rpm_max"][2].hide()
            self.forms_dict["rpm_num"][2].hide()
            self.forms_dict["tsr_min"][2].hide()
            self.forms_dict["tsr_max"][2].hide()
            self.forms_dict["tsr_num"][2].hide()
            self.forms_dict["J_min"][2].show()
            self.forms_dict["J_max"][2].show()
            self.forms_dict["J_num"][2].show()
            self.forms_dict["pitch_min"][2].show()
            self.forms_dict["pitch_max"][2].show()
            self.forms_dict["pitch_num"][2].show()
            pass

        current_index_constant_value = self.forms_dict["constant_selection"][0].currentIndex()
        if current_index != 3:
            if current_index_constant_value == 0:
                self.forms_dict["constant_speed"][2].show()
                self.forms_dict["constant_rpm"][2].hide()
            elif current_index_constant_value == 1:
                self.forms_dict["constant_speed"][2].hide()
                self.forms_dict["constant_rpm"][2].show()

        current_state_fix_reynolds = self.forms_dict["fix_reynolds"][0].isChecked()
        if current_state_fix_reynolds == True:
            self.forms_dict["reynolds"][2].show()
        else:
            self.forms_dict["reynolds"][2].hide()

    def refresh_methods(self):
        if self.turbine_type.currentIndex() == 0:
            for i in range(self.forms_dict["method"][0].count()):
                self.forms_dict["method"][0].model().item(i).setEnabled(True)
        else:
            for i in range(self.forms_dict["method"][0].count()):
                if i==0:
                    self.forms_dict["method"][0].model().item(i).setEnabled(True)
                else:
                    self.forms_dict["method"][0].model().item(i).setEnabled(False)
            self.forms_dict["method"][0].setCurrentIndex(0)


    def update_tsr_and_j(self):
        """

        """
        try:
            s = self.get_settings()
            R = float(self.main.wind_turbine_properties.R.text())
            tsr_min = 2 * np.pi * float(s["rpm_min"]) * R / 60 / float(s["v_max"])
            tsr_max = 2 * np.pi * float(s["rpm_max"]) * R / 60 / float(s["v_min"])
            self.tsr_string.setText("%.2f - %.2f" % (tsr_min, tsr_max))
            J_min = float(s["v_min"]) / (float(s["rpm_max"]) / 60 * 2 * R)
            J_max = float(s["v_max"]) / (float(s["rpm_min"]) / 60 * 2 * R)
            self.J_string.setText("%.2f - %.2f" % (J_min, J_max))
        except:
            print("couldnt update tsr min/max or J min/max")

    def check_forms(self):
        """

        :return:
        """
        out = ""
        for n, f, n_short in self.form_list:
            if isinstance(f, QLineEdit):
                state = self.validator.validate(f.text(), 0)[0]
                if state == QtGui.QValidator.Acceptable:
                    pass
                elif state == QtGui.QValidator.Intermediate:
                    out += "Form %s appears not to be valid.\n" % n
                else:
                    out += "Form %s is not of the valid type.\n" % n
        if out == "":
            return True
        return out

    def check_state(self, *args, **kwargs):
        """

        :param args:
        :param kwargs:
        """
        sender = self.sender()
        validator = sender.validator()
        state = validator.validate(sender.text(), 0)[0]
        if state == QtGui.QValidator.Acceptable:
            color = "#edf5e1"  # green
        elif state == QtGui.QValidator.Intermediate:
            color = "#fff79a"  # yellow
        else:
            color = "#f6989d"  # red
        sender.setStyleSheet("QLineEdit { color: %s;}" % color)

    def get_settings(self):
        """

        :return:
        """
        out_settings = {}
        for name, value, name_short in self.form_list:
            name = self.name_to_settings[name]
            if isinstance(value, QCheckBox):
                value = bool(value.checkState())
            if isinstance(value, QLineEdit):
                value = to_float(value.text())
            if isinstance(value, QComboBox):
                value = int(value.currentIndex())
            out_settings[name] = value
        return out_settings

    def set_settings(self, inp_dict):
        """

        :param inp_dict:
        """
        for name_long, item, name in self.form_list:
            if name in inp_dict:
                if isinstance(item, QComboBox):
                    _index = inp_dict[name]
                    index = _index if _index >= 0 and _index <= item.count()-1 else 0
                    item.setCurrentIndex(int(index))
                elif isinstance(item, QLineEdit):
                    item.setText(str(inp_dict[name]))
                elif isinstance(item, QCheckBox):
                    item.setChecked(inp_dict[name])

    def validate_inputs(self):
        """

        :return:
        """
        check = self.check_forms()
        if check != True:
            msg = MyMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText("Input validation error.")
            msg.setDetailedText(check)
            msg.exec_()
        return check

    def run(self):
        """

        :return:
        """
        self.clear()

        if not self.validate_inputs():
            return

        self.runner_input = self.main.get_input_params()

        if self.runner_input == None:
            print("No settings fetched... Exiting.")
            return

        self.main.emitter_add.connect(self.add_text)
        self.main.emitter_done.connect(self.terminate)

        if not self.main.running:
            self.main.set_process_running()
            self.main.getter.start()
            self.p = Process(target=calculate_power_3d, args=[self.runner_input, True])
            self.p.start()

    def add_text(self, string):
        """

        :param string:
        """
        if self.buttonEOF.checkState() == 2:
            self.right_output_text_area.moveCursor(QtGui.QTextCursor.End)
        self.right_output_text_area.insertPlainText(string)

    def clear(self):
        """

        """
        self.right_output_text_area.clear()

    def terminate(self):
        """

        """
        try:
            self.p.terminate()
        except:
            pass

        self.main.set_process_stopped()
        self.main.getter.quit()

        self.main.emitter_add.disconnect()
        self.main.emitter_done.disconnect()

        self.show_results()

    def show_results(self):
        """

        """
        if len(self.main.return_results) > 0:
            results = self.main.return_results[-1]
            if "v" in results:
                if len(results["v"]) > 0:
                    inp_params = self.runner_input
                    from results import ResultsWindow
                    self.r = ResultsWindow(None, self.main.screen_width, self.main.screen_width, results, inp_params, )
                else:
                    print("Not enough points to print results...")
            else:
                print("No results to print...")

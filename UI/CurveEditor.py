from PyQt5 import QtGui
from PyQt5.QtWidgets import QWidget, QGridLayout, QPushButton, QLineEdit, QLabel, QComboBox, QMessageBox

from UI.Curve import Curve
from utils import transpose, MyMessageBox
from UI.Table import Table


class CurveEditor(QWidget):
    """

    """
    def __init__(self, parent=None):
        super(CurveEditor, self).__init__(None)
        self.resize(1600, 768)
        self.parent = parent

        self.grid = QGridLayout()
        self.setLayout(self.grid)

        # self.tab_widget = TabWidget(self)
        # self.grid.addWidget(self.tab_widget, 2, 0)

        self.validator = QtGui.QDoubleValidator()

        self.table = Table()
        self.table.createTable([[0, 0, 0]])
        self.table.set_labels(["alpha [°]", "cL (lift coeff.)", "cD (drag coeff.)"])
        self.grid.addWidget(self.table, 2, 0)

        self.upper_widget = QWidget()
        self.upper_layout = QGridLayout()
        self.upper_widget.setLayout(self.upper_layout)
        self.grid.addWidget(self.upper_widget, 1, 0)

        self.button_remove_curve = QPushButton("Remove curve")
        self.button_remove_curve.clicked.connect(self.remove_curve)
        self.upper_layout.addWidget(self.button_remove_curve, 2, 3)
        self.button_save_curve = QPushButton("Update curve")
        self.button_save_curve.clicked.connect(self.save_curve)
        self.upper_layout.addWidget(self.button_save_curve, 1, 3)
        self.button_add_curve = QPushButton("Add curve")
        self.button_add_curve.clicked.connect(self.add_curve)
        self.upper_layout.addWidget(self.button_add_curve, 4, 3)

        self.ncrit_edit = QLineEdit("Insert ncrit")
        self.upper_layout.addWidget(self.ncrit_edit, 4, 1)
        self.ncrit_edit.setValidator(self.validator)
        self.ncrit_edit.textChanged.connect(self.check_state)
        self.re_edit = QLineEdit("Insert reynolds")
        self.upper_layout.addWidget(self.re_edit, 4, 2)
        self.re_edit.setValidator(self.validator)
        self.re_edit.textChanged.connect(self.check_state)

        # self.picker_mach_label = QLabel("Mach number:")
        # self.picker_mach = QComboBox()
        # self.picker_mach.setEditable(True)
        # self.upper_layout.addWidget(self.picker_mach_label,2,1)
        # self.upper_layout.addWidget(self.picker_mach,3,1)

        self.picker_ncrit_label = QLabel("NCrit:")
        self.picker_ncrit = QComboBox()
        self.picker_ncrit.setEditable(False)
        #
        self.upper_layout.addWidget(self.picker_ncrit_label, 1, 1)
        self.upper_layout.addWidget(self.picker_ncrit, 2, 1)

        self.picker_reynolds_label = QLabel("Reynolds:")
        self.picker_reynolds = QComboBox()
        self.picker_reynolds.setEditable(False)
        #
        self.upper_layout.addWidget(self.picker_reynolds_label, 1, 2)
        self.upper_layout.addWidget(self.picker_reynolds, 2, 2)

        # self.picker_mach.lineEdit().returnPressed.connect(self.refresh_dropdowns)
        # self.sig1 = self.picker_reynolds.lineEdit().returnPressed.connect(self.save_curve)
        # self.sig2 = self.picker_ncrit.lineEdit().returnPressed.connect(self.save_curve)

        self.connect()

        # self.load_curves()

    def save_curve(self):
        """

        """
        out_chosen_values = self.get_chosen_values_from_dropdowns()
        re_chosen, ncrit_chosen = out_chosen_values

        data_from_table = self.table.get_values()
        alpha_table, cl_table, cd_table = transpose(data_from_table)
        alpha, cl, cd = [], [], []

        for i in range(len(alpha_table)):
            if alpha_table[i] != "" and cl_table[i] != "" and cd_table[i] != "":
                alpha.append(float(alpha_table[i]))
                cl.append(float(cl_table[i]))
                cd.append(float(cd_table[i]))

        if self.parent.curves.get_curve(re_in=re_chosen, ncrit_in=ncrit_chosen) == None:
            msg = MyMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText("Curve with these Re - ncrit values does not exist yet. Did you mean to add a new curve?")
            msg.setDetailedText(
                "Curve with Re %s and ncrit %s values does not exist yet. Did you mean to add a new curve?" % (
                    re_chosen, ncrit_chosen))
            msg.exec_()
        else:
            self.current_curve = self.parent.curves.get_curve(re_in=re_chosen, ncrit_in=ncrit_chosen)
            self.current_curve.alpha = alpha
            self.current_curve.cl = cl
            self.current_curve.cd = cd
            self.load_curves()

    def add_curve(self):
        """

        :return:
        """
        self.disconnect()

        re_chosen, ncrit_chosen = self.re_edit.text(), self.ncrit_edit.text()
        try:
            re_chosen, ncrit_chosen = float(re_chosen), float(ncrit_chosen)
        except:
            msg = MyMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText("Values of Re and ncrit for new curve do not seem to be valid.")
            msg.setDetailedText(
                "Values of Re '%s' and ncrit '%s' could not be converted to float. Please double check the numbers." % (
                    re_chosen, ncrit_chosen))
            msg.exec_()
            return

        data_from_table = self.table.get_values()
        alpha_table, cl_table, cd_table = transpose(data_from_table)
        alpha, cl, cd = [], [], []

        for i in range(len(alpha_table)):
            if alpha_table[i] != "" and cl_table[i] != "" and cd_table[i] != "":
                alpha.append(float(alpha_table[i]))
                cl.append(float(cl_table[i]))
                cd.append(float(cd_table[i]))

        if self.parent.curves.get_curve(re_in=re_chosen, ncrit_in=ncrit_chosen) == None:
            print("item does not exist,creating new...")
            x, y = self.parent.get_x_y()
            self.current_curve = Curve()
            self.current_curve.create(x, y, re_chosen, ncrit_chosen, alpha, cl, cd)
            self.parent.curves.add(self.current_curve)
            print("self.parent.curves.curve_list", self.parent.curves.curve_list)
            self.load_curves()
        else:
            msg = MyMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText("Curve with this Re and ncrit already exists!")
            msg.setDetailedText(
                "Curve with Re %s and ncrit %s values already exists. Did you mean to update the existing curve?" % (
                    re_chosen, ncrit_chosen))
            msg.exec_()

        self.connect()

    def remove_curve(self):
        """

        :return:
        """
        if len(self.parent.curves.curve_list) == 0:
            return
        out_chosen_values = self.get_chosen_values_from_dropdowns()
        re_chosen, ncrit_chosen = out_chosen_values
        self.parent.curves.remove_curve(re_chosen, ncrit_chosen)
        self.load_curves()

    def get_chosen_values_from_dropdowns(self):
        """

        :return:
        """
        re_chosen = self.picker_reynolds.itemText(self.picker_reynolds.currentIndex())
        ncrit_chosen = self.picker_ncrit.itemText(self.picker_ncrit.currentIndex())
        return float(re_chosen), float(ncrit_chosen)

    def load_curves(self):
        """

        :return:
        """
        try:
            re_last, ncrit_last = self.get_chosen_values_from_dropdowns()
        except ValueError:
            re_last, ncrit_last = None, None

        if len(self.parent.curves.curve_list) == 0:
            self.disconnect()

            self.picker_reynolds.clear()
            self.picker_ncrit.clear()

            self.connect()
            return

        self.disconnect()

        self.picker_reynolds.clear()
        self.picker_ncrit.clear()

        ncrit_list = []
        for curve in self.parent.curves.curve_list:
            ncrit = curve.ncrit
            if not ncrit in ncrit_list:
                ncrit_list.append(ncrit)

        ncrit_list.sort()
        ncrit_list_str = [str(n) for n in ncrit_list]
        self.picker_ncrit.addItems(ncrit_list_str)

        if ncrit_last in ncrit_list:
            ncrit_index = ncrit_list.index(ncrit_last)
        else:
            ncrit_index = 0
            ncrit_last = ncrit_list[ncrit_index]

        self.picker_ncrit.setCurrentIndex(ncrit_index)

        re_list = []
        for curve in self.parent.curves.curve_list:
            re = curve.Re
            ncrit = curve.ncrit
            if not re in re_list:
                if ncrit == ncrit_last:
                    re_list.append(re)

        re_list.sort()
        re_list_str = [str(r) for r in re_list]
        self.picker_reynolds.addItems(re_list_str)

        if re_last in re_list:
            re_index = re_list.index(re_last)
        else:
            re_index = 0
            re_last = re_list[re_index]

        self.picker_reynolds.setCurrentIndex(re_index)

        self.connect()

        self.load_values_into_table()

    def disconnect(self):
        """

        """
        try:
            self.picker_reynolds.currentIndexChanged.disconnect()
            self.picker_ncrit.currentIndexChanged.disconnect()
        except TypeError:
            pass

    def connect(self):
        """

        """
        self.sig3 = self.picker_reynolds.currentIndexChanged.connect(self.load_curves)
        self.sig4 = self.picker_ncrit.currentIndexChanged.connect(self.load_curves)

    def load_values_into_table(self):
        """

        :return:
        """
        out_chosen_values = self.get_chosen_values_from_dropdowns()
        if out_chosen_values == None:
            return
        re_chosen, ncrit_chosen = out_chosen_values
        chosen_curve = self.parent.curves.get_curve(re_in=re_chosen, ncrit_in=ncrit_chosen)
        if chosen_curve != None:
            self.current_curve = chosen_curve
            array = [self.current_curve.alpha, self.current_curve.cl, self.current_curve.cd]
            array = transpose(array)
            self.table.clear_table()
            self.table.createTable(array)

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
        sender = self.sender()
        validator = sender.validator()
        state = validator.validate(sender.text(), 0)[0]
        if state == QtGui.QValidator.Acceptable:
            color = "#edf5e1"  # green
        elif state == QtGui.QValidator.Intermediate:
            color = "#fff79a"  # yellow
        else:
            color = "#f6989d"  # red

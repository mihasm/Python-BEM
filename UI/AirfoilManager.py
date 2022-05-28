from PyQt5.QtCore import pyqtSignal, QRect
from PyQt5.QtWidgets import QWidget, QGridLayout, QPushButton

from UI.helpers import PopupText, TabWidget
from UI.Airfoils import Airfoils


class AirfoilManager(QWidget):
    """

    """
    emitter = pyqtSignal(str)

    def __init__(self, parent=None):
        super(AirfoilManager, self).__init__(parent)

        self.main = self.parent()

        self.grid = QGridLayout()
        self.setLayout(self.grid)

        self.tab_widget = TabWidget(self)
        self.grid.addWidget(self.tab_widget, 2, 0)

        self.upper_widget = QWidget()
        self.upper_layout = QGridLayout()
        self.upper_widget.setLayout(self.upper_layout)
        self.grid.addWidget(self.upper_widget, 1, 0)

        self.button_add_foil = QPushButton("Add airfoil")
        self.button_add_foil.clicked.connect(self.add_foil_popup)
        self.button_remove_foil = QPushButton("Remove foil")
        self.button_remove_foil.clicked.connect(self.tab_widget.remove_current_tab)
        self.button_rename_foil = QPushButton("Rename foil")
        self.button_rename_foil.clicked.connect(self.rename_foil_popup)

        self.upper_layout.addWidget(self.button_add_foil, 0, 1)
        self.upper_layout.addWidget(self.button_remove_foil, 0, 2)
        self.upper_layout.addWidget(self.button_rename_foil, 0, 3)

    def add_foil_popup(self):
        """

        """
        self.emitter.connect(self.add_foil)
        self.p = PopupText("Add airfoil", "airfoil_name", self.emitter,"Airfoil name")
        self.p.show()

    def add_foil(self, string):
        """

        :param string:
        """
        c = Airfoils(string, self)
        self.tab_widget.add_tab(c, string)
        self.tab_widget.setCurrentIndex(len(self.tab_widget.tabs) - 1)

    def rename_foil_popup(self):
        """

        """
        self.emitter.connect(self.rename_foil)
        self.p = PopupText("Rename airfoil", self.tab_widget.current_tab_name(), self.emitter, "Rename airfoil")
        self.p.show()

    def rename_foil(self, string):
        """

        :param string:
        """
        self.tab_widget.rename_current_tab(string)  # self.tab_widget.tabs

    def get_settings(self):
        """

        :return:
        """
        out = {}
        i = 0

        # TODO Dont rely on name being set correctly in n!
        for w, n in self.tab_widget.tabs:
            out[n] = w.get_settings()
            i += 1

        return {"airfoils": out}

    def set_settings(self, dict_settings):
        """

        :param dict_settings:
        """
        self.tab_widget.remove_all_tabs()
        if "airfoils" in dict_settings:
            if len(dict_settings["airfoils"]) > 0:
                for c_name, c_dict in dict_settings["airfoils"].items():
                    curve_widget = Airfoils(c_name, self)
                    curve_widget.set_settings(c_dict)
                    self.tab_widget.add_tab(curve_widget, c_name)
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

from PyQt5.QtCore import pyqtSignal, QRect
from PyQt5.QtWidgets import QWidget, QGridLayout, QPushButton, QFileDialog

from UI.helpers import PopupText, TabWidget, PopupConfirmation
from UI.Airfoil import Airfoil
#from utils import fltr

import numpy as np
import json


class AirfoilManager(QWidget):
    """

    """
    emitter = pyqtSignal(str)
    emitter_yes = pyqtSignal(str)

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

        self.button_add_foil = QPushButton("Add foil")
        self.button_add_foil.clicked.connect(self.add_foil_popup)
        self.upper_layout.addWidget(self.button_add_foil, 0, 1)
        
        self.button_rename_foil = QPushButton("Rename foil")
        self.button_rename_foil.clicked.connect(self.rename_foil_popup)
        self.upper_layout.addWidget(self.button_rename_foil, 0, 2)
        
        self.button_duplicate_foil = QPushButton("Duplicate foil")
        self.button_duplicate_foil.clicked.connect(self.duplicate_foil)
        self.upper_layout.addWidget(self.button_duplicate_foil, 0, 3)
        
        self.button_remove_foil = QPushButton("Remove foil")
        self.button_remove_foil.clicked.connect(self.remove_foil_popup)
        self.upper_layout.addWidget(self.button_remove_foil, 0, 4)

        self.button_load_foil = QPushButton("Load foil")
        self.button_load_foil.clicked.connect(self.load_airfoil)
        self.upper_layout.addWidget(self.button_load_foil, 0, 5)

        self.button_save_foil = QPushButton("Save foil")
        self.button_save_foil.clicked.connect(self.save_airfoil)
        self.upper_layout.addWidget(self.button_save_foil, 0, 6)

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
        c = Airfoil(string, self)
        self.tab_widget.add_tab(c, string)
        self.tab_widget.setCurrentIndex(len(self.tab_widget.tabs) - 1)

    def remove_foil_popup(self):
        cur_widget, cur_name = self.tab_widget.tabs[self.tab_widget.currentIndex()]
        self.p = PopupConfirmation("Do you really want to delete "+cur_name+"?","Remove foil confirmation",self.emitter_yes)
        self.emitter_yes.connect(self.tab_widget.remove_current_tab)
        self.p.show()

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

    def duplicate_foil(self):
        cur_widget, cur_name = self.tab_widget.tabs[self.tab_widget.currentIndex()]
        names = [name for _,name in self.tab_widget.tabs]
        i=0
        while i<100:
            new_name = cur_name+"_"+str(i)
            if not new_name in names:
                break
            i+=1
        old_settings = cur_widget.get_settings()
        new_airfoil = Airfoil(new_name, self)
        new_airfoil.set_settings(old_settings)
        self.tab_widget.add_tab(new_airfoil,new_name)

    def load_airfoil(self):
        file_path = QFileDialog.getOpenFileName(self, "Load File", "", "BEM airfoil (*.bemfoil)")[0]
        if file_path != "":
            with open(file_path, "r") as fp:
                data = json.load(fp)
            airfoil_name,airfoil_settings = list(data.items())[0]
            new_airfoil = Airfoil(airfoil_name, self)
            new_airfoil.set_settings(airfoil_settings)
            self.tab_widget.add_tab(new_airfoil,airfoil_name)

    def save_airfoil(self):
        name = QFileDialog.getSaveFileName(self, 'Save File', "", "BEM airfoil (*.bemfoil)")[0]
        if name != "":
            cur_widget, cur_name = self.tab_widget.tabs[self.tab_widget.currentIndex()]
            cur_airfoil_settings = cur_widget.get_settings()
            d = {cur_name:cur_airfoil_settings}
            d_to_save = fltr(d, (float, int, list, str, bool, np.ndarray))
            json_d = json.dumps(d_to_save)
            file = open(name, 'w')
            file.write(json_d)
            file.close()


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
                    curve_widget = Airfoil(c_name, self)
                    curve_widget.set_settings(c_dict)
                    self.tab_widget.add_tab(curve_widget, c_name)
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

from PyQt5.QtWidgets import QWidget, QGridLayout, QPushButton, QScrollArea, QVBoxLayout, QLabel

from UI.CurveControl import CurveControl


class CurveExtrapolationEditor(QWidget):
    """

    """
    def __init__(self, parent=None):
        super(CurveExtrapolationEditor, self).__init__(None)
        self.resize(1600, 768)
        self.parent = parent

        self.setWindowTitle("Curve Extrapolation Editor")

        self.grid = QGridLayout()
        self.setLayout(self.grid)
        self.button = QPushButton("Close")
        self.grid.addWidget(self.button, 1, 1)
        self.button.clicked.connect(self.close)
        self.button_refresh = QPushButton("Refresh")
        self.grid.addWidget(self.button_refresh, 1, 2)
        self.button_refresh.clicked.connect(self.generate_views)

        self.bottom = QWidget()
        self.grid_curves = QGridLayout()
        self.bottom.setLayout(self.grid_curves)

        self.scroll_area = QScrollArea()
        self.scroll_widget = QWidget()
        self.scroll_widget_layout = QVBoxLayout()

        self.scroll_widget.setLayout(self.scroll_widget_layout)
        self.scroll_area.setWidget(self.bottom)
        self.scroll_area.setWidgetResizable(True)
        self.grid.addWidget(self.scroll_area, 2, 1, 2, 2)

        # self.generate_views()

    def generate_views(self):
        """

        """
        # delete stuff already here
        for i in reversed(range(self.grid_curves.count())):
            self.grid_curves.itemAt(i).widget().setParent(None)
        
        for curve in self.parent.curves.get_curves_sorted():
            label = QLabel("Re = " + str(round(curve.Re, 2)) + ", Ncrit = " + str(round(curve.ncrit, 2)))
            control = CurveControl(self, curve)
            control.update()
            self.grid_curves.addWidget(label)
            self.grid_curves.addWidget(control)
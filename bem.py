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

# coding=utf-8
__author__ = "Miha Smrekar"
__credits__ = ["Miha Smrekar"]
__license__ = "GPL"
__version__ = "0.4.7"
__maintainer__ = "Miha Smrekar"
__email__ = "xmiha.xsmrekar@gmail.com"
__status__ = "Development"

import ctypes
import multiprocessing
import os
import sys

import numpy as np
from PyQt5 import QtCore
from PyQt5 import QtGui
from PyQt5 import QtWidgets
from PyQt5.QtCore import QLocale
from PyQt5.QtWidgets import (QApplication)

from UI import MainWindow
from utils import (QDarkPalette)

np.set_printoptions(threshold=sys.maxsize)

TITLE_STR = "BEM analiza v%s" % __version__

# determine if application is a script file or frozen exe
if getattr(sys, 'frozen', False):
    application_path = os.path.dirname(sys.executable)
elif __file__:
    application_path = os.path.dirname(__file__)

# determine if application is a one-file exe
if hasattr(sys, "_MEIPASS"):
    # yes, resources are stored in temporary folder C:\TEMP or wherever it is
    data_path = sys._MEIPASS
else:
    # else, resources are stored in same folder as executable
    data_path = application_path

ICON_PATH = os.path.join(data_path, "icon_bem.ico")
DEFAULT_SETTINGS_PATH = os.path.join(data_path,"karlsen.bem")

QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)


def main(quick_results=False):
    """
    :param quick_results:
    """
    if sys.platform.startswith("win"):
        # On Windows calling this function is necessary for multiprocessing.
        multiprocessing.freeze_support()
        # To show icon in taskbar
        myappid = 'FEUM.BEM_Analiza.%s' % __version__  # arbitrary string
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    app = QApplication([])
    QLocale.setDefault(QLocale(QLocale.English))  # da je pika decimalno mesto
    app.setWindowIcon(QtGui.QIcon(ICON_PATH))
    app.setStyle("Fusion")
    if sys.platform.startswith("darwin") or True:
        # dark theme fix on OSX
        palette = QDarkPalette()
        palette.set_app(app)
        palette.set_stylesheet(app)
    screen = app.primaryScreen()
    size = screen.size()
    main_win = MainWindow.MainWindow(size.width(), size.height())
    font = main_win.font()
    font.setPointSize(7)
    app.instance().setFont(font)
    if quick_results:
        main_win.analysis.run()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()

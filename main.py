# coding=utf-8
__author__ = "Miha Smrekar"
__credits__ = ["Miha Smrekar"]
__license__ = "GPL"
__version__ = "0.4.5"
__maintainer__ = "Miha Smrekar"
__email__ = "miha.smrekar9@gmail.com"
__status__ = "Development"

import ctypes
import multiprocessing
import os
import sys

import numpy as np
from PyQt5 import QtGui
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
    app_icon = QtGui.QIcon("icon_bem.ico")
    app.setWindowIcon(app_icon)
    app.setStyle("Fusion")
    if sys.platform.startswith("darwin") or True:
        # dark theme fix on OSX
        palette = QDarkPalette()
        palette.set_app(app)
        palette.set_stylesheet(app)
    screen = app.primaryScreen()
    size = screen.size()
    main_win = MainWindow.MainWindow(size.width(), size.height())
    main_win.setWindowIcon(app_icon)
    if quick_results:
        main_win.analysis.run()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()

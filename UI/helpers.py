import sys

import pyqtgraph as pg
from PyQt5 import QtCore
from PyQt5.QtCore import QThread
from PyQt5.QtGui import QTextCursor
from PyQt5.QtWidgets import QMainWindow, QWidget, QGridLayout, QLabel, QLineEdit, QPushButton, QTabWidget, \
    QTextEdit, QFormLayout
from matplotlib import pyplot as plt
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

from scraping import scrape_data, get_x_y_from_link
from utils import generate_propeller_adkins
from xfoil import generate_polars_data


class ThreadGetter(QThread):
    """

    """
    def __init__(self, parent):
        super(ThreadGetter, self).__init__(parent)
        self.dataCollectionTimer = QtCore.QTimer()
        self.dataCollectionTimer.moveToThread(self)
        self.dataCollectionTimer.timeout.connect(self.updateInProc)

    def run(self):
        """

        """
        self.dataCollectionTimer.start(2)  # 0 causes freeze
        self.loop = QtCore.QEventLoop()
        self.loop.exec_()

    def updateInProc(self):
        """

        """
        if len(self.parent().return_print) > 0:
            t = self.parent().return_print.pop(0)
            self.parent().emitter_add.emit(str(t))
        if self.parent().end_of_file.value == True and len(self.parent().return_print) == 0:
            self.parent().emitter_done.emit()


class DataCaptureThread(QThread):
    """

    """
    def __init__(self, parent, *args, **kwargs):
        QThread.__init__(self, parent, *args, **kwargs)
        self.dataCollectionTimer = QtCore.QTimer()
        self.dataCollectionTimer.moveToThread(self)
        self.dataCollectionTimer.timeout.connect(self.updateInProc)

    def run(self):
        """

        """
        self.dataCollectionTimer.start(50)  # 0 causes freeze
        self.loop = QtCore.QEventLoop()
        self.loop.exec_()

    def updateInProc(self):
        """

        """
        if len(self.parent().parent.queue_pyqtgraph) > 0:
            item = self.parent().parent.queue_pyqtgraph[0]
            x = item[0]
            y = item[1]
            best_x = [item[2]]
            best_y = [item[3]]
            self.parent().curve.setData(x, y)
            self.parent().curve_red.setData(best_x, best_y)


class XFoilThread(QThread):
    """

    """
    startedSignal = QtCore.Signal(int)
    completeSignal = QtCore.Signal(str)

    def __init__(self, parent, *args, **kwargs):
        QThread.__init__(self, parent)
        self.parent = parent

    def set_params(self, dat_path,
                   alpha_from, alpha_to, alpha_num,
                   reynolds_from, reynolds_to, reynolds_num,
                   ncrit):
        """

        :param dat_path:
        :param alpha_from:
        :param alpha_to:
        :param alpha_num:
        :param reynolds_from:
        :param reynolds_to:
        :param reynolds_num:
        :param ncrit:
        """
        self.dat_path = dat_path
        self.alpha_from = alpha_from
        self.alpha_to = alpha_to
        self.alpha_num = alpha_num
        self.reynolds_from = reynolds_from
        self.reynolds_to = reynolds_to
        self.reynolds_num = reynolds_num
        self.ncrit = ncrit

    def run(self):
        """

        """
        self.startedSignal.emit("Started")
        out = generate_polars_data(
            self.dat_path,
            self.alpha_from,
            self.alpha_to,
            self.alpha_num,
            self.reynolds_from,
            self.reynolds_to,
            self.reynolds_num,
            self.ncrit
        )
        self.parent.xfoil_generated_data = out
        self.completeSignal.emit("Done")


class XfoilOptionsWindow(QWidget):
    """

    """
    def __init__(self, parent):
        super(XfoilOptionsWindow, self).__init__(None)
        self.parent = parent
        self.layout = QFormLayout()
        self.setLayout(self.layout)

        self.alpha_from_label = QLabel("Alpha from")
        self.alpha_from = QLineEdit("-10")
        self.layout.addRow(self.alpha_from_label, self.alpha_from)

        self.alpha_to_label = QLabel("Alpha to")
        self.alpha_to = QLineEdit("20")
        self.layout.addRow(self.alpha_to_label, self.alpha_to)

        self.alpha_num_label = QLabel("Alpha number")
        self.alpha_num = QLineEdit("31")
        self.layout.addRow(self.alpha_num_label, self.alpha_num)

        self.reynolds_from_label = QLabel("Reynolds from")
        self.reynolds_from = QLineEdit("50000")
        self.layout.addRow(self.reynolds_from_label, self.reynolds_from)

        self.reynolds_to_label = QLabel("Reynolds to")
        self.reynolds_to = QLineEdit("1000000")
        self.layout.addRow(self.reynolds_to_label, self.reynolds_to)

        self.reynolds_num_label = QLabel("Reynolds number")
        self.reynolds_num = QLineEdit("5")
        self.layout.addRow(self.reynolds_num_label, self.reynolds_num)

        self.ncrit_label = QLabel("Ncrit")
        self.ncrit = QLineEdit("9")
        self.layout.addRow(self.ncrit_label, self.ncrit)

        self.button_run_xfoil = QPushButton("Run")
        self.button_run_xfoil.clicked.connect(self.parent.generate_curves_xfoil)
        self.layout.addRow(self.button_run_xfoil)

        self.button_stop_xfoil = QPushButton("Stop")
        self.button_stop_xfoil.clicked.connect(self.parent.stop_xfoil)
        self.layout.addRow(self.button_stop_xfoil)
        self.button_stop_xfoil.setDisabled(True)

        self.show()

    def closeEvent(self, event):
        """

        :param event:
        """
        event.accept()  # let the window close


class AdkinsThread(QThread):
    """

    """
    progressSignal = QtCore.Signal(int)
    completeSignal = QtCore.Signal(str)

    def __init__(self, parent, *args, **kwargs):
        QThread.__init__(self, parent)
        self.parent = parent

    def set_params(self, inp):
        """

        :param inp:
        """
        self.inp = inp

    def run(self):
        """

        """
        out = generate_propeller_adkins(self.inp)
        self.parent.adkins_return = out
        self.completeSignal.emit("Done")


class ScrapeThread(QThread):
    """

    """
    progressSignal = QtCore.Signal(int)
    completeSignal = QtCore.Signal(str)

    def __init__(self, parent, *args, **kwargs):
        QThread.__init__(self, parent)
        self.parent = parent

    def set_params(self, link):
        """

        :param link:
        """
        self.link = link

    def run(self):
        """

        """
        print("Running scrape")
        data = scrape_data(self.link.text())
        x, y = get_x_y_from_link(self.link.text())
        out = [data, x, y]
        self.parent.scraping_generated_data = out
        self.completeSignal.emit("Done")
        print("Done scraping thread")


class PyQtGraphWindow(QMainWindow):
    """

    """
    def __init__(self, parent):
        super(PyQtGraphWindow, self).__init__(parent)
        self.obj = pg.PlotWidget()
        self.setCentralWidget(self.obj)
        self.curve = self.obj.plot(pen=None, symbol='o', symbolPen=None, symbolSize=4, symbolBrush='g')
        self.curve_red = self.obj.plot(pen=None, symbol='o', symbolPen=None, symbolSize=5, symbolBrush='r')
        self.obj.setLabel("left", "Optimization variable")
        self.obj.setLabel("bottom", "Theta [°]")
        self.parent = parent
        self.thread = DataCaptureThread(self)

    def start_update(self):
        """

        """
        self.thread.start()

    def stop_update(self):
        """

        """
        self.thread.quit()


class PopupText(QWidget):
    """

    """
    def __init__(self, message="message", default_str="", emitter=None):
        QWidget.__init__(self)

        self.emitter = emitter

        self.layout = QGridLayout()
        self.setLayout(self.layout)

        self.message = QLabel(message)

        self.layout.addWidget(self.message, 0, 0)

        self.inp = QLineEdit()
        self.inp.setText(default_str)
        self.layout.addWidget(self.inp, 1, 0)

        self.button = QPushButton("OK")
        self.button.clicked.connect(self.send_signal)
        self.layout.addWidget(self.button, 2, 0)

    def send_signal(self):
        """

        """
        if self.emitter != None:
            self.emitter.emit(self.inp.text())
        self.close()

    def closeEvent(self, event):
        """

        :param event:
        """
        self.emitter.disconnect()
        event.accept()


class MatplotlibWindow(QWidget):
    """

    """
    def __init__(self):
        super(MatplotlibWindow, self).__init__(None)
        self.layout = QGridLayout()
        self.setLayout(self.layout)
        self.figure = plt.figure(figsize=(10, 5))
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setMinimumSize(500, 500)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.layout.addWidget(self.canvas)
        self.layout.addWidget(self.toolbar)
        self.show()

    def closeEvent(self, event):
        """

        :param event:
        """
        self.figure.clear()
        plt.close(self.figure)
        event.accept()  # let the window close


class PrintoutWindow(QMainWindow):
    """

    """
    def __init__(self, parent):
        super(PrintoutWindow, self).__init__(parent)
        self.setWindowTitle("Progress")
        self.setGeometry(50, 50, 500, 300)
        self.parent = parent
        sys.stdout = Stream(newText=self.onUpdateText)
        sys.stderr = Stream(newText=self.onUpdateText)
        self.process = QTextEdit()
        self.setCentralWidget(self.process)
        self.show()

    def onUpdateText(self, text):
        """

        :param text:
        """
        cursor = self.process.textCursor()
        cursor.movePosition(QTextCursor.End)
        cursor.insertText(text)
        self.process.setTextCursor(cursor)
        self.process.ensureCursorVisible()

    def closeEvent(self, event):
        """

        :param event:
        """
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        QMainWindow.closeEvent(self, event)


class Stream(QtCore.QObject):
    newText = QtCore.pyqtSignal(str)

    def write(self, text):
        """

        :param text:
        """
        self.newText.emit(str(text))


class TabWidget(QTabWidget):
    """

    """
    def __init__(self, parent=None):
        super(TabWidget, self).__init__(parent)
        self.tabs = []

    def add_tab(self, widget, tab_name, tooltip=None):
        """

        :param widget:
        :param tab_name:
        :param tooltip:
        """
        self.tabs.append([widget, tab_name])
        self.addTab(widget, tab_name)
        if tooltip != None:
            self.setTabToolTip(len(self.tabs) - 1, tooltip)

    def remove_tab(self, index):
        """

        :param index:
        """
        self.removeTab(index)
        del self.tabs[index]

    def remove_all_tabs(self):
        """

        """
        while len(self.tabs) > 0:
            self.remove_tab(0)

    def remove_current_tab(self):
        """

        """
        self.remove_tab(self.currentIndex())

    def rename_current_tab(self, string):
        """

        :param string:
        """
        self.setTabText(self.currentIndex(), string)
        self.tabs[self.currentIndex()][1] = string

    def current_tab_name(self):
        """

        :return:
        """
        return self.tabText(self.currentIndex())

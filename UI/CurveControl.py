from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget, QGridLayout, QFormLayout, QLineEdit, QSlider, QPushButton
from matplotlib import pyplot as plt
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas


class CurveControl(QWidget):
    def __init__(self, parent=None, curve=None):
        super(CurveControl, self).__init__(parent)
        # self.setMinimumSize(300,400)
        self.parent = parent

        self.layout = QGridLayout()
        self.setLayout(self.layout)

        self.curve = curve

        self.right = QWidget()
        self.right_layout = QFormLayout()
        self.right.setLayout(self.right_layout)

        self.A = QLineEdit(str(self.curve.A))
        self.B = QLineEdit(str(self.curve.B))
        self.Am = QLineEdit(str(self.curve.Am))
        self.Bm = QLineEdit(str(self.curve.Bm))

        self.A = QSlider(Qt.Horizontal)
        self.A.setMinimum(-10)
        self.A.setMaximum(30)
        self.A.setValue(self.curve.A)
        self.A.setTickPosition(QSlider.TicksBelow)
        self.A.setTickInterval(1)
        self.A.valueChanged.connect(self.update)

        self.B = QSlider(Qt.Horizontal)
        self.B.setMinimum(1)
        self.B.setMaximum(100)
        self.B.setValue(self.curve.B)
        self.B.setTickPosition(QSlider.TicksBelow)
        self.B.setTickInterval(1)
        self.B.valueChanged.connect(self.update)

        self.Am = QSlider(Qt.Horizontal)
        self.Am.setMinimum(1)
        self.Am.setMaximum(80)
        self.Am.setValue(self.curve.Am)
        self.Am.setTickPosition(QSlider.TicksBelow)
        self.Am.setTickInterval(1)
        self.Am.valueChanged.connect(self.update)

        self.Bm = QSlider(Qt.Horizontal)
        self.Bm.setMinimum(1)
        self.Bm.setMaximum(70)
        self.Bm.setValue(self.curve.Bm)
        self.Bm.setTickPosition(QSlider.TicksBelow)
        self.Bm.setTickInterval(1)
        self.Bm.valueChanged.connect(self.update)

        self.m_CD90 = QLineEdit(str(self.curve.m_CD90))
        self.m_CD90.textChanged.connect(self.update)
        self.slope = QLineEdit(str(self.curve.slope))
        self.slope.textChanged.connect(self.update)

        self.right_layout.addRow("A", self.A)
        self.right_layout.addRow("B", self.B)
        self.right_layout.addRow("A-", self.Am)
        self.right_layout.addRow("B-", self.Bm)
        self.right_layout.addRow("CD@90°", self.m_CD90)
        self.right_layout.addRow("Slope", self.slope)

        self.layout.addWidget(self.right, 1, 2)

        self.left = QWidget()
        self.left_layout = QGridLayout()
        self.left.setLayout(self.left_layout)

        self.figure = plt.figure(figsize=(10, 5))
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setMinimumSize(500, 500)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.ax = self.figure.add_subplot(111)
        self.left_layout.addWidget(self.canvas)
        self.left_layout.addWidget(self.toolbar)

        self.layout.addWidget(self.left, 1, 1)

        self.button_update = QPushButton("update")
        self.button_update.clicked.connect(self.draw_extrapolation)
        self.left_layout.addWidget(self.button_update)

        # self.draw_base()

        self.show()

    def clear(self):
        self.ax.cla()

    def draw_base(self):
        self.ax.plot(self.curve.alpha, self.curve.cl)
        self.ax.plot(self.curve.alpha, self.curve.cd, "o-")
        self.canvas.draw()

    def draw_extrapolation(self):
        self.clear()

        self.draw_base()

        alpha, cl, cd = self.curve.get_extrapolated_curve()
        self.ax.plot(alpha, cl, "g.")
        self.ax.plot(alpha, cd, "r.")
        self.canvas.draw()

    def update(self):
        self.curve.A = int(self.A.value())
        self.curve.B = int(self.B.value())
        self.curve.Am = int(self.Am.value())
        self.curve.Bm = int(self.Bm.value())
        # print("A",self.curve.A,"B",self.curve.B,"A-",self.curve.Am,"B-",self.curve.Bm)
        try:
            self.curve.m_CD90 = float(self.m_CD90.text())
            self.curve.slope = float(self.slope.text())
            if self.curve.slope == 0:
                self.curve.slope = 1.0
        except:
            print("Error in slope or m_CD90")
        self.draw_extrapolation()
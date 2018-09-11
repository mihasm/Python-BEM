import sys
from PyQt5.QtWidgets import QMainWindow, QApplication, QWidget, QAction, QTableWidget, QTableWidgetItem, QVBoxLayout
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot
from PyQt5 import QtGui, QtCore
from utils import array_to_csv


# noinspection PyArgumentList
class Table(QWidget):

    def __init__(self, inp_ar):
        super().__init__()
        self.selected_array = []
        self.tableWidget = QTableWidget()
        self.layout = QVBoxLayout()
        self.title = 'Data'
        self.left = 300
        self.top = 300
        self.width = 1280
        self.height = 800
        self.array = inp_ar
        self.initUI()
        self.clip = QApplication.clipboard()

    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        self.createTable()

        # Add box layout, add table to box layout and add box layout to widget
        self.layout.addWidget(self.tableWidget)
        self.setLayout(self.layout)

        # Show widget
        self.show()

    def createTable(self):
        # Create table
        self.tableWidget.setRowCount(len(self.array))
        self.tableWidget.setColumnCount(len(self.array[0]))
        i = 0
        for r in self.array:
            j = 0
            for c in r:
                if not isinstance(c,str):
                    c=str(c)
                self.tableWidget.setItem(i, j, QTableWidgetItem(c))
                j += 1
            i += 1

        self.tableWidget.move(0, 0)

        # table selection change
        self.tableWidget.clicked.connect(self.on_click)

    @pyqtSlot()
    def on_click(self):
        self.get_selected()

    def get_selected(self):
        self.selected_array=[]

        rows_added = []
        columns_added = []
        for currentQTableWidgetItem in self.tableWidget.selectedItems():
            row = currentQTableWidgetItem.row()
            column = currentQTableWidgetItem.column()
            if row not in rows_added:
                rows_added.append(row)
            if column not in columns_added:
                columns_added.append(column)

        remap_row = {}
        remap_column = {}
        row_count = 0
        column_count = 0

        for r in range(len(rows_added)):
            self.selected_array.append([])
            for c in range(len(columns_added)):
                self.selected_array[r].append(None)
                if columns_added[c] not in remap_column:
                    remap_column[columns_added[c]] = column_count
                    column_count += 1
            remap_row[rows_added[r]] = row_count
            row_count += 1

        for currentQTableWidgetItem in self.tableWidget.selectedItems():
            row = currentQTableWidgetItem.row()
            column = currentQTableWidgetItem.column()
            text = currentQTableWidgetItem.text()
            self.selected_array[remap_row[row]][remap_column[column]] = text
        #print(self.selected_array)
        return self.selected_array

    def keyPressEvent(self, e):
        if e.modifiers() & QtCore.Qt.ControlModifier:
            if e.key() == QtCore.Qt.Key_C:  # copy
                s = array_to_csv(self.get_selected())
                self.clip.setText(s)


def draw_table(inp_ar):
    app = QApplication(sys.argv)
    ex = Table(inp_ar)
    sys.exit(app.exec_())

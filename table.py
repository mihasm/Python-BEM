__author__ = "Miha Smrekar"
__credits__ = ["Miha Smrekar"]
__license__ = "GPL"
__version__ = "0.3.0"
__maintainer__ = "Miha Smrekar"
__email__ = "miha.smrekar9@gmail.com"
__status__ = "Development"

import csv
import sys

from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import (
    QApplication, QWidget, QTableWidget, QTableWidgetItem, QVBoxLayout, QMenu, )

from utils import array_to_csv


# noinspection PyArgumentList
class Table(QWidget):
    def __init__(self):
        super().__init__()
        self.selected_array = []
        self.tableWidget = QTableWidget()
        self.tableWidget.setTabKeyNavigation(False)
        self.layout = QVBoxLayout()
        self.initUI()
        self.clip = QApplication.clipboard()
        self.set_headers()

    def set_headers(self):
        self.horizontal_headers = self.tableWidget.horizontalHeader()
        self.horizontal_headers.setContextMenuPolicy(
            QtCore.Qt.CustomContextMenu)
        self.horizontal_headers.customContextMenuRequested.connect(
            self.horizontal_header_popup)
        self.vertical_headers = self.tableWidget.verticalHeader()
        self.vertical_headers.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.vertical_headers.customContextMenuRequested.connect(
            self.vertical_header_popup)

    def set_labels(self, arr):
        self.tableWidget.setHorizontalHeaderLabels(arr)

    def initUI(self):
        self.createEmpty(4, 4)

        # Add box layout, add table to box layout and add box layout to widget
        self.layout.addWidget(self.tableWidget)
        self.setLayout(self.layout)

        # Show widget
        self.show()

    def createTable(self, array):
        # Create table
        self.tableWidget.setRowCount(len(array))
        self.tableWidget.setColumnCount(len(array[0]))
        i = 0
        for r in array:
            j = 0
            for c in r:
                if not isinstance(c, str):
                    c = str(c)
                self.tableWidget.setItem(i, j, QTableWidgetItem(c))
                j += 1
            i += 1

        self.tableWidget.move(0, 0)

        # table selection change
        self.tableWidget.clicked.connect(self.on_click)

    def createEmpty(self, x, y):
        # Create table
        self.tableWidget.setRowCount(y)
        self.tableWidget.setColumnCount(x)

        self.tableWidget.move(0, 0)

        # table selection change
        self.tableWidget.clicked.connect(self.on_click)

    @pyqtSlot()
    def on_click(self):
        self.get_selected()

    def get_selected(self):
        self.selected_array = []

        rows_added = sorted(set(index.row()
                                for index in self.tableWidget.selectedIndexes()))
        columns_added = sorted(set(index.column()
                                   for index in self.tableWidget.selectedIndexes()))

        delta_r = rows_added[0]
        delta_c = columns_added[0]

        for r in range(rows_added[-1] - rows_added[0] + 1):
            self.selected_array.append([])
            for c in range(columns_added[-1] - columns_added[0] + 1):
                self.selected_array[r].append(None)
        for currentQTableWidgetItem in self.tableWidget.selectedItems():
            row = currentQTableWidgetItem.row() - delta_r
            column = currentQTableWidgetItem.column() - delta_c
            text = currentQTableWidgetItem.text()
            self.selected_array[row][column] = text
        return self.selected_array

    def keyPressEvent(self, e):
        if e.modifiers() & QtCore.Qt.ControlModifier:
            if e.key() == QtCore.Qt.Key_C:  # copy
                s = array_to_csv(self.get_selected())
                self.clip.setText(s)
        if e.key() == QtCore.Qt.Key_Return or e.key() == QtCore.Qt.Key_Enter:
            self.select_next_row()
        if e.modifiers() & QtCore.Qt.ControlModifier:
            if e.key() == QtCore.Qt.Key_V:  # paste
                self.paste()

        if e.modifiers() & QtCore.Qt.ControlModifier:
            if e.key() == QtCore.Qt.Key_S:  # test
                self.get_values()

        if e.key() == QtCore.Qt.Key_Delete:
            self.delete_data()

    def paste(self):
        results = []
        text = self.clip.text()
        text = text.replace("   ", "\t")
        text = text.replace("  ", "\t")
        if len(text) > 0:
            # change contents to floats
            reader = csv.reader(text.splitlines(), delimiter="\t")
            for row in reader:  # each row is a list
                results.append(row)
            numrows = len(results)
            numcolumns = len(results[0])
            selected_row = sorted(
                set(index.row() for index in self.tableWidget.selectedIndexes()))[0]
            selected_column = sorted(
                set(index.column() for index in self.tableWidget.selectedIndexes()))[0]
            if selected_row + numrows >= self.tableWidget.rowCount():
                self.tableWidget.setRowCount(selected_row + numrows)
            if selected_column + numcolumns >= self.tableWidget.columnCount():
                self.tableWidget.setColumnCount(selected_column + numcolumns)
            currow = selected_row
            for r in results:
                curcolumn = selected_column
                for c in r:
                    self.tableWidget.setItem(
                        currow, curcolumn, QTableWidgetItem(c))
                    curcolumn += 1
                currow += 1
        return

    def delete_data(self):
        rows = sorted(set(index.row()
                          for index in self.tableWidget.selectedIndexes()))
        columns = sorted(set(index.column()
                             for index in self.tableWidget.selectedIndexes()))
        for r in rows:
            for c in columns:
                self.tableWidget.setItem(r, c, QTableWidgetItem(""))

    def select_next_row(self):
        rows = sorted(set(index.row()
                          for index in self.tableWidget.selectedIndexes()))
        columns = sorted(set(index.column()
                             for index in self.tableWidget.selectedIndexes()))
        last_selected_row = rows[-1]
        first_selected_column = columns[0]
        num_rows = self.tableWidget.rowCount()
        if last_selected_row + 1 >= num_rows:
            self.tableWidget.insertRow(num_rows)
        self.tableWidget.setCurrentCell(
            last_selected_row + 1, first_selected_column)

    def get_values(self):
        data = []
        for row in range(self.tableWidget.rowCount()):
            data.append([])
            for column in range(self.tableWidget.columnCount()):
                item = self.tableWidget.item(row, column)
                if item == None:
                    item = ""
                else:
                    item = item.text()
                data[row].append(item)
        return data

    def contextMenuEvent(self, event):
        menu = QMenu(self)
        item = self.tableWidget.itemAt(event.pos())
        if item != None:
            delete_row = menu.addAction("delete row(s)")
            delete_column = menu.addAction("delete column(s)")
        else:
            delete_row = False
            delete_column = False

        insert_row = menu.addAction("insert row")
        insert_column = menu.addAction("insert column")

        action = menu.exec_(self.mapToGlobal(event.pos()))

        rows = sorted(set(index.row()
                          for index in self.tableWidget.selectedIndexes()), reverse=True, )
        columns = sorted(set(index.column()
                             for index in self.tableWidget.selectedIndexes()), reverse=True, )

        if action == insert_row:
            if len(rows) == 0:
                rows.append(0)
            self.tableWidget.insertRow(rows[-1])
            for c in range(self.tableWidget.columnCount()):
                self.tableWidget.setItem(rows[-1], c, QTableWidgetItem(""))
        elif action == delete_row:
            for r in rows:
                self.tableWidget.removeRow(r)
        elif action == insert_column:
            if len(columns) == 0:
                columns.append(0)
            self.tableWidget.insertColumn(columns[-1])
            for r in range(self.tableWidget.rowCount()):
                self.tableWidget.setItem(r, columns[-1], QTableWidgetItem(""))
        elif action == delete_column:
            for c in columns:
                self.tableWidget.removeColumn(c)

    def horizontal_header_popup(self, position):
        pass

    def vertical_header_popup(self, position):
        pass


def draw_table():
    app = QApplication(sys.argv)
    ex = Table()
    ex.setGeometry(300, 300, 1280, 800)
    sys.exit(app.exec_())

# draw_table()

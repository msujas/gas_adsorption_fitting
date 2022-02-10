# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'adsorption_fit_gui.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets
import Ar_adsorption_fit_class

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(800, 380)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.inputfileLabel = QtWidgets.QLabel(self.centralwidget)
        self.inputfileLabel.setGeometry(QtCore.QRect(500, 30, 61, 21))
        self.inputfileLabel.setObjectName("inputfileLabel")
        self.InputFile = QtWidgets.QTextEdit(self.centralwidget)
        self.InputFile.setGeometry(QtCore.QRect(20, 30, 471, 21))
        self.InputFile.setObjectName("InputFile")
        self.InputType = QtWidgets.QComboBox(self.centralwidget)
        self.InputType.setGeometry(QtCore.QRect(50, 110, 171, 22))
        self.InputType.setEditable(True)
        self.InputType.setObjectName("InputType")
        self.InputType.addItem("")
        self.InputType.addItem("")
        self.InputTypeLabel = QtWidgets.QLabel(self.centralwidget)
        self.InputTypeLabel.setGeometry(QtCore.QRect(230, 100, 101, 31))
        self.InputTypeLabel.setObjectName("InputTypeLabel")
        self.MinMaxRefine = QtWidgets.QComboBox(self.centralwidget)
        self.MinMaxRefine.setGeometry(QtCore.QRect(50, 200, 161, 22))
        self.MinMaxRefine.setObjectName("MinMaxRefine")
        self.MinMaxRefine.addItem("")
        self.MinMaxRefine.addItem("")
        self.MinMaxLabel = QtWidgets.QLabel(self.centralwidget)
        self.MinMaxLabel.setGeometry(QtCore.QRect(220, 200, 151, 21))
        self.MinMaxLabel.setScaledContents(False)
        self.MinMaxLabel.setObjectName("MinMaxLabel")
        self.ModelType = QtWidgets.QComboBox(self.centralwidget)
        self.ModelType.setGeometry(QtCore.QRect(50, 280, 161, 22))
        self.ModelType.setObjectName("ModelType")
        self.ModelType.addItem("")
        self.ModelType.addItem("")
        self.ModelType.addItem("")
        self.ModelTypeLabel = QtWidgets.QLabel(self.centralwidget)
        self.ModelTypeLabel.setGeometry(QtCore.QRect(220, 280, 111, 21))
        self.ModelTypeLabel.setObjectName("ModelTypeLabel")
        self.optimiseButton = QtWidgets.QPushButton(self.centralwidget)
        self.optimiseButton.setGeometry(QtCore.QRect(430, 290, 91, 23))
        self.optimiseButton.setObjectName("optimiseButton")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 21))
        self.menubar.setObjectName("menubar")
        self.menufile = QtWidgets.QMenu(self.menubar)
        self.menufile.setObjectName("menufile")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionopen = QtWidgets.QAction(MainWindow)
        self.actionopen.setObjectName("actionopen")
        self.menufile.addAction(self.actionopen)
        self.menubar.addAction(self.menufile.menuAction())

        self.retranslateUi(MainWindow)
        self.InputType.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.inputfileLabel.setText(_translate("MainWindow", "Input File"))
        self.InputFile.setHtml(_translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">C:\\Users\\kenneth1a\\Documents\\beamline data\\Yaroslav_Iurii_data\\Ar_occupancy_1bar.txt</p></body></html>"))
        self.InputType.setCurrentText(_translate("MainWindow", "Two sites"))
        self.InputType.setItemText(0, _translate("MainWindow", "Single site"))
        self.InputType.setItemText(1, _translate("MainWindow", "Two sites"))
        self.InputTypeLabel.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:14pt;\">input type</span></p></body></html>"))
        self.MinMaxRefine.setItemText(0, _translate("MainWindow", "Off"))
        self.MinMaxRefine.setItemText(1, _translate("MainWindow", "On"))
        self.MinMaxLabel.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:14pt;\">min. max. refine</span></p></body></html>"))
        self.ModelType.setItemText(0, _translate("MainWindow", "Simple single site"))
        self.ModelType.setItemText(1, _translate("MainWindow", "Single site cooperative"))
        self.ModelType.setItemText(2, _translate("MainWindow", "Two sites, intra-pore interaction"))
        self.ModelTypeLabel.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:14pt;\">Model type</span></p></body></html>"))
        self.optimiseButton.setText(_translate("MainWindow", "Run Fit"))
        self.menufile.setTitle(_translate("MainWindow", "File"))
        self.actionopen.setText(_translate("MainWindow", "Open"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

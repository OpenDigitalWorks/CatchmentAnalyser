# -*- coding: utf-8 -*-
"""
/***************************************************************************
 CatchmentAnalyserDialog
                                 A QGIS plugin
 Network based catchment analysis
                             -------------------
        begin                : 2016-05-19
        git sha              : $Format:%H$
        copyright            : (C) 2016 by Laurens Versluis
        email                : l.versluis@spacesyntax.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

import os

from PyQt4 import QtGui, uic

FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'catchment_analyser_dialog_base.ui'))


class CatchmentAnalyserDialog(QtGui.QDialog, FORM_CLASS):
    def __init__(self, parent=None):
        """Constructor."""
        super(CatchmentAnalyserDialog, self).__init__(parent)
        # Set up the user interface from Designer.
        # After setupUI you can access any designer object by doing
        # self.<objectname>, and you can use autoconnect slots - see
        # http://qt-project.org/doc/qt-4.8/designer-using-a-ui-file.html
        # #widgets-and-dialogs-with-auto-connect
        self.setupUi(self)


        # Output internal GUI signals
        self.networkText.setPlaceholderText("Save as temporary layer...")
        self.networkSaveButton.clicked.connect(self.setNetworkOutput)
        self.polygonText.setPlaceholderText("Save as temporary layer...")
        self.polygonSaveButton.clicked.connect(self.setPolygonOutput)

        # setup the progress bar
        self.analysisProgress.setMinimum(0)
        self.analysisProgress.setMaximum(5)

    def setNetworkLayers(self, names):
        layers = ['-----']
        if names:
            layers = []
            layers.extend(names)
        self.networkCombo.clear()
        self.networkCombo.addItems(layers)


    def getNetwork(self):
        return self.networkCombo.currentText()


    def setCostFields(self, names):
        if self.costCheck.isChecked():
            self.costCombo.setEnabled(True)
            fields = ['-----']
            if names:
                fields = []
                fields.extend(names)
            self.costCombo.clear()
            self.costCombo.addItems(names)


    def getCostField(self):
        if self.costCheck.isChecked():
            cost_field = self.costCombo.currentText()
        else:
            cost_field = None
        return cost_field


    def setOriginLayers(self, names):
        layers = ['-----']
        if names:
            layers = []
            layers.extend(names)
        self.originsCombo.clear()
        self.originsCombo.addItems(layers)


    def getOrigins(self):
        return self.originsCombo.currentText()


    def setNameFields(self, names):
        if self.nameCheck.isChecked():
            self.nameCombo.setEnabled(True)
            fields = ['-----']
            if names:
                fields.extend(names)
            self.nameCombo.clear()
            self.nameCombo.addItems(names)


    def getName(self):
        return self.nameCombo.currentText()


    def getDistances(self):
        distances_text = self.distancesText.text().split(',')
        distances_integer = [int(i) for i in distances_text]
        return distances_integer


    def getNetworkTolerance(self):
        return self.networkTolSpin.value()


    def getPolygonTolerance(self):
        return self.polygonTolSpin.value()


    def setNetworkOutput(self):
        file_name = QtGui.QFileDialog.getSaveFileName(self, "Save output file ", "catchment_network", '*.shp')
        if file_name:
            self.networkText.setText(file_name)


    def getNetworkOutput(self):
        return self.networkText.text()


    def setPolygonOutput(self):
        file_name = QtGui.QFileDialog.getSaveFileName(self, "Save output file ", "catchment_polygon", '*.shp')
        if file_name:
            self.polygonText.setText(file_name)


    def getPolygonOutput(self):
        return self.polygonText.text()






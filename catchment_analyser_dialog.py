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

        # Network GUI signals
        self.choose_network.activated.connect(self.choose_network) # no repeating names
        self.check_cost.stateChanged.connect(self.browse_cost_fields)
        self.choose_cost.activated.connect(self.choose_cost)

        # Origin GUI signals
        self.choose_origins.activated.connect(self.choose_origins)
        self.check_name.stateChanged.connect(self.browse_name_fields)
        self.choose_name.activated.connect(self.choose_name_field)

        # Output GUI signals
        self.path_output_network.setPlaceholderText("Save as temporary layer...")
        self.browse_output_network.clicked.connect(self.browse_output_network)
        self.path_output_polygon.setPlaceholderText("Save as temporary layer...")
        self.browse_output_polygon.clicked.connect(self.browse_output_polygon)

        # connect the run button
        self.run_mca.clicked.connect(self.run_analysis)

        # setup the progress bar
        self.progress_mca.setMinimum(0)
        self.progress_mca.setMaximum(5)

    def choose_network(self):
        pass

    def browse_cost_fields(self):
        pass

    def choose_cost_fields(self):
        pass

    def choose_origin(self):
        pass

    def browse_name_fields(self):
        pass

    def choose_name_field(self):
        pass

    def browse_output_network(self):
        pass

    def browse_output_polygon(self):
        pass






    def updateLayers(self):
        pass

    def getLayer(self):
        pass




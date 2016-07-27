# -*- coding: utf-8 -*-
"""
/***************************************************************************
 CatchmentAnalyser
                             Catchment Analyser
 Network based catchment analysis
                              -------------------
        begin                : 2016-05-19
        author               : Laurens Versluis
        copyright            : (C) 2016 by Space Syntax Limited
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
from PyQt4.QtCore import *
from PyQt4.QtGui import *

# Initialize Qt resources from file resources.py
import resources

# Import the code for the dialog
from catchment_analyser_dialog import CatchmentAnalyserDialog
import os.path

# Import QGIS classes
from qgis.core import *
from qgis.gui import *
from qgis.utils import *

# Import tool classes
import catchment_tools

# Import utility tools
import utility_functions as uf

# Import the debug library
# set is_debug to False in release version
is_debug = False
try:
    import pydevd
    has_pydevd = False
except ImportError, e:
    has_pydevd = False
    is_debug = False

class CatchmentAnalyser:
    """QGIS Plugin Implementation."""

    def __init__(self, iface):
        """Constructor.

        :param iface: An interface instance that will be passed to this class
            which provides the hook by which you can manipulate the QGIS
            application at run time.
        :type iface: QgsInterface
        """
        # Save reference to the QGIS interface
        self.iface = iface
        # Initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)
        # # Initialize analysis
        # self.catchmentAnalysis = catchment_tools.catchmentAnalysis(self.iface)
        # initialize locale
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            'CatchmentAnalyser_{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)

            if qVersion() > '4.3.3':
                QCoreApplication.installTranslator(self.translator)

        # Create the dialog (after translation) and keep reference
        self.dlg = CatchmentAnalyserDialog()
        # Declare instance attributes
        self.actions = []
        self.menu = self.tr(u'&Space Syntax Toolkit')
        # TODO: We are going to let the user set this up in a future iteration
        self.toolbar = self.iface.addToolBar(u'CatchmentAnalyser')
        self.toolbar.setObjectName(u'CatchmentAnalyser')

        # Setup debugger
        if has_pydevd and is_debug:
            pydevd.settrace('localhost', port=53100, stdoutToServer=True, stderrToServer=True, suspend=True)

        # Setup GUI signals
        self.dlg.networkCombo.activated.connect(self.updateCost)
        self.dlg.originsCombo.activated.connect(self.updateName)
        self.dlg.costCheck.stateChanged.connect(self.updateCost)
        self.dlg.nameCheck.stateChanged.connect(self.updateName)
        self.dlg.analysisButton.clicked.connect(self.runAnalysis)



    # noinspection PyMethodMayBeStatic
    def tr(self, message):
        """Get the translation for a string using Qt translation API.

        We implement this ourselves since we do not inherit QObject.

        :param message: String for translation.
        :type message: str, QString

        :returns: Translated version of message.
        :rtype: QString
        """
        # noinspection PyTypeChecker,PyArgumentList,PyCallByClass
        return QCoreApplication.translate('CatchmentAnalyser', message)

    def add_action(
        self,
        icon_path,
        text,
        callback,
        enabled_flag=True,
        add_to_menu=True,
        add_to_toolbar=True,
        status_tip=None,
        whats_this=None,
        parent=None):
        """Add a toolbar icon to the toolbar.

        :param icon_path: Path to the icon for this action. Can be a resource
            path (e.g. ':/plugins/foo/bar.png') or a normal file system path.
        :type icon_path: str

        :param text: Text that should be shown in menu items for this action.
        :type text: str

        :param callback: Function to be called when the action is triggered.
        :type callback: function

        :param enabled_flag: A flag indicating if the action should be enabled
            by default. Defaults to True.
        :type enabled_flag: bool

        :param add_to_menu: Flag indicating whether the action should also
            be added to the menu. Defaults to True.
        :type add_to_menu: bool

        :param add_to_toolbar: Flag indicating whether the action should also
            be added to the toolbar. Defaults to True.
        :type add_to_toolbar: bool

        :param status_tip: Optional text to show in a popup when mouse pointer
            hovers over the action.
        :type status_tip: str

        :param parent: Parent widget for the new action. Defaults None.
        :type parent: QWidget

        :param whats_this: Optional text to show in the status bar when the
            mouse pointer hovers over the action.

        :returns: The action that was created. Note that the action is also
            added to self.actions list.
        :rtype: QAction
        """

        icon = QIcon(icon_path)
        action = QAction(icon, text, parent)
        action.triggered.connect(callback)
        action.setEnabled(enabled_flag)

        if status_tip is not None:
            action.setStatusTip(status_tip)

        if whats_this is not None:
            action.setWhatsThis(whats_this)

        if add_to_toolbar:
            self.toolbar.addAction(action)

        if add_to_menu:
            self.iface.addPluginToVectorMenu(
                self.menu,
                action)

        self.actions.append(action)

        return action

    def initGui(self):
        """Create the menu entries and toolbar icons inside the QGIS GUI."""

        icon_path = ':/plugins/CatchmentAnalyser/icon.png'
        self.add_action(
            icon_path,
            text=self.tr(u'Catchment Analyser'),
            callback=self.run,
            parent=self.iface.mainWindow())


    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""
        for action in self.actions:
            self.iface.removePluginVectorMenu(
                self.menu,
                action)
            self.iface.removeToolBarIcon(action)
        # remove the toolbar
        del self.toolbar


    def updateLayers(self):
        self.updateNetwork()
        self.updateOrigins()


    def updateNetwork(self):
        network_layers = uf.getLegendLayersNames(iface, geom=[1, ], provider='all')
        self.dlg.setNetworkLayers(network_layers)
        self.updateCost()


    def updateOrigins(self):
        origins_layers = uf.getLegendLayersNames(iface, geom=[0, ], provider='all')
        self.dlg.setOriginLayers(origins_layers)
        self.updateName()


    def updateCost(self):
        if self.dlg.costCheck.isChecked():
            network = self.getNetwork()
            self.dlg.setCostFields(uf.getNumericFieldNames(network))
        else:
            self.dlg.costCombo.clear()
            self.dlg.costCombo.setEnabled(False)


    def updateName(self):
        if self.dlg.nameCheck.isChecked():
            origins = self.getOrigins()
            self.dlg.setNameFields(uf.getFieldNames(origins))
        else:
            self.dlg.nameCombo.clear()
            self.dlg.nameCombo.setEnabled(False)


    def getNetwork(self):
        return uf.getLegendLayerByName(iface, self.dlg.getNetwork())


    def getOrigins(self):
        return uf.getLegendLayerByName(iface, self.dlg.getOrigins())


    def tempNetwork(self, epsg):
        if self.dlg.networkCheck.isChecked():
            output_network = uf.createTempLayer(
                'catchment_network',
                'LINESTRING',
                epsg,
                ['id',],
                [QVariant.Int,]
            )
            return output_network


    def tempPolygon(self, epsg):
        if self.dlg.polygonCheck.isChecked():
            output_polygon = uf.createTempLayer(
                'catchment_areas',
                'POLYGON',
                epsg,
                ['id', 'origin', 'distance'],
                [QVariant.Int, QVariant.Int, QVariant.Int]
            )
            return output_polygon

    def giveWarningMessage(self, message):
        # Gives warning according to message
        self.iface.messageBar().pushMessage(
            "Catchment Analyser: ",
            "%s" % (message),
            level=QgsMessageBar.WARNING,
            duration=5)


    def getAnalysisSettings(self):

        # Creating a combined settings dictionary
        settings = {}

        # Raise warnings
        if not self.getNetwork():
            self.giveWarningMessage("No network selected!")
        elif self.getNetwork().crs().geographicFlag() or self.getOrigins().crs().geographicFlag():
            self.giveWarningMessage("No layers with projection CRS!")
        elif not self.getOrigins():
            self.giveWarningMessage("Catchment Analyser: No origins selected!")
        elif not self.dlg.getDistances():
            self.giveWarningMessage("No distances defined!")
        else:
            try:
                distances = [int(i) for i in self.dlg.getDistances()]
            except ValueError:
                self.giveWarningMessage("No numerical distances!")
                return

            # Get settings from the dialog
            settings['network'] = self.getNetwork()
            settings['cost'] = self.dlg.getCostField()
            settings['origins'] = self.getOrigins()
            settings['name'] = self.dlg.getName()
            settings['distances'] = distances
            settings['network tolerance'] = self.dlg.getNetworkTolerance()
            settings['polygon tolerance'] = int(self.dlg.getPolygonTolerance())
            settings['crs'] = self.getNetwork().crs()
            settings['epsg'] = self.getNetwork().crs().authid()[5:]  # removing EPSG:
            settings['temp network'] = self.tempNetwork(settings['epsg'])
            settings['temp polygon'] = self.tempPolygon(settings['epsg'])
            settings['output network check'] = self.dlg.networkCheck.isChecked()
            settings['output network'] = self.dlg.getNetworkOutput()
            settings['output polygon check'] = self.dlg.polygonCheck.isChecked()
            settings['output polygon'] = self.dlg.getPolygonOutput()

            return settings


    def runAnalysis(self):
        self.dlg.analysisProgress.reset()
        # Create an analysis instance
        settings = self.getAnalysisSettings()
        analysis = catchment_tools.catchmentAnalysis(self.iface, settings)

        # Create new thread and move the analysis class to it
        analysis_thread = QThread()
        analysis.moveToThread(analysis_thread)

        # Setup signals
        analysis.finished.connect(self.analysisFinish)
        analysis.error.connect(self.analysisError)
        analysis.warning.connect(self.giveWarningMessage)
        analysis.progress.connect(self.dlg.analysisProgress.setValue)
        self.dlg.cancelButton.clicked.connect(analysis.kill_analysis)
        analysis.kill.connect(self.killAnalysis)

        # Start analysis
        analysis_thread.started.connect(analysis.analysis)
        analysis_thread.start()
        self.analysis_thread = analysis_thread
        self.analysis = analysis


    def analysisFinish(self, output):

        # Render output
        if output:
            output_network = output['output network']
            output_polygon = output['output polygon']
            distances = output['distances']
            if output_network:
                self.renderNetwork(output_network, distances)
            if output_polygon:
                self.renderPolygon(output_polygon)

        # Clean up thread and analysis
        self.analysis_thread.quit()
        self.analysis_thread.wait()
        self.analysis_thread.deleteLater()

        # Closing the dialog
        self.dlg.closeDialog()

    def renderNetwork(self, output_network, distances):

        # Settings
        catchment_threshold = int(max(distances))

        # settings for 10 color ranges depending on the radius
        color_ranges = (
            (0, (0.1 * catchment_threshold), '#ff0000'),
            ((0.1 * catchment_threshold), (0.2 * catchment_threshold), '#ff5100'),
            ((0.2 * catchment_threshold), (0.3 * catchment_threshold), '#ff9900'),
            ((0.3 * catchment_threshold), (0.4 * catchment_threshold), '#ffc800'),
            ((0.4 * catchment_threshold), (0.5 * catchment_threshold), '#ffee00'),
            ((0.5 * catchment_threshold), (0.6 * catchment_threshold), '#a2ff00'),
            ((0.6 * catchment_threshold), (0.7 * catchment_threshold), '#00ff91'),
            ((0.7 * catchment_threshold), (0.8 * catchment_threshold), '#00f3ff'),
            ((0.8 * catchment_threshold), (0.9 * catchment_threshold), '#0099ff'),
            ((0.9 * catchment_threshold), (1 * catchment_threshold), '#0033ff'))

        # list with all color ranges
        ranges = []

        # for each range create a symbol with its respective color
        for lower, upper, color in color_ranges:
            symbol = QgsSymbolV2.defaultSymbol(output_network.geometryType())
            symbol.setColor(QColor(color))
            symbol.setWidth(0.5)
            range = QgsRendererRangeV2(lower, upper, symbol, '')
            ranges.append(range)

        # create renderer based on ranges and apply to network
        renderer = QgsGraduatedSymbolRendererV2('min_dist', ranges)
        output_network.setRendererV2(renderer)

        # add network to the canvas
        QgsMapLayerRegistry.instance().addMapLayer(output_network)


    def renderPolygon(self, output_polygon):

        # create a black dotted outline symbol layer
        symbol = QgsFillSymbolV2().createSimple({'color': 'grey', 'outline_width': '0'})
        symbol.setAlpha(0.2)

        # create renderer and change the symbol layer in its symbol
        output_polygon.rendererV2().setSymbol(symbol)

        # add catchment to the canvas
        QgsMapLayerRegistry.instance().addMapLayer(output_polygon)


    def analysisError(self, e, exception_string):
        QgsMessageLog.logMessage(
            'Catchment Analyser raised an exception: %s' % exception_string,
            level=QgsMessageLog.CRITICAL)

        # Closing the dialog
        self.dlg.closeDialog()


    def killAnalysis(self):

        self.analysis_thread.quit()
        self.analysis_thread.wait()
        self.analysis_thread.deleteLater()
        self.analysis.deleteLater()

        # Closing the dialog
        self.dlg.closeDialog()


    def run(self):
        # Show the dialog
        self.dlg.show()

        # Update layers
        self.updateLayers()



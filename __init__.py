# -*- coding: utf-8 -*-
"""
/***************************************************************************
 CatchmentAnalyser
                                 A QGIS plugin
 Network based catchment analysis
                             -------------------
        begin                : 2016-05-19
        copyright            : (C) 2016 by Laurens Versluis
        email                : l.versluis@spacesyntax.com
        git sha              : $Format:%H$
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 This script initializes the plugin, making it known to QGIS.
"""


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load CatchmentAnalyser class from file CatchmentAnalyser.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """
    #
    from .catchment_analyser import CatchmentAnalyser
    return CatchmentAnalyser(iface)

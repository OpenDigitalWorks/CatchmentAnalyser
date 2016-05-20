from PyQt4.QtCore import *
from PyQt4.QtGui import *

from qgis.core import *
from qgis.gui import *
from qgis.networkanalysis import *
from qgis.utils import *

import math
from shapely.ops import cascaded_union, polygonize
from shapely.geometry import MultiPoint
from scipy.spatial import Delaunay



class customCost(QgsArcProperter):
    def __init__(self, costColumIndex, defaultValue):
        QgsArcProperter.__init__(self)
        self.cost_column_index = costColumIndex
        self.default_value = defaultValue

    def property(self, distance, Feature):
        cost = float(Feature.attributes()[self.cost_column_index])
        if cost <= 0.0:
            return QVariant(self.default_value)
        return cost

    def requiredAttributes(self):
        l = []
        l.append(self.cost_column_index);
        return l

class caTools():

    def __init__(self):
        self.network = None
        self.origins = None
        self.cost = None
        self.origin_name = None

    def graph_builder(self,network,cost_field, origin):
        # Settings
        crs = network.crs()
        epsg = crs.authid()
        otf = False
        network_fields = network_lines.pendingFields()
        custom_cost_index = network_fields.indexFromName(cost_field)
        # Setting up graph build director
        director = QgsLineVectorLayerDirector(network_lines, -1, '', '', '', 3)
        # Determining cost calculation
        if custom_cost == True:
            properter = customCost(custom_cost_index,0)
        else:
            properter = QgsDistanceArcProperter()
        # Building graph
        director.addProperter(properter)
        builder = QgsGraphBuilder(crs, otf, tolerance, epsg)
        # Reading origins and making list of coordinates
        origins = []
        origins_name = {}
        for i,f in enumerate(origin_points.getFeatures()):
            geom = f.geometry().asPoint()
            origins.append(geom)
            if origins_column:
                origins_name[i] = f[origins_column]
        # Connect origin points to the director and build graph
        tied_origins = director.makeGraph(builder, origins)
        graph = builder.graph()
        print origins_name
        return graph, tied_origins, origins_name

    def ca_polygon_builder(sdd):
        pass

    def ca_network_writer(self):
        pass

    def ca_polygon_writer(self):
        pass

    def ca_network_renderer(self):
        pass

    def ca_polygon_renderer(self):
        pass


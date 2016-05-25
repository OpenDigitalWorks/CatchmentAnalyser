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

class catchment_tools():

    def __init__(self):
        self.network = None
        self.origins = None
        self.cost = None
        self.origin_name = None

    def network_preparation(self, network_vector, topology_bool, stub_ratio, unlink_vector,):

        # Settings
        unlink_buffer = 5

        # Variables
        network = [] # Final list of network lines
        segment_index = QgsSpatialIndex() # Index of segment bounding boxes
        segment_dict = {} # Dictionary of segments indices and geometries
        unlinked_segments_ids = []
        origin_type = None # Origins can be point or polygon

        # Check network layer validity
        if not network_vector.isValid():
            self.warning_message("Invalid network layer!"),

        # Check if network layer contains lines
        elif not (network_vector.wkbType() == 2 or network_vector.wkbType() == 5):
            self.warning_message("Network layer contains no lines!")

        else :

            # Check origin layer validity
            if not unlink_vector.isValid():
                self.warning_message("Invalid origin layer!")

            # Check if origin layer contains lines
            elif not (unlink_vector.wkbType() == 1 or
            unlink_vector.wkbType() == 3 or
            unlink_vector.wkbType() == 4 or
            unlink_vector.wkbType() == 6):
                self.warning_message("Unlink layer contains no points or polygons!")

            else :

                # Check unlink geometry type
                if (unlink_vector.wkbType() == 1 or unlink_vector.wkbType() == 4):
                    origin_type = 'point'
                elif (unlink_vector.wkbType() == 1 or unlink_vector.wkbType() == 4):
                    origin_type = 'polygon'

                # If network is not topological start segmentation
                if topology_bool == True:

                    # Insert segments of network to the spatial index and dictionary
                    for segment in network_vector.getFeatures():
                        segment_index.insertFeature(segment)
                        segment_dict[segment.id()] = segment.geometry().asWkb()

                    # Loop through unlinks and list unlinked segments
                    for unlink in unlink_vector.getFeatures():

                        # Create unlink area when unlinks are points
                        if origin_type = 'point':
                            unlink_area = unlink.geometry().buffer(unlink_buffer,5)

                        # Create unlink area when unlinks are polygons
                        else:
                            unlink_area = unlink.geometry().boundingBox()

                        # Create list of id's of intersecting segments
                        nearest_segments = segment_index.intersects(unlink_area)

                        # Check number of intersecting segments
                        if nearest_segments > 2:
                            self.warning_message("Unlink layer references to many segments!")

                        # Add unlinked segments to the list
                        else:
                             for seg_id in nearest_segments:
                                 unlinked_segments_ids.append(seg_id)

                    # Insert segments of network to the spatial index
                    for segment in network_vector.getFeatures():

                        segment_id = segment.id()
                        segment_geom = segment.geometry()

                        # Add unlinked segments to the network
                        if segment_id in unlinked_segments_ids:
                            network.append(segment_geom)

                        # Split the remaining segments
                        else:

                            # Calculate segment length
                            segment_length = segment_geom.length()

                            # Identify intersecting segments
                            intersecting_segments = []
                            intersecting_segments_ids = segment_index.intersects(segment_geom)
                            inte

                            # Loop for segment parts excluding itself
                            for id in [i for i in nearest_segments if i != segment_id]:

                                # Loop through segment parts


                                # Loop for segment parts
                                if
                            intersecting_segments = QgsGeometry()
                                # Only append non-stubs based on stub ratio

                # Otherwise add all segments of the network layer
                else:
                    network.append(segment.geometry() for segment in network_vector.getFeatures())

        return network

    def origin_preparation(self, origin_layer):

        # Variables

		# Check origins validity
		
		# Check if origins are polygons
		
			# Get points from polygon

		# Combining points with their name

			# Use designated field
			
			# Use number
		
		return origin_points
		
		
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

        return graph, tied_origins

    def ca_analysis(self, graph, tied_origins):

        ca_network = []



    def ca_polygon_builder(self):
        pass

    def ca_network_writer(self):
        pass

    def ca_polygon_writer(self):
        pass

    def ca_network_renderer(self):
        pass

    def ca_polygon_renderer(self):
        pass

    def warning_message(self,message):

        # Gives warning according to message
        self.iface.messageBar().pushMessage(
            "Catchment Analyser: ",
            "%s" % (message),
            level=QgsMessageBar.WARNING,
            duration=5)
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from qgis.core import *
from qgis.gui import *
from qgis.networkanalysis import *
from qgis.utils import *

import math

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

class concaveHull():
    def clean_list(list_of_points):
        """
        Deletes duplicate points in list_of_points
        """
        return list(set(list_of_points))

    def length(vector):
        """
        Returns the number of elements in vector
        """
        return len(vector)

    def find_min_y_point(list_of_points):
        """
        Returns that point of *list_of_points* having minimal y-coordinate
        :param list_of_points: list of tuples
        :return: tuple (x, y)
        """
        min_y_pt = list_of_points[0]
        for point in list_of_points[1:]:
            if point[1] < min_y_pt[1] or (point[1] == min_y_pt[1] and point[0] < min_y_pt[0]):
                min_y_pt = point
        return min_y_pt

    def add_point(vector, element):
        """
        Returns vector with the given element append to the right
        """
        vector.append(element)
        return vector

    def remove_point(vector, element):
        """
        Returns a copy of vector without the given element
        """
        vector.pop(vector.index(element))
        return vector

    def euclidian_distance(point1, point2):
        """
        Returns the euclidian distance of the 2 given points.
        :param point1: tuple (x, y)
        :param point2: tuple (x, y)
        :return: float
        """
        return math.sqrt(math.pow(point1[0] - point2[0], 2) + math.pow(point1[1] - point2[1], 2))

    def nearest_points(list_of_points, point, k):
        # build a list of tuples of distances between point *point* and every point in *list_of_points*, and
        # their respective index of list *list_of_distances*
        list_of_distances = []
        for index in range(len(list_of_points)):
            list_of_distances.append((euclidian_distance(list_of_points[index], point), index))

        # sort distances in ascending order
        list_of_distances.sort()

        # get the k nearest neighbors of point
        nearest_list = []
        for index in range(min(k, len(list_of_points))):
            nearest_list.append((list_of_points[list_of_distances[index][1]]))
        return nearest_list

    def angle(from_point, to_point):
        """
        Returns the angle of the directed line segment, going from *from_point* to *to_point*, in radians. The angle is
        positive for segments with upward direction (north), otherwise negative (south). Values ranges from 0 at the
        right (east) to pi at the left side (west).
        :param from_point: tuple (x, y)
        :param to_point: tuple (x, y)
        :return: float
        """
        return math.atan2(to_point[1] - from_point[1], to_point[0] - from_point[0])

    def angle_difference(angle1, angle2):
        """
        Calculates the difference between the given angles in clockwise direction as radians.
        :param angle1: float
        :param angle2: float
        :return: float; between 0 and 2*Pi
        """
        if (angle1 > 0 and angle2 >= 0) and angle1 > angle2:
            return abs(angle1 - angle2)
        elif (angle1 >= 0 and angle2 > 0) and angle1 < angle2:
            return 2 * math.pi + angle1 - angle2
        elif (angle1 < 0 and angle2 <= 0) and angle1 < angle2:
            return 2 * math.pi + angle1 + abs(angle2)
        elif (angle1 <= 0 and angle2 < 0) and angle1 > angle2:
            return abs(angle1 - angle2)
        elif angle1 <= 0 < angle2:
            return 2 * math.pi + angle1 - angle2
        elif angle1 >= 0 >= angle2:
            return angle1 + abs(angle2)
        else:
            return 0

    def intersect(line1, line2):
        """
        Returns True if the two given line segments intersect each other, and False otherwise.
        :param line1: 2-tuple of tuple (x, y)
        :param line2: 2-tuple of tuple (x, y)
        :return: boolean
        """
        a1 = line1[1][1] - line1[0][1]
        b1 = line1[0][0] - line1[1][0]
        c1 = a1 * line1[0][0] + b1 * line1[0][1]
        a2 = line2[1][1] - line2[0][1]
        b2 = line2[0][0] - line2[1][0]
        c2 = a2 * line2[0][0] + b2 * line2[0][1]
        tmp = (a1 * b2 - a2 * b1)
        if tmp == 0:
            return False
        sx = (c1 * b2 - c2 * b1) / tmp
        if (sx > line1[0][0] and sx > line1[1][0]) or (sx > line2[0][0] and sx > line2[1][0]) or \
                (sx < line1[0][0] and sx < line1[1][0]) or (sx < line2[0][0] and sx < line2[1][0]):
            return False
        sy = (a1 * c2 - a2 * c1) / tmp
        if (sy > line1[0][1] and sy > line1[1][1]) or (sy > line2[0][1] and sy > line2[1][1]) or \
                (sy < line1[0][1] and sy < line1[1][1]) or (sy < line2[0][1] and sy < line2[1][1]):
            return False
        return True

    def point_in_polygon_q(point, list_of_points):
        """
        Return True if given point *point* is laying in the polygon described by the vertices *list_of_points*,
        otherwise False
        Based on the "Ray Casting Method" described by Joel Lawhead in this blog article:
        http://geospatialpython.com/2011/01/point-in-polygon.html
        """
        x = point[0]
        y = point[1]
        poly = [(pt[0], pt[1]) for pt in list_of_points]
        n = len(poly)
        inside = False

        p1x, p1y = poly[0]
        for i in range(n + 1):
            p2x, p2y = poly[i % n]
            if y > min(p1y, p2y):
                if y <= max(p1y, p2y):
                    if x <= max(p1x, p2x):
                        if p1y != p2y:
                            xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                        if p1x == p2x or x <= xints:
                            inside = not inside
            p1x, p1y = p2x, p2y

        return inside

    def write_wkt(point_list, file_name):
        """
        Writes the geometry described by *point_list* in Well Known Text format to file
        :param point_list: list of tuples (x, y)
        :return: None
        """
        if file_name is None:
            file_name = 'hull2.wkt'
        if os.path.isfile(file_name):
            outfile = open(file_name, 'a')
        else:
            outfile = open(file_name, 'w')
            outfile.write('%s\n' % 'WKT')
        wkt = 'POLYGON((' + str(point_list[0][0]) + ' ' + str(point_list[0][1])
        for p in point_list[1:]:
            wkt += ', ' + str(p[0]) + ' ' + str(p[1])
        wkt += '))'
        outfile.write('%s\n' % wkt)
        outfile.close()
        return None

    def as_wkt(point_list):
        """
        Returns the geometry described by *point_list* in Well Known Text format
        Example: hull = self.as_wkt(the_hull)
                 feature.setGeometry(QgsGeometry.fromWkt(hull))
        :param point_list: list of tuples (x, y)
        :return: polygon geometry as WTK
        """
        wkt = 'POLYGON((' + str(point_list[0][0]) + ' ' + str(point_list[0][1])
        for p in point_list[1:]:
            wkt += ', ' + str(p[0]) + ' ' + str(p[1])
        wkt += '))'
        return wkt

    def as_polygon(point_list):
        """
        Returns the geometry described by *point_list* in as QgsGeometry
        :param point_list: list of tuples (x, y)
        :return: QgsGeometry
        """
        # create a list of QgsPoint() from list of point coordinate strings in *point_list*
        points = [QgsPoint(point[0], point[1]) for point in point_list]
        # create the polygon geometry from list of point geometries
        poly = QgsGeometry.fromPolygon([points])
        return poly

    def enable_use_of_global_CRS():
        """
        Set new layers to use the project CRS.
        Code snipped taken from http://pyqgis.blogspot.co.nz/2012/10/basics-automatic-use-of-crs-for-new.html
        Example: old_behaviour = enable_use_of_global_CRS()
        :return: string
        """
        settings = QSettings()
        old_behaviour = settings.value('/Projections/defaultBehaviour')
        settings.setValue('/Projections/defaultBehaviour', 'useProject')
        return old_behaviour

    def disable_use_of_global_CRS(default_behaviour='prompt'):
        """
        Enables old settings again. If argument is missing then set behaviour to prompt.
        Example: disable_use_of_global_CRS(old_behaviour)
        :param default_behaviour:
        :return: None
        """
        settings = QSettings()
        settings.setValue('/Projections/defaultBehaviour', default_behaviour)
        return None

    def extract_points(geom):
        """
        Generate list of QgsPoints from QgsGeometry *geom* ( can be point, line, or polygon )
        Code taken from fTools plugin
        :param geom: an arbitrary geometry feature
        :return: list of points
        """
        multi_geom = QgsGeometry()
        temp_geom = []
        # point geometry
        if geom.type() == 0:
            if geom.isMultipart():
                temp_geom = geom.asMultiPoint()
            else:
                temp_geom.append(geom.asPoint())
        # line geometry
        if geom.type() == 1:
            # if multipart feature explode to single part
            if geom.isMultipart():
                multi_geom = geom.asMultiPolyline()
                for i in multi_geom:
                    temp_geom.extend(i)
            else:
                temp_geom = geom.asPolyline()
        # polygon geometry
        elif geom.type() == 2:
            # if multipart feature explode to single part
            if geom.isMultipart():
                multi_geom = geom.asMultiPolygon()
                # now single part polygons
                for i in multi_geom:
                    # explode to line segments
                    for j in i:
                        temp_geom.extend(j)
            else:
                multi_geom = geom.asPolygon()
                # explode to line segments
                for i in multi_geom:
                    temp_geom.extend(i)
        return temp_geom

    def sort_by_angle(list_of_points, last_point, last_angle):
        def getkey(item):
            return angle_difference(last_angle, angle(last_point, item))

        vertex_list = sorted(list_of_points, key=getkey, reverse=True)
        return vertex_list

    def concave_hull(points_list, k):
        """
        Calculates a valid concave hull polygon containing all given points. The algorithm searches for that
        point in the neighborhood of k nearest neighbors which maximizes the rotation angle in clockwise direction
        without intersecting any previous line segments.
        This is an implementation of the algorithm described by Adriano Moreira and Maribel Yasmina Santos:
        CONCAVE HULL: A K-NEAREST NEIGHBOURS APPROACH FOR THE COMPUTATION OF THE REGION OCCUPIED BY A SET OF POINTS.
        GRAPP 2007 - International Conference on Computer Graphics Theory and Applications; pp 61-68.
        :param points_list: list of tuples (x, y)
        :param k: integer
        :return: list of tuples (x, y)
        """
        # return an empty list if not enough points are given
        if k > len(points_list):
            return None

        # the number of nearest neighbors k must be greater than or equal to 3
        # kk = max(k, 3)
        kk = max(k, 2)

        # delete duplicate points
        point_set = clean_list(points_list)

        # if point_set has less then 3 points no polygon can be created and an empty list will be returned
        if len(point_set) < 3:
            return None

        # if point_set has 3 points then these are already vertices of the hull. Append the first point to
        # close the hull polygon
        if len(point_set) == 3:
            return add_point(point_set, point_set[0])

        # make sure that k neighbours can be found
        kk = min(kk, len(point_set))

        # start with the point having the smallest y-coordinate (most southern point)
        first_point = find_min_y_point(point_set)

        # add this points as the first vertex of the hull
        hull = [first_point]

        # make the first vertex of the hull to the current point
        current_point = first_point

        # remove the point from the point_set, to prevent him being among the nearest points
        point_set = remove_point(point_set, first_point)
        previous_angle = math.pi

        # step counts the number of segments
        step = 2

        # as long as point_set is not empty or search is returning to the starting point
        while (current_point != first_point) or (step == 2) and (len(point_set) > 0):

            # after 3 iterations add the first point to point_set again, otherwise a hull cannot be closed
            if step == 5:
                point_set = add_point(point_set, first_point)

            # search the k nearest neighbors of the current point
            k_nearest_points = nearest_points(point_set, current_point, kk)

            # sort the candidates (neighbors) in descending order of right-hand turn. This way the algorithm progresses
            # in clockwise direction through as many points as possible
            c_points = sort_by_angle(k_nearest_points, current_point, previous_angle)

            its = True
            i = -1

            # search for the nearest point to which the connecting line does not intersect any existing segment
            while its is True and (i < len(c_points) - 1):
                i += 1
                if c_points[i] == first_point:
                    last_point = 1
                else:
                    last_point = 0
                j = 2
                its = False

                while its is False and (j < len(hull) - last_point):
                    its = intersect((hull[step - 2], c_points[i]), (hull[step - 2 - j], hull[step - 1 - j]))
                    j += 1

            # there is no candidate to which the connecting line does not intersect any existing segment, so the
            # for the next candidate fails. The algorithm starts again with an increased number of neighbors
            if its is True:
                return concave_hull(points_list, kk + 1)

            # the first point which complies with the requirements is added to the hull and gets the current point
            current_point = c_points[i]
            hull = add_point(hull, current_point)

            # calculate the angle between the last vertex and his precursor, that is the last segment of the hull
            # in reversed direction
            previous_angle = angle(hull[step - 1], hull[step - 2])

            # remove current_point from point_set
            point_set = remove_point(point_set, current_point)

            # increment counter
            step += 1

        all_inside = True
        i = len(point_set) - 1

        # check if all points are within the created polygon
        while (all_inside is True) and (i >= 0):
            all_inside = point_in_polygon_q(point_set[i], hull)
            i -= 1

        # since at least one point is out of the computed polygon, try again with a higher number of neighbors
        if all_inside is False:
            return concave_hull(points_list, kk + 1)

        # a valid hull has been constructed
        return hull

class catchmentTools(customCost,concaveHull):

    def __init__(self):
        # self.network = None
        # self.origins = None
        # self.cost = None
        # self.origin_name = None

        # Add at the end of __init__ function of main plugin class
        pass

    def network_preparation(self, network_vector, unlink_vector, topology_bool, stub_ratio):

        # Settings
        unlink_buffer = 5

        # Variables
        network = []  # Final list of network lines
        segment_index = QgsSpatialIndex()  # Index of segment bounding boxes
        segment_dict = {}  # Dictionary of segments indices and geometries
        unlink_index = QgsSpatialIndex()


        if not network_vector:
            print "No network layer selected!"

        else:

            # Check network layer validity
            if not network_vector.isValid():
                print "Invalid network layer!"

            # Check if network layer contains lines
            elif not (network_vector.wkbType() == 2 or network_vector.wkbType() == 5):
                print "Network layer contains no lines!"

        # Check unlink layer
        if unlink_vector:

            # Check origin layer validity
            if not unlink_vector.isValid():
                print "Invalid origin layer!"

            # Check if origin layer contains lines
            elif not (unlink_vector.wkbType() == 1 or
                              unlink_vector.wkbType() == 3 or
                              unlink_vector.wkbType() == 4 or
                              unlink_vector.wkbType() == 6):
                print "Unlink layer contains no points or polygons!"

            # Check unlink geometry type
            if unlink_vector.wkbType() == 1 or unlink_vector.wkbType() == 4:
                origin_type = 'point'
            elif unlink_vector.wkbType() == 3 or unlink_vector.wkbType() == 6:
                origin_type = 'polygon'

        # If network is not topological start segmentation
        if topology_bool == False:

            # Insert segments of network to the spatial index and dictionary
            for segment in network_vector.getFeatures():
                segment_index.insertFeature(segment)
                segment_dict[segment.id()] = segment.geometryAndOwnership()

            # Create index of unlinks
            if unlink_vector:
                for unlink in unlink_vector.getFeatures():

                    # Create unlink area when unlinks are points
                    if origin_type == 'point':

                        # Create unlink area 5m around the point
                        unlink_geom = unlink.geometry().buffer(unlink_buffer, 5)
                        unlink_area = QgsFeature()
                        unlink_area.setGeometry(unlink_geom)

                    # Create unlink area when unlinks are polygons or lines
                    else:
                        unlink_area = unlink

                    # Add unlink to index and to dictionary
                    unlink_index.insertFeature(unlink_area)

            # Break each segment based on intersecting lines and unlinks
            for segment in network_vector.getFeatures():

                # Get id and geometry and length from segment
                segment_id = segment.id()
                segment_geom = segment.geometry()
                segment_length = segment_geom.length()

                # Get points from original segment
                seg_start_point = segment_geom.asPolyline()[0]
                seg_end_point = segment_geom.asPolyline()[-1]

                # Create list of id's of intersecting segments
                nearest_segment_ids = segment_index.intersects(segment_geom.boundingBox())

                # List of break points for the new segments
                break_points = []

                # Identify intersecting segments
                intersecting_segments_ids = segment_index.intersects(segment_geom.boundingBox())

                # Loop for intersecting segments excluding itself
                for id in intersecting_segments_ids:

                    # Skip if segment is itself
                    if id == segment_id:
                        continue

                    # Break segment according to remaining intersecting segment
                    else:

                        # Get geometry of intersecting segment
                        int_seg_geom = segment_dict[id]

                        # Identify the construction point of the new segment
                        if segment_geom.crosses(int_seg_geom) or segment_geom.touches(int_seg_geom):

                            # Create point where lines cross
                            point_geom = segment_geom.intersection(int_seg_geom)

                            # Create polygon of inters
                            point_buffer_geom = point_geom.buffer(1, 1).boundingBox()

                            # Check if cross point is an unlink
                            if not unlink_index.intersects(point_buffer_geom):
                                # Break points of intersecting lines
                                break_points.append(point_geom.asPoint())

                # Sort break_points according to distance to start point
                break_points.sort(key=lambda x: QgsDistanceArea().measureLine(seg_start_point, x))

                # Create segments using break points
                for i in range(0, len(break_points) - 1):
                    # Set end points
                    start_geom = QgsPoint(break_points[i])
                    end_geom = QgsPoint(break_points[i + 1])

                    # Create line and add to network
                    network.append(QgsGeometry.fromPolyline([start_geom, end_geom]))

                # Check if first segment is a potential stub
                for point in break_points:

                    if point != seg_start_point:

                        # Calculate distance between point and start point
                        distance_nearest_break = QgsDistanceArea().measureLine(seg_start_point, break_points[0])

                        # Only add first segment if it is a dead end
                        if distance_nearest_break > (stub_ratio * segment_length):
                            network.append(QgsGeometry.fromPolyline([seg_start_point, break_points[0]]))

                    # Check if last segment is a potential stub
                    elif point != seg_end_point:

                        # Calculate distance between point and end point
                        distance_nearest_break = QgsDistanceArea().measureLine(seg_end_point, break_points[-1])

                        # Only add last segment if it is a dead end
                        if distance_nearest_break > (stub_ratio * segment_length):
                            network.append(QgsGeometry.fromPolyline([seg_end_point, break_points[-1]]))

        # If topological network add all segments of the network layer straight away
        else:

            # Loop through features and add them to network
            for segment in network_vector.getFeatures():
                # Get geometry from features
                geom = segment.geometryAndOwnership()

                # Append geometry to network list
                network.append(geom)

        return network

    def origin_preparation(self, origin_vector, origin_name_field):

        # Create a list of origin point dictionaries containing name and geometry
        origins = []

        # Check origin layer validity
        if not origin_vector.isValid():
            self.warning_message("Invalid origin layer!")

        else:

            # Check origin layer geometry
            if origin_vector.wkbType() == 7:
                self.warning_message("Invalid origin geometry!")

            # Loop through origin and get or create points
            for i,f in enumerate(origin_vector.getFeatures()):

                # Create origin dictionary
                origin = {}

                # If origin name field is given get name
                if origin_name_field:
                    origin_name = f[origin_name_field]

                # Otherwise use index as name
                else:
                    origin_name = i

                # Depending on type of origin create and append points
                if f.vectorType() == Qgis.Point:
                    origin[i] = {origin_name : f.geometry()}

                elif f.vectorType() == Qgis.Polygon or f.vectorType() == Qgis.Line:
                    origin[i] = {"name" : origin_name, "geom" : f.geometry().centroid()}

                # Append origin names and geometry to origin points list
                origins.append(origin)

            return origins

    def graph_builder(self,network, cost_field, origins, tolerance):

        # Settings
        network_crs = network.crs()
        network_epsg = network_crs.authid()
        otf = False
        network_fields = network.pendingFields()
        custom_cost_index = network_fields.indexFromName(cost_field)

        # Setting up graph build director
        director = QgsLineVectorLayerDirector(network, -1, '', '', '', 3)

        # Determining cost calculation
        if cost_field == True:
            properter = customCost(custom_cost_index,0)
        else:
            properter = QgsDistanceArcProperter()

        # Creating graph builder
        director.addProperter(properter)
        builder = QgsGraphBuilder(network_crs, otf, tolerance, network_epsg)

        # Reading origins and making list of coordinates
        graph_origin_points = []

        # Loop through the origin points and add graph vertex indices
        for index,origin in enumerate(origins):
            graph_origin_points.append(origin[index][name])

        # Get origin graph vertex index
        tied_origin_vertices = director.makeGraph(builder, graph_origin_points)

        # Build the graph
        graph = builder.graph()

        # Create dictionary of origin names and tied origins
        tied_origins = {}

        # Combine origin names and tied point vertices
        for index, tied_origin in enumerate(tied_origin_vertices):
            tied_origins[index] = {"name" : origin[index][name], "vertex" : tied_origins[index]}

        return graph, tied_origins

    def ca_arc_analysis(self, graph, tied_origins, radii):

        # Settings
        ca_threshold = max(radii)

        # Variables
        catchment_network = {}
        catchment_points = {}

        # Loop through graph edges, update line dictionary and append
        for i in graph.arcCount():

            # Get geometry of edge from graph
            arc_geom = graph.arc(i)

            # Update catchment network
            catchment_network[i] = {'geom': arc_geom}

        # Loop through tied origins
        for origin in tied_origins:

            origin_name = origin["name"]

            # Find origin vertex id
            origin_vertex_id = graph.findVertex(tied_origins[i][vertex])

            # Run dijkstra and get tree and cost
            (tree, cost) = QgsGraphAnalyzer.dijkstra(graph, origin_vertex_id, 0)

            # Loop through graph arcs
            for index in graph.ArcCount():

                # Define the arc end point
                arc_outer_vertex_id = graph.arc(index).inVertex()
                arc_outer_vertex_geom = graph.vertex(arc_outer_vertex_id).point()

                # If arc is the origin set zero cost
                if arc_outer_vertex_id == origin_vertex_id:
                    catchment_network[index]['cost'][origin_name].append(0)

                # If arc is within connected and within the maximum radius set cost accordingly
                elif cost[arc_outer_vertex_id] < ca_threshold and cost[arc_outer_vertex_id] != -1:
                    arc_cost = cost[arc_outer_vertex_id]

                    # Update cost in catchment network dictionary
                    catchment_network[index]['cost'][origin_name].append(arc_cost)

                    # Add catchment points for each given radius
                    for radius in radii:
                        catchment_points[origin_name]["radius"] = radius
                        catchment_points[origin_name]["geom"].append(arc_outer_vertex_geom)

        return catchment_network, catchment_points

    def ca_network_writer(self, origins, catchment_network, output_network):

        # Variables
        arc_length_list = []
        arc_cost_list = []

        # Setup output network id column
        output_network.dataProvider().addAttributes([QgsField("id", QVariant.Int)])

        # Setup all unique origin columns
        unique_origin_list = []
        for origin in origins:
            if not origin in unique_origin_list:
                output_network.dataProvider().addAttributes([QgsField("%s" % origin, QVariant.Int)])
                unique_origin_list.append(origin)

        # Setup minimum origin distance column
        output_network.dataProvider().addAttributes(QgsField["min_dis", QVariant.Int])

        # Loop through arcs in catchment network and write geometry and costs
        for arc in catchment_network:

            # Get arc geometry
            arc_geom = QgsGeometry.fromPolyline(arc['geom'])
            arc_length = arc_geom.length()

            # Ignore arc if already processed
            if arc_length in arc_length_list:
                pass

            # Ignore arc is not connected to origins or outside catchment
            elif not arc['cost']:
                pass

            else:

                # Create feature and write id and geom
                f = QgsFeature(output_network.pendingFields())
                f.setAttribute("id", arc)
                f.setGeometry(arc_geom)

                # Read the list of costs and write them to output network
                for cost_origin in arc['cost'].keys():

                    # Get cost
                    cost = arc['cost'][cost_origin]

                    # Current cost entry
                    current_cost = f[cost_origin]

                    # If no entry set cost
                    if not current_cost:
                        f.setAttribute("%s" % cost_origin, cost)

                    else:

                        # Replace current cost when less than cost
                        if current_cost < cost:
                            f.setAttribute("%s" % cost_origin, cost)

                    # Add to list of costs
                    arc_cost_list.append(cost)

                # Calculate and write distance to nearest origin
                if arc_cost_list > 0:
                    f.setAttribute('min_dist', int(min(cost_list)))

                # Write feature to output network layer
                output_network.dataProvider().addFeatures([f])

                # Add the length of arc to length list in order to ignore duplicates;d
                arc_length_list.append(arc_length)

        return output_network

    def ca_polygon_writer(self, catchment_points, radii, output_polygon, alpha):

        # Variables
        unique_origin_list = []
        polygon_points = {}

        # Setup output polygon table
        output_polygon.dataProvider().addAttributes([
            QgsField("id", QVariant.Int),
            QgsField("origin", QVariant.Int),
            QgsField("radius", QVariant.Int)])
        output_network.updateFields()

        # Loop through origins and create list of aggregate polygon points
        for name in catchment_points:

            # Check if origin is not duplicated
            if name not in unique_origin_list:

                # Append origin to unique origin list
                unique_origin_list.append(name)

                # Copy origin points for all radii
                polygon_points[name] = catchment_points[name]

            # Otherwise append point to existing entry
            else:

                # Loop through radii and append points
                for radius in radii:
                    points = catchment_points[origin][radius]
                    polygon_points[name][radius].append(points)

        # Create index for the id column
        index = 1

        # Loop through polygon points and create their concave hull
        for name,points in polygon_points:

            # Loop through radii
            for radius in points[radius]:

                # Create polygon feature
                p = QgsFeature()
                p.setAttribute("id", index)
                p.setAttribute("origin", name)
                p.setAttribute("radius", radius)
                polygon_geom = QgsGeometry.fromWkt((concave_hull(points, alpha)).wkt())
                p.setGeometry(polygon_geom)
                output_polygon.dataProvider().addFeatures([p])

                # Next feature index
                index += 1

        return output_polygon

    def ca_network_renderer(self, output_network, radii):

        # Settings
        ca_threshold = max(int(radii))

        # settings for 10 color ranges depending on the radius
        color_ranges = (
            (0, (0.1 * ca_threshold), '#ff0000'),
            ((0.1 * ca_threshold), (0.2 * ca_threshold), '#ff5100'),
            ((0.2 * ca_threshold), (0.3 * ca_threshold), '#ff9900'),
            ((0.3 * ca_threshold), (0.4 * ca_threshold), '#ffc800'),
            ((0.4 * ca_threshold), (0.5 * ca_threshold), '#ffee00'),
            ((0.5 * ca_threshold), (0.6 * ca_threshold), '#a2ff00'),
            ((0.6 * ca_threshold), (0.7 * ca_threshold), '#00ff91'),
            ((0.7 * ca_threshold), (0.8 * ca_threshold), '#00f3ff'),
            ((0.8 * ca_threshold), (0.9 * ca_threshold), '#0099ff'),
            ((0.9 * ca_threshold), (1 * ca_threshold), '#0033ff'))

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

    def ca_polygon_renderer(self, output_polygon):

        # create a black dotted outline symbol layer
        symbol_layer = QgsMarkerLineSymbolLayerV2()
        symbol_layer.setColor(QColor('black'))
        symbol_layer.setWidth(1)

        # create renderer and change the symbol layer in its symbol
        renderer = output_polygon.rendererV2()
        renderer.symbols()[0].changeSymbolLayer(0, symbol_layer)
        output_polygon.setRendererV2(renderer)

        # add catchment to the canvas
        QgsMapLayerRegistry.instance().addMapLayer(output_polygon)

    def warning_message(self,message):
        # Gives warning according to message
        self.iface.messageBar().pushMessage(
            "Catchment Analyser: ",
            "%s" % (message),
            level=QgsMessageBar.WARNING,
            duration=5)
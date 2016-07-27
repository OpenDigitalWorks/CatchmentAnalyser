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

from qgis.core import *
from qgis.gui import *
from qgis.networkanalysis import *
from qgis.utils import *

import math

import utility_functions as uf

class customCost(QgsArcProperter):
    def __init__(self, costColumIndex, defaultValue):
        QgsArcProperter.__init__(self)
        self.cost_column_index = costColumIndex
        self.default_value = defaultValue

    def property(self, distance, feature):
        if not feature.attributes()[self.cost_column_index]:
            return self.default_value
        else:
            cost = float(feature.attributes()[self.cost_column_index])
            if cost <= 0.0:
                return self.default_value
        return cost

    def requiredAttributes(self):
        l = []
        l.append(self.cost_column_index);
        return l

class concaveHull():
    def clean_list(self, list_of_points):
        """
        Deletes duplicate points in list_of_points
        """
        return list(set(list_of_points))

    def length(self, vector):
        """
        Returns the number of elements in vector
        """
        return len(vector)

    def find_min_y_point(self, list_of_points):
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

    def add_point(self, vector, element):
        """
        Returns vector with the given element append to the right
        """
        vector.append(element)
        return vector

    def remove_point(self, vector, element):
        """
        Returns a copy of vector without the given element
        """
        vector.pop(vector.index(element))
        return vector

    def euclidian_distance(self, point1, point2):
        """
        Returns the euclidian distance of the 2 given points.
        :param point1: tuple (x, y)
        :param point2: tuple (x, y)
        :return: float
        """
        return math.sqrt(math.pow(point1[0] - point2[0], 2) + math.pow(point1[1] - point2[1], 2))

    def nearest_points(self, list_of_points, point, k):
        # build a list of tuples of distances between point *point* and every point in *list_of_points*, and
        # their respective index of list *list_of_distances*
        list_of_distances = []
        for index in range(len(list_of_points)):
            list_of_distances.append((self.euclidian_distance(list_of_points[index], point), index))

        # sort distances in ascending order
        list_of_distances.sort()

        # get the k nearest neighbors of point
        nearest_list = []
        for index in range(min(k, len(list_of_points))):
            nearest_list.append((list_of_points[list_of_distances[index][1]]))
        return nearest_list

    def angle(self, from_point, to_point):
        """
        Returns the angle of the directed line segment, going from *from_point* to *to_point*, in radians. The angle is
        positive for segments with upward direction (north), otherwise negative (south). Values ranges from 0 at the
        right (east) to pi at the left side (west).
        :param from_point: tuple (x, y)
        :param to_point: tuple (x, y)
        :return: float
        """
        return math.atan2(to_point[1] - from_point[1], to_point[0] - from_point[0])

    def angle_difference(self, angle1, angle2):
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

    def intersect(self, line1, line2):
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

    def point_in_polygon_q(self, point, list_of_points):
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

    def write_wkt(self, point_list, file_name):
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

    def as_wkt(self, point_list):
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

    def as_polygon(self, point_list):
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

    def enable_use_of_global_CRS(self):
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

    def disable_use_of_global_CRS(self, default_behaviour='prompt'):
        """
        Enables old settings again. If argument is missing then set behaviour to prompt.
        Example: disable_use_of_global_CRS(old_behaviour)
        :param default_behaviour:
        :return: None
        """
        settings = QSettings()
        settings.setValue('/Projections/defaultBehaviour', default_behaviour)
        return None

    def extract_points(self, geom):
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

    def sort_by_angle(self, list_of_points, last_point, last_angle):
        def getkey(item):
            return self.angle_difference(last_angle, self.angle(last_point, item))

        vertex_list = sorted(list_of_points, key=getkey, reverse=True)
        return vertex_list

    def concave_hull(self, points_list, k):
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
        point_set = self.clean_list(points_list)

        # if point_set has less then 3 points no polygon can be created and an empty list will be returned
        if len(point_set) < 3:
            return None

        # if point_set has 3 points then these are already vertices of the hull. Append the first point to
        # close the hull polygon
        if len(point_set) == 3:
            return self.add_point(point_set, point_set[0])

        # make sure that k neighbours can be found
        kk = min(kk, len(point_set))

        # start with the point having the smallest y-coordinate (most southern point)
        first_point = self.find_min_y_point(point_set)

        # add this points as the first vertex of the hull
        hull = [first_point]

        # make the first vertex of the hull to the current point
        current_point = first_point

        # remove the point from the point_set, to prevent him being among the nearest points
        point_set = self.remove_point(point_set, first_point)
        previous_angle = math.pi

        # step counts the number of segments
        step = 2

        # as long as point_set is not empty or search is returning to the starting point
        while (current_point != first_point) or (step == 2) and (len(point_set) > 0):

            # after 3 iterations add the first point to point_set again, otherwise a hull cannot be closed
            if step == 5:
                point_set = self.add_point(point_set, first_point)

            # search the k nearest neighbors of the current point
            k_nearest_points = self.nearest_points(point_set, current_point, kk)

            # sort the candidates (neighbors) in descending order of right-hand turn. This way the algorithm progresses
            # in clockwise direction through as many points as possible
            c_points = self.sort_by_angle(k_nearest_points, current_point, previous_angle)

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
                    its = self.intersect((hull[step - 2], c_points[i]), (hull[step - 2 - j], hull[step - 1 - j]))
                    j += 1

            # there is no candidate to which the connecting line does not intersect any existing segment, so the
            # for the next candidate fails. The algorithm starts again with an increased number of neighbors
            if its is True:
                return self.concave_hull(points_list, kk + 1)

            # the first point which complies with the requirements is added to the hull and gets the current point
            current_point = c_points[i]
            hull = self.add_point(hull, current_point)

            # calculate the angle between the last vertex and his precursor, that is the last segment of the hull
            # in reversed direction
            previous_angle = self.angle(hull[step - 1], hull[step - 2])

            # remove current_point from point_set
            point_set = self.remove_point(point_set, current_point)

            # increment counter
            step += 1

        all_inside = True
        i = len(point_set) - 1

        # check if all points are within the created polygon
        while (all_inside is True) and (i >= 0):
            all_inside = self.point_in_polygon_q(point_set[i], hull)
            i -= 1

        # since at least one point is out of the computed polygon, try again with a higher number of neighbors
        if all_inside is False:
            return self.concave_hull(points_list, kk + 1)

        # a valid hull has been constructed
        return hull

class catchmentAnalysis(QObject):

    # Setup signals
    finished = pyqtSignal(object)
    error = pyqtSignal(Exception, basestring)
    kill = pyqtSignal(bool)
    progress = pyqtSignal(float)
    warning = pyqtSignal(str)

    def __init__(self, iface, settings):
        QObject.__init__(self)
        self.concave_hull = concaveHull()
        self.iface = iface
        self.settings = settings
        self.killed = False

    def origin_preparation(self, origin_vector, origin_name_field):

        # Create a dictionary of origin point dictionaries containing name and geometry
        origins = []

        # Loop through origin and get or create points
        for i, f in enumerate(origin_vector.getFeatures()):

            # Get origin name
            if origin_name_field:
                origin_name = f[origin_name_field]
            else:
                origin_name = "origin_" + "%s" % (i+1)

            origins.append({'name': origin_name, 'geom': f.geometry().centroid()})

        return origins


    def graph_builder(self, network, cost_field, origins, tolerance, crs, epsg):

        # Settings
        otf = False

        # Get index of cost field
        network_fields = network.pendingFields()
        network_cost_index = network_fields.indexFromName(cost_field)

        # Setting up graph build director
        director = QgsLineVectorLayerDirector(network, -1, '', '', '', 3)

        # Determining cost calculation
        if cost_field:
            properter = customCost(network_cost_index, 0)
        else:
            properter = QgsDistanceArcProperter()

        # Creating graph builder
        director.addProperter(properter)
        builder = QgsGraphBuilder(crs, otf, tolerance, epsg)

        # Reading origins and making list of coordinates
        graph_origin_points = []

        # Loop through the origin points and add graph vertex indices
        for index, origin in enumerate(origins):
            graph_origin_points.append(origins[index]['geom'].asPoint())

        # Get origin graph vertex index
        tied_origin_vertices = director.makeGraph(builder, graph_origin_points)

        # Build the graph
        graph = builder.graph()

        # Create dictionary of origin names and tied origins
        tied_origins = {}

        # Combine origin names and tied point vertices
        for index, tied_origin in enumerate(tied_origin_vertices):
            tied_origins[index] = {'name': origins[index]['name'], 'vertex': tied_origin}

        return graph, tied_origins


    def graph_analysis(self, graph, tied_origins, distances):

        # Settings
        catchment_threshold = max(distances)

        # Variables
        catchment_network = {}
        catchment_points = {}

        # Loop through graph and get geometry and write to catchment network
        for index in range(0, graph.arcCount()):
            inVertexId = graph.arc(index).inVertex()
            outVertexId = graph.arc(index).outVertex()
            inVertexGeom = graph.vertex(inVertexId).point()
            outVertexGeom = graph.vertex(outVertexId).point()
            arcGeom = QgsGeometry.fromPolyline([inVertexGeom, outVertexGeom])
            catchment_network[index] = {'geom': arcGeom, 'cost': {}}

        # Loop through tied origins and write origin names
        for tied_point, origin in enumerate(tied_origins):
            # origin_name = tied_origins[i]['name']
            catchment_points[tied_point] = {distance: [] for distance in distances}

        # Loop through tied origins and write costs and polygon points
        for tied_point, origin in enumerate(tied_origins):

            # Kill check
            if self.killed == True:
                self.kill.emit(True)
                break

            origin_name = tied_origins[tied_point]['name']
            originVertexId = graph.findVertex(tied_origins[tied_point]['vertex'])

            # Run dijkstra and get tree and cost
            (tree, cost) = QgsGraphAnalyzer.dijkstra(graph, originVertexId, 0)

            # Loop through graph arcs
            for index in range(0, graph.arcCount()):

                # Define the arc properties
                inVertexId = graph.arc(index).inVertex()
                outVertexId = graph.arc(index).outVertex()
                inVertexGeom = graph.vertex(inVertexId).point()
                # outVertexGeom = graph.vertex(outVertexId).point()
                arcCost = max(cost[outVertexId], cost[inVertexId])

                # If arc is the origin set cost to 0
                if outVertexId == originVertexId or inVertexId == originVertexId:
                    catchment_network[index]['cost'][origin_name] = 0

                # If arc is connected and within the maximum radius set cost
                elif arcCost < catchment_threshold and tree[inVertexId] != -1:
                    if origin_name in catchment_network[index]['cost']:
                        if catchment_network[index]['cost'][origin_name] > int(arcCost):
                            catchment_network[index]['cost'][origin_name] = int(arcCost)
                    else:
                        catchment_network[index]['cost'][origin_name] = int(arcCost)

                # Add catchment points for each given radius
                for distance in distances:
                    if arcCost < distance:
                        catchment_points[tied_point][distance].extend([inVertexGeom, ])

        return catchment_network, catchment_points


    def network_writer(self, origins, catchment_network, output_network):

        # Variables
        arc_length_list = []

        # Setup all unique origin columns and minimum origin distance column
        unique_origin_list = []
        for origin in origins:
            name = origin['name']
            if not name in unique_origin_list:
                output_network.dataProvider().addAttributes([QgsField("%s" % name, QVariant.Int)])
                unique_origin_list.append(name)
        output_network.dataProvider().addAttributes([QgsField('min_dist', QVariant.Int)])
        output_network.updateFields()

        # Loop through arcs in catchment network and write geometry and costs
        for index in catchment_network:

            # Kill check
            if self.killed == True:
                self.kill.emit(True)
                break

            # Get arc properties
            arc_geom = catchment_network[index]['geom']
            arc_length = arc_geom.length()
            arc_cost_dict = catchment_network[index]['cost']
            arc_cost_list = []

            # Ignore arc if already processed
            if arc_length in arc_length_list:
                pass

            # Ignore arc is not connected to origins or outside catchment
            elif not arc_cost_dict:
                pass

            else:
                # Create feature and write id and geom
                f = QgsFeature(output_network.pendingFields())
                f.setAttribute("id", index)
                f.setGeometry(arc_geom)

                # Read the list of costs and write them to output network
                for name in arc_cost_dict:
                    cost = arc_cost_dict[name]
                    arc_cost_list.append(cost)
                    f.setAttribute(name, cost)
                    f.setAttribute("%s" % name, cost)

                # Set minimum cost
                if arc_cost_list > 0:
                    f.setAttribute('min_dist', min(arc_cost_list))

                # Write feature to output network layer
                output_network.dataProvider().addFeatures([f])

                # Add the length of arc to length list in order to ignore duplicates
                arc_length_list.append(arc_length)

        return output_network


    def polygon_writer(self, catchment_points, tied_origins, distances, output_polygon, polygon_tolerance):

        # Setup unique origin dictionary containing all distances
        polygon_dict = {}
        for tied_point in tied_origins:
            name = tied_origins[tied_point]['name']
            polygon_dict[name] = {distance: [] for distance in distances}

        # Create list of hulls for each distance and each unique origin
        index = 1
        hull_validity = True
        for tied_point in catchment_points:
            if hull_validity:
                name = tied_origins[tied_point]['name']
                for distance in distances:
                    points = catchment_points[tied_point][distance]
                    if len(points) > 2:  # Only three points can create a polygon
                        hull = self.concave_hull.concave_hull(points, polygon_tolerance)
                        polygon_dict[name][distance].append(hull)
                        print len(polygon_dict[name][distance])
                        if len(polygon_dict[name][distance]) < 2: # Write if only one hull
                            # Check if hull is a actual polygon
                            try:
                                p = QgsFeature(output_polygon.pendingFields())
                                p.setAttribute('id', index)
                                p.setAttribute('origin', name)
                                p.setAttribute('distance', distance)
                                polygon_geom = QgsGeometry.fromPolygon([hull, ])
                                p.setGeometry(polygon_geom)
                                output_polygon.dataProvider().addFeatures([p])
                                index += 1
                            except TypeError:
                                hull_validity = False
                                self.warning.emit('Polygon tolerance too high for cost band')
                                self.kill = True
                                break
                        else:
                            # Merge polygons when they intersect
                            polygon = None
                            for hull in polygon_dict[name][distance]:
                                geom = QgsGeometry.fromPolygon([hull, ])
                                if polygon:
                                    if polygon.intersects(geom):
                                        polygon = polygon.combine(geom)
                                else:
                                    polygon = geom

                            # Check if hull is a actual polygon
                            try:
                                p = QgsFeature(output_polygon.pendingFields())
                                p.setAttribute('id', index)
                                p.setAttribute('origin', name)
                                p.setAttribute('distance', distance)
                                p.setGeometry(polygon)
                                output_polygon.dataProvider().addFeatures([p])
                                index += 1
                            except TypeError:
                                hull_validity = False
                                self.warning.emit('Polygon tolerance too high for cost band')
                                self.kill = True
                                break

        return output_polygon


    def analysis(self):
        if self.settings:
            output = None
            try:
                # Prepare the origins
                origins = self.origin_preparation(
                    self.settings['origins'],
                    self.settings['name']
                )
                self.progress.emit(10)
                # Kill check
                if self.killed == True:
                    self.kill.emit(True)
                # Build the graph
                graph, tied_origins = self.graph_builder(
                    self.settings['network'],
                    self.settings['cost'],
                    origins,
                    self.settings['network tolerance'],
                    self.settings['crs'],
                    self.settings['epsg']
                )
                self.progress.emit(30)
                # Kill check
                if self.killed == True:
                    self.kill.emit(True)
                    # Run the analysis
                catchment_network, catchment_points = self.graph_analysis(
                    graph,
                    tied_origins,
                    self.settings['distances']
                )
                self.progress.emit(50)
                # Kill check
                if self.killed == True:
                    self.kill.emit(True)
                    # Write and render the catchment polygons
                if self.settings['output polygon check']:
                    output_polygon = self.polygon_writer(
                        catchment_points,
                        tied_origins,
                        self.settings['distances'],
                        self.settings['temp polygon'],
                        self.settings['polygon tolerance']
                    )
                    if self.settings['output polygon']:
                        uf.createShapeFile(output_polygon, self.settings['output polygon'], self.settings['crs'])
                        output_polygon = QgsVectorLayer(self.settings['output polygon'], 'catchment_areas', 'ogr')
                self.progress.emit(80)
                # Kill check
                if self.killed == True:
                    self.kill.emit(True)
                    # Write and render the catchment network
                if self.settings['output network check']:
                    output_network = self.network_writer(
                        origins,
                        catchment_network,
                        self.settings['temp network']
                    )
                    if self.settings['output network']:
                        uf.createShapeFile(output_network, self.settings['output network'], self.settings['crs'])
                        output_network = QgsVectorLayer(self.settings['output network'], 'catchment_network', 'ogr')

                if self.killed == False:
                    self.kill.emit(True)
                    self.progress.emit(100)
                    output = {'output network': output_network,
                              'output polygon': output_polygon,
                              'distances': self.settings['distances']}

            except Exception, e:
                self.error.emit(e, traceback.format_exc())
            self.finished.emit(output)


    def kill_analysis(self):
        self.killed = True



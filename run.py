
execfile(u'/Users/i.kolovou/Documents/Github/CatchmentAnalyser/catchment_analysis.py'.encode('utf-8'))
execfile(u'/Users/i.kolovou/Documents/Github/CatchmentAnalyser/utility_functions.py'.encode('utf-8'))



network = getLayerByName('test')
origins_vl = getLayerByName('test_u')
crs = network.crs()
epsg = network.crs().authid()[5:]
otf = False
tolerance = 0.01

network_fields = network.pendingFields()
# Setting up graph build director
director = QgsLineVectorLayerDirector(network, -1, '', '', '', 3)

# Determining cost calculation
#if cost_field:
#    properter = ct.CustomCost(network_cost_index, 0.01)
#else:
properter = QgsDistanceArcProperter()

# Creating graph builder
director.addProperter(properter)
builder = QgsGraphBuilder(crs, otf, tolerance, epsg)


origins = []
origin_name_field = None
# Loop through origin and get or create points
for i, f in enumerate(origins_vl.getFeatures()):
    # Get origin name
    if origin_name_field:
        origin_name = f[origin_name_field]
    else:
        origin_name = "origin_" + "%s" % (i+1)
    origins.append({'name': origin_name, 'geom': f.geometry().centroid()})


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

''''''''''''''''''''''''''''''



catchment_threshold = max(1200)

# Variables
catchment_network = {}
catchment_points = {}

# Loop through graph and get geometry and write to catchment network

# TODO: index is same as f.id()

for index in range(graph.arcCount()):
    inVertexId = graph.arc(index).inVertex()
    outVertexId = graph.arc(index).outVertex()
    inVertexGeom = graph.vertex(inVertexId).point()
    outVertexGeom = graph.vertex(outVertexId).point()
    # only include one of the two possible arcs
    if inVertexId < outVertexId:
        arcGeom = QgsGeometry.fromPolyline([inVertexGeom, outVertexGeom])
        catchment_network[index] = {'geom': arcGeom, 'start':inVertexId, 'end':outVertexId, 'cost': {}}

# Loop through tied origins and write origin names
for tied_point, origin in enumerate(tied_origins):
    origin_name = tied_origins[tied_point]['name']
    catchment_points[tied_point] = {'name': origin_name}
    catchment_points[tied_point].update({distance: [] for distance in distances})

# Loop through tied origins and write costs and polygon points
i = 1
for tied_point, origin in enumerate(tied_origins):
    self.progress.emit(20 + int(20 * i / len(tied_origins)))
    origin_name = tied_origins[tied_point]['name']
    originVertexId = graph.findVertex(tied_origins[tied_point]['vertex'])

    # Run dijkstra and get tree and cost
    (tree, cost) = QgsGraphAnalyzer.dijkstra(graph, originVertexId, 0)

    # Loop through graph arcs
    for index in catchment_network.iterkeys():
        if self.killed == True: break
        # Define the arc properties
        inVertexId = catchment_network[index]['start']
        outVertexId = catchment_network[index]['end']
        inVertexCost = cost[inVertexId]
        outVertexCost = cost[outVertexId]
        # this is the permissive option gives cost to the arc based on the closest point,
        # it just needs to be reach by one node
        arcCost = min(inVertexCost, outVertexCost)
        # this is the restrictive option, gives cost to the arc based on the furtherst point,
        # it needs to be entirely within distance
        #arcCost = max(inVertexCost, outVertexCost)

        # If arc is the origin set cost to 0
        if outVertexId == originVertexId or inVertexId == originVertexId:
            catchment_network[index]['cost'][origin_name] = 0

        # If arc is connected and within the maximum radius set cost
        elif arcCost <= catchment_threshold and tree[inVertexId] != -1:
            if origin_name in catchment_network[index]['cost']:
                if catchment_network[index]['cost'][origin_name] > int(arcCost):
                    catchment_network[index]['cost'][origin_name] = int(arcCost)
            else:
                catchment_network[index]['cost'][origin_name] = int(arcCost)

            # Add catchment points for each given radius
            inVertexGeom = graph.vertex(inVertexId).point()
            outVertexGeom = graph.vertex(outVertexId).point()
            seg_length = catchment_network[index]['geom'].length() #  math.sqrt(inVertexGeom.sqrDist(outVertexGeom))
            target_dist = 0
            for distance in distances:
                if self.killed == True: break
                # this option includes both nodes as long as arc is within distance
                # the polygon is the same as the network output
                #if arcCost <= distance:
                #    catchment_points[tied_point][distance].extend([inVertexGeom, outVertexGeom])
                # this option only includes nodes within distance
                # it does linear interpolation for extra points
                if inVertexCost <= distance:
                    catchment_points[tied_point][distance].append(inVertexGeom)
                    # add an extra point with linear referencing
                    if outVertexCost > distance:
                        target_dist = distance - inVertexCost
                        midVertexGeom = catchment_network[index]['geom'].interpolate(target_dist).asPoint()
                        catchment_points[tied_point][distance].append(midVertexGeom)
                if outVertexCost <= distance:
                    catchment_points[tied_point][distance].append(outVertexGeom)
                    # add an extra point with linear referencing
                    if inVertexCost > distance:
                        target_dist = distance - outVertexCost
                        midVertexGeom = catchment_network[index]['geom'].interpolate(seg_length-target_dist).asPoint()
                        catchment_points[tied_point][distance].append(midVertexGeom)

    i += 1

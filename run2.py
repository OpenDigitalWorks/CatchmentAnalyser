from qgis.core import *
from qgis.gui import *
from qgis.networkanalysis import *
from qgis.PyQt.QtCore import *
from qgis.PyQt.QtGui import *


# if geometry is null fails

execfile(u'/Users/i.kolovou/Documents/Github/CatchmentAnalyser/catchment_analysis.py'.encode('utf-8'))
execfile(u'/Users/i.kolovou/Documents/Github/CatchmentAnalyser/utility_functions.py'.encode('utf-8'))

vl = getLayerByName('test_cl')
#vl = getLayerByName('spm_cl_seg_small')
#vl = getLayerByName('grid_100_seg')

origins_vl = getLayerByName('test_u')
indices = {}#{ f.id(): [] for f in vl.getFeatures()}

class CustomCost(QgsArcProperter):
    def __init__(self, costColumIndex, defaultValue):
        QgsArcProperter.__init__(self)
        self.cost_column_index = costColumIndex
        self.default_value = defaultValue
    def property(self, distance, feature):
        cost = float(feature.attributes()[self.cost_column_index])
        if not cost or cost <= 0.0:
            return self.default_value
        else:
            return cost
    def requiredAttributes(self):
        return [0, self.cost_column_index]


director = QgsLineVectorLayerDirector(vl, -1, '', '', '', 3)
properter = CustomCost(13, 0.01)
properter = QgsDistanceArcProperter()
director.addProperter(properter)
crs = vl.crs()

builder = QgsGraphBuilder(crs)

pStart = QgsPoint(475869, 4203007)

tiedPoints = director.makeGraph(builder, [pStart])
graph = builder.graph()
tStart = tiedPoints[0]

idStart = graph.findVertex(tStart)

graph.arc(0).properties()
dir(graph.arc(0))

graph.arc(0).property(0)



indices = {}
idx = 0
for f in vl.getFeatures():
    if f.geometry().wkbType() == 2:
        for i in f.geometry().asPolyline()[1:]:
            indices[idx] = f.id()
            indices[idx + 1] = f.id()
            idx += 2
    elif f.geometry().wkbType() == 5:
        for pl in f.geometry().asMultiPolyline():
            for p in pl[1:]:
                indices[idx] = f.id()
                indices[idx + 1] = f.id()
                idx += 2

# write lines to see the indices
from PyQt4.QtCore import  QVariant
arc_features = []
flds = QgsFields()
flds.append(QgsField('arc_index', QVariant.Int))
flds.append(QgsField('feat_id', QVariant.Int))
for index in range(graph.arcCount()):
    new_feat = QgsFeature()
    inVertexId = graph.arc(index).inVertex()
    outVertexId = graph.arc(index).outVertex()
    inVertexGeom = graph.vertex(inVertexId).point()
    outVertexGeom = graph.vertex(outVertexId).point()
    # only include one of the two possible arcs
    if inVertexId < outVertexId:
        arcGeom = QgsGeometry.fromPolyline([inVertexGeom, outVertexGeom])
        new_feat.setGeometry(arcGeom)
        new_feat.setFields(flds)
        try:
            new_feat.setAttributes([index, indices[index]])
            arc_features.append(new_feat)
        except KeyError:
            print 'error', index


arc_layer = to_layer(arc_features, vl.crs(), vl.dataProvider().encoding(), 'Linestring', 'memory', None, 'arc_layer')
QgsMapLayerRegistry.instance().addMapLayer(arc_layer)













def to_layer(features, crs, encoding, geom_type, layer_type, path, name):
    first_feat = features[0]
    fields = first_feat.fields()
    layer = None
    if layer_type == 'memory':
        layer = QgsVectorLayer(geom_type + '?crs=' + crs.authid(), name, "memory")
        pr = layer.dataProvider()
        pr.addAttributes(fields.toList())
        layer.updateFields()
        layer.startEditing()
        pr.addFeatures(features)
        layer.commitChanges()
    elif layer_type == 'shapefile':
        wkbTypes = { 'Point': QGis.WKBPoint, 'Linestring': QGis.WKBLineString, 'Polygon': QGis.WKBPolygon }
        file_writer = QgsVectorFileWriter(path, encoding, fields, wkbTypes[geom_type], crs, "ESRI Shapefile")
        if file_writer.hasError() != QgsVectorFileWriter.NoError:
            print "Error when creating shapefile: ", file_writer.errorMessage()
        del file_writer
        layer = QgsVectorLayer(path, name, "ogr")
        pr = layer.dataProvider()
        layer.startEditing()
        pr.addFeatures(features)
        layer.commitChanges()
    return layer
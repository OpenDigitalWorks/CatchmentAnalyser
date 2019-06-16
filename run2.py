from qgis.core import *
from qgis.gui import *
from qgis.networkanalysis import *
from qgis.PyQt.QtCore import *
from qgis.PyQt.QtGui import *


# if geometry is null fails

execfile(u'/Users/i.kolovou/Documents/Github/CatchmentAnalyser/catchment_analysis.py'.encode('utf-8'))
execfile(u'/Users/i.kolovou/Documents/Github/CatchmentAnalyser/utility_functions.py'.encode('utf-8'))

vl = getLayerByName('test_cl')

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
#properter = QgsDistanceArcProperter()
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
            idx += 1
            indices[idx] = f.id()
            idx += 1
    elif f.geometry().wkbType() == 5:
        for pl in f.geometry().asMultiPolyline():
            for p in pl[1:]:
                indices[idx] = f.id()
                idx += 1
                indices[idx] = f.id()
                idx += 1












(tree, cost) = QgsGraphAnalyzer.dijkstra(graph, idStart, 0)

upperBound = []
r = 2000.0
i = 0
while i < len(cost):
  if cost[i] > r and tree[i] != -1:
    outVertexId = graph.arc(tree [i]).outVertex()
    if cost[outVertexId] < r:
      upperBound.append(i)
  i = i + 1

for i in upperBound:
  centerPoint = graph.vertex(i).point()
  rb = QgsRubberBand(qgis.utils.iface.mapCanvas(), True)
  rb.setColor(Qt.red)
  rb.addPoint(QgsPoint(centerPoint.x() - delta, centerPoint.y() - delta))
  rb.addPoint(QgsPoint(centerPoint.x() + delta, centerPoint.y() - delta))
  rb.addPoint(QgsPoint(centerPoint.x() + delta, centerPoint.y() + delta))
  rb.addPoint(QgsPoint(centerPoint.x() - delta, centerPoint.y() + delta))
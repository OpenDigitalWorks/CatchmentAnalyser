# Catchment Analyser
### Calculate catchments from multiple origins at multiple distances

## About
This plugin takes a line-based network and point-based origin layer and calculated the distance from each segment within a given distance or list of distances to each of the origins. The tool outputs the catchment as lines and as a concave hull polygon layer. Credit for the concave hull functionality goes to the algorithm described by Adriano Moreira and Maribel Yasmina Santos.

## Installation
Currently the plugin is not available through the QGIS plugins repository. To install you need to download the latest Plugin.zip file from the releases page (link). Unzip and copy the entire folder into the QGIS plugins directory. This directory can be found here:

* MS Windows: C:\Users\[your user name]\.qgis2\python\plugins\
* Mac OSX: Users/[your user name]/.qgis2/python/plugins/
* Linux: home/[your user name]/.qgis2/python/plugins/

This directory is usually hidden and you must make hidden files visible.Under Mac OSX you can open the folder in Finder by selecting 'Go > Go To Folder...' and typing '~/.qgis2/python/plugins/'. If you haven’t installed any QGIS plugins so far, you need to create the ‘plugins’ directory in the ‘.qgis2/python/’ directory. 

After copying the plugin, it will be available in the plugin manager window once you (re)start QGIS. Check the box next to the plugin to load it.

## How to
**Network layer**
Choose the line-based vector layer that comprises a topological network you want to analyse. The vector layer needs to have a projected CRS. Vector layers with a geographic CRS are ignored.

**Custom cost** 
By default the tool will generate catchments based on the length of the network segments. If you want to apply a custom cost for time-based catchments for example please check the 'Custom cost' checkbox. If checked, choose the field of the network layer which contains the cost information from the drop-down menu. Costs need to be numerical.

**Origin layer**
Choose the point-based vector layer containing the origins from which catchment will be calculated. 

**Cost bands**:     
By default the origin names will be based on the the feature ids. If you have specific names of the origins please check 'Origin names' checkbox. If checked, choose the field of the origin layer that contains the origin names information. The output network will be aggregated according the origin name.

**Network tolerance**: 
The network layer might contain un-connected lines. Gabs up to the network tolerance will be read as connected.

**Polygon tolerance**: 
The Catchment Analyser tool creates concave hull polygons describing the catchment from a specific origin. The polygon tolerance defines the level of 'concaveness' of the catchment. The lower the value, the more concave the catchment. The higher the value the more convex the catchment. 

**Catchment network**: 
The tool provides a catchment network output based on the original network layer with cost information on every origin. If checked the tool will generate the network as a temporary layer or as a shapefile using the browse button. By default the output network is renderer using the cost to the nearest origin.

**Catchment polygon**: 
The tool provides a catchment polygon output for each origin and each specified cost. If checked the tool will generate the polygons as a temporary layer or as a shapefile using the browse button. By default the output polygons are rendered as grey opaque areas. 

**Run**: 
Pressing the run button will activate the analysis for all the current settings.

**Cancel**: 
Pressing cancel will close and terminate the Catchment Analyser.

## Software Requirements
QGIS (2.0 or above) - http://www.qgis.org/en/site/


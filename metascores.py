## Calculate ecosystem services metascores for individual ES rasters
## For each input layer sum up all pixel values outside of clump polygon extents
## Developed for Richard Morris
## C. Doscher October 2022

import arcpy
# from arcpy.ia import RasterCalculator
# import math
from arcpy.sa import *

arcpy.env.overwriteOutput = True

# Get Clump polygon layer from user
inputFC = arcpy.GetParameterAsText(0)

# Get input ES raster from user
inputRas = arcpy.GetParameterAsText(1)

# Get workspace from user
arcpy.env.workspace = arcpy.GetParameterAsText(2)

# Get extent from user
arcpy.env.extent = arcpy.GetParameterAsText(3)

# Get boundary zome FC from user
boundary = arcpy.GetParameterAsText(4)

# Get cell size from user
# get cellsize from user
cellSize = arcpy.GetParameterAsText(5)

# get output file name from user
outName = arcpy.GetParameterAsText(6)

# convert clump polygon layer to raster called clumpras - use later in Zonal Stats tool
arcpy.PolygonToRaster_conversion(inputFC, "FID", "clumpras", "CELL_CENTER", "", 5)

# use Con to change NoData values to 0 and existing values to NoData
mask = Con((IsNull(Raster("clumpras"))), 0)
mask.save("mask.tif")

# Sum up pixel values outside of clumps using mask
arcpy.env.mask = mask
outZStat = arcpy.sa.ZonalStatisticsAsTable(boundary, "FID", inputRas, outName, "DATA", "SUM")
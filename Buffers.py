# Development script for FourInOneV2.py.  Create masks for nonlinear ES calculations from SPUs
# 27 May 2023 - CD

import arcpy
import math
from arcpy.sa import *


arcpy.env.overwriteOutput = True

# get input layer from user
inputFC = arcpy.GetParameterAsText(0)

# Get workspace from user
ws = arcpy.GetParameterAsText(1)
arcpy.env.workspace = ws

# Get boundary zone FC from user
boundary = arcpy.GetParameterAsText(2)
arcpy.env.extent = boundary

# get cellsize from user
cellSize = arcpy.GetParameterAsText(3)

# Ndist = 7.0
# dcalcBB = 500.0
# dcalcFT = 100.0

arcpy.CheckOutExtension("Spatial")

# raster lists here
listN = []
listFT = []
listBB = []

# Create buffers and raster grids for use in ES code
with arcpy.da.SearchCursor(inputFC, ['OBJECTID', 'Shape@', 'Shape_Area', 'CC', 'd', 'NEAR_DIST']) as cursor:
    for row in cursor:
        fid = row[0]
        arcpy.env.workspace = ws

        buffOutN = "buf_" + str(fid) + "_N.shp"
        arcpy.Buffer_analysis(row[1], buffOutN, "7 meters", "OUTSIDE_ONLY", "", "NONE", "", "PLANAR")
        buffrasN = "brasN_" + str(fid) + ".tif"
        arcpy.PolygonToRaster_conversion(buffOutN, "FID", buffrasN, "", "", cellSize)
        # brasreclassN = arcpy.Reclassify(buffOutN, "Value", RemapValue([1,1], ["NOTDATA", 0]))
        brasreclassN = Con(IsNull(Raster(buffrasN)), 0, 1)
        brasreclassN.save("brasrcN_" + str(fid) + ".tif")
        listN.append("brasrcN_" + str(fid) + ".tif")

        buffOutFT = "buf_" + str(fid) + "_FT.shp"
        arcpy.Buffer_analysis(row[1], buffOutFT, "100 meters", "OUTSIDE_ONLY", "", "NONE", "", "PLANAR")
        buffrasFT = "brasFT_" + str(fid) + ".tif"
        arcpy.PolygonToRaster_conversion(buffOutFT, "FID", buffrasFT, "", "", cellSize)
        # brasreclassFT = arcpy.Reclassify(buffOutN, "Value", RemapValue([1, 1], ["NOTDATA", 0]))
        brasreclassFT = Con(IsNull(Raster(buffrasFT)), 0, 1)
        brasreclassFT.save("brasrcFT_" + str(fid) + ".tif")
        listFT.append("brasrcFT_" + str(fid) + ".tif")

        buffOutBB = "buf_" + str(fid) + "_BB.shp"
        arcpy.Buffer_analysis(row[1], buffOutBB, "500 meters", "OUTSIDE_ONLY", "", "NONE", "", "PLANAR")
        buffrasBB = "brasBB_" + str(fid) + ".tif"
        arcpy.PolygonToRaster_conversion(buffOutBB, "FID", buffrasBB, "", "", cellSize)
        # brasreclassBB = arcpy.Reclassify(buffOutN, "Value", RemapValue([1, 1], ["NOTDATA", 0]))
        brasreclassBB = Con(IsNull(Raster(buffrasBB)), 0, 1)
        brasreclassBB.save("brasrcBB_" + str(fid) + ".tif")
        listBB.append("brasrcBB_" + str(fid) + ".tif")

# Add up N grids and create masks
addN = CellStatistics(listN, "SUM", "NODATA")
addNras = Raster(addN)
maxAddN = int(addNras.maximum)
addN.save("addN.tif")
# output mask for N
maskN = Reclassify(Raster("addN.tif"), "VALUE", RemapRange([[0, 0, 0], [0, maxAddN, 1]]))
maskN.save("maskN.tif")
# output nonlinear mask for N
nonlinmaskN = Reclassify(Raster("addN.tif"), "VALUE", RemapRange([[0, 1, 0], [1, maxAddN, 1]]))
nonlinmaskN.save("maskNnonlin.tif")

# Add up FT grids and create masks
addFT = CellStatistics(listFT, "SUM", "NODATA")
addFTras = Raster(addFT)
maxAddFT = int(addFTras.maximum)
addFT.save("addFT.tif")
# output mask for FT
maskFT = Reclassify(Raster("addFT.tif"), "VALUE", RemapRange([[0, 0, 0], [0, maxAddFT, 1]]))
maskFT.save("maskFT.tif")
# output nonlinear mask for N
nonlinmaskFT = Reclassify(Raster("addFT.tif"), "VALUE", RemapRange([[0, 1, 0], [1, maxAddFT, 1]]))
nonlinmaskFT.save("maskFTnonlin.tif")

# Add up BB grids and create masks
addBB = CellStatistics(listBB, "SUM", "NODATA")
addBBras = Raster(addFT)
maxAddBB = int(addBBras.maximum)
addBB.save("addBB.tif")
# output mask for BB
maskBB = Reclassify(Raster("addBB.tif"), "VALUE", RemapRange([[0, 0, 0], [0, maxAddBB, 1]]))
maskBB.save("maskBB.tif")
# output nonlinear mask for BB
nonlinmaskBB = Reclassify(Raster("addBB.tif"), "VALUE", RemapRange([[0, 1, 0], [1, maxAddBB, 1]]))
nonlinmaskBB.save("maskBBnonlin.tif")

# clear lists for next iteration
listN.clear()
listFT.clear()
listBB.clear()

# delete all except final masks
for ras in arcpy.ListRasters("*", "TIF"):
        if not ras.startswith("mask"):
                arcpy.Delete_management(ras)

for shp in arcpy.ListFeatureClasses():
    arcpy.Delete_management(shp)
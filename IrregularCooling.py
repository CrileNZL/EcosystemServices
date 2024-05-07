# Script to calculate cooling algorithm for Richard Morris
# C. Doscher - August 2022.  Updated May 2024 for irregular polygons

import arcpy
# from arcpy.ia import RasterCalculator
import math
from arcpy.sa import *

# arcpy.env.workspace = r"D:\OLW Nitrates\OLWNitrates.gdb" # memory?
arcpy.env.overwriteOutput = True

# get input layer from user
# parameters needed in input layer:
# * CC - cooling coefficient (float)
# * area to be calculated from Shape_Area (float)
# * d - distance parameter (float)

# inputFC = arcpy.GetParameterAsText(0)
inputFC = arcpy.GetParameterAsText(0)
# get workspace from user
ws = arcpy.GetParameterAsText(1)
# arcpy.env.workspace = arcpy.GetParameterAsText(1)
# arcpy.env.

# get extent from user
arcpy.env.extent = arcpy.GetParameterAsText(2)

# get cellsize from user
cellSize = arcpy.GetParameterAsText(3)

# get output file name from user
outName = arcpy.GetParameterAsText(4)

arcpy.CheckOutExtension("Spatial")
arcpy.CheckOutExtension("Foundation")

# while loop to iterate through each feature in input layer - cursor
# cur = arcpy.da.SearchCursor(input, #, "FID")
# row = cur.Next()\

# Setup irregular polygon for analysis: convert to line, generate points, create centreline, calculate distance to points
# create centerline
outCenterline = "centerline"
arcpy.PolygonToCenterline_topographic(inputFC, outCenterline, "")
# polygon to line
outLine = "SPULine"
arcpy.PolygonToLine_management(inputFC, outLine, "")
# autogenerate points
bdyp = "BDYPoints"
arcpy.GeneratePointsAlongLines_management(outLine, bdyp, 'DISTANCE', Distance="100 meters")
# Near bondary points to centerline
arcpy.Near_analysis(bdyp, outCenterline)
with arcpy.da.SearchCursor(bdyp, ['FID', 'Shape_Area', 'NEAR_DIST']) as cursor:
    for row in cursor:
        fid = row[0]
        d = row[2]
        distOut = "dist_" + str(fid) + ".tif"
        dist = arcpy.sa.EucDistance(row[0], d, 5)
        dist.save(distOut)

        # calculate cooling effect
        arcpy.env.workspace = ws
        distIn = Raster(distOut)
        ha = row[1] / 10000
        cc = 10
        # coolOut = (cc/(ha * Exp(Raster(distIn)/d)))
        coolOut = ha * cc * Exp(Raster(-1 * distIn) / d)
        coolOut.save("cool_" + str(fid) + ".tif")

del row
del cursor

# Combine all grids using Mosaic to New Raster
rasters = arcpy.ListRasters("cool*", "TIF")
output = outName + ".tif"
proj = arcpy.SpatialReference(2193)
arcpy.MosaicToNewRaster_management(rasters, ws, output, proj, "32_BIT_FLOAT", cellSize, "1", "MAXIMUM")

arcpy.CheckInExtension("Spatial")

# with points as cursor
#  grab NEAR_DIST = d
#  Euclidean Distance, max dist = d
#  Calculate cooling effect
# Combine all rasters into one Mosaic to New Raster sum? max? min? custom?

#
# with arcpy.da.SearchCursor(inputFC, ['FID', 'Shape@', 'Shape_Area', 'CC', 'd']) as cursor:
#     for row in cursor:
#
# # # Set up loop
# # while row: # don't know how many features the input will have
#     # for each feature, derive distance grid - Euclidean Distance
#     # get FID here for iterating
#         if(row[2]>=15000):
#                 fid = row[0]
#                 distOut = "dist_" + str(fid) + ".tif"
#                 arcpy.env.workspace = ws
#                 #maxDist = '#' # can use ""?
#                 dist = arcpy.sa.EucDistance(row[1], "", cellSize)
#                 dist.save(distOut)
#
#                 # get parameter values from input
#                 cc = float(row[3])
#                 d = float(row[4])
#                 ha = float(row[2])/10000
#                 # expression = "%cc%/(%ha% * exp(distOut/%d%))"
#
#                 # Calculate cooling raster
#                 arcpy.env.workspace = ws
#                 distIn = Raster(distOut)
#                 #coolName = "cool_" + str(fid) + ".tif"
#                 # expExp = Raster(math.exp(distIn/d))
#                 # nextExp = Raster(ha * expExp)
#                 # coolOut = Raster(cc/nextExp)
#                 # coolOut = float(%cc%)/( float(%ha%) * Exp(distOut/ float(%d%)
#
#                 # coolOut = (cc/(ha * Exp(Raster(distIn)/d)))
#                 coolOut = ha * cc * Exp(Raster(-1 * distIn)/d)
#                 coolOut.save("cool_" + str(fid) + ".tif")
#                 #RasterCalculator(distOut, "distOut", expression, coolName) # workspace?
#
# del row
# del cursor
#
# rasters = arcpy.ListRasters("cool*", "TIF")
# output = outName + ".tif"
# proj = arcpy.SpatialReference(2193)
# arcpy.MosaicToNewRaster_management(rasters, ws, output, proj, "32_BIT_FLOAT", cellSize, "1", "MAXIMUM")
#
# arcpy.CheckInExtension("Spatial")

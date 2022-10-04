# Script to implement an algorithm for Richard Morris based on Laca, 2021
# Multi-scape interventions to match spatial scale of demand and supply of ecosystem services
# Frontiers in Sustainable Food Systems, vol 4
# C. Doscher - September 2022

import arcpy
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

# get value of d from user
d = arcpy.GetParameterAsText(4)
dcalc = float(d)

# get output file name from user
outName = arcpy.GetParameterAsText(5)

arcpy.CheckOutExtension("Spatial")

# while loop to iterate through each feature in input layer - cursor
# cur = arcpy.da.SearchCursor(input, #, "FID")
# row = cur.Next()

with arcpy.da.SearchCursor(inputFC, ['FID', 'Shape@', 'Shape_Area']) as cursor:
    for row in cursor:

# Set up loop
# while row: # don't know how many features the input will have
    # for each feature, derive distance grid - Euclidean Distance
    # get FID here for iterating
        fid = row[0]
        distOut = "dist_" + str(fid) + ".tif"
        arcpy.env.workspace = ws
        #maxDist = '#' # can use ""?
        dist = arcpy.sa.EucDistance(row[1], "", cellSize)
        dist.save(distOut)

        # Calculate supply raster
        arcpy.env.workspace = ws
        distIn = Raster(distOut)


        # whereClause = "VALUE <= dcalc"
        lacaOut = Con(distIn <= dcalc, ((1/dcalc)*1.094*(1 - (1/dcalc)**2*Raster(distIn)**2)**3), 0)
        lacaOut.save("laca_" + str(fid) + ".tif")
        #RasterCalculator(distOut, "distOut", expression, coolName) # workspace?

del row
del cursor

rasters = arcpy.ListRasters("laca*", "TIF")
output = outName + ".tif"
proj = arcpy.SpatialReference(2193)
#proj = "clump_export.prj"
# wkt = PROJCS["NZGD2000 / New Zealand Transverse Mercator 2000",GEOGCS["NZGD2000",DATUM["D_NZGD_2000",SPHEROID["GRS_1980",6378137,298.257222101]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",173],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",1600000],PARAMETER["false_northing",10000000],UNIT["Meter",1]]
arcpy.MosaicToNewRaster_management(rasters, ws, output, proj, "32_BIT_FLOAT", cellSize, "1", "SUM")

arcpy.CheckInExtension("Spatial")
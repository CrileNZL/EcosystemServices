## Calculate ecosystem services layers and metascores for individual ES rasters
## Combines Cooling.py, Nitrogen.py, Laca.py and metascores.py into one script
## Developed for Richard Morris
## C. Doscher October 2022 - Updated 7 Nov 2022

import arcpy
# from arcpy.ia import RasterCalculator
import math

from arcpy import Raster
from arcpy.sa import *

arcpy.env.overwriteOutput = True

# Get Clump polygon layer from user
inputFC = arcpy.GetParameterAsText(0)

# Get input ES raster from user - raster names will be internally generated
# inputRas = arcpy.GetParameterAsText(1)

# Get workspace from user
ws = arcpy.GetParameterAsText(1)
arcpy.env.workspace = ws

# Get extent from user
# arcpy.env.extent = arcpy.GetParameterAsText(3)

# get value of d from user
species = arcpy.GetParameterAsText(2)
# if species == "Fantail":
    # dcalc = 50.0
# elif species == "Bellbird":
    # dcalc = 500.0

# Get boundary zone FC from user
boundary = arcpy.GetParameterAsText(3)
arcpy.env.extent = boundary

# Get cell size from user
# get cellsize from user
cellSize = arcpy.GetParameterAsText(4)

# get output file name from user
outName = arcpy.GetParameterAsText(5)

arcpy.CheckOutExtension("Spatial")

# Run NEAR on inputFC for FT habitat calc
arcpy.Near_analysis(inputFC, inputFC)

## Create distance grids for use in ES code
# can I do this in one script given minimmum clump size in Cooling.py?
with arcpy.da.SearchCursor(inputFC, ['OBJECTID', 'Shape@', 'Shape_Area', 'CC', 'd', 'NEAR_DIST']) as cursor:
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

        # get parameter values from input
        # cc = float(row[3])
        # d = float(row[4])
        # ha = float(row[2])/10000
        # expression = "%cc%/(%ha% * exp(distOut/%d%))"

        # Calculate cooling raster
        distIn = Raster(distOut)

## Calculate individual ESs
# Cooling and Bellbird habitat here
        if(row[2]>=15000):
                cc = float(row[3])
                d = float((row[2]/math.pi)**0.5)
                ha = float(row[2])/10000
                coolOut = ha * cc * Exp(Raster(-1 * distIn)/d)
                coolOut.save("cool_" + str(fid) + ".tif")

                dcalcBB = 500.00
                lacaBBOut = Con(distIn <= dcalcBB,
                                ((1 / dcalcBB) * 1.094 * (1 - (1 / dcalcBB) ** 2 * Raster(distIn) ** 2) ** 3), 0)
                lacaBBOut.save("lacaBB_" + str(fid) + ".tif")

# Nitrogen here
        nOut = Con(distIn <= 7, (-3.9 * (distIn)**3 + 89.1 * (distIn)**2 - (814.5 * distIn) + 2968), 0)
        nOut.save("N_" + str(fid) + ".tif")


# Fantail Habitat here - include SPUs greater than 1.5 ha and within
        if(row[2] < 15000 and row[5] < 150.0) or row[2] >= 15000:
                dcalcFT = 50.0
                lacaFTOut = Con(distIn <= dcalcFT, ((1/dcalcFT)*1.094*(1 - (1/dcalcFT)**2*Raster(distIn)**2)**3), 0)
                lacaFTOut.save("lacaFT_" + str(fid) + ".tif")

del row
del cursor

## Calculate final layer for each ES
# Cooling
rasterList = arcpy.ListRasters("cool*", "TIF")
outputC = outName + "_cool.tif"
proj = arcpy.SpatialReference(2193)
if len(rasterList) > 0:
        # outputC = outName + "_cool.tif"
        # proj = arcpy.SpatialReference(2193)
        arcpy.MosaicToNewRaster_management(rasterList, ws, "midcool", proj, "32_BIT_FLOAT", cellSize, "1", "MAXIMUM")
        inConstant = 0.75
        outTimes = Times("midcool", inConstant)
        outTimes.save(outputC)

# Nitrogen MS raster
rasters = arcpy.ListRasters("N_*", "TIF")
outputN = outName + "_nitrogen.tif"
proj = arcpy.SpatialReference(2193)
arcpy.MosaicToNewRaster_management(rasters, ws, outputN, proj, "32_BIT_FLOAT", cellSize, "1", "SUM")


# Bellbird Habitat MS raster
rasters = arcpy.ListRasters("lacaBB*", "TIF")
outputHBB = outName + "_habitatBB.tif"
proj = arcpy.SpatialReference(2193)
arcpy.MosaicToNewRaster_management(rasters, ws, outputHBB, proj, "32_BIT_FLOAT", cellSize, "1", "MAXIMUM")

# Fantail  Habitat MS raster
rasters = arcpy.ListRasters("lacaFT*", "TIF")
outputHFT = outName + "_habitatFT.tif"
proj = arcpy.SpatialReference(2193)
arcpy.MosaicToNewRaster_management(rasters, ws, outputHFT, proj, "32_BIT_FLOAT", cellSize, "1", "MAXIMUM")

## Calculate metascores for each ES
# Need this to loop through each ES output raster
# convert clump polygon layer to raster called clumpras - use later in Zonal Stats tool
arcpy.PolygonToRaster_conversion(inputFC, "OBJECTID", "clumpras.tif", "CELL_CENTER", "", 5)

# use Con to change NoData values to 0 and existing values to NoData
mask = Con((IsNull("clumpras.tif")), 0)
mask.save("mask.tif")

# Sum up pixel values outside of clumps using mask
arcpy.env.mask = mask

# Calculate metascores.  Skip C if it doesn't exist
arcpy.sa.ZonalStatisticsAsTable(boundary, "OBJECTID", Raster(outputN), outName + "_MSNitrogen.dbf", "DATA", "SUM")
arcpy.sa.ZonalStatisticsAsTable(boundary, "OBJECTID", Raster(outputHBB), outName + "_MSBBHabitat.dbf", "DATA", "SUM")
arcpy.sa.ZonalStatisticsAsTable(boundary, "OBJECTID", Raster(outputHFT), outName + "_MSFTHabitat.dbf", "DATA", "SUM")
if arcpy.Exists(Raster(outputC)):
    arcpy.sa.ZonalStatisticsAsTable(boundary, "OBJECTID", Raster(outputC), outName + "_MSCool.dbf", "DATA", "SUM")

## Clean up
# Delete distance and individual ES grids
for ras in arcpy.ListRasters("*", "TIF"):
        if not ras.startswith(outName):
                arcpy.Delete_management(ras)
# oldRasters = arcpy.ListRasters("", "TIF")
# for ras in oldRasters:
#         arcpy.Delete_management(ras)

arcpy.CheckInExtension("Spatial")
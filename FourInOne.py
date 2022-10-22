## Calculate ecosystem services layers and metascores for individual ES rasters
## Combines Cooling.py, Nitrogen.py, Laca.py and metascores.py into one script
## Developed for Richard Morris
## C. Doscher October 2022

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
if species == "Fantail":
    dcalc = 100.0
elif species == "Bellbird":
    dcalc = 500.0

# Get boundary zone FC from user
boundary = arcpy.GetParameterAsText(3)
arcpy.env.extent = boundary

# Get cell size from user
# get cellsize from user
cellSize = arcpy.GetParameterAsText(4)

# get output file name from user
outName = arcpy.GetParameterAsText(5)

arcpy.CheckOutExtension("Spatial")

## Create distance grids for use in ES code
# can I do this in one script given minimmum clump size in Cooling.py?
with arcpy.da.SearchCursor(inputFC, ['FID', 'Shape@', 'Shape_Area', 'CC', 'd']) as cursor:
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
# Cooling here
        if(row[2]>=15000):
                cc = float(row[3])
                d = float((row[2]/math.pi)**0.5)
                ha = float(row[2])/10000
                coolOut = ha * cc * Exp(Raster(-1 * distIn)/d)
                coolOut.save("cool_" + str(fid) + ".tif")

# Nitrogen here
        nOut = Con(distIn <= 7, (-3.9 * (distIn)**3 + 89.1 * (distIn)**2 - (814.5 * distIn) + 2968), 0)
        nOut.save("N_" + str(fid) + ".tif")

# Habitat here
        lacaOut = Con(distIn <= dcalc, ((1/dcalc)*1.094*(1 - (1/dcalc)**2*Raster(distIn)**2)**3), 0)
        lacaOut.save("laca_" + str(fid) + ".tif")

del row
del cursor

## Calculate final layer for each ES
# Cooling
rasters = arcpy.ListRasters("cool*", "TIF")
outputC = outName + "_cool.tif"
proj = arcpy.SpatialReference(2193)
if len(rasters) > 0:
        # outputC = outName + "_cool.tif"
        # proj = arcpy.SpatialReference(2193)
        arcpy.MosaicToNewRaster_management(rasters, ws, outputC, proj, "32_BIT_FLOAT", cellSize, "1", "MAXIMUM")

# Nitrogen
rasters = arcpy.ListRasters("N_*", "TIF")
outputN = outName + "_nitrogen.tif"
proj = arcpy.SpatialReference(2193)
arcpy.MosaicToNewRaster_management(rasters, ws, outputN, proj, "32_BIT_FLOAT", cellSize, "1", "SUM")


# Habitat
rasters = arcpy.ListRasters("laca*", "TIF")
outputH = outName + "_habitat.tif"
proj = arcpy.SpatialReference(2193)
#proj = "clump_export.prj"
# wkt = PROJCS["NZGD2000 / New Zealand Transverse Mercator 2000",GEOGCS["NZGD2000",DATUM["D_NZGD_2000",SPHEROID["GRS_1980",6378137,298.257222101]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",173],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",1600000],PARAMETER["false_northing",10000000],UNIT["Meter",1]]
arcpy.MosaicToNewRaster_management(rasters, ws, outputH, proj, "32_BIT_FLOAT", cellSize, "1", "SUM")


## Calculate metascores for each ES
# Need this to loop through each ES output raster
# convert clump polygon layer to raster called clumpras - use later in Zonal Stats tool
arcpy.PolygonToRaster_conversion(inputFC, "FID", "clumpras.tif", "CELL_CENTER", "", 5)

# use Con to change NoData values to 0 and existing values to NoData
mask = Con((IsNull("clumpras.tif")), 0)
mask.save("mask.tif")

# Sum up pixel values outside of clumps using mask
arcpy.env.mask = mask

# Calculate metascores.  Skip C if it doesn't exist
arcpy.sa.ZonalStatisticsAsTable(boundary, "FID", Raster(outputN), outName + "_MSNitrogen.dbf", "DATA", "SUM")
arcpy.sa.ZonalStatisticsAsTable(boundary, "FID", Raster(outputH), outName + "_MSHabitat.dbf", "DATA", "SUM")
if arcpy.Exists('Raster(outputC)'):
    arcpy.sa.ZonalStatisticsAsTable(boundary, "FID", Raster(outputC), outName + "_MSCool.dbf", "DATA", "SUM")

## Clean up
# Delete distance and individual ES grids
for ras in arcpy.ListRasters("*", "TIF"):
        if not ras.startswith(outName):
                arcpy.Delete_management(ras)
# oldRasters = arcpy.ListRasters("", "TIF")
# for ras in oldRasters:
#         arcpy.Delete_management(ras)

arcpy.CheckInExtension("Spatial")
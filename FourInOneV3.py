# Calculate ecosystem services layers and metascores for individual ES rasters
# Combines Cooling.py, Nitrogen.py, Laca.py and metascores.py into one script
# Developed for Richard Morris
# C. Doscher October 2022 - Updated 19 May 2023
# New version of 24 April script
# N calculation updated on 26 May 2023
# new BB and error catching for FT on 2 June 2023
# Added code to create masks for final analysis
# Updated 4 June - added error catching code for nonlinear mask creation - if max value = 0

import arcpy

import math

from arcpy.sa import *

arcpy.env.overwriteOutput = True

# Get Clump polygon layer from user
inputFC = arcpy.GetParameterAsText(0)

# Get workspace from user
ws = arcpy.GetParameterAsText(1)
arcpy.env.workspace = ws

# EL model parameters for cooling
asymC = float(arcpy.GetParameterAsText(2))
midC = float(arcpy.GetParameterAsText(3))
kC = float(arcpy.GetParameterAsText(4))

# EL model parameters for nitrogen
asymN = float(arcpy.GetParameterAsText(5))
midN = float(arcpy.GetParameterAsText(6))
kN = float(arcpy.GetParameterAsText(7))

# EL model parametersf for Bellbirds
asymBB = float(arcpy.GetParameterAsText(8))
midBB = float(arcpy.GetParameterAsText(9))
kBB = float(arcpy.GetParameterAsText(10))

# EL parameters for Fantails
asymFT = float(arcpy.GetParameterAsText(11))
midFT = float(arcpy.GetParameterAsText(12))
kFT = float(arcpy.GetParameterAsText(13))

# Get boundary zone FC from user
boundary = arcpy.GetParameterAsText(14)
arcpy.env.extent = boundary

# Get cell size from user
# get cellsize from user
cellSize = arcpy.GetParameterAsText(15)
# cellSizeHa = (float(cellSize)**2)/10000.0

# get output file name from user
outName = arcpy.GetParameterAsText(16)

dcalcBB = 500.0
dcalcFT = 100.0

arcpy.CheckOutExtension("Spatial")

# Create masks to limit area of application for nonlinear equation

listN = []
listFT = []
listBB = []

# Create buffers and raster grids for use in ES code
with arcpy.da.SearchCursor(inputFC, ['OBJECTID', 'Shape@']) as cursor:
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

del row
del cursor

# Add up N grids and create masks
addN = CellStatistics(listN, "SUM", "NODATA")
addNras = Raster(addN)
maxAddN = int(addNras.maximum)
addN.save("addN.tif")

# output mask for N
if maxAddN <= 0:
    pass
else:
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
if maxAddFT <= 0:
    pass
else:
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
if maxAddBB <= 0:
    pass
else:
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

# Create mask to set SPU area to 0 for Zonal Statistics as Table
arcpy.PolygonToRaster_conversion(inputFC, "Shape_Area", "Mask.tif", "CELL_CENTER", "", 5)
# use Con to change NoData values to 0 and existing values to NoData
nCon = Con(IsNull("Mask.tif"), 1, 0)
nCon.save("SPUMask.tif")

# Run NEAR on inputFC for FT habitat calc
arcpy.Near_analysis(inputFC, inputFC)

# Create distance grids for use in ES code
with arcpy.da.SearchCursor(inputFC, ['OBJECTID', 'Shape@', 'Shape_Area', 'CC', 'd', 'NEAR_DIST']) as cursor:
    for row in cursor:

        # Set up loop
        # for each feature, derive distance grid - Euclidean Distance
        # get FID here for iterating
        fid = row[0]
        distOut = "dist_" + str(fid) + ".tif"
        arcpy.env.workspace = ws
        # maxDist = '#' # can use ""?
        dist = arcpy.sa.EucDistance(row[1], "", cellSize)
        dist.save(distOut)

        # get parameter values from input
        # cc = float(row[3])
        # d = float(row[4])
        # ha = float(row[2])/10000
        # expression = "%cc%/(%ha% * exp(distOut/%d%))"

        # Calculate cooling raster
        distIn = Raster(distOut)

# Calculate individual ESs
# Cooling and Bellbird habitat here
        # select all SPUs great than 4900 m2 in area and those smaller that are within 2*R of another SPU
        if row[2] >= 4900 or (row[2] < 4900 and row[5] <= float(2 * ((row[2] / math.pi)**0.5))):
                cc = float(row[3])
                radius = float((row[2] / math.pi)**0.5)
                if row[2] > 4000.00 and row[2] < 5000.00:
                    d = 2 * radius/1.75
                else:
                    d = 2 * radius

                ha = float(row[2])/10000
                coolOut = ha * cc * Exp(Raster(-1 * distIn)/d)
                coolOut.save("cool_" + str(fid) + ".tif")

                # dcalcBB = 500.00
        if row[2] >= 15000 or (row[2] > 950 and row[5] < 10.0):
                lacaBBOut = Con(distIn <= dcalcBB,
                                ((1 / dcalcBB) * 1.094 * (1 - (1 / dcalcBB) ** 2 * Raster(distIn) ** 2) ** 3), 0)
                lacaBBOut.save("lacaBB_" + str(fid) + ".tif")

# Nitrogen here
        # nOut = Con(distIn <= 7, (-3.9 * distIn ** 3 + 89.1 * distIn ** 2 - (814.5 * distIn) + 2968), 0)
        # Updated by RM 26 May 23
        nOut = Con(distIn <= 7, (-235 * distIn + 2065), 0)
        nOut.save("N_" + str(fid) + ".tif")


# Fantail Habitat here - include SPUs greater than 1.5 ha and within 150 m of another SPU of any size
        # if row[2] >= 15000 or (row[2] < 15000 and row[5] < 150.0):
        # dcalcFT = 50.0
        lacaFTOut = Con(distIn <= dcalcFT, ((1/dcalcFT)*1.094*(1 - (1/dcalcFT)**2*Raster(distIn)**2)**3), 0)
        lacaFTOut.save("lacaFT_" + str(fid) + ".tif")

del row
del cursor

# Calculate final layer for each ES
# Cooling
rasterC = arcpy.ListRasters("cool*", "TIF")
outputC = outName + "_cool.tif"
proj = arcpy.SpatialReference(2193)
if len(rasterC) > 0:
        # outputC = outName + "_cool.tif"
        # proj = arcpy.SpatialReference(2193)
        arcpy.MosaicToNewRaster_management(rasterC, ws, "midcool.tif", proj, "32_BIT_FLOAT", cellSize, "1", "SUM")
        # multiply output by 0.75 and then 100/122 to standardise 0 - 100
        stdCoolValue = 0.6148
        # inConstant = 0.75
        outTimes = Times(Raster("midcool.tif"), stdCoolValue)
        # outTimes = Times(Raster("midcool.tif"), inConstant)
        outTimes.save(outName + "_baseCooling.tif")
        # coolOverlap = 30 / (1 + Exp(4.365 - Raster("timescool.tif"))
        # user inputs values for ASYM, MID and k
        coolOverlap = asymC / (1 + Exp((midC - Raster(outTimes))/kC))
        coolOverlap.save("nonlinCool.tif")
        finalCool = Times(Raster("nonlinCool.tif"), Raster("SPUMask.tif"))
        finalCool.save(outputC)

# Nitrogen MS raster
rastersN = arcpy.ListRasters("N_*", "TIF")
stdNit = outName + "_baseNitrogen.tif"
outputN = outName + "_nitrogen.tif"
proj = arcpy.SpatialReference(2193)
arcpy.MosaicToNewRaster_management(rastersN, ws, "tempNit.tif", proj, "32_BIT_FLOAT", cellSize, "1", "SUM")
# standardise base nitrogen by multiplying by 100/7603
nitStdValue = 0.01315
stdNitras = Times(Raster("tempNit.tif"), nitStdValue)
stdNitras.save(stdNit)
# Use EL model for nonlinear effects - parameters set by user
if arcpy.Exists("maskNnonlin.tif"):
    ELNitrogen = Con(Raster("maskNnonlin.tif") == 1, (asymN / (1 + Exp((midN - Raster(stdNit))/kN))), Raster(stdNit))
    ELNitrogen.save(outputN)
    # finalN = Times(Raster("nonlinNit.tif"), Raster("SPUMask.tif"))
    # finalN.save(outputN)
else:
    stdNitras.save(outputN)

# Bellbird Habitat MS raster
rastersBB = arcpy.ListRasters("lacaBB*", "TIF")
stdHBBras = outName + "_baseBellbirdHabitat.tif"
outputHBB = outName + "_BellbirdHabitat.tif"
proj = arcpy.SpatialReference(2193)
if len(rastersBB) > 0:
    arcpy.MosaicToNewRaster_management(rastersBB, ws, "tempBB.tif", proj, "32_BIT_FLOAT", cellSize, "1", "SUM")
    # standardise base Bellbird by multiplying by 100/0.016
    bbStdValue = 6250
    stdBB = Times(Raster("tempBB.tif"), bbStdValue)
    stdBB.save(stdHBBras)
    # Use EL Model for nonlinear effects - user supplies parameters
    if arcpy.Exists("maskBBnonlin.tif"):
        ELHBB = Con(Raster("maskBBnonlin.tif") == 1, (asymBB / (1 + Exp((midBB - Raster(stdBB))/kBB))), Raster(stdBB))
        ELHBB.save("nonlinBB.tif")
        finalBB = Times(Raster("nonlinBB.tif"), Raster("SPUMask.tif"))
        finalBB.save(outputHBB)
    else:
        stdBB.save(outputHBB)


# Fantail Habitat MS raster
rastersFT = arcpy.ListRasters("lacaFT*", "TIF")
midFTras = outName + "_baseFantailHabitat.tif"
outputHFT = outName + "_FantailHabitat.tif"
proj = arcpy.SpatialReference(2193)
if len(rastersFT) > 0:
    arcpy.MosaicToNewRaster_management(rastersFT, ws, "tempFT.tif", proj, "32_BIT_FLOAT", cellSize, "1", "SUM")
    # standardise FT by multiplying by 100/0.24
    stdFTValue = 416.7
    stdFT = Times(Raster("tempFT.tif"), stdFTValue)
    stdFT.save(midFTras)

    if arcpy.Exists("maskFTnonlin.tif"):
        ELHFT = Con(Raster("maskFTnonlin.tif") == 1, (asymFT / (1 + Exp((midFT - Raster(stdFT))/kFT))), Raster(stdFT))
        ELHFT.save("nonlinFT.tif")
        finalFT = Times(Raster("nonlinFT.tif"), Raster("SPUMask.tif"))
        finalFT.save(outputHFT)
    else:
        stdFT.save(outputHFT)

# use Zonal Statistics as Table to calculate ES metascore
arcpy.sa.ZonalStatisticsAsTable(boundary, "OBJECTID", Raster(outputN), outName + "_MSNitrogen.dbf", "DATA", "SUM")
# arcpy.sa.ZonalStatisticsAsTable(boundary, "OBJECTID", Raster(outName + "_Ncontrol.tif"), outName + "_MSNControl.dbf", "DATA", "SUM")

# Calculate Cooling and Bellbird  control metascore DBFs
if arcpy.Exists(outputHBB):
    arcpy.sa.ZonalStatisticsAsTable(boundary, "OBJECTID", Raster(outputHBB), outName + "_MSBBHabitat.dbf", "DATA", "SUM")
    # arcpy.sa.ZonalStatisticsAsTable(boundary, "OBJECTID", Raster(outName + "_BBcontrol.tif"), outName + "_MSBBControl.dbf",
                                        # "DATA", "SUM")
if arcpy.Exists(outputC):
    arcpy.sa.ZonalStatisticsAsTable(boundary, "OBJECTID", Raster(outputC), outName + "_MSCool.dbf", "DATA", "SUM")
    # arcpy.sa.ZonalStatisticsAsTable(boundary, "OBJECTID", Raster(outName + "_CoolControl.tif"), outName + "_MSCControl.dbf", "DATA",
                                    # "SUM")

# Calculate control metascores.  Skip C if it doesn't exist
# need a custom mask for each ES - Nitrogen can use the already existing mask from above
# arcpy.env.mask = hFTCM
if (arcpy.Exists(outputHFT)):
    arcpy.sa.ZonalStatisticsAsTable(boundary, "OBJECTID", Raster(outputHFT), outName + "_MSFTHabitat.dbf", "DATA", "SUM")
    # arcpy.sa.ZonalStatisticsAsTable(boundary, "OBJECTID", Raster(outName + "_FTcontrol.tif"), outName + "_MSFTControl.dbf", "DATA", "SUM")

# Clean up
# Delete distance and individual ES grids
for ras in arcpy.ListRasters("*", "TIF"):
        if not ras.startswith(outName):
                arcpy.Delete_management(ras)


arcpy.CheckInExtension("Spatial")

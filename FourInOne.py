# Calculate ecosystem services layers and metascores for individual ES rasters
# Combines Cooling.py, Nitrogen.py, Laca.py and metascores.py into one script
# Developed for Richard Morris
# C. Doscher October 2022 - Updated 26 March 2023

import arcpy

import math

from arcpy.sa import *

arcpy.env.overwriteOutput = True

# Get Clump polygon layer from user
inputFC = arcpy.GetParameterAsText(0)
# get attribute values for inputFC
# with arcpy.da.SearchCursor(inputFC, ['Shape_Area', 'CC']) as cursor:
#     for row in cursor:
#         cc = float(row[1])
#         d = 2 * math.sqrt(float((row[0] / math.pi)))
#         ha = float(row[0]) / 10000
#
# del row
# del cursor

# Get input ES raster from user - raster names will be internally generated
# inputRas = arcpy.GetParameterAsText(1)

# Get workspace from user
ws = arcpy.GetParameterAsText(1)
arcpy.env.workspace = ws

# Get extent from user
# arcpy.env.extent = arcpy.GetParameterAsText(3)

# get value of d from user
# species = arcpy.GetParameterAsText(2)
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
dcalcFT = 50.0

arcpy.CheckOutExtension("Spatial")

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
        # select all SPUs great than 4900 m2 in area and those smaller that are within 2*D of another SPU
        if row[2] >= 4900 or (row[2] < 4900 and row[5] <= float(4 * ((row[2] / math.pi)**0.5))):
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
                lacaBBOut = Con(distIn <= dcalcBB,
                                ((1 / dcalcBB) * 1.094 * (1 - (1 / dcalcBB) ** 2 * Raster(distIn) ** 2) ** 3), 0)
                lacaBBOut.save("lacaBB_" + str(fid) + ".tif")

# Nitrogen here
        nOut = Con(distIn <= 7, (-3.9 * distIn ** 3 + 89.1 * distIn ** 2 - (814.5 * distIn) + 2968), 0)
        nOut.save("N_" + str(fid) + ".tif")


# Fantail Habitat here - include SPUs greater than 1.5 ha and within 150 m of another SPU of any size
        if(row[2] < 15000 and row[5] < 150.0) or row[2] >= 15000:
                # dcalcFT = 50.0
                lacaFTOut = Con(distIn <= dcalcFT, ((1/dcalcFT)*1.094*(1 - (1/dcalcFT)**2*Raster(distIn)**2)**3), 0)
                lacaFTOut.save("lacaFT_" + str(fid) + ".tif")

del row
del cursor

# Calculate final layer for each ES
# Cooling
rasterList = arcpy.ListRasters("cool*", "TIF")
outputC = outName + "_cool.tif"
proj = arcpy.SpatialReference(2193)
if len(rasterList) > 0:
        # outputC = outName + "_cool.tif"
        # proj = arcpy.SpatialReference(2193)
        arcpy.MosaicToNewRaster_management(rasterList, ws, "midcool.tif", proj, "32_BIT_FLOAT", cellSize, "1", "SUM")
        inConstant = 0.75
        outTimes = Times(Raster("midcool.tif"), inConstant)
        outTimes.save(outName + "_baseCooling.tif")
        # coolOverlap = 30 / (1 + Exp(4.365 - Raster("timescool.tif"))
        # user inputs values for ASYM, MID and k
        coolOverlap = asymC / (1 + Exp((midC - Raster(outTimes))/kC))
        coolOverlap.save(outputC)

# Nitrogen MS raster
rasters = arcpy.ListRasters("N_*", "TIF")
midN = outName + "_baseNitrogen.tif"
outputN = outName + "_nitrogen.tif"
proj = arcpy.SpatialReference(2193)
arcpy.MosaicToNewRaster_management(rasters, ws, midN, proj, "32_BIT_FLOAT", cellSize, "1", "SUM")
# Use EL model for nonlinear effects - parameters set by user
ELNitrogen = asymN / (1 + Exp((midN - Raster(midN))/kN))
ELNitrogen.save(outputN)


# Bellbird Habitat MS raster
rasters = arcpy.ListRasters("lacaBB*", "TIF")
midHBB = outName + "_baseBellbirdHabitat.tif"
outputHBB = outName + "_BellbirdHabitat.tif"
proj = arcpy.SpatialReference(2193)
if len(rasterList) > 0:
    arcpy.MosaicToNewRaster_management(rasters, ws, midHBB, proj, "32_BIT_FLOAT", cellSize, "1", "SUM")
# Use EL Model for nonlinear effects - user supplies parameters
ELHBB = asymBB / (1 + Exp((midBB - Raster(midHBB))/kBB))
ELHBB.save(outputHBB)


# Fantail  Habitat MS raster
rasters = arcpy.ListRasters("lacaFT*", "TIF")
midFT = outName + "_baseFantailHabitat.tif"
outputHFT = outName + "_FantailHabitat.tif"
proj = arcpy.SpatialReference(2193)
arcpy.MosaicToNewRaster_management(rasters, ws, midFT, proj, "32_BIT_FLOAT", cellSize, "1", "SUM")
ELHFT = asymFT / (1 + Exp((midFT - Raster(midHBB))/kFT))
ELHFT.save(outputHFT)

# Create Control rasters
# Nitrogen control - uses all SPUs
# arcpy.PolygonToRaster_conversion(inputFC, "Shape_Area", "NC.tif", "CELL_CENTER", "", cellSize)
# use Con to change NoData values to 0 and existing values to NoData
# nCon = Con((IsNull("NC.tif")), 0, Raster("NC.tif"))
# nCon.save("NControl.tif")

# Nitrogen control and ES calculation - uses all SPUs
# ncDist = arcpy.sa.EucDistance(inputFC, "", cellSize)
# ncCalc = Con(ncDist <= 0, (-3.9 * ncDist ** 3 + 89.1 * ncDist ** 2 - (814.5 * ncDist) + 2968), 0)
# ncCalc.save(outName + "_Ncontrol.tif")
#
# # Cooling and Bellbird Habitat control ES calculation - uses subset of SPUs
# whereClause1 = "Shape_Area >= 4900"
# arcpy.SelectLayerByAttribute_management(inputFC, "NEW_SELECTION", whereClause1)
# # ccDist = arcpy.sa.EucDistance(inputFC, "", cellSize)
# arcpy.PolygonToRaster_conversion(inputFC, "Shape_Area", "ccPoly.tif", "CELL_CENTER", "", cellSize)
#
# ccCalc = Con(IsNull(Raster("ccPoly.tif")), 0, (10 * Raster("ccPoly.tif")))
# ccCalc.save(outName + "_CoolControl.tif")
#
# bbCalc = Con(IsNull(Raster("ccPoly.tif")), 0, (1 / dcalcBB) * 1.094)
# bbCalc.save(outName + "_BBcontrol.tif")
#
# # Fantail Habitat control and ES calculation
# whereClause2 = "Shape_Area >= 15000 or (Shape_Area < 15000 and (NEAR_DIST <= 150 and NEAR_DIST > 0))"
# arcpy.SelectLayerByAttribute_management(inputFC, "NEW_SELECTION", whereClause2)
# ftDist = arcpy.sa.EucDistance(inputFC, "", cellSize)
# ftCalc = Con(ftDist <= 0, ((1/dcalcFT)*1.094*(1 - (1/dcalcFT)**2*Raster(ftDist)**2)**3), 0)
# ftCalc.save(outName + "_FTcontrol.tif")
#
# arcpy.SelectLayerByAttribute_management(inputFC, "CLEAR_SELECTION")

# Calculate Nitrogen and control metascore DBFs using mask
# arcpy.env.mask = mask
arcpy.sa.ZonalStatisticsAsTable(boundary, "OBJECTID", Raster(outputN), outName + "_MSNitrogen.dbf", "DATA", "SUM")
# arcpy.sa.ZonalStatisticsAsTable(boundary, "OBJECTID", Raster(outName + "_Ncontrol.tif"), outName + "_MSNControl.dbf", "DATA", "SUM")

# Calculate Cooling and Bellbird  control metascore DBFs
if arcpy.Exists(outputC):
    arcpy.sa.ZonalStatisticsAsTable(boundary, "OBJECTID", Raster(outputHBB), outName + "_MSBBHabitat.dbf", "DATA", "SUM")
    # arcpy.sa.ZonalStatisticsAsTable(boundary, "OBJECTID", Raster(outName + "_BBcontrol.tif"), outName + "_MSBBControl.dbf",
                                        # "DATA", "SUM")
    arcpy.sa.ZonalStatisticsAsTable(boundary, "OBJECTID", Raster(outputC), outName + "_MSCool.dbf", "DATA", "SUM")
    # arcpy.sa.ZonalStatisticsAsTable(boundary, "OBJECTID", Raster(outName + "_CoolControl.tif"), outName + "_MSCControl.dbf", "DATA",
                                    # "SUM")

# Calculate control metascores.  Skip C if it doesn't exist
# need a custom mask for each ES - Nitrogen can use the already existing mask from above
# arcpy.env.mask = hFTCM
arcpy.sa.ZonalStatisticsAsTable(boundary, "OBJECTID", Raster(outputHFT), outName + "_MSFTHabitat.dbf", "DATA", "SUM")
# arcpy.sa.ZonalStatisticsAsTable(boundary, "OBJECTID", Raster(outName + "_FTcontrol.tif"), outName + "_MSFTControl.dbf", "DATA", "SUM")

# Clean up
# Delete distance and individual ES grids
for ras in arcpy.ListRasters("*", "TIF"):
        if not ras.startswith(outName):
                arcpy.Delete_management(ras)


arcpy.CheckInExtension("Spatial")

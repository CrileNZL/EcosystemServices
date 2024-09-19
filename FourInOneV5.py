# Calculate ecosystem services layers and metascores for individual ES rasters
# Combines Cooling.py, Nitrogen.py, Laca.py and metascores.py into one script
# Developed for Richard Morris
# C. Doscher October 2022 - Updated 19 May 2023
# New version of 24 April script
# N calculation updated on 26 May 2023
# new BB and error catching for FT on 2 June 2023
# Added code to create masks for final analysis
# Updated 4 June - added error catching code for nonlinear mask creation - if max value = 0
# Updated July 2024 - new code to handle irregular polygons and cooling
# Updating September 2024 - new code to set minimum cooling distance internally

# Version 5 - Sept 2024 - replace nonlinear adjustment with Sigmoid function
# rescale final values to a 1 - 10 scale

import arcpy

import math

from arcpy.sa import *

arcpy.env.overwriteOutput = True

# Get Clump polygon layer from user
inputFC = arcpy.GetParameterAsText(0)

# Get centerline shapefile from user
cline = arcpy.GetParameterAsText(1)

# Get workspace from user
ws = arcpy.GetParameterAsText(2)
arcpy.env.workspace = ws

dcool = float(arcpy.GetParameterAsText(3))


# Get boundary zone FC from user
boundary = arcpy.GetParameterAsText(4)
arcpy.env.extent = boundary

# Get cell size from user
# get cellsize from user
cellSize = arcpy.GetParameterAsText(5)
# cellSizeHa = (float(cellSize)**2)/10000.0

# get output file name from user
outName = arcpy.GetParameterAsText(6)

dcalcBB = 25.0
dcalcFT = 100.0

arcpy.CheckOutExtension("Spatial")

# Create masks to limit area of application for nonlinear equation

listN = []
listFT = []
listBB = []

# Create buffers and raster grids for N, BB and FT to use in ES code later
with arcpy.da.SearchCursor(inputFC, ['FID', 'Shape@']) as cursor:
    for row in cursor:
        fid = row[0]
        arcpy.env.workspace = ws

        arcpy.AddMessage("Creating buffers")

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

# del row
# del cursor
arcpy.AddMessage("Finished buffers")

# Add up N grids and create masks
addN = CellStatistics(listN, "SUM", "NODATA")
addNras = Raster(addN)
maxAddN = int(addNras.maximum)
addN.save("addN.tif")

# create nonlinear adjustment layer for N
nonlinNadj = Con(Raster("addN.tif") > 0, 1 / (1 + Exp(-1 * Raster("addN.tif") + Exp(0.1))), 0)
nonlinNadj.save("masknonlinNadj.tif")

# output mask for N
if maxAddN <= 0:
    pass
else:
    maskN = Reclassify(Raster("addN.tif"), "VALUE", RemapRange([[0, 0, 0], [0, maxAddN, 1]]))
    maskN.save("maskN.tif")
    # output nonlinear mask for N
    nonlinmaskN = Reclassify(Raster("addN.tif"), "VALUE", RemapRange([[0, 1, 0], [1, maxAddN, 1]]))
    nonlinmaskN.save("maskNnonlin.tif")

arcpy.AddMessage("N masks done")

# Add up FT grids and create masks
addFT = CellStatistics(listFT, "SUM", "NODATA")
addFTras = Raster(addFT)
maxAddFT = int(addFTras.maximum)
addFT.save("addFT.tif")

# create nonlinear adjustment layer for FT
nonlinFTadj = Con(Raster("addFT.tif") > 0, 1 / (1 + Exp(-1 * Raster("addFT.tif") + Exp(0.1))), 0)
nonlinFTadj.save("masknonlinFTadj.tif")

# output mask for FT
if maxAddFT <= 0:
    pass
else:
    maskFT = Reclassify(Raster("addFT.tif"), "VALUE", RemapRange([[0, 0, 0], [0, maxAddFT, 1]]))
    maskFT.save("maskFT.tif")
    # output nonlinear mask for N
    nonlinmaskFT = Reclassify(Raster("addFT.tif"), "VALUE", RemapRange([[0, 1, 0], [1, maxAddFT, 1]]))
    nonlinmaskFT.save("maskFTnonlin.tif")

arcpy.AddMessage("Fantail masks done")

# Add up BB grids and create masks
addBB = CellStatistics(listBB, "SUM", "NODATA")
addBBras = Raster(addBB)
maxAddBB = int(addBBras.maximum)
addBB.save("addBB.tif")

# create nonlinear adjustment layer for BB
nonlinBBadj = Con(Raster("addBB.tif") > 0, 1 / (1 + Exp(-1 * Raster("addBB.tif") + Exp(0.1))), 0)
nonlinBBadj.save("masknonlinBBadj.tif")

# output mask for BB
if maxAddBB <= 0:
    pass
else:
    maskBB = Reclassify(Raster("addBB.tif"), "VALUE", RemapRange([[0, 0, 0], [0, maxAddBB, 1]]))
    maskBB.save("maskBB.tif")
    # output nonlinear mask for BB
    nonlinmaskBB = Reclassify(Raster("addBB.tif"), "VALUE", RemapRange([[0, 1, 0], [1, maxAddBB, 1]]))
    nonlinmaskBB.save("maskBBnonlin.tif")

arcpy.AddMessage("Bellbird masks done")

# clear lists for next iteration
listN.clear()
listFT.clear()
listBB.clear()

# delete all except final masks
for ras in arcpy.ListRasters("*", ""):
    if not ras.startswith("mask"):
        arcpy.Delete_management(ras)
# Delete buffer shapefiles
for shp in arcpy.ListFeatureClasses():
    if shp.startswith("buf"):
        arcpy.Delete_management(shp)

# Create mask to set SPU area to 0 for Zonal Statistics as Table
arcpy.PolygonToRaster_conversion(inputFC, "Shape_Area", "Mask.tif", "CELL_CENTER", "", cellSize)
# use Con to change NoData values to 0 and existing values to NoData
nCon = Con(IsNull("Mask.tif"), 1, 0)
nCon.save("SPUMask.tif")

# Run NEAR on inputFC for FT habitat calc
arcpy.Near_analysis(inputFC, inputFC)

# Create distance grids for use in ES code
with arcpy.da.SearchCursor(inputFC, ['FID', 'Shape@', 'Shape_Area', 'CC', 'd', 'NEAR_DIST']) as cursor:
    for row in cursor:

        # Set up loop
        # for each feature, derive distance grid - Euclidean Distance
        # distances grids used by Nitrogen, Bellbird and Fantail ESs only
        # Cooling uses custom distance grids
        # get FID here for iterating
        fid = row[0]
        distOut = "dist_" + str(fid) + ".tif"
        arcpy.env.workspace = ws
        # maxDist = '#' # can use ""?
        dist = arcpy.sa.EucDistance(row[1], "", cellSize)
        dist.save(distOut)

        # Calculate cooling raster
        distIn = Raster(distOut)

        # Calculate individual ESs
        # Cooling and Bellbird habitat here
        # select all SPUs greater than 4900 m2 in area and those smaller that are within 2*R of another SPU
        # new Cooling code here
        # For each SPU, create centerline and points around boundary
        # For each point, calculate distance to centerline
        # Use distance as Max Distance to calculate distance then calculate cooling effect and later mosaic to
        # new raster.
        arcpy.AddMessage("Starting cooling calculations")
        listCool = []
        if row[2] >= 4900.00 or (row[2] < 4900.00 and row[5] <= float(2 * ((row[2] / math.pi) ** 0.5))):

            # cc = float(row[3])
            cc = 10

            # new Cooling code here
            # polygon to line
            arcpy.AddMessage("Creating SPU Line")
            outLine = "SPULine_" + str(fid) + ".shp"
            arcpy.AddMessage("Creating points on SPU Line")
            arcpy.PolygonToLine_management(row[1], outLine, "")
            # autogenerate points
            bdyp = "BDYPoints_" + str(fid) + ".shp"
            arcpy.GeneratePointsAlongLines_management(outLine, bdyp, 'DISTANCE', Distance="5 meters", Include_End_Points='END_POINTS')
            # Near boundary points to centerline
            if int(arcpy.GetCount_management(bdyp)["row_count"]) > 0:
                arcpy.Near_analysis(bdyp, cline)
                minDistC = (max([cur[0] for cur in arcpy.da.SearchCursor(bdyp, "NEAR_DIST")]) / dcool)
                arcpy.AddMessage("minDistC = " + str(minDistC))

                # Loop through BDYPoints and do cooling calcs
                with arcpy.da.SearchCursor(bdyp, ['FID', 'SHAPE@', 'NEAR_DIST']) as cursorp:
                    for rowp in cursorp:
                        pfid = rowp[0]
                        if rowp[2] <= minDistC:
                            d = minDistC
                        else:
                            d = rowp[2]
                        distOutC = "cooldist_" + str(fid) + "_" + str(pfid) + ".tif"
                        dist = arcpy.sa.EucDistance(rowp[1], d, cellSize)
                        dist.save(distOutC)

                        # calculate cooling effect
                        distInC = Raster(distOutC)
                        ha = float(row[2])/10000

                        coolOut = ha * cc * Exp(-5 * (distInC / d))
                        coolOut.save("coolcalc_" + str(fid) + "_" + str(pfid) + ".tif")
                        listCool.append(Raster("coolcalc_" + str(fid) + "_" + str(pfid) + ".tif"))

            # del rowp
            # del cursorp

            # make cooling overlaps
            if len(listCool) > 0:
                overlapC = CellStatistics(listCool, "SUM", "NODATA")
                overlapCreclass = Con(IsNull(Raster(overlapC)), 0, 1)
                overlapCreclass.save("OverlapC_" + str(fid) + ".tif")

            listCool.clear()

            arcpy.AddMessage("Cooling calcs done.")

        if row[2] >= 15000 or (row[2] > 950 and row[5] < 10.0):
            lacaBBOut = Con(distIn <= dcalcBB,
                            ((1 / dcalcBB) * 1.094 * (1 - (1 / dcalcBB) ** 2 * Raster(distIn) ** 2) ** 3), 0)
            lacaBBOut.save("lacaBB_" + str(fid) + ".tif")

        arcpy.AddMessage("Bellbird calcs done")

# Nitrogen here
        # nOut = Con(distIn <= 7, (-3.9 * distIn ** 3 + 89.1 * distIn ** 2 - (814.5 * distIn) + 2968), 0)
        # Updated by RM 26 May 23
        nOut = Con(distIn <= 7, (-235 * distIn + 2065), 0)
        nOut.save("N_" + str(fid) + ".tif")


# Fantail Habitat here - include SPUs greater than 1.5 ha and within 150 m of another SPU of any size
        lacaFTOut = Con(distIn <= dcalcFT, ((1/dcalcFT)*1.094*(1 - (1/dcalcFT)**2*Raster(distIn)**2)**3), 0)
        lacaFTOut.save("lacaFT_" + str(fid) + ".tif")

# del row
# del cursor

# Calculate final layer for each ES
# Cooling
rasterC = arcpy.ListRasters("coolcalc*", "")
outputC = outName + "_cool" + ".tif"
proj = arcpy.SpatialReference(2193)
if len(rasterC) > 0:

    arcpy.MosaicToNewRaster_management(rasterC, ws, "midcool.tif", proj, "32_BIT_FLOAT", cellSize, "1", "SUM")

    # Bring cooling overlaps together
    overlaps = arcpy.ListRasters("OverlapC*", "")
    finalCoverlaps = CellStatistics(overlaps, "SUM", "NODATA")
    finalCoverlaps.save("AllCoverlaps.tif")

    # nonlinear cooling adjustment here
    nonlinCadj = Con(Raster("AllCoverlaps.tif") > 0, 1 / (1 + Exp(-1 * Raster("AllCoverlaps.tif") + Exp(0.1))), 0)
    nonlinCadj.save("masknonlinCadj.tif")

    nlC = Raster("midcool.tif") + (Raster("midcool.tif") * Raster("masknonlinCadj.tif"))
    nlC.save("Cadjusted.tif")

    # Rescale N to 1 - 10 scale
    nlCrescale = ((9 * Raster("Cadjusted.tif")) / (Raster("Cadjusted.tif").maximum - Raster("Cadjusted.tif").minimum)) + \
                 (10 - ((9 * Raster("Cadjusted.tif").maximum) / (Raster("Cadjusted.tif").maximum - Raster("Cadjusted.tif").minimum)))
    finalCool = Times(nlCrescale, Raster("SPUMask.tif"))
    finalCool.save(outputC)

    arcpy.AddMessage("Final Cooling layer done.")

# Nitrogen MS raster
rastersN = arcpy.ListRasters("N_*", "")
outputN = outName + "_nitrogen.tif"
proj = arcpy.SpatialReference(2193)
if len(rastersN) > 0:
    arcpy.MosaicToNewRaster_management(rastersN, ws, "tempNit.tif", proj, "32_BIT_FLOAT", cellSize, "1", "SUM")

    # nonlinear adjustment of Nitrogen raw values
    if arcpy.Exists("maskNnonlin.tif"):
        nlN = Raster("tempNit.tif") + (Raster("tempNit.tif") * Raster("masknonlinNadj.tif"))
        nlN.save("Nadjusted.tif")

        # Rescale N to 1 - 10 scale
        nlNrescale = ((9 * Raster("Nadjusted.tif")) / (Raster("Nadjusted.tif").maximum - Raster("Nadjusted.tif").minimum)) + \
                     (10 - ((9 * Raster("Nadjusted.tif").maximum) / (Raster("Nadjusted.tif").maximum - Raster("Nadjusted.tif").minimum)))
        finalNnl = Times(nlNrescale, Raster("SPUMask.tif"))
        finalNnl.save(outputN)

        arcpy.AddMessage("Final Nitrogen calcs done (with NL adjustment.")
    else:
        Nrescale = ((9 * Raster("tempNit.tif")) / (Raster("tempNit.tif").maximum - Raster("tempNit.tif").minimum)) + \
                   (10 - ((9 * Raster("tempNit.tif").maximum) / (Raster("tempNit.tif").maximum - Raster("tempNit.tif").minimum)))
        finalN = Times(Nrescale, Raster("SPUMask.tif"))
        finalN.save(outputN)

        arcpy.AddMessage("Final Nitrogen calcs done.")

# Bellbird Habitat MS raster
rastersBB = arcpy.ListRasters("lacaBB*", "")
stdHBBras = outName + "_baseBellbirdHabitat.tif"
outputHBB = outName + "_BellbirdHabitat.tif"
proj = arcpy.SpatialReference(2193)
if len(rastersBB) > 0:
    arcpy.MosaicToNewRaster_management(rastersBB, ws, "tempBB.tif", proj, "32_BIT_FLOAT", cellSize, "1", "SUM")

    if arcpy.Exists("maskBBnonlin.tif"):
        # nolinear BB adjustment
        nlBB = Raster("tempBB.tif") + (Raster("tempBB.tif") * Raster("masknonlinBBadj.tif"))
        nlBB.save("BBadjusted.tif")

        # Rescale N to 1 - 10 scale
        nlBBrescale = ((9 * Raster("BBadjusted.tif")) / (Raster("BBadjusted.tif").maximum - Raster("BBadjusted.tif").minimum)) + \
                      (10 - ((9 * Raster("BBadjusted.tif").maximum) / (Raster("BBadjusted.tif").maximum - Raster("BBadjusted.tif").minimum)))
        finalBBnl = Times(nlBBrescale, Raster("SPUMask.tif"))
        nlBBrescale.save(outputHBB)

        arcpy.AddMessage("Final Bellbird calcs done (with NL adjustment.")
    else:
        BBrescale = ((9 * Raster("tempBB.tif")) / (Raster("tempBB.tif").maximum - Raster("tempBB.tif").minimum)) + \
                    (10 - ((9 * Raster("tempBB.tif").maximum) / (Raster("tempBB.tif").maximum - Raster("tempBB.tif").minimum)))
        finalBB = Times(BBrescale, Raster("SPUMask.tif"))
        finalBB.save(outputHBB)

        arcpy.AddMessage("Final Bellbird calcs done.")


# Fantail Habitat MS raster
rastersFT = arcpy.ListRasters("lacaFT*", "")
midFTras = outName + "_baseFantailHabitat.tif"
outputHFT = outName + "_FantailHabitat.tif"
proj = arcpy.SpatialReference(2193)
if len(rastersFT) > 0:
    arcpy.MosaicToNewRaster_management(rastersFT, ws, "tempFT.tif", proj, "32_BIT_FLOAT", cellSize, "1", "SUM")

    if arcpy.Exists("maskFTnonlin.tif"):
        # nolinear FT adjustment
        nlFT = Raster("tempFT.tif") + (Raster("tempFT.tif") * Raster("masknonlinFTadj.tif"))
        nlFT.save("FTadjusted.tif")

        # Rescale N to 1 - 10 scale
        nlFTrescale = ((9 * Raster("FTadjusted.tif")) / (Raster("FTadjusted.tif").maximum - Raster("FTadjusted.tif").minimum)) + \
                      (10 - ((9 * Raster("FTadjusted.tif").maximum) / (Raster("FTadjusted.tif").maximum - Raster("FTadjusted.tif").minimum)))
        finalFTnl = Times(nlFTrescale, Raster("SPUMask.tif"))
        finalFTnl.save(outputHFT)

        arcpy.AddMessage("Final Fantail calcs done (with NL adjustment).")
    else:
        FTrescale = ((9 * Raster("tempFT.tif")) / (Raster("tempFT.tif").maximum - Raster("tempFT.tif").minimum)) + \
                    (10 - ((9 * Raster("tempFT.tif").maximum) / (Raster("tempFT.tif").maximum - Raster("tempFT.tif").minimum)))
        finalFT = Times(FTrescale, Raster("SPUMask.tif"))
        finalFT.save(outputHFT)
        # Raster("tempFT.tif").save(outputHFT)
        arcpy.AddMessage("Final Fantail calcs done.")

arcpy.AddMessage("Starting ES metascores.")

# use Zonal Statistics as Table to calculate ES metascore
arcpy.sa.ZonalStatisticsAsTable(boundary, "FID", Raster(outputN), outName + "_MSNitrogen.dbf", "DATA", "SUM")

if arcpy.Exists(outputHBB):
    arcpy.sa.ZonalStatisticsAsTable(boundary, "FID", Raster(outputHBB), outName + "_MSBBHabitat.dbf", "DATA", "SUM")

if arcpy.Exists(outputC):
    arcpy.sa.ZonalStatisticsAsTable(boundary, "FID", Raster(outputC), outName + "_MSCool.dbf", "DATA", "SUM")

# Calculate control metascores.  Skip C if it doesn't exist
# need a custom mask for each ES - Nitrogen can use the already existing mask from above
if arcpy.Exists(outputHFT):
    arcpy.sa.ZonalStatisticsAsTable(boundary, "FID", Raster(outputHFT), outName + "_MSFTHabitat.dbf", "DATA", "SUM")

# Clean up
# Delete distance and individual ES grids
for ras in arcpy.ListRasters("*", ""):
    if not ras.startswith(outName):
        arcpy.Delete_management(ras)

arcpy.AddMessage("Script complete.")


arcpy.CheckInExtension("Spatial")

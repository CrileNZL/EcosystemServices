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

# March 2025 - adding Select by ATtribute so each SPU's correct centerline gets selected before running Near

# October 2025 - adapting script for vineyard application

import arcpy

import math

from arcpy.sa import *

arcpy.env.overwriteOutput = True

# Get Clump polygon layer from user
inputFC = arcpy.GetParameterAsText(0)
proj = arcpy.Describe(inputFC).spatialReference

# Get centerline shapefile from user
cline = arcpy.GetParameterAsText(1)

# Get workspace from user
ws = arcpy.GetParameterAsText(2)
arcpy.env.workspace = ws

# get minimum distance divisor from user - used on line 236
# dcool = float(arcpy.GetParameterAsText(3))
dcool = 10

# get CC value from user - used on line 255
# cc = float(arcpy.GetParameterAsText(4))
cc = 10

# Get boundary zone FC from user
boundary = arcpy.GetParameterAsText(3)
arcpy.env.extent = boundary

# Get cell size from user
# get cellsize from user
cellSize = arcpy.GetParameterAsText(4)
# cellSizeHa = (float(cellSize)**2)/10000.0

# get output file name from user
outName = arcpy.GetParameterAsText(5)

dcalcBB = 500.0
dcalcFT = 100.0

az = int(arcpy.GetParameterAsText(6))

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

# del row
# del cursor
arcpy.AddMessage("Finished N and FT buffers")

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
with arcpy.da.SearchCursor(inputFC, ['FID', 'Shape@', 'Shape_Area', 'CC', 'd', 'NEAR_DIST', 'SPUDist', 'MinWidth']) as cursor:
    for row in cursor:

        fid = row[0]

        distOut = "dist_" + str(fid) + ".tif"
        arcpy.env.workspace = ws
        # maxDist = '#' # can use ""?
        dist = arcpy.sa.EucDistance(row[1], "", cellSize)
        dist.save(distOut)

        # Calculate cooling raster
        distIn = Raster(distOut)

        listCool = []
        listBB = []

        # check if right size/proximity for cooling
        if row[2] >= 35.00 or (row[2] < 35.00 and row[5] <= row[6]):
            print("SPU FID " + str(fid) + " is okay for cooling")
            # cc = float(row[3]) - now set by user

            # new Cooling code here
            # polygon to line
            # If the right conditions, make points
            arcpy.AddMessage("Creating SPU Line")
            # print("Creating SPU Line")
            outLine = "SPULine_" + str(fid) + ".shp"
            arcpy.AddMessage("Creating points on SPU Line")
            # print("Creaing points on SPU Line")
            arcpy.PolygonToLine_management(row[1], outLine, "")
            # autogenerate points
            bdyp = "BDYPoints_" + str(fid) + ".shp"
            arcpy.GeneratePointsAlongLines_management(outLine, bdyp, 'DISTANCE', Distance="5 meters",
                                                      Include_End_Points='END_POINTS')
            # Near boundary points to centerline
            if int(arcpy.GetCount_management(bdyp)["row_count"]) > 0:
                # select centerline for that SPU
                arcpy.management.SelectLayerByAttribute(cline, "NEW_SELECTION", '"FID" =' + str(fid))
                arcpy.Near_analysis(bdyp, cline)
                arcpy.management.SelectLayerByAttribute(cline, "CLEAR_SELECTION")
                minDistC = (max([cur[0] for cur in arcpy.da.SearchCursor(bdyp, "NEAR_DIST")]) / dcool)
                arcpy.AddMessage("minDistC = " + str(minDistC))
                # print("minDistC = " + str(minDistC))

                # find max NEAR_DIST to centerline
                maxWidth = 2 * max([cur[0] for cur in arcpy.da.SearchCursor(bdyp, "NEAR_DIST")])
                print("maxWidth = " + str(maxWidth))

                # if maxWidth > MinWidth, loop through BDYPoints and do cooling calcs
                if maxWidth >= row[7]:
                    with arcpy.da.SearchCursor(bdyp, ['FID', 'SHAPE@', 'NEAR_DIST']) as cursorp:
                        for rowp in cursorp:
                            pfid = rowp[0]
                            if rowp[2] <= minDistC:
                                d = minDistC
                            else:
                                d = rowp[2]
                            # if 2 * max.NEAR_DIST > 5, create disance grid
                            distOutC = "cooldist_" + str(fid) + "_" + str(pfid) + ".tif"
                            dist = arcpy.sa.EucDistance(rowp[1], "", cellSize)
                            dist.save(distOutC)

                            # calculate cooling effect
                            distInC = Raster(distOutC)
                            ha = float(row[2]) / 10000

                            print("Calculating " + str(pfid) + " cooling layer")
                            if 2 * rowp[2] >= 1:
                                # calc cooling using CON
                                coolOut = Con(distInC < d, ha * cc * Exp(-5 * (distInC / d)), 0)
                                coolOut.save("coolcalc_" + str(fid) + "_" + str(pfid) + ".tif")
                                # add outputs to list
                                listCool.append(Raster("coolcalc_" + str(fid) + "_" + str(pfid) + ".tif"))
                else:
                    print("SPU FID " + str(fid) + " isn't okay for cooling")
        arcpy.AddMessage("Cooling calcs done")
        # print("Cooling steps done, checking BB calcs (print)")

        SPUoverlaps = arcpy.ListRasters("coolcalc_" + str(fid) + "*")
        # proj = arcpy.SpatialReference(2193)
        if len(SPUoverlaps) > 0:
            SPUoverlapAdd = arcpy.MosaicToNewRaster_management(SPUoverlaps, ws, "thisOverlapC_" + str(fid) + ".tif", proj, "32_BIT_FLOAT", cellSize, "1", "SUM")
            SPUoverlap = Con(Raster("thisOverlapC_" + str(fid) + ".tif") > 0, 1, 0)
            SPUoverlap.save("OverlapC_" + str(fid) + ".tif")

        # check if right size/proximity for Bellbird Habitat
        if row[2] >= 15000 or (row[2] > 950 and row[5] < 10.0):
            # if right conditions - check if bdy points layer exists
            if arcpy.Exists(bdyp):
                if maxWidth >= 25:
                    with arcpy.da.SearchCursor(bdyp, ['FID', 'SHAPE@', 'NEAR_DIST']) as cursorbb:
                        for rowbb in cursorbb:
                            pfid = rowbb[0]
                            if 2 * rowbb[2] >= 25:
                                if arcpy.Exists("cooldist_" + str(fid) + "_" + str(pfid) + ".tif"):
                                    print("SPU FID " + str(fid) + "_" + str(pfid) + " is all good for BB calc")

                                    # if yes, calc BB using CON
                                    bbcalc = "bbcalc_" + str(fid) + "_" + str(pfid) + ".tif"
                                    lacaBBOut = Con(Raster("cooldist_" + str(fid) + "_" + str(pfid) + ".tif") <= dcalcBB,
                                                    ((1 / dcalcBB) * 1.094 * (1 - (1 / dcalcBB) ** 2 * Raster("cooldist_" + str(fid) + "_" + str(pfid) + ".tif") ** 2) ** 3), 0)
                                    lacaBBOut.save(bbcalc)
                                    listBB.append(bbcalc)
                                    print("Cooldist exists - BB calculated")
                                # open search cursor to iterate through each point
                                else:

                                    distOutBB = "bbdist_" + str(fid) + "_" + str(pfid) + ".tif"
                                    dist = arcpy.sa.EucDistance(rowp[1], "", cellSize)
                                    dist.save(distOutBB)

                                    # calc BB using CON
                                    lacaBBOut = Con(dist <= dcalcBB, ((1 / dcalcBB) * 1.094 * (1 - (1 / dcalcBB) ** 2 * Raster("bbdist_" + str(fid) + "_" + str(pfid) + ".tif") ** 2) ** 3), 0)
                                    lacaBBOut.save(bbcalc)
                                    listBB.append(bbcalc)

                                    # print("DD calc with existing layers - cooldist didn't exist so I made one")
                    # cursor.reset()

            else:
                # If the right Bellbird conditions but no BDYPoints, make points
                # print("No cooldist exists - creating now")
                arcpy.AddMessage("No cooldist exists - creating BB SPU Line")
                outLine = "SPULine_" + str(fid) + ".shp"
                arcpy.AddMessage("Creating points on BB SPU Line")
                arcpy.PolygonToLine_management(row[1], outLine, "")
                # autogenerate points
                bdyp = "BDYPoints_" + str(fid) + ".shp"
                arcpy.GeneratePointsAlongLines_management(outLine, bdyp, 'DISTANCE', Distance="5 meters",
                                                          Include_End_Points='END_POINTS')
                # Near boundary points to centerline
                if int(arcpy.GetCount_management(bdyp)["row_count"]) > 0:
                    arcpy.Near_analysis(bdyp, cline)
                    minDistBB = (max([cur[0] for cur in arcpy.da.SearchCursor(bdyp, "NEAR_DIST")]) / dcool)
                    arcpy.AddMessage("minDistBB = " + str(minDistBB))
                    # print("minDistBB = " + str(minDistBB))

                    # find max NEAR_DIST to centerline
                    maxWidth = 2 * max([cur[0] for cur in arcpy.da.SearchCursor(bdyp, "NEAR_DIST")])

                    # if maxWidth > 25, loop through BDYPoints and do cooling calcs
                    if maxWidth > 25:
                        with arcpy.da.SearchCursor(bdyp, ['FID', 'SHAPE@', 'NEAR_DIST']) as cursorp:
                            for rowp in cursorp:
                                pfid = rowp[0]

                                # if 2 * max.NEAR_DIST > 5, create disance grid
                                if 2 * rowp[2] >= 25:
                                    distOutBB = "bbdist_" + str(fid) + "_" + str(pfid) + ".tif"
                                    dist = arcpy.sa.EucDistance(rowp[1], "", cellSize)
                                    dist.save(distOutBB)

                                    # calc BB using CON
                                    lacaBBOut = Con(dist <= dcalcBB,
                                                    ((1 / dcalcBB) * 1.094 * (
                                                                1 - (1 / dcalcBB) ** 2 * Raster("bbdist_" + str(fid) + "_" + str(pfid) + ".tif") ** 2) ** 3), 0)
                                    lacaBBOut.save(bbcalc)
                                    listBB.append(bbcalc)
        else:
            print("SPU FID " + str(fid) + " isn't right for BB calc")
        arcpy.AddMessage("Bellbird calcs done")
        # print("Bellbird steps done")

        SPUoverlapsBB = arcpy.ListRasters("bbcalc_" + str(fid) + "*")
        # proj = arcpy.SpatialReference(2193)
        if len(SPUoverlapsBB) > 0:
            SPUoverlapBBAdd = arcpy.MosaicToNewRaster_management(SPUoverlaps, ws, "thisOverlapBB_" + str(fid) + ".tif", proj, "32_BIT_FLOAT", cellSize, "1", "SUM")
            SPUoverlapBB = Con(Raster("thisOverlapBB_" + str(fid) + ".tif") > 0, 1, 0)
            SPUoverlapBB.save("OverlapBB_" + str(fid) + ".tif")

# Nitrogen here
        # nOut = Con(distIn <= 7, (-3.9 * distIn ** 3 + 89.1 * distIn ** 2 - (814.5 * distIn) + 2968), 0)
        # Updated by RM 26 May 23
        # nOut = Con(distIn <= 7, (-235 * distIn + 2065), 0)
        # nOut.save("N_" + str(fid) + ".tif")

# Exit loop

# Nitrogen here
# calc distance grid from SPUs - used for N and for FT
nDist = arcpy.sa.EucDistance(inputFC, "", cellSize)
nDist.save("nDist.tif")

# calc N = 2065 - 235 * distIn within 7 m of SPU
nOut = Con(nDist <= 7, (-235 * nDist + 2065), 0)
nOut.save("NStart.tif")

# Hillshade, azimuth = 315, altitude = 45, no shadows
nHS = Hillshade(nOut, az, 45, "SHADOWS", 0.01)
nHS.save("nHS.tif")

# Reclassify: 0 = 1, > 0 = NODATA
hsReclass = Reclassify(Raster("nHS.tif"), "Value", RemapValue([[0,0,1]]), "NODATA")
hsReclass.save("hsReclass.tif")

# HSReclass * SPUDistance, get max distance
hsDist = hsReclass * nDist
hsDist.save("hsDist.tif")
nDistMax = Raster("hsDist.tif").maximum

# Calc ESn as 2065 - 1645/maxDist
nHSOut = Con(Raster("hsDist.tif") <= nDistMax, (2065 - (1645/nDistMax)*Raster("hsDist.tif")), 0)
nHSOut.save("nHSOut.tif")

# Combine using Con
nOut = Con(IsNull(Raster("nHSOut.tif")), "NStart.tif", "nHSOut.tif")
nOut.save("nOut.tif")

# Fantail Habitat here - include SPUs greater than 1.5 ha and within 150 m of another SPU of any size
tempFT = Con(nDist <= dcalcFT, ((1/dcalcFT)*1.094*(1 - (1/dcalcFT)**2*Raster(nDist)**2)**3), 0)
tempFT.save("tempFT.tif")

# Calculate final layer for each ES
# Cooling
# Combine all cooling rasters together by adding them all up with Mosaic to New Raster
# First add all cooling rasters to a list
rasterC = arcpy.ListRasters("coolcalc*", "")
outputC = outName + "_cool" + ".tif"
# proj = arcpy.SpatialReference(2193)
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
# rastersN = arcpy.ListRasters("N_*", "")
outputN = outName + "_nitrogen.tif"
# proj = arcpy.SpatialReference(2193)
# if len(rastersN) > 0:
#     arcpy.MosaicToNewRaster_management(rastersN, ws, "tempNit.tif", proj, "32_BIT_FLOAT", cellSize, "1", "SUM")

# nonlinear adjustment of Nitrogen raw values
if arcpy.Exists("maskNnonlin.tif"):
    nlN = Raster("nOut.tif") + (Raster("nOut.tif") * Raster("masknonlinNadj.tif"))
    nlN.save("Nadjusted.tif")

    # Rescale N to 1 - 10 scale
    nlNrescale = ((9 * Raster("Nadjusted.tif")) / (Raster("Nadjusted.tif").maximum - Raster("Nadjusted.tif").minimum)) + \
                 (10 - ((9 * Raster("Nadjusted.tif").maximum) / (Raster("Nadjusted.tif").maximum - Raster("Nadjusted.tif").minimum)))
    finalNnl = Times(nlNrescale, Raster("SPUMask.tif"))
    finalNnl.save(outputN)

    arcpy.AddMessage("Final Nitrogen calcs done (with NL adjustment.")
else:
    Nrescale = ((9 * Raster("nOut.tif")) / (Raster("nOut.tif").maximum - Raster("nOut.tif").minimum)) + \
               (10 - ((9 * Raster("nOut.tif").maximum) / (Raster("nOut.tif").maximum - Raster("nOut.tif").minimum)))
    finalN = Times(Nrescale, Raster("SPUMask.tif"))
    finalN.save(outputN)

    arcpy.AddMessage("Final Nitrogen calcs done.")

# Bellbird Habitat MS raster
rastersBB = arcpy.ListRasters("thisOverlapBB*", "")
stdHBBras = outName + "_baseBellbirdHabitat.tif"
outputHBB = outName + "_BellbirdHabitat.tif"
# proj = arcpy.SpatialReference(2193)
if len(rastersBB) > 0:
    arcpy.MosaicToNewRaster_management(rastersBB, ws, "midBB.tif", proj, "32_BIT_FLOAT", cellSize, "1", "SUM")

    # Bring Bellbird overlaps together
    overlapsBB = arcpy.ListRasters("OverlapBB*", "")
    finalCoverlaps = CellStatistics(overlapsBB, "SUM", "NODATA")
    finalCoverlaps.save("AllBBoverlaps.tif")

    # nonlinear cooling adjustment here
    nonlinBBadj = Con(Raster("AllBBoverlaps.tif") > 0, 1 / (1 + Exp(-1 * Raster("AllBBoverlaps.tif") + Exp(0.1))), 0)
    nonlinBBadj.save("masknonlinBBadj.tif")

    nlC = Raster("midBB.tif") + (Raster("midBB.tif") * Raster("masknonlinBBadj.tif"))
    nlC.save("BBadjusted.tif")

    if arcpy.Exists("maskBBnonlin.tif"):
        # nolinear BB adjustment
        nlBB = Raster("midBB.tif") + (Raster("midBB.tif") * Raster("masknonlinBBadj.tif"))
        nlBB.save("BBadjusted.tif")

        # Rescale N to 1 - 10 scale
        nlBBrescale = ((9 * Raster("BBadjusted.tif")) / (Raster("BBadjusted.tif").maximum - Raster("BBadjusted.tif").minimum)) + \
                      (10 - ((9 * Raster("BBadjusted.tif").maximum) / (Raster("BBadjusted.tif").maximum - Raster("BBadjusted.tif").minimum)))
        finalBBnl = Times(nlBBrescale, Raster("SPUMask.tif"))
        nlBBrescale.save(outputHBB)

        arcpy.AddMessage("Final Bellbird calcs done (with NL adjustment.")
    else:
        BBrescale = ((9 * Raster("midBB.tif")) / (Raster("midBB.tif").maximum - Raster("midBB.tif").minimum)) + \
                    (10 - ((9 * Raster("midBB.tif").maximum) / (Raster("midBB.tif").maximum - Raster("midBB.tif").minimum)))
        finalBB = Times(BBrescale, Raster("SPUMask.tif"))
        finalBB.save(outputHBB)

        arcpy.AddMessage("Final Bellbird calcs done.")


# Fantail Habitat MS raster
# rastersFT = arcpy.ListRasters("lacaFT*", "")
# midFTras = outName + "_baseFantailHabitat.tif"
outputHFT = outName + "_FantailHabitat.tif"
# proj = arcpy.SpatialReference(2193)
# if len(rastersFT) > 0:
#     arcpy.MosaicToNewRaster_management(rastersFT, ws, "tempFT.tif", proj, "32_BIT_FLOAT", cellSize, "1", "SUM")

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
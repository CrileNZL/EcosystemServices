import arcpy
# from arcpy.ia import RasterCalculator
import math
from arcpy.sa import *

# arcpy.env.workspace = r"D:\OLW Nitrates\OLWNitrates.gdb" # memory?
arcpy.env.overwriteOutput = True

ws = arcpy.env.workspace = r"E:\IRPs\Test"

inputFC = "Polys.shp"

cline = "cline.shp"

bdy = "BDY.shp"
arcpy.env.extent = bdy

cellSize = 5
cc = 10
dcalcBB = 500

# Run NEAR on inputFC for FT habitat calc
arcpy.Near_analysis(inputFC, inputFC)

# Set up search cursor on SPUs
with arcpy.da.SearchCursor(inputFC, ['FID', 'Shape@', 'Shape_Area', 'CC', 'd', 'NEAR_DIST']) as cursor:
    for row in cursor:

        fid = row[0]
        dcool = 5

        listCool = []
        listBB = []

        # check if right size/proximity for cooling
        if row[2] >= 4900.00 or (row[2] < 4900.00 and row[5] <= float(2 * ((row[2] / math.pi) ** 0.5))):
            print("SPU FID " + str(fid) + " is okay for cooling")
            # cc = float(row[3]) - now set by user

            # new Cooling code here
            # polygon to line
            # If the right conditions, make points
            # arcpy.AddMessage("Creating SPU Line")
            print("Creating SPU Line")
            outLine = "SPULine_" + str(fid) + ".shp"
            # arcpy.AddMessage("Creating points on SPU Line")
            print("Creaing points on SPU Line")
            arcpy.PolygonToLine_management(row[1], outLine, "")
            # autogenerate points
            bdyp = "BDYPoints_" + str(fid) + ".shp"
            arcpy.GeneratePointsAlongLines_management(outLine, bdyp, 'DISTANCE', Distance="5 meters",
                                                      Include_End_Points='END_POINTS')
            # Near boundary points to centerline
            if int(arcpy.GetCount_management(bdyp)["row_count"]) > 0:
                arcpy.Near_analysis(bdyp, cline)
                minDistC = (max([cur[0] for cur in arcpy.da.SearchCursor(bdyp, "NEAR_DIST")]) / dcool)
                # arcpy.AddMessage("minDistC = " + str(minDistC))
                print("minDistC = " + str(minDistC))

                # find max NEAR_DIST to centerline
                maxWidth = 2 * max([cur[0] for cur in arcpy.da.SearchCursor(bdyp, "NEAR_DIST")])
                print("maxWidth = " + str(maxWidth))

                # if maxWidth . 5, loop through BDYPoints and do cooling calcs
                if maxWidth >= 5:
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
                            if d >= 5:
                                # calc cooling using CON
                                coolOut = Con(distInC < d, ha * cc * Exp(-5 * (distInC / d)), 0)
                                coolOut.save("coolcalc_" + str(fid) + "_" + str(pfid) + ".tif")
                                # add outputs to list
                                listCool.append(Raster("coolcalc_" + str(fid) + "_" + str(pfid) + ".tif"))
                else:
                    print("SPU FID " + str(fid) + " isn't okay for cooling")
        # arcpy.AddMessage("Cooling calcs done, starting BB calcs")
        print("Cooling steps done, checking BB calcs (print)")

        # check if right size/proximity for BB
        if row[2] >= 15000 or (row[2] > 950 and row[5] < 10.0):
            # if right conditions - check if bdy points layer exists
            if arcpy.Exists(bdyp):
                if maxWidth >= 25:
                    with arcpy.da.SearchCursor(bdyp, ['FID', 'SHAPE@', 'NEAR_DIST']) as cursorbb:
                        for rowbb in enumerate(cursorbb, start=1):
                            pfid = rowbb[0]
                            if rowbb[2] >= 25:
                                if arcpy.Exists("cooldist_" + str(fid) + "_" + str(pfid) + ".tif"):
                                    print("SPU FID " + str(fid) + "_" + str(pfid) + " is all good for BB calc")

                                    # if yes, calc BB using CON
                                    bbcalc = "bbcalc_" + str(fid) + "_" + str(pfid) + ".tif"
                                    lacaBBOut = Con(Raster("cooldist_" + str(fid) + "_" + str(pfid) + ".tif") <= dcalcBB,
                                                    ((1 / dcalcBB) * 1.094 * (1 - (1 / dcalcBB) ** 2 * Raster(dist) ** 2) ** 3), 0)
                                    lacaBBOut.save(bbcalc)
                                    listBB.append(bbcalc)
                                    print("Cooldist exists - BB calculated")
                                # open search cursor to iterate through each point
                                else:
                                    with arcpy.da.SearchCursor(bdyp, ['FID', 'SHAPE@', 'NEAR_DIST']) as cursorp2:
                                        print("Inside BB search cursor")
                                        for rowp in enumerate(cursorp2, start=1):
                                            pfid = rowp[0]
                                            if rowp[2] <= minDistC:
                                                d = minDistC
                                            else:
                                                d = rowp[2]
                                            print("d = " + str(d))
                                            # if 2 * max.NEAR_DIST > 5, create disance grid
                                            distOutBB = "bbdist_" + str(fid) + "_" + str(pfid) + ".tif"
                                            dist = arcpy.sa.EucDistance(rowp[1], "", cellSize)
                                            dist.save(distOutBB)

                                            # calc BB using CON
                                            lacaBBOut = Con(dist <= dcalcBB,
                                                            ((1 / dcalcBB) * 1.094 * (
                                                                    1 - (1 / dcalcBB) ** 2 * Raster(dist) ** 2) ** 3), 0)
                                            lacaBBOut.save(bbcalc)
                                            listBB.append(bbcalc)

                                            print("DD calc with existing layers - cooldist didn't exist so I made one")
                                        cursor.reset()

            else:
                # If the right conditions, make points
                print("No cooldist exists - creating now")
                # arcpy.AddMessage("Creating BB SPU Line")
                outLine = "SPULine_" + str(fid) + ".shp"
                # arcpy.AddMessage("Creating points on BB SPU Line")
                arcpy.PolygonToLine_management(row[1], outLine, "")
                # autogenerate points
                bdyp = "BDYPoints_" + str(fid) + ".shp"
                arcpy.GeneratePointsAlongLines_management(outLine, bdyp, 'DISTANCE', Distance="5 meters",
                                                          Include_End_Points='END_POINTS')
                # Near boundary points to centerline
                if int(arcpy.GetCount_management(bdyp)["row_count"]) > 0:
                    arcpy.Near_analysis(bdyp, cline)
                    minDistBB = (max([cur[0] for cur in arcpy.da.SearchCursor(bdyp, "NEAR_DIST")]) / dcool)
                    # arcpy.AddMessage("minDistBB = " + str(minDistBB))
                    print("minDistBB = " + str(minDistBB))

                    # find max NEAR_DIST to centerline
                    maxWidth = 2 * max([cur[0] for cur in arcpy.da.SearchCursor(bdyp, "NEAR_DIST")])

                    # if maxWidth > 25, loop through BDYPoints and do cooling calcs
                    if maxWidth > 25:
                        with arcpy.da.SearchCursor(bdyp, ['FID', 'SHAPE@', 'NEAR_DIST']) as cursorp:
                            for rowp in enumerate(cursorp, start=1):
                                pfid = rowp[0]
                                if rowp[2] <= minDistC:
                                    d = minDistC
                                else:
                                    d = rowp[2]
                                # if 2 * max.NEAR_DIST > 5, create disance grid
                                if d >= 25:
                                    distOutBB = "bbdist_" + str(fid) + "_" + str(pfid) + ".tif"
                                    dist = arcpy.sa.EucDistance(rowp[1], "", cellSize)
                                    dist.save(distOutBB)

                                    # calc BB using CON
                                    lacaBBOut = Con(dist <= dcalcBB,
                                                    ((1 / dcalcBB) * 1.094 * (
                                                                1 - (1 / dcalcBB) ** 2 * Raster(dist) ** 2) ** 3), 0)
                                    lacaBBOut.save(distOutBB)
                                    listBB.append(distOutBB)
        else:
            print("SPU FID " + str(fid) + " isn't right for BB calc")
        # arcpy.AddMessage("Bellbird calcs done")
        print("Bellbird steps done")

# arcpy.AddMessage("Done here")
print("Done here")

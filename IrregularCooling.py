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
# inputFC = arcpy.GetParameterAsText(0)
ws = arcpy.env.workspace = 'E:\IRPs'
inputFC = 'Polys.shp'
# get workspace from user
# ws = arcpy.GetParameterAsText(1)

# arcpy.env.

# get extent from user
arcpy.env.extent = 'BDY.shp'

cline = "cline.shp"
dcool = 5

# get cellsize from user
# cellSize = arcpy.GetParameterAsText(3)

# get output file name from user
outName = "IrCool"

arcpy.CheckOutExtension("Spatial")
arcpy.CheckOutExtension("Foundation")

# while loop to iterate through each feature in input layer - cursor
# cur = arcpy.da.SearchCursor(input, #, "FID")
# row = cur.Next()\

# Setup irregular polygon for analysis: convert to line, generate points, create centreline, calculate distance to points
# create centerline

with arcpy.da.SearchCursor(inputFC, ['FID', 'Shape@', 'Shape_Area', 'CC', 'd', 'NEAR_DIST']) as cursorfc:
    for rowfc in cursorfc:
        fid = rowfc[0]
        # inputFCFL = arcpy.MakeFeatureLayer_management(inputFC, "FCFL" + str(fid))
        # if rowfc[2] < 4900:
        if rowfc[2] >= 4900.00 or (rowfc[2] < 4900.00 and rowfc[5] <= float(2 * ((rowfc[2] / math.pi) ** 0.5))):

            cc = float(rowfc[3])

            # new Cooling code here
            # tCenterline1 = gdb + "\cline_" + str(fid)
            # arcpy.PolygonToCenterline_topographic(inputFC, outCenterline1)
            # outCenterline = arcpy.conversion.FeatureClassToShapefile(outCenterline1, ws)
            # polygon to line
            arcpy.AddMessage("Creating SPU Line")
            outLine = "SPULine_" + str(fid) + ".shp"
            arcpy.AddMessage("Creating points on SPU Line")
            arcpy.PolygonToLine_management(rowfc[1], outLine, "")
            # autogenerate points
            bdyp = "BDYPoints_" + str(fid) + ".shp"
            arcpy.GeneratePointsAlongLines_management(outLine, bdyp, 'DISTANCE', Distance="20 meters")
            # Near boundary points to centerline
            # count = int(arcpy.GetCount_management(bdyp)["row_count"])
            # print(type(count))
            # print("using row_count: " + str(count["row_count"]) + " " + str(type(count["row_count"])))
            # print("for feature Count is: " + str(count))
            if int(arcpy.GetCount_management(bdyp)["row_count"]) > 0:
                print(arcpy.GetCount_management(bdyp))
                # print(str(fid) + " is too small to get points")
                # pass
            # else:
                print(str(fid) + " is big enough")
                arcpy.Near_analysis(bdyp, cline)
                minDistC = (max([cur[0] for cur in arcpy.da.SearchCursor(bdyp, "NEAR_DIST")]) / dcool)
                print("minDistC = " + str(minDistC))
                arcpy.AddMessage("minDistC = " + str(minDistC))


            # del rowp
            #del cursorp

            print("Done")
            arcpy.AddMessage("Cooling calcs done.")

            # radius = float((row[2] / math.pi)**0.5)
            # if row[2] > 4000.00 and row[2] < 5000.00:
            #     d = 2 * radius/1.75
            # else:
            #     d = 2 * radius

            # ha = float(row[2])/10000
            # coolOut = ha * cc * Exp(Raster(-1 * distIn)/d)
            # coolOut.save("cool_" + str(fid) + ".tif")

            # dcalcBB = 500.00

del rowfc
del cursorfc

# Combine all grids using Mosaic to New Raster
# rasters = arcpy.ListRasters("cool*", "")
# output = outName + ".tif"
# proj = arcpy.SpatialReference(2193)
# arcpy.MosaicToNewRaster_management(rasters, ws, output, proj, "32_BIT_FLOAT", 5, "1", "MAXIMUM")

arcpy.CheckInExtension("Spatial")
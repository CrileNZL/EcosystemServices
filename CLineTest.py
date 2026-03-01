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

import arcpy

import math

from arcpy.sa import *

arcpy.env.overwriteOutput = True

# Get Clump polygon layer from user
inputFC = r'J:\Current_Projects\RichardM\IRPs\IR.gdb\Polys'

# Get workspace from user
ws = arcpy.env.workspace = 'J:\Current_Projects\RichardM\IRPs'
# arcpy.env.workspace = ws

gdbs = arcpy.ListWorkspaces("*", "FileGDB")
arcpy.AddMessage(gdbs)
print(gdbs)
gdb = gdbs[0]
arcpy.AddMessage(gdb)
print(gdb)


# Run NEAR on inputFC for FT habitat calc
arcpy.Near_analysis(inputFC, inputFC)

# Create distance grids for use in ES code
with arcpy.da.SearchCursor(inputFC, ['OBJECTID', 'Shape@', 'Shape_Area', 'CC', 'd', 'NEAR_DIST']) as cursor:
    for row in cursor:

        # Set up loop
        # for each feature, derive distance grid - Euclidean Distance
        # distances grids used by Nitrogen, Bellbird and Fantail ESs only
        # Cooling uses custom distance grids
        # get FID here for iterating
        # fid = row[0]
        # distOut = "dist_" + str(fid) + ".tif"
        # arcpy.env.workspace = ws
        # maxDist = '#' # can use ""?
        # dist = arcpy.sa.EucDistance(row[1], "", cellSize)
        # dist.save(distOut)

        # Calculate cooling raster
        # distIn = Raster(distOut)

        # Calculate individual ESs
        # Cooling and Bellbird habitat here
        # select all SPUs great than 4900 m2 in area and those smaller that are within 2*R of another SPU
        # new Cooling code here
        # For each SPU, create centerline and points around boundary
        # For each point, calculate distance to centerline
        # Use distance as Max Distance to calculate distance then calculate cooling effect and later mosaic to
        # new raster.
        # if row[2] >= 4900.00 or (row[2] < 4900.00 and row[5] <= float(2 * ((row[2] / math.pi) ** 0.5))):
        if row[2] >= 400.00:

            # cc = float(row[3])

            # new Cooling code here
            fid = row[0]
            outCenterline = gdb + "\clinenew"
            # outCenterline = r"cline_" + str(fid) + ".shp"
            arcpy.PolygonToCenterline_topographic(inputFC, outCenterline)
del row
del cursor
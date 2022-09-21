# Script to calculate cooling algorithm for Richard Morris
# C. Doscher - August 2022

import arcpy
# from arcpy.ia import RasterCalculator
import math
from arcpy.sa import *


arcpy.env.overwriteOutput = True

# dcalc = float(d)

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
        arcpy.env.workspace = ws
        distIn = Raster(distOut)
        #coolName = "cool_" + str(fid) + ".tif"
        # expExp = Raster(math.exp(distIn/d))
        # nextExp = Raster(ha * expExp)
        # coolOut = Raster(cc/nextExp)
        # coolOut = float(%cc%)/( float(%ha%) * Exp(distOut/ float(%d%)

# model: y = -3.9x3 + 89.1x2 - 814.5x + 2968
nOut = (-3.9 * distIn)**3 + (89.1 * distIn)**2 - (814.5 * distIn) + 2968
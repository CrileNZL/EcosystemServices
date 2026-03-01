# 2 March 2026
# Script to transform LWS metric rasters for use in ESMAX

import arcpy
import math
from arcpy.sa import *
from arcpy import env

arcpy.env.overwriteOutput = True

# set working directory
ws = arcpy.env.workspace = r"G:\LandscapeMAX\Clint\""

# Read in runoff raster
runoff = arcpy.Raster(r"G:\LandscapeMAX\Clint\TIFFs\runoff.tif")

# Calculate L+ raster
ROPos = runoff * 5
ROPos.save(r"G:\LandscapeMAX\Clint\LMax\ROPos.tif")

# Calculate L- raster
ROInv = (1 - runoff) * 5
ROInv.save(r"G:\LandscapeMAX\Clint\LMax\ROInv.tif")

import arcpy
# from arcpy.ia import RasterCalculator
import math
from arcpy.sa import *

# arcpy.env.workspace = r"D:\OLW Nitrates\OLWNitrates.gdb" # memory?
arcpy.env.overwriteOutput = True

ws = arcpy.env.workspace = r"M:\DataRM"

inputFC = "SPU_2.shp"

bdy = "R_bndy_ESMAX.shp"
arcpy.env.extent = bdy

cellSize = 1
cc = 10

# calc distance grid from SPUs - used for N and for FT
nDist = arcpy.sa.EucDistance(inputFC, "", cellSize)
nDist.save("nDist.tif")

# calc N = 2065 - 235 * distIn within 7 m of SPU
nOut = Con(nDist <= 7, (-235 * nDist + 2065), 0)
nOut.save("NStart.tif")

# Hillshade, azimuth = 315, altitude = 45, no shadows
nHS = Hillshade(nOut, 315, 45, "SHADOWS", 0.01)
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
# 2 March 2026
# Script to transform LWS metric rasters for use in ESMAX

import arcpy
import math
from arcpy.sa import *
from arcpy import env

arcpy.env.overwriteOutput = True

bdy = r"G:\LandscapeMAX\Clint\BDY.gdb\Catchment_Export"
extent = arcpy.Describe(bdy).extent
sr = arcpy.SpatialReference(2193)

# set working directory
ws = arcpy.env.workspace = r"G:\LandscapeMAX\Clint\Test"

# generate tesselations
tesselation = arcpy.management.GenerateTessellation(r"G:\LandscapeMAX\Clint\LMax\tiles.shp", extent, "TRANSVERSE_HEXAGON", "2 Hectares", sr)
# whereclause1 = "'tesselation', 'INTERSECT', 'bdy'"
print("Hexagons done"
      "")
select = arcpy.management.SelectLayerByLocation(tesselation,'INTERSECT', bdy)
arcpy.conversion.ExportFeatures(select, r"G:\LandscapeMAX\Clint\LMAX\Hex.shp")
print("New Hexes done")

# Take LWS inputs and transform to Pos, Inv or Log outputs
# ee = extreme event, dr = drainage, fa = flow accumulation, ro = runoff propensity, wp = wetness persistence
nameList = ["ee", "dr", "fa", "ro", "wp"]
rasList = arcpy.ListRasters()

for raster in rasList:
      ras = Raster(raster)
      for name in nameList:
            pos = ras * 5
            pos.save(name + r"Pos.tif")
            inv = (1 - ras) * 5
            inv.save(name + "Inv.tif")
            log = ROLog = 5 / (1 + (49 * Exp(-10 * ras)))
            log.save(name + "Log.tif")

# calculate LMAX cases for trees and wetlands
# get LMAX layers
lmaxList = arcpy.ListRasters("*_*")

treeh = (3*Raster("ro_Pos.tif")) + (3*Raster("wp_Pos.tif")) + (2*Raster("ee_Pos.tif")) + (3*Raster("fa_Pos.tif")) + (3*Raster("dr_Pos.tif"))
treeh.save("treeh.tif")

treenh = ((2*Raster("ro_Pos.tif")) + Raster("ro_Inv.tif")) + ((2*Raster("wp_Inv.tif")) + Raster("wp_Pos.tif")) + (2*Raster("ee_Pos.tif")) + ((2*Raster("fa_Pos.tif")) + Raster("fa_Inv.tif")) + ((2*Raster("dr_Pos.tif")) + Raster("dr_Inv.tif"))
treenh.save("treenh.tif")

wetland = ((2*Raster("ro_Pos.tif")) + Raster("ro_Inv.tif")) + ((2*Raster("wp_Pos.tif")) + Raster("wp_Log.tif")) + (Raster("ee_Log.tif") + Raster("ee_Inv.tif")) + ((2*Raster("fa_Pos.tif")) + Raster("fa_Inv.tif")) + (3*Raster("dr_Pos.tif"))
wetland.save(wetland.tif)

# Zonal Statistics as Table
# arcpy.sa.ZonalStatisticsAsTable(r"G:\LandscapeMAX\Clint\LMAX\Hex.shp", "GRID_ID", r"G:\LandscapeMAX\Clint\LMAX\ROPos.tif", r"G:\LandscapeMAX\Clint\LMAX\ROPosOut.dbf", "DATA", "MEAN", "", "", "", "", "", "ROPosJoin")
# print("ZSaT done")

# arcpy.management.CopyFeatures("ROPosJoin", r"G:\LandscapeMAX\Clint\LMAX\ROPos.shp")

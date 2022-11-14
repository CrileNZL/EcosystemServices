import arcpy

from arcpy.sa import *

arcpy.env.overwriteOutput = True

# Get Clump polygon layer from user
inputFC = arcpy.GetParameterAsText(0)

boundary = arcpy.GetParameterAsText(1)
arcpy.env.extent = boundary

ws = arcpy.GetParameterAsText(2)
arcpy.env.workspace = ws

outName = "Test"

# Create control grids here

# Nitrogen  and control
arcpy.PolygonToRaster_conversion(inputFC, "Shape_Area", "NC.tif", "CELL_CENTER", "", 5)
# use Con to change NoData values to 0 and existing values to NoData
nCon = Con(IsNull("NC.tif"), 0, Raster("NC.tif"))
nCon.save("NControl.tif")

# Cooling and Bellbird Habitat and controls
whereClause1 = "Shape_Area >= 15000"
arcpy.SelectLayerByAttribute_management(inputFC, "NEW_SELECTION", whereClause1)
arcpy.PolygonToRaster_conversion(inputFC, "Shape_Area", "CC.tif", "CELL_CENTER", "", 5)
# use Con to change NoData values to 0 and existing values to NoData
cC = Con(IsNull("CC.tif"), 0, Raster("CC.tif"))
cC.save("CControl.tif")
cC.save("HBBControl.tif")

whereClause2 = "Shape_Area < 15000"
arcpy.SelectLayerByAttribute_management(inputFC, "NEW_SELECTION", whereClause2)
arcpy.PolygonToRaster_conversion(inputFC, "Shape_Area", "CCM.tif", "CELL_CENTER", "", 5)
cCM = Con(IsNull("CCM.tif"), 0)
cCM.save("CCmask.tif")

# Fantail Habitat and control
whereClause3 = "Shape_Area >= 15000 or (Shape_Area < 15000 and (NEAR_DIST <= 150 and NEAR_DIST > 0))"
arcpy.SelectLayerByAttribute_management(inputFC, "NEW_SELECTION", whereClause3)
arcpy.PolygonToRaster_conversion(inputFC, "Shape_Area", "HFTC.tif", "CELL_CENTER", "", 5)
# use Con to change NoData values to 0 and existing values to NoData
hFTC = Con((IsNull("HFTC.tif")), 0, Raster("HFTC.tif"))
hFTC.save("HFTControl.tif")


whereClause4 = "Shape_Area < 15000 and NEAR_DIST > 150"
arcpy.SelectLayerByAttribute_management(inputFC, "NEW_SELECTION", whereClause4)
arcpy.PolygonToRaster_conversion(inputFC, "Shape_Area", "HFTCM.tif", "CELL_CENTER", "", 5)
hFTCM = Con(IsNull("HFTCM.tif"), 0)
hFTCM.save("HFTCmask.tif")

arcpy.SelectLayerByAttribute_management(inputFC, "CLEAR_SELECTION")

# Create metascores here
arcpy.sa.ZonalStatisticsAsTable(boundary, "OBJECTID", Raster("NControl.tif"), outName + "_MSNControl.dbf", "DATA", "SUM")
arcpy.sa.ZonalStatisticsAsTable(boundary, "OBJECTID", Raster("HBBControl.tif"), outName + "_MSBBControl.dbf", "DATA", "SUM")
arcpy.sa.ZonalStatisticsAsTable(boundary, "OBJECTID", Raster("HFTControl.tif"), outName + "_MSFTControl.dbf", "DATA", "SUM")
if arcpy.Exists(Raster("CControl.tif")):
    arcpy.sa.ZonalStatisticsAsTable(boundary, "OBJECTID", Raster("CControl.tif"), outName + "_MSCControl.dbf", "DATA", "SUM")

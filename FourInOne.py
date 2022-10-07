## Calculate ecosystem services layers and metascores for individual ES rasters
## Combines Cooling.py, Nitrogen.py, Laca.py and metascores.py into one script
## Developed for Richard Morris
## C. Doscher October 2022

import arcpy
# from arcpy.ia import RasterCalculator
import math
from arcpy.sa import *

arcpy.env.overwriteOutput = True

# Get Clump polygon layer from user
inputFC = arcpy.GetParameterAsText(0)

# Get input ES raster from user
inputRas = arcpy.GetParameterAsText(1)

# Get workspace from user
arcpy.env.workspace = arcpy.GetParameterAsText(2)

# Get extent from user
arcpy.env.extent = arcpy.GetParameterAsText(3)

# Get boundary zome FC from user
boundary = arcpy.GetParameterAsText(4)

# Get cell size from user
# get cellsize from user
cellSize = arcpy.GetParameterAsText(5)

# get output file name from user
outName = arcpy.GetParameterAsText(6)

## Create distance grids for use in ES code
# can I do this in one script given minimmum clump size in Cooling.py?


## Calculate individual ESs
# Cooling here


# Nitrogen here


# Habitat here


## Calculate metascores for each ES



## Clean up
# Delete distance and individual ES grids
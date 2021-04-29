#!/bin/sh
#===================================================================================#
# Convert .e00 archives to shapefiles using avcimport and ogr2ogr. It accepts file
# names or absolut paths. No relative paths.
#-----------------------------------------------------------------------------------#
#
# usage :		$ e00toshp.sh input.e00 output
#
# The second argument is used as the path and  file name of the output shp, shx,
# dbf, and prj files. 
#

# Separate path and filename
IPATH=`dirname $1`

# Convert from e00      to coverage
cd $IPATH
avcimport $1 COV 

# Convert from coverage to shapefiles
ogr2ogr -f "ESRI Shapefile" SHP COV

# Move where appropriate
mv SHP/PAL.shp $2.shp
mv SHP/PAL.shx $2.shx
mv SHP/PAL.dbf $2.dbf
mv SHP/PAL.prj $2.prj

# Clean up
rm -rf COV info SHP

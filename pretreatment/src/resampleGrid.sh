#!/bin/bash

DIRECTORY="$1"
FILE_IN="$2"
FILE_OUT="$3"
VAR_NAME="$4"
xres="$5"
yres="$6"
fillval="$7"
longMini="$8"
latiMini="$9"
longMax="${10}"
latiMax="${11}"

# Translate the desired variable to GeoTIFF
gdal_translate -ot Float64 NETCDF:$DIRECTORY/$FILE_IN:$VAR_NAME -unscale $DIRECTORY/file.tiff

#resampling
gdalwarp $DIRECTORY/file.tiff -tr $xres $yres -r near -srcnodata $fillval -te $longMini $latiMini $longMax $latiMax $DIRECTORY/resampled_file.tiff -t_srs EPSG:4326

# Translate back to NetCDF
gdal_translate -of NetCDF $DIRECTORY/resampled_file.tiff $DIRECTORY/$FILE_OUT
ncatted -O -a valid_max,time,d,, $DIRECTORY/$FILE_OUT

# Remove unnecessary files
rm $DIRECTORY/file.tiff
rm $DIRECTORY/file.tiff.aux.xml
rm $DIRECTORY/resampled_file.tiff
rm $DIRECTORY/resampled_file.tiff.aux.xml
#!/bin/bash

DIRECTORY="$1"
FILE_IN="$2"
FILE_OUT="$3"
VAR_NAME="$4"

# Translate the desired variable to GeoTIFF
gdal_translate -ot Float64 NETCDF:$DIRECTORY/$FILE_IN:$VAR_NAME -unscale $DIRECTORY/file_Stereo.tiff

# resampling
gdalwarp -s_srs "+proj=stere +lat_0=90 +lat_ts=90 +lon_0=-45 +k=0.00001 +x_0=0 +y_0=0 +a=6378273 +b=6378273  +units=m +no_defs" $DIRECTORY/file_Stereo.tiff -t_srs EPSG:4326 $DIRECTORY/file_EPSG.tiff -overwrite

# Translate back to NetCDF
gdal_translate -of NetCDF $DIRECTORY/file_EPSG.tiff $DIRECTORY/$FILE_OUT

# Remove unnecessary files
rm $DIRECTORY/file_Stereo.tiff
rm $DIRECTORY/file_Stereo.tiff.aux.xml
rm $DIRECTORY/file_EPSG.tiff
rm $DIRECTORY/file_EPSG.tiff.aux.xml
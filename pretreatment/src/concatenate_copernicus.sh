#!/bin/bash

ZONE="$1"
PARAM="$2"

DIR_IN="/media/share/data/$ZONE/$PARAM"
DIR_OUT="/media/share/data_merged/$ZONE/$PARAM"
mkdir -p $DIR_OUT

if [[ ! -d $DIR_IN ]]; then
    echo "Directory $DIR_IN does not exist"
fi

# Get the name of the first file, without the directories
first_file_name=$(ls $DIR_IN/* | head -n 1 | xargs -n 1 basename)
echo "first_file_name: $first_file_name"

echo "On first file, make 'time' into the record dimension:"
# On first file, make 'time' into the record dimension:
ncks -O --mk_rec_dim time $DIR_IN/$first_file_name $DIR_OUT/$first_file_name

echo "Then concatenate all files (the modified first file, then all the unmodified files):"
# Then concatenate all files (the modified first file, then all the unmodified files):
shopt -s extglob
ncrcat -h $DIR_OUT/$first_file_name $DIR_IN/!($first_file_name) $DIR_OUT/${ZONE}_${PARAM}_merged.nc

echo "Finally change 'valid_max' attribute to a large value:"
# Finally change 'valid_max' attribute to a large value:
ncatted -O -a valid_max,time,o,d,999999. $DIR_OUT/${ZONE}_${PARAM}_merged.nc

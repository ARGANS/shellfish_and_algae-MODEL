#!/bin/bash
DIR_IN="$1"
DIR_OUT="$2"
FILE_NAME="$3"
mkdir -p $DIR_OUT
if [[ ! -d $DIR_IN ]]; then
    echo "Directory $DIR_IN does not exist"
fi

# Get the name of the first file, without the directories
first_file_name=$(ls $DIR_IN/* | head -n 1 | xargs -n 1 basename)
echo "first_file_name: $first_file_name"

# On first file, make 'time' into the record dimension and store to a temporary file:
ncks --mk_rec_dim time $DIR_IN/$first_file_name $DIR_OUT/tmp.nc

# Then concatenate all files (the modified first file, then all the unmodified files):
shopt -s extglob
ncrcat -h $DIR_OUT/tmp.nc $DIR_IN/!($first_file_name) $DIR_OUT/$FILE_NAME

# remove the temporary file
rm -f $DIR_OUT/tmp.nc

# Finally change 'valid_max' attribute to a large value:
#ncatted -O -a valid_max,time,o,d,999999. $DIR_OUT/$FILE_NAME
# Finally delete the 'valid_max' attribute
ncatted -O -a valid_max,time,d,, $DIR_OUT/$FILE_NAME
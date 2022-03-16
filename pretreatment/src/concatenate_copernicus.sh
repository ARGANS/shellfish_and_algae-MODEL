ZONE="$1"
PARAM="$2"

DIR_IN="/media/share/data/$ZONE/$PARAM"
DIR_OUT="/media/share/data_merged/$ZONE/$PARAM"

# Get the name of the first file, without the directories
first_file=$(ls -t $DIR_IN/* | head -n 1)
first_file_name=$($first_file#$DIR_IN/)

# On first file, make 'time' into the record dimension:
ncks -O --mk_rec_dim time $DIR_IN/$first_file_name $DIR_OUT/$first_file_name

# Then concatenate all files (the modified first file, then all the unmodified files):
ncrcat -h $DIR_OUT/$first_file_name $DIR_IN/!($first_file_name) $DIR_OUT/$ZONE_$PARAM_merged.nc

# Finally change 'valid_max' attribute to a large value:
ncatted -O -a valid_max,time,o,d,999999. $DIR_OUT/$ZONE_$PARAM_merged.nc
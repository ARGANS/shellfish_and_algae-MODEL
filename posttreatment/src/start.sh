#!/bin/bash
input_path="$SOURCE_DIR"
destination=$input_path/posttreatment
tmp_path=/tmp/$(basename $input_path)


mkdir -p $destination
mkdir -p $tmp_path
rm -rf $tmp_path/*


echo -n $(date +%s) > $destination/start.mark 
cp $input_path/parameters.json $destination/parameters.json

#  make_interest_vars.R mutates the inputed data
cp $input_path/concat.nc $tmp_path/concat.nc
./make_interest_vars.R $tmp_path/concat.nc $input_path/parameters.json


for variable in 'DW' 'DW_line' 'DW_PUA' 'FW' 'FW_line' 'FW_PUA' 'kcal_PUA' 'protein_PUA' 'Biomass_CO2' 'CO2_uptake_PUA'; do
    gdal_translate NETCDF:"$tmp_path/concat.nc":$variable $destination/$variable.tif
done 


echo -n $(date +%s) > $destination/end.mark 
rm -rf $tmp_path

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

# ??? catch error log
error_log=$destination/error.txt
print_log=$destination/print.txt

# Store the type of simulation
sim_type=`cat $destination/parameters.json | jq -r ".type"`

if [ $sim_type == "Algae" ]; then

    python ./make_interest_vars.py -i $tmp_path/concat.nc -j $input_path/parameters.json 1>$print_log 2>$error_log
    if [ $? -eq 0 ]; then
        echo "Success"
    else
        echo "Failure:"
        echo 'Finished '$(date "+%d/%m/%Y %H:%M:%S") >> $error_log
        cat $error_log
    fi

    for variable in 'NO3' 'NH4' 'DW' 'DW_line' 'DW_PUA' 'FW' 'FW_line' 'FW_PUA' 'kcal_PUA' 'protein_PUA' 'Biomass_CO2' 'CO2_uptake_PUA' 'NO3field' 'NH4field' 'D' 'N_f' 'N_s' 'avNO3' 'avNH4' 'cNO3' 'cNH4'; do
        gdal_translate NETCDF:"$tmp_path/concat.nc":$variable $destination/$variable.tif 1>$print_log 2>$error_log
    done
    #nbrIntermediateFiles = ls $input_path/concat?*.nc| wc -l
    for fileName in $input_path/concat?*.nc; do
        base_name=$(basename ${fileName%.nc})
        for variable in 'cNO3' 'cNH4'; do
          gdal_translate NETCDF:"$fileName":$variable $destination/$base_name.tif 1>$print_log 2>$error_log
        done
    done

elif [ $sim_type == "Shellfish" ]; then
    for variable in 'DSTW' 'STE' 'FW' 'DWW' 'SHL' 'NH4_production' 'CO2_production'; do
        gdal_translate NETCDF:"$tmp_path/concat.nc":$variable $destination/$variable.tif 1>$print_log 2>$error_log
    done
fi

# Finishing steps
echo '-------------------- Last 100 lines printed to stdout --------------------' >> $error_log
tail -n 100 $print_log >> $error_log
echo $input_path/parameters.json >> $error_log

echo -n $(date +%s) > $destination/end.mark 
rm -rf $tmp_path

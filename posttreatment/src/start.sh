#!/bin/bash
input_path="$SOURCE_DIR"
destination=$input_path/posttreatment
tmp_path=/tmp/$(basename $input_path)

mkdir -p $destination
mkdir -p $tmp_path
rm -rf $tmp_path/*


echo -n $(date +%s) > $destination/start.mark 
cp $input_path/parameters.json $destination/parameters.json

# ??? catch error log
error_log=$destination/error.txt
print_log=$destination/print.txt

# Store the type of simulation and the zone
sim_type=`cat $destination/parameters.json | jq -r ".type"`
sim_zone=`cat $destination/parameters.json | jq -r ".metadata.zone"`


function posttreatment_Algae {
    concat_name="$1"

    mkdir $destination/$concat_name

    #  make_interest_vars.py mutates the inputed data
    cp $input_path/$concat_name.nc $tmp_path/$concat_name.nc

    # First add the interest variables to concat.nc
    python ./make_interest_vars.py -i $tmp_path/$concat_name.nc -j $input_path/parameters.json

    for variable in 'CMEMS_NO3' 'CMEMS_NH4' 'DW' 'DW_line' 'DW_PUA' 'FW' 'FW_line' 'FW_PUA' 'kcal_PUA' 'protein_PUA' 'Biomass_CO2' 'CO2_uptake_PUA' 'D' 'N_f' 'N_s' 'avNO3' 'avNH4' 'cNO3' 'cNH4'; do
        gdal_translate NETCDF:"$tmp_path/$concat_name.nc":$variable $destination/$concat_name/$variable.tif
    done

    # Translate the intermediate files.
    for fileName in $input_path/$concat_name??.nc; do
        base_name=$(basename ${fileName%.nc})
        mkdir $destination/$concat_name/$base_name
        for variable in 'cNO3' 'cNH4' 'CMEMS_NO3' 'CMEMS_NH4'; do
            gdal_translate NETCDF:"$fileName":$variable $destination/$concat_name/$base_name/$variable.tif
        done
    done
}


function posttreatment_Shellfish {
    concat_name="$1"

    for variable in 'DSTW' 'STE' 'FW' 'DWW' 'SHL' 'NH4_production' 'CO2_production'; do
        gdal_translate NETCDF:"$input_path/$concat_name.nc":$variable $destination/$variable.tif
    done
}


if [ $sim_zone == "Europe" ]; then
    for zone_name in 'IBI' 'NWS' 'MED' 'Baltic' 'BS' 'Arctic'; do
        posttreatment_$sim_type "concat_$zone_name" 1>$print_log 2>$error_log
    done
    # TODO: then merge all areas
else
    posttreatment_$sim_type "concat" 1>$print_log 2>$error_log
    # When only one area, the full map is this area's map
    cp $destination/concat/* $destination/.
fi


# Finishing steps
echo '-------------------- Last 100 lines printed to stdout --------------------' >> $error_log
tail -n 100 $print_log >> $error_log
echo $input_path/parameters.json >> $error_log

echo -n $(date +%s) > $destination/end.mark 
rm -rf $tmp_path

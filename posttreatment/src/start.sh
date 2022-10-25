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
    mkdir $destination/$concat_name/params

    #  make_interest_vars.py mutates the inputed data
    cp $input_path/$concat_name.nc $tmp_path/$concat_name.nc

    # First add the interest variables to concat.nc
    python ./make_interest_vars.py -i $tmp_path/$concat_name.nc -j $input_path/parameters.json

    for variable in 'CMEMS_NO3' 'CMEMS_NH4' 'DW' 'DW_line' 'DW_PUA' 'FW' 'FW_line' 'FW_PUA' 'kcal_PUA' 'protein_PUA' 'Biomass_CO2' 'CO2_uptake_PUA' 'D' 'N_f' 'N_s' 'avNO3' 'avNH4' 'cNO3' 'cNH4'; do
        gdal_translate NETCDF:"$tmp_path/$concat_name.nc":$variable $destination/$concat_name/params/$variable.tif
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

    mkdir $destination/$concat_name
    mkdir $destination/$concat_name/params

    for variable in 'DSTW' 'STE' 'FW' 'DWW' 'SHL' 'NH4_production' 'CO2_production'; do
        gdal_translate NETCDF:"$input_path/$concat_name.nc":$variable $destination/$concat_name/params/$variable.tif
    done
}


function read_zee_params {
    zee_info=`gdalinfo -json /media/global/zee_europe.tif`

    echo "zee info: $zee_info"

    lon_min=`echo $zee_info | jq -r ".cornerCoordinates.lowerLeft[0]"`
    lat_min=`echo $zee_info | jq -r ".cornerCoordinates.lowerLeft[1]"`
    lon_max=`echo $zee_info | jq -r ".cornerCoordinates.upperRight[0]"`
    lat_max=`echo $zee_info | jq -r ".cornerCoordinates.upperRight[1]"`

    lon_step=`echo $zee_info | jq -r ".geoTransform[1]"`
    lat_step=`echo $zee_info | jq -r ".geoTransform[5]"`

    # absolute values of step
    #lon_step=${lon_step#-}
    #lat_step=${lat_step#-}
}


function resample_to_ZEE {
    # Resample a tif file to the ZEE grid

    local file_in="$1"
    local file_out="$2"

    local tmpfile="$destination/tmp.tif"

    cmd="gdalwarp $file_in $tmpfile -tr $lon_step $lat_step -te $lon_min $lat_min $lon_max $lat_max -t_srs EPSG:4326 -r cubicspline"
    echo "COMMAND: $cmd"
    $cmd

    cmd="gdal_translate -co compress=deflate $tmpfile $file_out"
    echo "COMMAND: $cmd"
    $cmd

    rm $tmpfile
}


function fusion_zones {
    # Fusion all zones into one European map

    # Superimpose zones in the correct order (left is bottom)
    ordered_zones=(Arctic Baltic BS NWS IBI MED)

    for param_name in 'CMEMS_NO3' 'CMEMS_NH4' 'DW' 'DW_line' 'DW_PUA' 'FW' 'FW_line' 'FW_PUA' 'kcal_PUA' 'protein_PUA' 'Biomass_CO2' 'CO2_uptake_PUA' 'D' 'N_f' 'N_s' 'avNO3' 'avNH4' 'cNO3' 'cNH4'; do

        local file_out=$destination/$param_name.tif
        local tmpfile="$destination/tmp.tif"

        # File names in the correct order
        local file_in_names=(${ordered_zones[@]/#/"$destination/concat_"})
        file_in_names=(${file_in_names[@]/%/"/params_1km/${param_name}_1km.tif"})

        #cmd="gdalwarp -overwrite -srcnodata \"9.96921e+36\" "$destination/*/params_1km/${param_name}_1km.tif" $tmpfile -tr $lon_step $lat_step -te $lon_min $lat_min $lon_max $lat_max -t_srs EPSG:4326"
        cmd="gdalwarp -overwrite "${file_in_names[*]}" $tmpfile -tr $lon_step $lat_step -te $lon_min $lat_min $lon_max $lat_max -t_srs EPSG:4326"
        echo "COMMAND: $cmd"
        $cmd

        cmd="gdal_translate -co compress=deflate $tmpfile $file_out"
        echo "COMMAND: $cmd"
        $cmd

        rm $tmpfile
    done
}


if [ $sim_zone == "Europe" ]; then
    for zone_name in 'IBI' 'NWS' 'MED' 'Baltic' 'BS' 'Arctic'; do
        posttreatment_$sim_type "concat_$zone_name" 1>>$print_log 2>>$error_log
    done
    # read zee tif parameters
    read_zee_params 1>>$print_log 2>>$error_log

    # Resample all parameters to 1km
    for zone_name in 'IBI' 'NWS' 'MED' 'Baltic' 'BS' 'Arctic'; do
        mkdir $destination/concat_$zone_name/params_1km
        for param_file in $destination/concat_$zone_name/params/*; do
            param_file_1km=$destination/concat_$zone_name/params_1km/$(basename ${param_file%.tif})_1km.tif
            resample_to_ZEE "$param_file" "$param_file_1km" 1>>$print_log 2>>$error_log
        done
    done

    # Fusion all zones for each parameter
    fusion_zones 1>>$print_log 2>>$error_log
else
    posttreatment_$sim_type "concat" 1>>$print_log 2>>$error_log
    # When only one area, the full map is this area's map
    cp $destination/concat/params/* $destination/.
fi


# Finishing steps
echo '-------------------- Last 100 lines printed to stdout --------------------' >> $error_log
tail -n 100 $print_log >> $error_log
echo $input_path/parameters.json >> $error_log

echo -n $(date +%s) > $destination/end.mark 
rm -rf $tmp_path

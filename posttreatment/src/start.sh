#!/bin/bash
input_path="$SOURCE_DIR"
destination="$INPUT_DESTINATION"
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

params_algae=('CMEMS_NO3' 'CMEMS_NH4' 'DW' 'DW_line' 'DW_PUA' 'FW' 'FW_line' 'FW_PUA' 'kcal_PUA' 'protein_PUA' 'Biomass_CO2' 'CO2_uptake_PUA' 'D' 'N_f' 'N_s' 'avNO3' 'avNH4' 'cNO3' 'cNH4')
params_shellfish=('DSTW' 'STE' 'FW' 'DWW' 'SHL' 'NH4_production' 'CO2_production')
continuous_param="CMEMS_NO3 CMEMS_NH4 D cNO3 cNH4 avNO3 avNH4"
function posttreatment_Algae {
    concat_name="$1"

    mkdir $destination/$concat_name
    mkdir $destination/$concat_name/params

    #  make_interest_vars.py mutates the inputed data
    cp $input_path/$concat_name.nc $tmp_path/$concat_name.nc

    # First add the interest variables to concat.nc
    python ./make_interest_vars.py -i $tmp_path/$concat_name.nc -j $input_path/parameters.json

    for variable in ${params_algae[*]}; do
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

    for variable in ${params_shellfish[*]}; do
        gdal_translate NETCDF:"$input_path/$concat_name.nc":$variable $destination/$concat_name/params/$variable.tif
    done
}

function exists_in_list() {
    LIST=$1
    DELIMITER=$2
    VALUE=$3
    LIST_WHITESPACES=`echo $LIST | tr "$DELIMITER" " "`
    for x in $LIST_WHITESPACES; do
        echo $x
        echo $VALUE
        if [ $x = $VALUE ]; then
            return 0
        fi
    done
    return 1
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
    local var_name=$(basename ${file_in%.tif})

    local out_basename="${destination##*/}"

    echo $out_basename

    local tmpfile="$destination/tmp.tif"
    local tmpfile2="$destination/tmp2.tif"

    if [ $out_basename = '_posttreatment' ]; then
        if exists_in_list "$continuous_param" " " $var_name; then
            interpol='cubicspline'
        else
            interpol='near'
            cmd="gdalwarp $file_in $tmpfile2 -tr $lon_step $lat_step -srcnodata 0 -dstnodata 9.969209968386869e+36 -te $lon_min $lat_min $lon_max $lat_max -t_srs EPSG:4326 -r $interpol"
            $cmd
            interpol='cubicspline'
            local file_in="$destination/tmp2.tif"
        fi
    else
        if exists_in_list "$continuous_param" " " $var_name; then
            interpol='cubicspline'
        else
            interpol='near'
        fi
    fi
    echo $var_name
    echo $interpol
    cmd="gdalwarp $file_in $tmpfile -tr $lon_step $lat_step -te $lon_min $lat_min $lon_max $lat_max -t_srs EPSG:4326 -r $interpol"
    echo "COMMAND: $cmd"
    $cmd

    cmd="gdal_translate -co compress=deflate $tmpfile $file_out"
    echo "COMMAND: $cmd"
    $cmd

    rm $tmpfile
    rm $tmpfile2
}


function fusion_zones {
    # Fusion all zones into one European map

    simulation_type="$1"
    if [ $simulation_type = "Algae" ]; then
        param_name_table=${params_algae[*]}
    else
        param_name_table=${params_shellfish[*]}
    fi

    # Superimpose zones in the correct order (left is bottom)
    ordered_zones=(Arctic Baltic BS NWS IBI MED)

    for param_name in ${param_name_table[*]}; do

        local file_out=$destination/$param_name.tif
        local ncfile_out=$destination/$param_name.nc
        local tmpfile="$destination/tmp.tif"

        # File names in the correct order
        local file_in_names=(${ordered_zones[@]/#/"$destination/concat_"})
        file_in_names=(${file_in_names[@]/%/"/params_1km/${param_name}_1km.tif"})

        cmd="gdalwarp -overwrite "${file_in_names[*]}" $tmpfile -tr $lon_step $lat_step -te $lon_min $lat_min $lon_max $lat_max -t_srs EPSG:4326"
        echo "COMMAND: $cmd"
        $cmd

        cmd="gdal_translate -co compress=deflate $tmpfile $file_out"
        echo "COMMAND: $cmd"
        $cmd

        cmd="gdal_translate -of NetCDF $file_out $ncfile_out"
        echo "COMMAND: $cmd"
        $cmd
        rm $tmpfile

        python ./add_metadata.py -i $ncfile_out -j $input_path/parameters.json -n ${param_name} 1>>$error_log 2>>$error_log
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
    fusion_zones $sim_type 1>>$print_log 2>>$error_log
else
    posttreatment_$sim_type "concat" 1>>$print_log 2>>$error_log
    # When only one area, the full map is this area's map
    cp $destination/concat/params/* $destination/.

    # read zee tif parameters
    read_zee_params 1>>$print_log 2>>$error_log

    # Resample to ZEE for farm optimization
    mkdir $destination/concat/params_1km
    for param_file in $destination/concat/params/*; do
        param_file_1km=$destination/concat/params_1km/$(basename ${param_file%.tif})_1km.tif
        ncfile_out=$destination/$(basename ${param_file%.tif})_1km.nc
        resample_to_ZEE "$param_file" "$param_file_1km" 1>>$print_log 2>>$error_log
        gdal_translate -of NetCDF "$param_file_1km" $ncfile_out
        python ./add_metadata.py -i $ncfile_out -j $input_path/parameters.json -n $(basename ${param_file%.tif})
    done
fi


# Finishing steps
echo '-------------------- Last 100 lines printed to stdout --------------------' >> $error_log
tail -n 100 $print_log >> $error_log
echo $input_path/parameters.json >> $error_log

echo -n $(date +%s) > $destination/end.mark 
rm -rf $tmp_path

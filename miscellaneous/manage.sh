#!/bin/bash

DIR=$(dirname $0)
echo $DIR
source $DIR/_common.sh
source $DIR/_dataread.sh
source $DIR/_dataread_b.sh
source $DIR/_dataimport.sh
source $DIR/_pretreatment.sh
source $DIR/_posttreatment.sh


function action_bash {
    local image_tag="ac-dataread/runtime"

    docker run \
        --rm \
        --volume "$SHARED_VOLUME_NAME:/media/share" \
        --volume "${HOME}:/media/home" \
        --volume "$GLOBAL_VOLUME_NAME":/media/global \
        --workdir=/media/share \
        -e PYTHONDONTWRITEBYTECODE=1 \
        --entrypoint=/bin/bash \
        -it $image_tag:latest
}


# DATASET PROPERTIES: year depth-min depth-max datasets (from metadata)
datasetProperties_json=`cat $DIR/__dataimport_input_parameters.json` 
modelProperties_json=`cat $DIR/__user1_alaria_IBI_26-07-2022.json`

datasetProperties_hash=`echo -n "$datasetProperties_json" | md5sum | head -c 32`
modelProperties_hash=`echo -n "$modelProperties_json" | md5sum | head -c 32`

dataimport_destination="/media/share/data/$datasetProperties_hash"

pretreatment_source="$dataimport_destination"
pretreatment_destination="$dataimport_destination/_pretreated"

dataread_source="$pretreatment_destination"
dataread_destination="$dataimport_destination/_dataread/$modelProperties_hash"

posttreatment_source="$dataread_destination"


scenario_b_json='{"parameters":{"species":{"alaria":{"options":{},"parameters":{"mu":0.1,"V_NH4":60,"V_NO3":25,"K_NH4":1,"K_NO3":5,"Q_max":70,"Q_min":14,"N_to_P":12,"K_c":7,"T_O":5,"T_min":1,"T_max":16,"I_s":60,"a_cs":0.00036,"d_m":0.003,"h_MA":0.4,"w_MA":0.2,"r_L":0.2,"r_N":0.1,"CN_MA":21,"kcal_MA":2.29,"prot_MA":0.08,"DF_MA":0.2}}},"farm":{"default":{"options":{},"parameters":{"y_farm":1000,"density_MA":0.4,"x_farm":1000,"z":2}}},"harvest":{"Summer_growth":{"options":{},"parameters":{"deployment_Nf":1000,"deployment_month":3,"harvesting_month":9}}},"run":{"default":{"options":{},"parameters":{"K_d490":0.1,"Von_Karman":0.4,"Detritus":0.1,"min_lon":-180,"max_lon":180,"min_lat":-90,"max_lat":90,"n_cores":10}}}},"metadata":{"name":"user1_alaria_IBI_11-08-2022","zone":"IBI","_suggested":{"login":"user1","species":"alaria","zone":"IBI","date":"11-08-2022"},"scenario":"A"},"dataset_parameters":{"depth-min":0,"depth-max":20,"year":2021,"datasets":{"Temperature":{"Name":"Copernicus Analysis & Forecast","type":"model","Parameter":"Temperature","Place":"IBI","motu":"https://nrt.cmems-du.eu/motu-web/Motu","service-id":"IBI_ANALYSISFORECAST_PHY_005_001-TDS","product-id":"cmems_mod_ibi_phy_anfc_0.027deg-3D_P1D-m","longitude-min":-19,"longitude-max":"5","latitude-min":26,"latitude-max":56,"resolution":"2.8*2.8","level":"4","first obs":"16/11/2019","last obs":"present","source":"marineCopernicus","variable":"thetao","link":"https://resources.marine.copernicus.eu/product-detail/IBI_ANALYSISFORECAST_PHY_005_001/INFORMATION","mask":"0","depth-min":0.494,"depth-max":5727.917,"fileType":"NetCDF","time":"","latName":"latitude","longName":"longitude","depthName":"depth","timeName":"time","frequency":"daily","unitFactor":"","timeOrigin":"1950-01-01T00:00:00Z","timeUnit":"hours","dataimport":"Copernicus","pretreatment":"Copernicus"},"Nitrate":{"Name":"Copernicus Analysis & Forecast","type":"model","Parameter":"Nitrate","Place":"IBI","motu":"https://nrt.cmems-du.eu/motu-web/Motu","service-id":"IBI_ANALYSISFORECAST_BGC_005_004-TDS","product-id":"cmems_mod_ibi_bgc_anfc_0.027deg-3D_P1M-m","longitude-min":-19,"longitude-max":"5","latitude-min":26,"latitude-max":56,"resolution":"2.8*2.8","level":"4","first obs":"16/12/2018","last obs":"present","source":"marineCopernicus","variable":"no3","link":"https://resources.marine.copernicus.eu/product-detail/IBI_ANALYSISFORECAST_BGC_005_004/INFORMATION","mask":"0","depth-min":0.494,"depth-max":5727.917,"fileType":"NetCDF","time":"","latName":"latitude","longName":"longitude","depthName":"depth","timeName":"time","frequency":"monthly","unitFactor":14,"timeOrigin":"1950-01-01T00:00:00Z","timeUnit":"hours","dataimport":"Copernicus","pretreatment":"Copernicus"},"Phosphate":{"Name":"Copernicus Analysis & Forecast","type":"model","Parameter":"Phosphate","Place":"IBI","motu":"https://nrt.cmems-du.eu/motu-web/Motu","service-id":"IBI_ANALYSISFORECAST_BGC_005_004-TDS","product-id":"cmems_mod_ibi_bgc_anfc_0.027deg-3D_P1M-m","longitude-min":-19,"longitude-max":"5","latitude-min":26,"latitude-max":56,"resolution":"2.8*2.8","level":"4","first obs":"16/12/2018","last obs":"present","source":"marineCopernicus","variable":"po4","link":"https://resources.marine.copernicus.eu/product-detail/IBI_ANALYSISFORECAST_BGC_005_004/INFORMATION","mask":"0","depth-min":0.494,"depth-max":5727.917,"fileType":"NetCDF","time":"","latName":"latitude","longName":"longitude","depthName":"depth","timeName":"time","frequency":"monthly","unitFactor":"","timeOrigin":"1950-01-01T00:00:00Z","timeUnit":"hours","dataimport":"Copernicus","pretreatment":"Copernicus"},"Ammonium":{"Name":"Copernicus Analysis & Forecast","type":"model","Parameter":"Ammonium","Place":"IBI","motu":"https://nrt.cmems-du.eu/motu-web/Motu","service-id":"IBI_ANALYSISFORECAST_BGC_005_004-TDS","product-id":"cmems_mod_ibi_bgc_anfc_0.027deg-3D_P1M-m","longitude-min":-19,"longitude-max":"5","latitude-min":26,"latitude-max":56,"resolution":"2.8*2.8","level":"4","first obs":"16/12/2018","last obs":"present","source":"marineCopernicus","variable":"nh4","link":"https://resources.marine.copernicus.eu/product-detail/IBI_ANALYSISFORECAST_BGC_005_004/INFORMATION","mask":"0","depth-min":0.494,"depth-max":5727.917,"fileType":"NetCDF","time":"","latName":"latitude","longName":"longitude","depthName":"depth","timeName":"time","frequency":"monthly","unitFactor":14,"timeOrigin":"1950-01-01T00:00:00Z","timeUnit":"hours","dataimport":"Copernicus","pretreatment":"Copernicus"},"eastward_Water_current":{"Name":"Copernicus Analysis & Forecast","type":"model","Parameter":"eastward_Water_current","Place":"IBI","motu":"https://nrt.cmems-du.eu/motu-web/Motu","service-id":"IBI_ANALYSISFORECAST_PHY_005_001-TDS","product-id":"cmems_mod_ibi_phy_anfc_0.027deg-3D_P1D-m","longitude-min":-19,"longitude-max":"5","latitude-min":26,"latitude-max":56,"resolution":"2.8*2.8","level":"4","first obs":"16/11/2019","last obs":"present","source":"marineCopernicus","variable":"uo","link":"https://resources.marine.copernicus.eu/product-detail/IBI_ANALYSISFORECAST_PHY_005_001/INFORMATION","mask":"0","depth-min":0.494,"depth-max":5727.917,"fileType":"NetCDF","time":"","latName":"latitude","longName":"longitude","depthName":"depth","timeName":"time","frequency":"daily","unitFactor":86400,"timeOrigin":"1950-01-01T00:00:00Z","timeUnit":"hours","dataimport":"Copernicus","pretreatment":"Copernicus"},"northward_Water_current":{"Name":"Copernicus Analysis & Forecast","type":"model","Parameter":"northward_Water_current","Place":"IBI","motu":"https://nrt.cmems-du.eu/motu-web/Motu","service-id":"IBI_ANALYSISFORECAST_PHY_005_001-TDS","product-id":"cmems_mod_ibi_phy_anfc_0.027deg-3D_P1M-m","longitude-min":-19,"longitude-max":"5","latitude-min":26,"latitude-max":56,"resolution":"2.8*2.8","level":"4","first obs":"19/12/2018","last obs":"present","source":"marineCopernicus","variable":"vo","link":"https://resources.marine.copernicus.eu/product-detail/IBI_ANALYSISFORECAST_PHY_005_001/INFORMATION","mask":"0","depth-min":0.494,"depth-max":5727.917,"fileType":"NetCDF","time":"","latName":"latitude","longName":"longitude","depthName":"depth","timeName":"time","frequency":"monthly","unitFactor":86400,"timeOrigin":"1950-01-01T00:00:00Z","timeUnit":"hours","dataimport":"Copernicus","pretreatment":"Copernicus"},"par":{"Name":"NASA Ocean Colour 2020, filled","type":"sat","Parameter":"par","Place":"IBI","motu":"","service-id":"","product-id":"","longitude-min":"","longitude-max":"","latitude-min":"","latitude-max":"","resolution":"","level":"","first obs":"","last obs":"","source":"","variable":"","link":"","mask":"","depth-min":"","depth-max":"","fileType":"NetCDF","time":"","latName":"lat","longName":"lon","depthName":"depth","timeName":"time","frequency":"monthly","unitFactor":11.574,"timeOrigin":"2020-01-01T00:00:00Z","timeUnit":"days","dataimport":"","pretreatment":""}}},"type":"Algae"}'
scenario_b_hash=`echo -n "$scenario_b_json" | md5sum | head -c 32`
scenario_b_source="/media/share/data/16eae479a19122c3b866dee78a3a442a/_pretreated"
scenario_b_destination="/media/share/data/16eae479a19122c3b866dee78a3a442a/_dataread/${scenario_b_hash}_b"

function handle_arguments {
    local command="$1"

    case $command in
        'build_dataimport')
            build_images_for_dataimport 
            ;;
        'execute_dataimport')
            run_container_for_dataimport "$datasetProperties_json" "$dataimport_destination"
            ;;
        'run_dataimport')
            run_in_interactive_mode "$datasetProperties_json" "$dataimport_destination"
            ;;


        'build_pretreatment')
            build_images_for_pretreatment 
            ;;
        'execute_pretreatment')
            run_pretreatment "$pretreatment_source" "$pretreatment_destination"
            ;;
        'run_pretreatment')
            run_pretreatment_in_interactive_mode "$pretreatment_source" "$pretreatment_destination"
            ;;


        
        'build_dataread')
            build_dataread_image
            ;;
        'execute_dataread')
            run_dataread "$dataread_source" "$dataread_destination" "$modelProperties_json"
            ;;
        'run_dataread')
            run_dataread_in_interactive_mode "$dataread_source" "$dataread_destination" "$modelProperties_json"
            ;;
        

        'build_datareadb')
            build_datareadb_image
            ;;
        'execute_datareadb')
            run_datareadb "$scenario_b_source" "$scenario_b_destination" "$scenario_b_json"
            ;;
        'run_datareadb')
            run_datareadb_in_interactive_mode "$scenario_b_source" "$scenario_b_destination" "$scenario_b_json"
            ;;


        'ls')
            sudo ls -al /var/lib/docker/volumes/$SHARED_VOLUME_NAME/_data/$2
            ;;
        'ls2')
            docker run --rm -i -v=$SHARED_VOLUME_NAME:/media/volume busybox sh
            ;;
        

        'build_posttreatment')
            build_posttreatment_action
            ;;
        'execute_posttreatment')
            run_posttreatment_action "$posttreatment_source"
            ;;


        'create_volumes')
            create_volumes
            ;;
        

        'bash')
            action_bash
            ;;
        *)
            echo 'commands:'
            echo 'build_dataread, execute_dataread'
            echo 'ls, ls2'
            echo 'build_posttreatment'
            echo 'bash'
            ;;
    esac
}

handle_arguments $@

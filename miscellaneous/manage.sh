#!/bin/bash

DIR=$(dirname $0)
echo $DIR
source $DIR/_common.sh
source $DIR/_dataread.sh
source $DIR/_dataread_b.sh
source $DIR/_dataread_shellfish.sh
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

dataread_shellfish_source="$pretreatment_destination"
dataread_shellfish_destination="$dataimport_destination/_dataread_shellfish/$modelProperties_hash"

posttreatment_source="$dataread_destination"


# scenario_b_json=`cat $DIR/__user1_alaria_NWS_11-08-2022.json`
# dataset_id='3a881fd5c4aff18e57c35523c289fe88'

scenario_b_json=`cat $DIR/__user1_alaria_MED_06-09-2022.json`
dataset_id='d066c005901af8a84dd90448ccef90c0'

scenario_b_hash=`echo -n "$scenario_b_json" | md5sum | head -c 32`
scenario_b_source="/media/share/data/$dataset_id/_pretreated"
scenario_b_destination="/media/share/data/$dataset_id/_dataread-b/${scenario_b_hash}_b"


function action_update {
    local container_name=${2:-'ac-dataread/runtime'}
    local container_id=$( docker ps -q -f name=$container_name)
    echo "Update $container_id"
    docker cp global/dataCmd.csv ${container_id}:/opt/dataCmd.csv
}



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
            # TODO deprecated
            build_datareadb_image
            ;;
        'execute_datareadb')
            run_datareadb "$scenario_b_source" "$scenario_b_destination" "$scenario_b_json"
            ;;
        'run_datareadb')
            run_datareadb_in_interactive_mode "$scenario_b_source" "$scenario_b_destination" "$scenario_b_json"
            ;;


        'build_datareadshellfish')
            build_datareadb_image
            ;;
        'execute_datareadshellfish')
            run_datareadb "$scenario_b_source" "$scenario_b_destination" "$scenario_b_json"
            ;;
        'run_datareadshellfish')
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
        
        'run_pretreatment')
            run_posttreatment_in_interactive_mode "$posttreatment_source" "$posttreatment_destination"
            ;;


        'create_volumes')
            create_volumes
            ;;

        'update')
            action_update $@
            ;;
        
        'build_world')
            build_images_for_dataimport
            build_images_for_pretreatment
            build_dataread_image
            build_posttreatment_action
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
            echo 'build_world'
            ;;
    esac
}

handle_arguments $@

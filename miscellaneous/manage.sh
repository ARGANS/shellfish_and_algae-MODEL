#!/bin/bash

DIR=$(dirname $0)
echo $DIR
source $DIR/_common.sh
source $DIR/_dataread.sh
source $DIR/_dataimport.sh

## Properties of scripts used to load datasets:



function build_images_for_posttreatment {
    local dir="./posttreatment"
    local base_image_tag="ac-posttreatment/base"
    local runtime_image_tag="ac-posttreatment/runtime"
    local base_image_dockerfile="./miscellaneous/posttreatmentBase.Dockerfile"
    local runtime_image_dockerfile="./miscellaneous/posttreatmentRuntime.Dockerfile"

    docker build \
        --network host \
        -t $base_image_tag:v1 -t $base_image_tag:latest \
        -f $base_image_dockerfile \
        $dir && \
    docker build \
        --network host \
        --build-arg BASE_IMAGE="$base_image_tag" \
        -t $runtime_image_tag:v1 -t $runtime_image_tag:latest \
        -f $runtime_image_dockerfile \
        $dir
}


function run_container_for_posttreatment {
    local container_name="$1"
    local image_tag="$2"

    #--volume "$(pwd)/posttreatment/src:/opt" \
    # TODO
    # prepare_runtime $container_name $image_tag
    docker run \
        --rm \
        --name $container_name \
        --volume "$SHARED_VOLUME_NAME:/media/share" \
        -e PYTHONDONTWRITEBYTECODE=1 \
        -e SOURCE_DIR='/media/share/results/IBI-2021-0-20/521cac05234a4b7afb01d3a624924d53' \
        -it $image_tag:latest
}

function action_bash {
    local image_tag="ac-posttreatment/runtime"

    docker run \
        --rm \
        --volume "$SHARED_VOLUME_NAME:/media/share" \
        --volume "${HOME}:/media/home" \
        --workdir=/media/share \
        -e PYTHONDONTWRITEBYTECODE=1 \
        --entrypoint=/bin/bash \
        -it $image_tag:latest
}

function handle_arguments {
    local command="$1"

    case $command in
        'build_dataimport')
            build_images_for_dataimport 
            ;;
        'execute_dataimport')
            data=`cat $DIR/__dataimport_input_parameters.json` 
            run_container_for_dataimport "ac-dataimport_run" "ac-import/runtime" "$data"
            ;;
        'run_dataimport')
            data=`cat $DIR/__dataimport_input_parameters.json` 
            run_in_interactive_mode "ac-dataimport_run" "ac-import/runtime" "$data"
            ;;
        'build_dataread')
            build_images_for_model_execution
            ;;
        'execute_dataread')
            run_container_for_model_execution "ac-model_run" "ac-processing/runtime"
            ;;
        'ls')
            sudo ls -al /var/lib/docker/volumes/$SHARED_VOLUME_NAME/_data/$2
            ;;
        'ls2')
            docker run --rm -i -v=$SHARED_VOLUME_NAME:/media/volume busybox sh
            ;;
        'build_posttreatment')
            build_images_for_posttreatment
            ;;
        'execute_posttreatment')
            run_container_for_posttreatment "ac-posttreatment_run" "ac-posttreatment/runtime"
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

handle_arguments "$1" "$2"

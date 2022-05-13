#!/bin/bash
SHARED_VOLUME_NAME='ac_share'

## Properties of scripts used to load datasets:

function build_images_for_model_execution {
    local dir="./dataread"
    local base_image_tag="ac-processing/base"
    local runtime_image_tag="ac-processing/runtime"
    local base_image_dockerfile="./miscellaneous/pythonBase.Dockerfile"
    local runtime_image_dockerfile="./miscellaneous/pythonRuntime.Dockerfile"

    ### for GDAL: --build-arg WITH_GDAL="true" \
    docker build \
        --network host \
        --build-arg REQUIREMENTS_PATH="./requirements.txt" \
        --build-arg WITH_R="true" \
        --build-arg WITH_NETCDF="true" \
        -t $base_image_tag:v1 -t $base_image_tag:latest \
        -f $base_image_dockerfile \
        $dir && \
    docker build \
        --network host \
        --build-arg BASE_IMAGE="$base_image_tag" \
        -t $runtime_image_tag:v2 -t $runtime_image_tag:latest \
        -f $runtime_image_dockerfile \
        $dir
}

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

function prepare_runtime {
    local container_name="$1"
    local image_tag="$2"

    local container_id=$( docker ps -q -f name=$container_name )
    local image_id=$( docker images -q $image_tag )

    if [[ ! -z "$container_id" ]]; then
        docker stop $container_name || true
    fi

    if [[ -z "$image_id" ]]; then
        build_images_for_model_execution
    fi

    docker volume create --name $SHARED_VOLUME_NAME
}

function run_container_for_model_execution {
    local container_name="$1"
    local image_tag="$2"

    prepare_runtime $container_name $image_tag

    # All model properties received from the application
    data=`cat macroalgae/macroalgae_model_parameters_input.json`

    # use -d to start a container in detached mode
    # use --entrypoint=/bin/bash \ to override the command
    docker run \
        --rm \
        --name $container_name \
        --volume "$SHARED_VOLUME_NAME:/media/share" \
        -e DATASET_ID="IBI-2021-3-18" \
        -e TASK_ID="abcdef" \
        -e PARAMETERS_JSON="$data" \
        -e PYTHONDONTWRITEBYTECODE=1 \
        -it $image_tag:latest
}

function run_container_for_posttreatment {
    local container_name="$1"
    local image_tag="$2"

    #--volume "$(pwd)/posttreatment/src:/opt" \
        
    prepare_runtime $container_name $image_tag
    docker run \
        --rm \
        --name $container_name \
        --volume "$SHARED_VOLUME_NAME:/media/share" \
        -e PYTHONDONTWRITEBYTECODE=1 \
        -e SOURCE_DIR='/media/share/results/IBI-2021-0-20/521cac05234a4b7afb01d3a624924d53' \
        -it $image_tag:latest
}

function handle_arguments {
    local command="$1"

    case $command in
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
        *)
            echo 'commands:'
            echo 'build_dataread, execute_dataread'
            echo 'ls, ls2'
            echo 'build_posttreatment'
            ;;
    esac
}

handle_arguments "$1" "$2"

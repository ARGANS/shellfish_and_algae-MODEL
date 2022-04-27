#!/bin/bash
SHARED_VOLUME_NAME='ac_share'
SHARED_VOLUME_NAME='ac_share_year'

## Properties of scripts used to load datasets:

function build_images_for_model_execution {
    local dir="./dataread"
    local base_image_tag="ac-processing/base"
    local runtime_image_tag="ac-processing/runtime"
    local base_image_dockerfile="$dir/pythonBase.Dockerfile"
    local runtime_image_dockerfile="$dir/pythonRuntime.Dockerfile"

    ### for GDAL: --build-arg WITH_GDAL="true" \
    docker build \
        --network host \
        --build-arg REQUIREMENTS_PATH="./requirements.txt" \
        --build-arg WITH_R="true" \
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

function run_container_for_model_execution {
    local container_name="$1"
    local container_id=$( docker ps -q -f name=$container_name )
            
    if [[ ! -z "$container_id" ]]; then
        docker stop $container_name || true
    fi
    local image_tag="ac-processing/runtime"
    local image_id=$( docker images -q $image_tag )

    if [[ -z "$image_id" ]]; then
        build_images_for_model_execution
    fi

    docker volume create --name $SHARED_VOLUME_NAME

    # use -d to start a container in detached mode
    # use --entrypoint=/bin/bash \ to override the command
    docker run \
        --rm \
        --name $container_name \
        --volume "$SHARED_VOLUME_NAME":/media/share \
        -e PYTHONDONTWRITEBYTECODE=1 \
        --entrypoint=/bin/bash \
        -it $image_tag:latest
}

function handle_arguments {
    local command="$1"

    case $command in
        'build')
            build_images_for_model_execution
            ;;
        'execute')
            run_container_for_model_execution "ac-model_run"
            ;;
        'ls')
            sudo ls -al /var/lib/docker/volumes/$SHARED_VOLUME_NAME/_data/$2
            ;;
        'ls2')
            docker run --rm -i -v=$SHARED_VOLUME_NAME:/media/volume busybox sh
            ;;
        *)
            echo 'todo';;
    esac
}

handle_arguments "$1" "$2"

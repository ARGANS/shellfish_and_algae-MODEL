SHARED_VOLUME_NAME='ac_share'

## Properties of scripts used to load datasets:
zone='IBI'
year=2021
deepthmin=0
deepthmax=20

output_dir="/media/share/data/$zone-$year-$deepthmin-$deepthmax/"

function build_images_for_data_import {
    local dir="./dataimport"
    local base_image_tag="ac-import/base"
    local runtime_image_tag="ac-import/runtime"
    local base_image_dockerfile="./miscellaneous/pythonBase.Dockerfile"
    local runtime_image_dockerfile="./miscellaneous/dataimportRuntime.Dockerfile"

    docker build \
        --network host \
        --build-arg REQUIREMENTS_PATH="./src/requirements.txt" \
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

function run_container_for_data_import {
    local container_name="$1"
    local image_tag="$2"
    local data="$3"
    prepare_runtime "$container_name" "$image_tag"

    # use -d to start a container in detached mode
    # use --entrypoint=/bin/bash \ to override the command
    docker run \
        --rm \
        --name $container_name \
        --volume "$SHARED_VOLUME_NAME":/media/share \
        -e AC_OUTPUT_DIR="$output_dir" \
        -e parameters_json="$data" \
        -e PYTHONDONTWRITEBYTECODE=1 \
        -it $image_tag:latest
}

function run_in_interactive_mode {
    local container_name="$1"
    local image_tag="$2"
    prepare_runtime "$container_name" "$image_tag"

    docker run \
        --rm \
        --name $container_name \
        --volume "$SHARED_VOLUME_NAME":/media/share \
        -e AC_OUTPUT_DIR="$output_dir" \
        -e PYTHONDONTWRITEBYTECODE=1 \
        --entrypoint=/bin/bash \
        -it $image_tag:latest
}

function run {
    local command="$1"

    case $command in
        'build')
            build_images_for_data_import 
            ;;
        'execute')
            run_container_for_data_import "ac-dataimport_run" "ac-import/runtime" "{\"zone\":\"$zone\",\"depth_min\":$deepthmin,\"depth_max\":$deepthmax,\"year\":$year}"
            ;;
        'run')
            run_in_interactive_mode "ac-dataimport_run" "ac-import/runtime"
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

run "$1" "$2"

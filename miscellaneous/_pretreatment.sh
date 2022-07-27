#!/bin/bash
PRETREATMENT_IMAGE='ac-pretreatment/runtime'
PRETREATMENT_CONTAINER='ac-pretreatment_run'

function build_images_for_pretreatment {
    local dir="./pretreatment/src"
    local base_image_tag="ac-pretreatment/base"
    local base_image_dockerfile="./miscellaneous/pythonBase.Dockerfile"
    local runtime_image_dockerfile="./miscellaneous/pythonRuntime.Dockerfile"

    docker build \
        --network host \
        --build-arg WITH_NETCDF="true" \
        -t $base_image_tag:v1 -t $base_image_tag:latest \
        -f $base_image_dockerfile \
        $dir && \
    docker build \
        --network host \
        --build-arg BASE_IMAGE="$base_image_tag" \
        -t $PRETREATMENT_IMAGE:v1 -t $PRETREATMENT_IMAGE:latest \
        -f $runtime_image_dockerfile \
        $dir
}

function run_pretreatment {
    stop_existed_container $PRETREATMENT_CONTAINER
    create_volumes
    
    # use -d to start a container in detached mode
    # use --entrypoint=/bin/bash \ to override the command
    docker run \
        --rm \
        --name $PRETREATMENT_CONTAINER \
        --volume "$SHARED_VOLUME_NAME":/media/share \
        --volume $(pwd)/global:/media/global \
        -e INPUT_SOURCE="$1" \
        -e INPUT_DESTINATION="$2" \
        -e PYTHONDONTWRITEBYTECODE=1 \
        -it $PRETREATMENT_IMAGE:latest
}

function run_pretreatment_in_interactive_mode {
    stop_existed_container $$PRETREATMENT_CONTAINER
    create_volumes

    docker run \
        --rm \
        --name $PRETREATMENT_CONTAINER \
        --volume "$SHARED_VOLUME_NAME":/media/share \
        --volume $(pwd)/global:/media/global \
        --volume $(pwd)/pretreatment/src/:/opt/ \
        -e INPUT_SOURCE="$1" \
        -e INPUT_DESTINATION="$2" \
        -e PYTHONDONTWRITEBYTECODE=1 \
        --entrypoint=/bin/bash \
        -it $PRETREATMENT_IMAGE:latest
}

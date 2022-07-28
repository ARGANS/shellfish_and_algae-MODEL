#!/bin/bash
DATAREAD_IMAGE='ac-dataread/runtime'
DATAREAD_CONTAINER='ac-dataread_run'

function build_dataread_image {
    local dir="./dataread"
    local base_image_tag="ac-dataread/base"
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
        -t $DATAREAD_IMAGE:v1 -t $DATAREAD_IMAGE:latest \
        -f $runtime_image_dockerfile \
        $dir
}

function run_dataread {
    stop_existed_container DATAREAD_CONTAINER
    create_volumes

    # use -d to start a container in detached mode
    # use --entrypoint=/bin/bash \ to override the command
    docker run \
        --rm \
        --name $DATAREAD_CONTAINER \
        --volume "$SHARED_VOLUME_NAME:/media/share" \
        --volume $(pwd)/global:/media/global \
        -e INPUT_SOURCE="$1" \
        -e INPUT_DESTINATION="$2" \
        -e INPUT_MODEL_PROPERTIES_JSON="$3" \
        -e PYTHONDONTWRITEBYTECODE=1 \
        -it $DATAREAD_IMAGE:latest
}

function run_dataread_in_interactive_mode {
    docker run \
        --rm \
        --name $DATAIMPORT_CONTAINER \
        --volume "$SHARED_VOLUME_NAME":/media/share \
        --volume $(pwd)/global:/media/global \
        --volume $(pwd)/dataread:/opt/ \
        -e INPUT_SOURCE="$1" \
        -e INPUT_DESTINATION="$2" \
        -e INPUT_MODEL_PROPERTIES_JSON="$3" \
        -e PYTHONDONTWRITEBYTECODE=1 \
        --entrypoint=/bin/bash \
        -it $DATAREADB_IMAGE:latest
}

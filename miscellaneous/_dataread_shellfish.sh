#!/bin/bash
DATAREAD_SHELLFISH_IMAGE='ac-dataread_shellfish/runtime'
DATAREAD_SHELLFISH_CONTAINER='ac-dataread_shellfish_run'

function build_dataread_shellfish_image {
    local dir="./dataread_shellfish"
    local base_image_tag="ac-dataread_shellfish/base"
    local base_image_dockerfile="./miscellaneous/pythonBase.Dockerfile"
    local runtime_image_dockerfile="./miscellaneous/pythonRuntime.Dockerfile"

    ### for GDAL: --build-arg WITH_GDAL="true" \
    docker build \
        --network host \
        --build-arg REQUIREMENTS_PATH="./dataread_shellfish/requirements.txt" \
        --build-arg WITH_R="true" \
        --build-arg WITH_NETCDF="true" \
        -t $base_image_tag:v1 -t $base_image_tag:latest \
        -f $base_image_dockerfile \
        $dir && \
    docker build \
        --network host \
        --build-arg BASE_IMAGE="$base_image_tag" \
        -t $DATAREAD_SHELLFISH_IMAGE:v1 -t $DATAREAD_SHELLFISH_IMAGE:latest \
        -f $runtime_image_dockerfile \
        $dir
}

function run_dataread_shellfish {
    stop_existed_container DATAREAD_SHELLFISH_CONTAINER
    create_volumes

    # use -d to start a container in detached mode
    # use --entrypoint=/bin/bash \ to override the command
    docker run \
        --rm \
        --name $DATAREAD_SHELLFISH_CONTAINER \
        --volume "$SHARED_VOLUME_NAME:/media/share" \
        --volume "$GLOBAL_VOLUME_NAME":/media/global \
        -e INPUT_SOURCE="$1" \
        -e INPUT_DESTINATION="$2" \
        -e INPUT_MODEL_PROPERTIES_JSON="$3" \
        -e PYTHONDONTWRITEBYTECODE=1 \
        -it $DATAREAD_IMAGE:latest
}

function run_dataread_shellfish_in_interactive_mode {
    docker run \
        --rm \
        --name $DATAIMPORT_CONTAINER \
        --volume "$SHARED_VOLUME_NAME":/media/share \
        --volume "$GLOBAL_VOLUME_NAME":/media/global \
        --volume $(pwd)/dataread:/opt/ \
        -e INPUT_SOURCE="$1" \
        -e INPUT_DESTINATION="$2" \
        -e INPUT_MODEL_PROPERTIES_JSON="$3" \
        -e PYTHONDONTWRITEBYTECODE=1 \
        --entrypoint=/bin/bash \
        -it $DATAREADB_IMAGE:latest
}

#!/bin/bash
DATAREADB_IMAGE='ac-datareadb/runtime'
DATAREADB_CONTAINER='ac-datareadb_run'

function build_datareadb_image {
    local dir="./"
    local base_image_tag="ac-datareadb/base"
    local base_image_dockerfile="./miscellaneous/pythonBase.Dockerfile"
    local runtime_image_dockerfile="./miscellaneous/scenario_b.Dockerfile"

    ### for GDAL: --build-arg WITH_GDAL="true" \
            # --build-arg WITH_NETCDF="true" \
    docker build \
        --network host \
        --build-arg REQUIREMENTS_PATH="./scenario_b/requirements.txt" \
        --build-arg WITH_R="true" \
        --build-arg WITH_GDAL="true" \
        -t $base_image_tag:v1 -t $base_image_tag:latest \
        -f $base_image_dockerfile \
        $dir && \
    docker build \
        --network host \
        --build-arg BASE_IMAGE="$base_image_tag" \
        -t $DATAREADB_IMAGE:v1 -t $DATAREADB_IMAGE:latest \
        -f $runtime_image_dockerfile \
        $dir
}

function run_datareadb {
    stop_existed_container DATAREADB_CONTAINER
    create_volumes

    # use -d to start a container in detached mode
    # use --entrypoint=/bin/bash \ to override the command
    docker run \
        --rm \
        --name $DATAREADB_CONTAINER \
        --volume "$SHARED_VOLUME_NAME:/media/share" \
        --volume $(pwd)/global:/media/global \
        -e INPUT_SOURCE="$1" \
        -e INPUT_DESTINATION="$2" \
        -e INPUT_MODEL_PROPERTIES_JSON="$3" \
        -e PYTHONDONTWRITEBYTECODE=1 \
        -it $DATAREADB_IMAGE:latest
}

function run_datareadb_in_interactive_mode {
    docker run \
        --rm \
        --name $DATAIMPORT_CONTAINER \
        --volume "$SHARED_VOLUME_NAME":/media/share \
        --volume $(pwd)/global:/media/global \
        --volume $(pwd)/dataread/src:/opt/dataread \
        --volume $(pwd)/advectionPrototype:/opt/advectionPrototype \
        --volume $(pwd)/scenario_b/start.py:/opt/start.py \
        --volume $(pwd)/scenario_b/start.sh:/opt/start.sh \
        -e INPUT_SOURCE="$1" \
        -e INPUT_DESTINATION="$2" \
        -e INPUT_MODEL_PROPERTIES_JSON="$3" \
        -e PYTHONDONTWRITEBYTECODE=1 \
        --entrypoint=/bin/bash \
        -it $DATAREADB_IMAGE:latest
}

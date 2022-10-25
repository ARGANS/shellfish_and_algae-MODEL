#!/bin/bash
DATAREADB_IMAGE='ac-farmrepartition/runtime'
DATAREADB_CONTAINER='ac-farmrepartition_run'

function build_farmrepartition_image {
    local dir="./"
    local base_image_tag="ac-datareadb/base"
    local base_image_dockerfile="./miscellaneous/pythonBase.Dockerfile"
    local runtime_image_dockerfile="./miscellaneous/farmrepartition.Dockerfile"

    docker build \
        --network host \
        --build-arg REQUIREMENTS_PATH="./farm_repartition/requirements.txt" \
        --build-arg WITH_NETCDF="true" \
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

function run_farmrepartition {
    stop_existed_container DATAREADB_CONTAINER
    create_volumes

    # use -d to start a container in detached mode
    # use --entrypoint=/bin/bash \ to override the command
    docker run \
        --rm \
        --name $DATAREADB_CONTAINER \
        --volume "$SHARED_VOLUME_NAME:/media/share" \
        --volume "$GLOBAL_VOLUME_NAME":/media/global \
        -e INPUT_SOURCE="$1" \
        -e INPUT_DESTINATION="$2" \
        -e INPUT_MODEL_PROPERTIES_JSON="$3" \
        -e PYTHONDONTWRITEBYTECODE=1 \
        -it $DATAREADB_IMAGE:latest
}

function run_farmrepartition_in_interactive_mode {
    docker run \
        --rm \
        --name $DATAIMPORT_CONTAINER \
        --volume "$SHARED_VOLUME_NAME":/media/share \
        --volume "$GLOBAL_VOLUME_NAME":/media/global \
        --volume $(pwd)/farm_repartition:/opt \
        --memory=4g \
		--memory-swap=8g \
        --label 'task.model.id=012' \
        -e INPUT_SOURCE="$1" \
        -e INPUT_DESTINATION="$2" \
        -e INPUT_MODEL_PROPERTIES_JSON="$3" \
        -e PYTHONDONTWRITEBYTECODE=1 \
        --entrypoint=/bin/bash \
        -it $DATAREADB_IMAGE:latest
}

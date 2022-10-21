DATAIMPORT_IMAGE='ac-import/runtime'
DATAIMPORT_CONTAINER='ac-dataimport_run'


function build_images_for_dataimport {
    local dir="./dataimport"
    local base_image_tag="ac-import/base"
    local runtime_image_tag="ac-import/runtime"
    local base_image_dockerfile="./miscellaneous/pythonBase.Dockerfile"
    local runtime_image_dockerfile="./miscellaneous/dataimportRuntime.Dockerfile"

    docker build \
        --network host \
        --build-arg REQUIREMENTS_PATH="./src/requirements.txt" \
        --build-arg WITH_JQ="true" \
        -t $base_image_tag:v1 -t $base_image_tag:latest \
        -f $base_image_dockerfile \
        $dir && \
    docker build \
        --network host \
        --build-arg BASE_IMAGE="$base_image_tag" \
        -t $DATAIMPORT_IMAGE:v1 -t $DATAIMPORT_IMAGE:latest \
        -f $runtime_image_dockerfile \
        $dir
}

function run_container_for_dataimport {
    stop_existed_container $DATAIMPORT_CONTAINER
    create_volumes

    # use -d to start a container in detached mode
    # use --entrypoint=/bin/bash \ to override the command
    docker run \
        --rm \
        --name $DATAIMPORT_CONTAINER \
        --volume "$SHARED_VOLUME_NAME":/media/share \
        --volume "$GLOBAL_VOLUME_NAME":/media/global \
        -e INPUT_DESTINATION="$2" \
        -e INPUT_PARAMETERS="$1" \
        -e MOTU_LOGIN="mjaouen" \
        -e MOTU_PASSWORD="Azerty123456" \
        -e PYTHONDONTWRITEBYTECODE=1 \
        -it $DATAIMPORT_IMAGE:latest
}

function run_in_interactive_mode {
    stop_existed_container $DATAIMPORT_CONTAINER
    create_volumes

    docker run \
        --rm \
        --name $DATAIMPORT_CONTAINER \
        --volume "$SHARED_VOLUME_NAME":/media/share \
        --volume "$GLOBAL_VOLUME_NAME":/media/global \
        --volume $(pwd)/dataimport/src/start.py:/opt/start.py \
        --volume $(pwd)/dataimport/src/general.py:/opt/general.py \
        -e INPUT_DESTINATION="$2" \
        -e INPUT_PARAMETERS="$1" \
        -e PYTHONDONTWRITEBYTECODE=1 \
        -e MOTU_LOGIN="mjaouen" \
        -e MOTU_PASSWORD="Azerty123456" \
        --entrypoint=/bin/bash \
        -it $DATAIMPORT_IMAGE:latest
}

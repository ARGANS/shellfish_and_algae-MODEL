POSTTREATMENT_IMAGE='ac-posttreatment/runtime'
POSTTREATMENT_CONTAINER='ac-posttreatment_run'


function build_posttreatment_action {
    local dir="./posttreatment"
    local base_image_tag="ac-posttreatment/base"
    local base_image_dockerfile="./miscellaneous/pythonBase.Dockerfile"
    local runtime_image_dockerfile="./miscellaneous/posttreatmentRuntime.Dockerfile"

    docker build \
        --network host \
        --build-arg REQUIREMENTS_PATH="./src/requirements.txt" \
        --build-arg WITH_GDAL="true" \
        --build-arg WITH_JQ="true" \
        -t $base_image_tag:v1 -t $base_image_tag:latest \
        -f $base_image_dockerfile \
        $dir && \
    docker build \
        --network host \
        --build-arg BASE_IMAGE="$base_image_tag" \
        -t $POSTTREATMENT_IMAGE:v1 -t $POSTTREATMENT_IMAGE:latest \
        -f $runtime_image_dockerfile \
        $dir
}


function run_posttreatment_action {
    docker run \
        --rm \
        --name POSTTREATMENT_CONTAINER \
        --volume "$SHARED_VOLUME_NAME:/media/share" \
        -e PYTHONDONTWRITEBYTECODE=1 \
        -e SOURCE_DIR="$1" \
        -it $POSTTREATMENT_IMAGE:latest
}

function run_posttreatment_in_interactive_mode {
    stop_existed_container $$POSTTREATMENT_CONTAINER
    create_volumes

    docker run \
        --rm \
        --name $POSTTREATMENT_CONTAINER \
        --volume "$SHARED_VOLUME_NAME":/media/share \
        --volume "$GLOBAL_VOLUME_NAME":/media/global \
        --volume $(pwd)/posttreatment/src/:/opt/ \
        -e INPUT_SOURCE="$1" \
        -e INPUT_DESTINATION="$2" \
        -e PYTHONDONTWRITEBYTECODE=1 \
        --entrypoint=/bin/bash \
        -it $POSTTREATMENT_IMAGE:latest
}
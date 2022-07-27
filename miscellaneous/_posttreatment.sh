POSTTREATMENT_IMAGE='ac-posttreatment/runtime'
POSTTREATMENT_CONTAINER='ac-posttreatment_run'


function build_posttreatment_action {
    local dir="./posttreatment"
    local base_image_tag="ac-posttreatment/base"
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

## Properties of scripts used to load datasets:
# zone='IBI'
# year=2021
# deepthmin=0
# deepthmax=20

# 

function build_images_for_dataimport {
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

function run_container_for_dataimport {
    local container_name="$1"
    local image_tag="$2"
    local data="$3"
    stop_existed_container $container_name
    create_volumes
    local hash=`echo -n "$data" | md5sum | head -c 32`
    local output_dir="/media/share/data/$hash/"

    # use -d to start a container in detached mode
    # use --entrypoint=/bin/bash \ to override the command
    docker run \
        --rm \
        --name $container_name \
        --volume "$SHARED_VOLUME_NAME":/media/share \
        --volume $(pwd)/global:/media/global \
        -e INPUT_DESTINATION="$output_dir" \
        -e INPUT_PARAMETERS="$data" \
        -e PYTHONDONTWRITEBYTECODE=1 \
        -it $image_tag:latest
}

function run_in_interactive_mode {
    local container_name="$1"
    local image_tag="$2"
    local data="$3"
    stop_existed_container $container_name
    create_volumes
    local hash=`echo -n "$data" | md5sum | head -c 32`
    local output_dir="/media/share/data/$hash/"


    docker run \
        --rm \
        --name $container_name \
        --volume "$SHARED_VOLUME_NAME":/media/share \
        --volume $(pwd)/global:/media/global \
        -e INPUT_DESTINATION="$output_dir" \
        -e INPUT_PARAMETERS="$data" \
        -e PYTHONDONTWRITEBYTECODE=1 \
        --entrypoint=/bin/bash \
        -it $image_tag:latest
}
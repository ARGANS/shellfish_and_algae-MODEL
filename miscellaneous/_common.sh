SHARED_VOLUME_NAME='ac_share'

function stop_existed_container {
    local container_name="$1"
    local container_id=$( docker ps -q -f name=$container_name )
    

    if [[ ! -z "$container_id" ]]; then
        docker stop $container_name || true
    fi
}

function create_volumes {
    docker volume create --name $SHARED_VOLUME_NAME
}

# DEPRECATED
# function prepare_runtime {
#     local image_tag="$1"
#     local image_id=$( docker images -q $image_tag )
#     if [[ -z "$image_id" ]]; then
#         build_images_for_model_execution
#     fi
# }
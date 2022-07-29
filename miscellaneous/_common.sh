SHARED_VOLUME_NAME='ac_share'
GLOBAL_VOLUME_NAME='ac_global'

function stop_existed_container {
    local container_name="$1"
    local container_id=$( docker ps -q -f name=$container_name )
    

    if [[ ! -z "$container_id" ]]; then
        docker stop $container_name || true
    fi
}

function create_volumes {
    echo "::" $(pwd)/global
    docker volume create --name $SHARED_VOLUME_NAME
    docker volume create \
        --driver local \
        --opt type=none \
        --opt device=$(pwd)/global \
        --opt o=bind \
        $GLOBAL_VOLUME_NAME
}

#  docker volume create -d local-persist -o mountpoint=/mnt/ --name=extra-addons


# DEPRECATED
# function prepare_runtime {
#     local image_tag="$1"
#     local image_id=$( docker images -q $image_tag )
#     if [[ -z "$image_id" ]]; then
#         build_images_for_model_execution
#     fi
# }

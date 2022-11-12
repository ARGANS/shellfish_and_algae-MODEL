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
    local destination=$(pwd)/global/maps
    mkdir -p $destination

    local host='https://213.166.43.12'
    local login='bot'
    local password='bot_ma5ter'
    local cookie_file='cookies.txt'

    local redirectURL=$(curl -sik -o /dev/null -F "username=$login" -F "password=$password" -H 'Cache-Control: no-cache'  -X POST --write-out '%{redirect_url}' -c $cookie_file $host/api/v1/auth/login)
    # curl -k -b $cookie_file $host/api/v1/auth/whoami
    curl -sk -b $cookie_file $host/api/v2/file?path=/media/global/maps/Bathy.TIF -o $destination/Bathy.TIF
    curl -sk -b $cookie_file $host/api/v2/file?path=/media/global/maps/zee_europe.tif -o $destination/zee_europe.tif

    docker volume create --name $SHARED_VOLUME_NAME
    docker volume create \
        --driver local \
        --opt type=none \
        --opt device=$(pwd)/global \
        --opt o=bind \
        $GLOBAL_VOLUME_NAME
}

#  docker volume create -d local-persist -o mountpoint=/mnt/ --name=extra-addons

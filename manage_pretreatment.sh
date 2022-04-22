DIR='pretreatment'
DD_TAG='aquaculture/pretreatment'
DOCKERFILE='pretreatment.Dockerfile'
CONTAINER_NAME='NCOcontainer'
# SHARED_VOLUME_NAME='/profils/mjaouen/share'
SHARED_VOLUME_NAME='ac_share'

function run {
    local command="$1"

    case $command in
        'build')
            docker build \
                --network host \
                -t $DD_TAG-base:v1 -t $DD_TAG-base:latest \
                -f $DIR/base.$DOCKERFILE \
                $DIR
            docker build \
                --network host \
                -t $DD_TAG:v1 -t $DD_TAG:latest \
                -f $DIR/runtime.$DOCKERFILE \
                $DIR
            ;;
        'run')
            container_id=$( docker ps -q -f name=$CONTAINER_NAME )
            
            if [[ ! -z "$container_id" ]]; then
                docker stop $CONTAINER_NAME || true
            fi

            image_id=$( docker images -q $DD_TAG )

            if [[ -z "$image_id" ]]; then
                run "build"
            fi

            # --volume "/profils/mjaouen/share/data":/media/share/data \
            # --volume "$HOME/share/data_merged":/media/share/data_merged \
            # --volume "$HOME/share/results":/media/share/results \
            docker run \
                --rm \
                --name $CONTAINER_NAME \
                --volume "$SHARED_VOLUME_NAME":/media/share \
                --entrypoint=/bin/bash \
                -it $DD_TAG:latest
            ;;
        *)
            echo 'todo';;
    esac
}

run "$1"

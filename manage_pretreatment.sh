DIR='pretreatment'
DD_TAG='aquaculture/pretreatment'
DOCKERFILE="./$DIR/pretreatment.Dockerfile"
CONTAINER_NAME='NCOcontainer'

function run {
    local command="$1"

    case $command in
        'build')
            docker build \
                --network host \
                -t $DD_TAG:v1 -t $DD_TAG:latest \
                -f $DOCKERFILE ./$DIR
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

            docker run \
                --rm \
                --name $CONTAINER_NAME \
                --volume "/profils/mjaouen/share/data":/media/share/data \
                --volume "$HOME/share/data_merged":/media/share/data_merged \
                --volume "$HOME/share/results":/media/share/results \
                --entrypoint '/usr/local/bin/python' \
                -dit $DD_TAG:latest

            docker exec -it $CONTAINER_NAME bash 
            ;;
        *)
            echo 'todo';;
    esac
}

run "$1"

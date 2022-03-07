DIR='dataimport'
DD_TAG='aquaculture/datadownload'
DOCKERFILE="./$DIR/dataimport.Dockerfile"
CONTAINER_NAME='aquacontainer'

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
            # docker stop $CONTAINER_NAME || true
            if [[ -z "$container_id" ]]; then
                image_id=$( docker images -q $DD_TAG )

                if [[ -z "$image_id" ]]; then
                    run "build"
                fi

                docker run \
                    --rm \
                    --name $CONTAINER_NAME \
                    --volume "$PWD/share":/media/share \
                    --volume "$PWD/$DIR/src/general.py":/opt/general.py \
                    --volume "$PWD/$DIR/src/main.py":/opt/main.py \
                    --entrypoint '/usr/local/bin/python' \
                    -dit $DD_TAG:latest
            fi
            
            docker exec -it $CONTAINER_NAME bash 
            ;;
        *)
            echo 'todo';;
    esac
}

run "$1"

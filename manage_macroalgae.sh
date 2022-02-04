DIR='macroalgae'
MA_TAG='aquaculture/macroalgae_model1'
DOCKERFILE="./miscellaneous/macroalgae.dev.Dockerfile"
CONTAINER_NAME='macroalgae_container1'

function run {
    local command="$1"

    case $command in
        'build')
            docker build \
                -t $MA_TAG:v1 -t $MA_TAG:latest \
                -f $DOCKERFILE ./$DIR
            ;;
        'run')
            # we can get the container id if the container is running
            container_id=$( docker ps -q -f name=$CONTAINER_NAME )

            if [[ -z "$container_id" ]]; then
                image_id=$( docker images -q $MA_TAG )

                if [[ -z "$image_id" ]]; then
                    run "build"
                fi

                docker run \
                    --rm \
                    --name $CONTAINER_NAME \
                    --volume "$PWD/share":/media/share \
                    --entrypoint '/bin/bash' \
                    -dit $MA_TAG:latest
            fi
            
            docker exec -it $CONTAINER_NAME bash 
            ;;
        *)
            echo 'todo';;
    esac
}

run "$1"

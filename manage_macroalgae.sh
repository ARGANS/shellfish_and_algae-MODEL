DIR='macroalgae'
MA_TAG='aquaculture/macroalgae_model1'
DOCKERFILE="./miscellaneous/macroalgae.dev.Dockerfile"
CONTAINER_NAME='macroalgae_container1'

function start_container {
    # we can get the container id if the container is running
    local container_id=$( docker ps -q -f name=$CONTAINER_NAME )

    if [[ -z "$container_id" ]]; then
        local image_id=$( docker images -q $MA_TAG )

        if [[ -z "$image_id" ]]; then
            run "build"
        fi

        docker run \
            --rm \
            --name $CONTAINER_NAME \
            --volume "$HOME/share/data_merged":/media/share/data_merged \
            --volume "$HOME/share/results":/media/share/results \
            --entrypoint '/bin/bash' \
            -dit $MA_TAG:latest

        run "update"
    fi
}

function run {
    local command="$1"

    case $command in
        'build')
            docker build \
                --network host \
                -t $MA_TAG:v1 -t $MA_TAG:latest \
                -f $DOCKERFILE ./
                #-f $DOCKERFILE ./$DIR
            ;;
        'shell')
            start_container
            docker exec -it $CONTAINER_NAME bash 
            ;;
        'execute')
            start_container
            docker exec -itd $CONTAINER_NAME bash -c "python make_runs.py"
            #docker cp $CONTAINER_NAME:/opt/Rplots.pdf ~/Documents
            ;;
        'stop')
            local container_id=$( docker ps -q -f name=$CONTAINER_NAME )
            if [[ ! -z "$container_id" ]]; then
                docker stop $CONTAINER_NAME || true
            fi
            ;;
        'update')
            local container_id=$( docker ps -q -f name=$CONTAINER_NAME )
            docker cp dataread/. ${container_id}:/opt/
            docker cp macroalgae/. ${container_id}:/opt/
            docker cp dataimport/src/dataCmd.csv ${container_id}:/opt/dataCmd.csv
            ;;
        *)
            echo 'todo help';;
    esac
}

run "$1"

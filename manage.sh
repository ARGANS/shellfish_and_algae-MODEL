DIR='dataimport'
TAG='aquaculture/datadownload'
DOCKERFILE="./$DIR/dataimport.Dockerfile"
CONTAINERNAME='aquacontainer'

# 

function run {
    local command="$1"

    echo $TAG
    if [[ "$command" == "build" ]]; then
        docker build \
            -t $TAG:v1 -t $TAG:latest \
            -f $DOCKERFILE ./$DIR
    elif [[ "$command" == "run" ]]; then
        docker stop $CONTAINERNAME || true
        docker run \
            --rm \
            --name $CONTAINERNAME \
            --volume "$PWD/share":/media/share \
            --volume "$PWD/$DIR/src/general.py":/opt/general.py \
            --entrypoint '/usr/local/bin/python' \
            -dit $TAG:latest
        docker exec -it $CONTAINERNAME bash 
    fi
    echo ""
}

run "$1"

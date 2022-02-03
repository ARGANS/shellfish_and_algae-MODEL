DIR='dataimport'
TAG='aquaculture/datadownload'
DOCKERFILE="./$DIR/dataimport.Dockerfile"
CONTAINERNAME='aquacontainer'

function run {
    local command="$1"

    case $command in
        'build')
            docker build \
                -t $TAG:v1 -t $TAG:latest \
                -f $DOCKERFILE ./$DIR
            ;;
        'run')
            docker stop $CONTAINERNAME || true
            docker run \
                --rm \
                --name $CONTAINERNAME \
                --volume "$PWD/share":/media/share \
                --volume "$PWD/$DIR/src/general.py":/opt/general.py \
                --volume "$PWD/$DIR/src/main.py":/opt/main.py \
                --entrypoint '/usr/local/bin/python' \
                -dit $TAG:latest
            docker exec -it $CONTAINERNAME bash 
            ;;
        *)
            echo 'todo';;
    esac
}

run "$1"

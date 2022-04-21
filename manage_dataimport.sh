DIR='dataimport'
DD_TAG='aquaculture/datadownload'
DOCKERFILE="./$DIR/dataimport.Dockerfile"
CONTAINER_NAME='aquacontainer'

# output_dir='I:/work-he/apps/safi/data/IBI/'
output_dir='/media/share/data/IBI/'

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
                --add-host=ftp.hermes.acri.fr:5.252.148.37 \
                --volume "$PWD/share":/media/share \
                --volume "$PWD/$DIR/src/.":"/opt/." \
                --entrypoint '/usr/local/bin/python' \
                -e OUTPUT_DIR="$output_dir" \
                -dit $DD_TAG:latest

            docker exec -it $CONTAINER_NAME bash 
            ;;
        *)
            echo 'todo';;
    esac
}

run "$1"

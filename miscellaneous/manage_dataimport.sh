DIR='./dataimport'
DD_TAG='aquaculture/dataimport'
CONTAINER_NAME='ac-dataimport'

SHARED_VOLUME_NAME='ac_share'

zone='IBI'
year=2020
# output_dir='I:/work-he/apps/safi/data/IBI/'
output_dir='/media/share/data/IBI/'
deepthmin=0
deepthmax=20

function run {
    local command="$1"

    case $command in
        'build')
            docker build \
                --network host \
                -t aquaculture/base:v1 -t aquaculture/base:latest \
                -f $DIR/base.Dockerfile \
                $DIR
            docker build \
                --network host \
                -t $DD_TAG:v1 -t $DD_TAG:latest \
                -f $DIR/runtime.Dockerfile \
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

            docker volume create --name $SHARED_VOLUME_NAME

            # use -d to start a container in detached mode
            # use --entrypoint=/bin/bash \ to override the command
            docker run \
                --rm \
                --name $CONTAINER_NAME \
                --add-host=ftp.hermes.acri.fr:5.252.148.37 \
                --volume "$SHARED_VOLUME_NAME":/media/share \
                -e AC_OUTPUT_DIR="$output_dir" \
                -e AC_YEAR="$year" \
                -e AC_ZONE="$zone" \
                -e AC_DEPTHMIN="$deepthmin" \
                -e AC_DEPTHMAX="$deepthmax" \
                -it $DD_TAG:latest
            ;;
        
        'ls')
            sudo ls -al /var/lib/docker/volumes/ac_share/_data/$2
            ;;
        *)
            echo 'todo';;
    esac
}

run "$1" "$2"

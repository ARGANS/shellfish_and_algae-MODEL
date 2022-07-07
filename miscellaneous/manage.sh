#!/bin/bash
SHARED_VOLUME_NAME='ac_share'

## Properties of scripts used to load datasets:

function build_images_for_model_execution {
    local dir="./dataread"
    local base_image_tag="ac-processing/base"
    local runtime_image_tag="ac-processing/runtime"
    local base_image_dockerfile="./miscellaneous/pythonBase.Dockerfile"
    local runtime_image_dockerfile="./miscellaneous/pythonRuntime.Dockerfile"

    ### for GDAL: --build-arg WITH_GDAL="true" \
    docker build \
        --network host \
        --build-arg REQUIREMENTS_PATH="./requirements.txt" \
        --build-arg WITH_R="true" \
        --build-arg WITH_NETCDF="true" \
        -t $base_image_tag:v1 -t $base_image_tag:latest \
        -f $base_image_dockerfile \
        $dir && \
    docker build \
        --network host \
        --build-arg BASE_IMAGE="$base_image_tag" \
        -t $runtime_image_tag:v2 -t $runtime_image_tag:latest \
        -f $runtime_image_dockerfile \
        $dir
}

function build_images_for_posttreatment {
    local dir="./posttreatment"
    local base_image_tag="ac-posttreatment/base"
    local runtime_image_tag="ac-posttreatment/runtime"
    local base_image_dockerfile="./miscellaneous/posttreatmentBase.Dockerfile"
    local runtime_image_dockerfile="./miscellaneous/posttreatmentRuntime.Dockerfile"

    docker build \
        --network host \
        -t $base_image_tag:v1 -t $base_image_tag:latest \
        -f $base_image_dockerfile \
        $dir && \
    docker build \
        --network host \
        --build-arg BASE_IMAGE="$base_image_tag" \
        -t $runtime_image_tag:v1 -t $runtime_image_tag:latest \
        -f $runtime_image_dockerfile \
        $dir
}

function prepare_runtime {
    local container_name="$1"
    local image_tag="$2"

    local container_id=$( docker ps -q -f name=$container_name )
    local image_id=$( docker images -q $image_tag )

    if [[ ! -z "$container_id" ]]; then
        docker stop $container_name || true
    fi

    if [[ -z "$image_id" ]]; then
        build_images_for_model_execution
    fi

    docker volume create --name $SHARED_VOLUME_NAME
}

function run_container_for_model_execution {
    local container_name="$1"
    local image_tag="$2"

    prepare_runtime $container_name $image_tag

    # All model properties received from the application
    # data=`cat macroalgae/macroalgae_model_parameters_input.json`

    zone='IBI'
    year=2021
    depthmin=0
    depthmax=20
    
    dataset_id="$zone-$year-$depthmin-$depthmax"

    # Form properties (ALARIA)
    json='{"parameters":{"species":{"alaria":{"options":{},"parameters":{"mu":0.33,"V_NH4":60,"V_NO3":25,"K_NH4":700,"K_NO3":100,"Q_max":70,"Q_min":14,"N_to_P":12,"K_c":7,"T_O":12,"T_min":1,"T_max":25,"I_s":277,"a_cs":0.00036,"d_m":0.003,"h_MA":0.4,"w_MA":0.2,"r_L":0.2,"r_N":0.1}}},"farm":{"default":{"options":{},"parameters":{"y_farm":1000,"density_MA":0.4,"x_farm":1000,"z":2}}},"harvest":{"CCA":{"options":{},"parameters":{"deployment_day":2,"harvest_first":65,"harvest_freq":50,"harvest_fraction":0.75,"deployment_Nf":10000}}},"run":{"default":{"options":{"harvest_method":0,"light_scheme":3},"parameters":{}}}},"metadata":{"name":"user1_Alaria_IBI_04-05-2022","zone":"IBI","_suggested":{"login":"user1","species":"Alaria","zone":"'$zone'","date":"04-05-2022"},"depth_min":'$depthmin',"depth_max":"'$depthmax'","year":"'$year'","data_import_container":"8496183ebd"}}'
    task_id="abcdef_alaria"
    
    # Form properties (Saccharina)
    # json="{\"parameters\":{\"species\":{\"saccharina\":{\"options\":{},\"parameters\":{\"mu\":0.18,\"V_NH4\":100,\"V_NO3\":200,\"K_NH4\":11,\"K_NO3\":200,\"Q_max\":22,\"Q_min\":10,\"N_to_P\":12,\"K_c\":8,\"T_O\":12.5,\"T_min\":0,\"T_max\":20,\"I_s\":90,\"a_cs\":0.00036,\"d_m\":0.0003,\"h_MA\":2,\"w_MA\":0.3,\"r_L\":0.1,\"r_N\":0.1,\"CN_MA\":21,\"kcal_MA\":2.29,\"prot_MA\":0.08,\"DF_MA\":0.113}}},\"farm\":{\"default\":{\"options\":{},\"parameters\":{\"y_farm\":1000,\"density_MA\":0.4,\"x_farm\":1000,\"z\":2}}},\"harvest\":{\"CCA\":{\"options\":{},\"parameters\":{\"deployment_day\":2,\"harvest_first\":65,\"harvest_freq\":50,\"harvest_fraction\":0.75,\"deployment_Nf\":10000}}},\"run\":{\"default\":{\"options\":{\"harvest_method\":0,\"light_scheme\":3},\"parameters\":{}}}},\"metadata\":{\"name\":\"user1_saccharina_IBI_16-05-2022\",\"zone\":\"$zone\",\"_suggested\":{\"login\":\"user1\",\"species\":\"saccharina\",\"zone\":\"IBI\",\"date\":\"16-05-2022\"},\"depth_min\":$depthmin,\"depth_max\":\"$depthmax\",\"year\":\"$year\",\"data_read_container\":\"5a4c48c2dc\"}}"
    # task_id="abcdef_saccharina"

    # use -d to start a container in detached mode
    # use --entrypoint=/bin/bash \ to override the command
    docker run \
        --rm \
        --name $container_name \
        --volume "$SHARED_VOLUME_NAME:/media/share" \
        -e DATASET_ID="$dataset_id" \
        -e TASK_ID="$task_id" \
        -e PARAMETERS_JSON="$json" \
        -e PYTHONDONTWRITEBYTECODE=1 \
        -it $image_tag:latest
}

function run_container_for_posttreatment {
    local container_name="$1"
    local image_tag="$2"

    #--volume "$(pwd)/posttreatment/src:/opt" \
        
    prepare_runtime $container_name $image_tag
    docker run \
        --rm \
        --name $container_name \
        --volume "$SHARED_VOLUME_NAME:/media/share" \
        -e PYTHONDONTWRITEBYTECODE=1 \
        -e SOURCE_DIR='/media/share/results/IBI-2021-0-20/521cac05234a4b7afb01d3a624924d53' \
        -it $image_tag:latest
}

function action_bash {
    local image_tag="ac-posttreatment/runtime"

    docker run \
        --rm \
        --volume "$SHARED_VOLUME_NAME:/media/share" \
        --volume "${HOME}:/media/home" \
        --workdir=/media/share \
        -e PYTHONDONTWRITEBYTECODE=1 \
        --entrypoint=/bin/bash \
        -it $image_tag:latest
}

function handle_arguments {
    local command="$1"

    case $command in
        'build_dataread')
            build_images_for_model_execution
            ;;
        'execute_dataread')
            run_container_for_model_execution "ac-model_run" "ac-processing/runtime"
            ;;
        'ls')
            sudo ls -al /var/lib/docker/volumes/$SHARED_VOLUME_NAME/_data/$2
            ;;
        'ls2')
            docker run --rm -i -v=$SHARED_VOLUME_NAME:/media/volume busybox sh
            ;;
        'build_posttreatment')
            build_images_for_posttreatment
            ;;
        'execute_posttreatment')
            run_container_for_posttreatment "ac-posttreatment_run" "ac-posttreatment/runtime"
            ;;
        'bash')
            action_bash
            ;;
        *)
            echo 'commands:'
            echo 'build_dataread, execute_dataread'
            echo 'ls, ls2'
            echo 'build_posttreatment'
            echo 'bash'
            ;;
    esac
}

handle_arguments "$1" "$2"

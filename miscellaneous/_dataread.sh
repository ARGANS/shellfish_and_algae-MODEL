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

function run_container_for_model_execution {
    local container_name="$1"
    local image_tag="$2"

    stop_existed_container $container_name
    create_volumes

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

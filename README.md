# shellfish_and_algae-MODEL

## Data downloading subtask 
- command to build docker image: `./miscellaneous/manage.sh build_dataimport`
- command to run container from the image created by the command above: `./miscellaneous/manage.sh run_dataimport`


## Commands to deploy
```
./miscellaneous/manage.sh build_dataimport
./miscellaneous/manage.sh build_dataread
./miscellaneous/manage.sh build_posttreatment
```

## /miscellaneous/manage.sh bash

Use one of the following commands to access the volume:
- `source miscellaneous/manage.sh bash`
- `. miscellaneous/manage.sh bash`


## Description of tasks in terms of manual execution step by step in order to debug the pipeline:

1) Data Import
1.1) Imagine that you need a new dataset, so you need to tweek the `miscellaneous/manage.sh` script.
1.2) In future tasks, this data set will be identified as a collection of its properties: IBI-2021-0-20. 
1.3) Then you need to build the image with the current code and run the container:
- `./miscellaneous/manage.sh build_dataimport`
- `./miscellaneous/manage.sh execute_dataimport`
1.4) When the command completes, you should be able to find the new files in the docker volume directory:
`sudo ls /var/lib/docker/volumes/ac_share/_data/data/IBI-2021-0-20`
	
2) Data read
2.1) You have made some changes in the code so you need to rebuild the image:
`./miscellaneous/manage.sh build_dataread`
2.2) You need to specify form properties and dataset id
``` SH
data='{"parameters":{"species":{"alaria":{"options":{},"parameters":{"mu":0.33,"V_NH4":60,"V_NO3":25,"K_NH4":700,"K_NO3":100,"Q_max":70,"Q_min":14,"N_to_P":12,"K_c":7,"T_O":12,"T_min":1,"T_max":25,"I_s":277,"a_cs":0.00036,"d_m":0.003,"h_MA":0.4,"w_MA":0.2,"r_L":0.2,"r_N":0.1}}},"farm":{"default":{"options":{},"parameters":{"y_farm":1000,"density_MA":0.4,"x_farm":1000,"z":2}}},"harvest":{"CCA":{"options":{},"parameters":{"deployment_day":2,"harvest_first":65,"harvest_freq":50,"harvest_fraction":0.75,"deployment_Nf":10000}}},"run":{"default":{"options":{"harvest_method":0,"light_scheme":3},"parameters":{}}}},"metadata":{"name":"user1_Alaria_IBI_04-05-2022","zone":"IBI","_suggested":{"login":"user1","species":"Alaria","zone":"IBI","date":"04-05-2022"},"depth_min":0,"depth_max":"20","year":"2021","data_import_container":"8496183ebd"}}'

docker run \
    --rm \
    --name $container_name \
    --volume "$SHARED_VOLUME_NAME:/media/share" \
    -e DATASET_ID="IBI-2021-0-20" \ # <- the dataset id specified at 1.1)
    -e TASK_ID="abcdef" \ # <- hash of the task, which will be needed to track the results of this task
    -e PARAMETERS_JSON="$data" \ # <- form properties
    -e PYTHONDONTWRITEBYTECODE=1 \
    -it $image_tag:latest
```	
2.3) Then You can use the `./miscellaneous/manage.sh execute_dataread` command to run a new container with code


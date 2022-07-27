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

### step 0: build images
0.1. `./miscellaneous/manage.sh build_dataimport`
0.2. `./miscellaneous/manage.sh build_pretreatment`
0.3. `./miscellaneous/manage.sh build_dataread`
0.4. 
0.5. `./miscellaneous/manage.sh build_posttreatment`

### step 1: dataimport
1.0. `./miscellaneous/manage.sh execute_dataimport`

### step 2: pretreatment
2.0. `./miscellaneous/manage.sh execute_pretreatment`

### step 3: dataread
3.1. `./miscellaneous/manage.sh execute_dataread`

### step 4: B

### step 5: posttreatment
5.1. `./miscellaneous/manage.sh execute_posttreatment`

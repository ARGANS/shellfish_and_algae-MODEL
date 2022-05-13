# shellfish_and_algae-MODEL

## Data downloading subtask 
- command to build docker image: `./manage_dataimport.sh build`
- command to run container from the image created by the command above: `./manage_dataimport.sh run`


## Commands to deploy
```
./manage_dataimport.sh build
./miscellaneous/manage.sh build_dataread
./miscellaneous/manage.sh build_posttreatment
```
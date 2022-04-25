#!/bin/bash

# https://docs.docker.com/engine/reference/commandline/system_events/

# Filter options:
# --filter 'type=container'
# --filter 'event=stop', start stop die destroy create ...

# --format '{{json .}}'
# --format 'Type={{.Type}}  Status={{.Status}}  ID={{.ID}}'
CONTAINER='ac-dataimport'

opts="--filter 'container=$CONTAINER' --filter 'event=start' --filter 'event=die'"
docker events $opts |  while read line; 
do  
    if [[ ${line} = *"start"* ]]; then 
        echo "Container started ${line}" ;
    elif [[ ${line} = *"die"* ]]; then 
        echo "Container terminated ${line}" ;
    fi
done

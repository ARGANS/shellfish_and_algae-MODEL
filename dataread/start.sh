#!/bin/bash
workdir=/media/share/results/$TASK_ID
mkdir -p $workdir
echo -n $PARAMETERS_JSON > $workdir/parameters.json
python start.py
echo "Start concatenation"
. concatenate_longitude.sh $workdir $workdir/concat.nc
now=$(date +%s)
echo -n "$now" > $workdir/task.mark

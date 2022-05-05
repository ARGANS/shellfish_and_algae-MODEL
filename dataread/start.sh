#!/bin/bash
workdir=/media/share/results/$DATASET_ID/$TASK_ID
mkdir -p $workdir
echo -n $PARAMETERS_JSON > $workdir/parameters.json
python start.py || rm $workdir/parameters.json
echo "Start concatenation"
. concatenate_longitude.sh $workdir $workdir/concat.nc
now=$(date +%s)
echo -n "$now" > $workdir/task.mark

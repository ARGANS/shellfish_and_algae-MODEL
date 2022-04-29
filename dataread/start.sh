#!/bin/bash
mkdir -p /media/share/results/$TASK_ID
python start.py
echo "Start concatenation"
. concatenate_longitude.sh /media/share/results/$TASK_ID /media/share/results/$TASK_ID/concat.nc

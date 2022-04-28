#!/bin/bash

mkdir -p /media/share/
### rm -rf /media/share/*
# DEPRECATED
# mkdir -p /media/share/data/IBI/{eastward_Water_current,northward_Water_current,Salinity,Phosphate,Ammonium,Nitrate,Temperature}
mkdir -p $AC_OUTPUT_DIR/{eastward_Water_current,northward_Water_current,Salinity,Phosphate,Ammonium,Nitrate,Temperature}

python start.py 
# TODO
now=$(date +%s)
echo -n "$now" > $AC_OUTPUT_DIR/task.mark
# DEPRECATED
# chmod -R 777 /media/share/data/IBI
chmod -R 777 $AC_OUTPUT_DIR

#!/bin/bash

mkdir -p /media/share/
rm -rf $AC_OUTPUT_DIR/*

mkdir -p $AC_OUTPUT_DIR/{eastward_Water_current,northward_Water_current,Salinity,Phosphate,Ammonium,Nitrate,Temperature,pCO2,disolved_inorganic_carbon,primary_producer_POC}
echo -n $parameters_json > $AC_OUTPUT_DIR/parameters.json

python start.py || rm $AC_OUTPUT_DIR/parameters.json
now=$(date +%s)
echo -n "$now" > $AC_OUTPUT_DIR/task.mark

chmod -R 777 $AC_OUTPUT_DIR

#!/bin/bash

destination="$INPUT_DESTINATION"
mkdir -p $destination
rm -rf $destination/*

echo -n $(date +%s) > $destination/start.mark
echo -n $INPUT_PARAMETERS > $destination/parameters.json
error_log=$destination/error.txt

# Detect if the zone is 'Europe'
sim_zone=`cat $destination/parameters.json | jq -r ".zone"`
if [ $sim_zone == "Europe" ]; then
    echo 'Will not download anything when "Europe" is selected, proceed to next step.' >> $error_log
    echo -n $(date +%s) > $destination/end.mark
    exit 0
fi

mkdir -p $destination/{eastward_Water_current,northward_Water_current,Phosphate,Ammonium,Nitrate,Temperature,par,Chlorophyll-a}

python start.py 2>$error_log
if [ $? -eq 0 ]
then
  echo "Success"
else
  echo "Failure:"
  echo 'Finished '$(date "+%d/%m/%Y %H:%M:%S") >> $error_log
  cat $error_log
fi

echo -n $(date +%s) > $destination/end.mark

chmod -R 777 $destination

#!/bin/bash

destination="$INPUT_DESTINATION"
mkdir -p $destination
rm -rf $destination/*

echo -n $(date +%s) > $destination/start.mark
mkdir -p $destination/{eastward_Water_current,northward_Water_current,Phosphate,Ammonium,Nitrate,Temperature,par,Chlorophyll-a}

echo -n $INPUT_PARAMETERS > $destination/parameters.json
error_log=$destination/error.txt

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

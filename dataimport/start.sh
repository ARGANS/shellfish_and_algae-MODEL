#!/bin/bash

destination="$INPUT_DESTINATION"
rm -rf $destination/*
mkdir -p $destination/{eastward_Water_current,northward_Water_current,Salinity,Phosphate,Ammonium,Nitrate,Temperature,pCO2,disolved_inorganic_carbon,primary_producer_POC}

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

echo -n $(date +%s) > $destination/task.mark

chmod -R 777 $destination
#!/bin/bash

destination=$INPUT_DESTINATION
mkdir -p $destination
rm -rf $destination/*

echo -n $(date +%s) > $destination/start.mark
mkdir -p $destination/{eastward_Water_current,northward_Water_current,Phosphate,Ammonium,Nitrate,Temperature,pCO2,disolved_inorganic_carbon,primary_producer_POC,par}

error_log=$destination/error.txt
print_log=$destination/print.txt

python start.py 1>$print_log 2>$error_log

if [ $? -eq 0 ]
then
  echo "Success"
else
  echo "Failure:"
  echo 'Finished '$(date "+%d/%m/%Y %H:%M:%S") >> $error_log
  cat $error_log
fi
echo '-------------------- Last 100 lines printed to stdout --------------------' >> $error_log
tail -n 100 $print_log >> $error_log

echo -n $(date +%s) > $destination/end.mark

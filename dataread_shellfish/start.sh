#!/bin/bash

workdir=$INPUT_DESTINATION
mkdir -p $workdir
rm -rf $workdir/*

echo -n $(date +%s) > $workdir/start.mark
echo -n $INPUT_MODEL_PROPERTIES_JSON > $workdir/parameters.json
error_log=$workdir/error.txt
print_log=$workdir/print.txt

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


#echo "Start concatenation"
#. concatenate_longitude.sh $workdir $workdir/concat.nc 2>>$error_log
#if [ ! $? -eq 0 ]; then
#    cat $error_log
#fi

echo -n $(date +%s) > $workdir/end.mark

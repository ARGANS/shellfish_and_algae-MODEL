#!/bin/bash

INPUT_FOLDER="$1"
OUTPUT_FILE="$2"

for filein in ${INPUT_FOLDER}/*; do

    # make longitude the first dimension
    ncpdq -O -a longitude,time $filein $filein

    # make longitude the record dimension
    ncks -O --mk_rec_dim longitude $filein $filein
done

#concatenate all
ncrcat -h ${INPUT_FOLDER}/* $OUTPUT_FILE

# revert time and longitude
ncpdq -O -a time,longitude $OUTPUT_FILE $OUTPUT_FILE
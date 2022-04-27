#!/bin/bash
mkdir -p /media/share/results/simulations/monthly
python start.py
echo "Start concatenation"
. concatenate_longitude.sh /media/share/results/simulations/monthly /media/share/results/simulations/concat.nc
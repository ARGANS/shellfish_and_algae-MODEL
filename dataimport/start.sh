#!/bin/bash

mkdir -p /media/share
mkdir -p /media/share/data/IBI/{eastward_Water_current,northward_Water_current,Salinity,Phosphate,Ammonium,Nitrate,Temperature}

python start.py

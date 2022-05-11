FROM python:3.10.1-slim-bullseye

LABEL NM <nmaltsev@argans.eu>
ARG WORK_DIR="/opt"


WORKDIR $WORK_DIR
COPY ./src ./
RUN chmod u+x concatenate_longitude.sh

RUN python -m pip install --upgrade pip

RUN apt update
RUN apt -y install nco

# install GDAL
RUN apt-get update && apt-get install -y build-essential binutils libproj-dev gdal-bin libgdal-dev

# Install R
RUN apt-get install -y r-base
RUN R -e "install.packages(c('ncdf4', 'rjson', 'seacarb'))"

RUN mkdir -p /media/share && \
    mkdir -p /media/share/data/IBI/{eastward_Water_current,northward_Water_current,Salinity,Phosphate,Ammonium,Nitrate,Temperature} && \
    mkdir -p /media/share/data_merged/IBI/{eastward_Water_current,northward_Water_current,Salinity,Phosphate,Ammonium,Nitrate,Temperature} &&\
    mkdir -p /media/share/results

# Install utilities for convenience
RUN apt-get install -y vim
RUN apt-get install -y netcdf-bin

ENTRYPOINT ["/usr/bin/bash"]
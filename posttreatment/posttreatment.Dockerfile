FROM python:3.10.1-slim-bullseye

LABEL NM <nmaltsev@argans.eu>
ARG WORK_DIR="/opt"


WORKDIR $WORK_DIR
COPY ./src ./
RUN chmod u+x concatenate_longitude.sh

RUN python -m pip install --upgrade pip

RUN apt update
RUN apt -y install nco

RUN mkdir -p /media/share && \
    mkdir -p /media/share/data/IBI/{eastward_Water_current,northward_Water_current,Salinity,Phosphate,Ammonium,Nitrate,Temperature} && \
    mkdir -p /media/share/data_merged/IBI/{eastward_Water_current,northward_Water_current,Salinity,Phosphate,Ammonium,Nitrate,Temperature} &&\
    mkdir -p /media/share/results

# Install utilities for convenience
RUN apt-get install -y vim
RUN apt-get install -y netcdf-bin

ENTRYPOINT ["/usr/bin/bash"]
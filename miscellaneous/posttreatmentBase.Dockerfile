FROM python:3.10.1-slim-bullseye

RUN apt update && \
    apt-get install -y nco

# install GDAL
RUN apt-get install -y build-essential binutils libproj-dev gdal-bin libgdal-dev

# Install R
RUN apt-get install -y r-base
RUN R -e "install.packages(c('ncdf4', 'rjson', 'seacarb'))"

# Install utilities for convenience
RUN apt-get install -y vim netcdf-bin

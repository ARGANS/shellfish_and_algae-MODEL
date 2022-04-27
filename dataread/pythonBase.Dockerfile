FROM python:3.10.1-slim-bullseye
WORKDIR /opt

RUN apt-get update && \
    python -m pip install --upgrade pip

# Install GDAL
ARG WITH_GDAL
RUN if [ ! -z "$WITH_GDAL" ] ; then \
    apt-get install -y build-essential binutils libproj-dev gdal-bin libgdal-dev python3-gdal python-dev nano && \
    pip install GDAL==$(gdal-config --version); \
    else echo "Without GDAL"; \
    fi

ARG WITH_R
RUN if [ ! -z "$WITH_R" ] ; then \
    apt-get install -y r-base libcurl4-openssl-dev libssl-dev libxml2-dev; \
    else echo "Without R"; \
    fi

ARG WITH_NETCDF
RUN if [ ! -z "$WITH_NETCDF" ] ; then \
    apt-get install -y nco netcdf-bin; \
    else echo "Without NetCDF"; \
    fi

ARG REQUIREMENTS_PATH
RUN echo ":: $REQUIREMENTS_PATH"
COPY "$REQUIREMENTS_PATH" .
RUN python -m pip --no-cache-dir install -r requirements.txt

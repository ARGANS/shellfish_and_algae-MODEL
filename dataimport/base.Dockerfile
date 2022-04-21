FROM python:3.10.1-slim-bullseye
WORKDIR /opt

RUN python -m pip install --upgrade pip
# install GDAL
RUN apt-get update && \
    apt-get install -y build-essential binutils libproj-dev gdal-bin libgdal-dev python3-gdal python-dev nano && \
    pip install GDAL==$(gdal-config --version)

COPY ./src/requirements.txt .
RUN python -m pip --no-cache-dir install -r requirements.txt

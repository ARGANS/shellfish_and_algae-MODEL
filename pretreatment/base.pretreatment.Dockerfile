FROM python:3.10.1-slim-bullseye

RUN apt update && \
    apt -y install nco netcdf-bin vim

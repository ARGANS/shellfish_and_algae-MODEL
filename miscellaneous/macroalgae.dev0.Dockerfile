FROM python:3.10.1-slim-bullseye

ARG WORK_DIR="/opt"

WORKDIR $WORK_DIR
RUN mkdir -p /media/share
COPY ./ ./


RUN python -m pip install --upgrade pip

# R v4.0.4
RUN apt-get update && apt-get install -y r-base
RUN apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev
RUN R -e "install.packages(c('deSolve', 'rootSolve', 'bvpSolve', 'deTestSet', 'ReacTran', 'simecol', 'tidyverse', 'GA'), dependencies = TRUE)"

# https://cran.r-project.org/#debian-bullseye
# https://stackoverflow.com/questions/69971991/r-backports-and-invalid-signature
# ?? apt-get --allow-unauthenticated update

ENTRYPOINT ["/bin/bash"]

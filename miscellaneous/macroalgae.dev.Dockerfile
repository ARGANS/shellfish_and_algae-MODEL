FROM python:3.10.1-slim-bullseye

ARG WORK_DIR="/opt"

WORKDIR $WORK_DIR
RUN mkdir -p /media/share && \
    mkdir -p /media/share/data_merged/IBI/{eastward_Water_current,northward_Water_current,Salinity,Phosphate,Ammonium,Nitrate,Temperature}
COPY dataread ./
COPY macroalgae ./
COPY dataimport/src/dataCmd.csv ./dataCmd.csv

RUN apt-get update && apt-get install -y gnupg apt-transport-https ca-certificates software-properties-common libcurl4-openssl-dev libssl-dev libxml2-dev
RUN echo "deb [trusted=yes] http://cloud.r-project.org/bin/linux/debian bullseye-cran40/" >> /etc/apt/sources.list.d/backports.list
RUN apt-get update
#RUN apt-get -y upgrade
## to show the list of available versions, use the command 'apt list -a r-base`
#RUN apt-get install -y r-base=4.1.2-1~bullseyecran.0
RUN apt-get install -y r-base
RUN apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev
#RUN R -e "install.packages(c('deSolve', 'rootSolve', 'bvpSolve', 'deTestSet', 'ReacTran', 'simecol', 'tidyverse', 'GA'))"
RUN R -e "install.packages(c('deSolve', 'rootSolve', 'tidyverse'))"

RUN python -m pip install --upgrade pip
RUN python -m pip --no-cache-dir install -r requirements.txt

ENTRYPOINT ["/bin/bash"]

FROM python:3.10.1-slim-bullseye

LABEL NM <nmaltsev@argans.eu>
ARG WORK_DIR="/opt"
ARG CDSAPIRC_DIR="/root"



WORKDIR $WORK_DIR
COPY ./src ./
RUN echo "HOME ${CDSAPIRC_DIR} :: ${HOME}"
COPY ./secret/main.cdsapirc $CDSAPIRC_DIR/.cdsapirc

RUN python -m pip install --upgrade pip
# install GDAL
RUN apt-get update && apt-get install -y build-essential binutils libproj-dev gdal-bin libgdal-dev python3-gdal python-dev && \
    pip install GDAL==$(gdal-config --version)
RUN python -m pip --no-cache-dir install -r requirements.txt

RUN mkdir -p /media/share


# RUN chmod +x ./main.sh
# CMD ["./main.sh"]
# CMD ["bash"]
ENTRYPOINT ["/usr/bin/bash"]

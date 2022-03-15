FROM python:3.10.1-slim-bullseye

LABEL NM <nmaltsev@argans.eu>
ARG WORK_DIR="/opt"



WORKDIR $WORK_DIR
COPY ./src ./

RUN python -m pip install --upgrade pip
# install GDAL
#RUN apt-get update && apt-get install -y build-essential binutils libproj-dev gdal-bin libgdal-dev python3-gdal python-dev nano && \
    pip install GDAL==$(gdal-config --version)
#RUN python -m pip --no-cache-dir install -r requirements.txt

RUN aptitude install nco

RUN mkdir -p /media/share && \
    mkdir -p /media/share/data/IBI/{eastward_Water_current,northward_Water_current,Salinity,Phosphate,Ammonium,Nitrate,Temperature}


# RUN chmod +x ./main.sh
# CMD ["./main.sh"]
# CMD ["bash"]
ENTRYPOINT ["/usr/bin/bash"]

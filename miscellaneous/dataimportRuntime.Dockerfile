ARG BASE_IMAGE
FROM $BASE_IMAGE

LABEL NM <nmaltsev@argans.eu>
ARG CDSAPIRC_DIR="/root"

WORKDIR /opt
COPY ./src ./
COPY ./secret/main.cdsapirc $CDSAPIRC_DIR/.cdsapirc
COPY ./start.sh .
RUN chmod +x *.sh
RUN echo "HOME ${CDSAPIRC_DIR} :: ${HOME}"

CMD ["/opt/start.sh"]

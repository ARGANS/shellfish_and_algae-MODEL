FROM aquaculture/base

LABEL NM <nmaltsev@argans.eu>
ARG CDSAPIRC_DIR="/root"

WORKDIR /opt
COPY ./src ./
COPY ./start.sh .
RUN chmod +x start.sh
RUN echo "HOME ${CDSAPIRC_DIR} :: ${HOME}"
COPY ./secret/main.cdsapirc $CDSAPIRC_DIR/.cdsapirc

CMD ["/opt/start.sh"]
# CMD ["bash"]


ARG BASE_IMAGE
FROM $BASE_IMAGE

LABEL NM <nmaltsev@argans.eu>

WORKDIR /opt
RUN mkdir -p /opt/dataread && \
    mkdir -p /opt/advectionPrototype

COPY ./dataread/src dataread
COPY ./advectionPrototype advectionPrototype
COPY ./scenario_b/. .

RUN chmod +x start.sh

CMD ["/opt/start.sh"]

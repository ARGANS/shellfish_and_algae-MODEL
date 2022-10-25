ARG BASE_IMAGE
FROM $BASE_IMAGE


WORKDIR /opt
RUN mkdir -p /opt

COPY ./farm_repartition/. .

RUN chmod +x start.sh

CMD ["/opt/start.sh"]

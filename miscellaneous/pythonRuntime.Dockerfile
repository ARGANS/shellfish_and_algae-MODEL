ARG BASE_IMAGE
FROM $BASE_IMAGE

LABEL NM <nmaltsev@argans.eu>

WORKDIR /opt
COPY . .

COPY ./start.sh .
RUN chmod +x *.sh

CMD ["/opt/start.sh"]
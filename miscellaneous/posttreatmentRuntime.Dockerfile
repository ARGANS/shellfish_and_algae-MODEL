ARG BASE_IMAGE
FROM $BASE_IMAGE

WORKDIR /opt
COPY ./src/ .
RUN chmod +x *.sh && \
    chmod +x make_interest_vars.R

CMD ["/opt/start.sh"]

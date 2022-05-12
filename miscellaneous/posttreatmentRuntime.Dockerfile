ARG BASE_IMAGE
FROM $BASE_IMAGE

WORKDIR /opt
COPY ./src/ .


# COPY ./start.sh .
# RUN chmod +x *.sh
# CMD ["/opt/start.sh"]
CMD ["/bin/bash"]
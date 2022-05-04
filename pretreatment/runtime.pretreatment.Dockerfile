FROM aquaculture/pretreatment-base

WORKDIR /opt
COPY ./src ./
RUN chmod u+x concatenate_copernicus.sh
RUN chmod u+x start.sh

CMD ["/opt/start.sh"]

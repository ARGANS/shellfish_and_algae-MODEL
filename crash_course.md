## Introduction to Docker terminology

You need to know about 3 terms that are often used - Containers, Images and Docker files:
- Containers are processes that are isolated from other processes on the host system.
- Images are archives that store the artifacts and dependencies needed to run processes in a container.
- Docker files are text files with instructions describing how to build an image.

So, a typical workflow consists of these steps
- creating docker files
- building images
- starting containers

* As a result, the running container looks like an instance of a virtual machine;
* Containers can communicate with the host system via the http protocol, shared volumes.

In the repository, I create a [bash script](https://github.com/ARGANS/shellfish_and_algae-MODEL/blob/main/manage_dataimport.sh) with commands that facilitate management. From this script, the code below uses the docker CLI utility to build a Docker image:

```
docker build \
	-t aquaculture/datadownload:v1 \
	-t aquaculture/datadownload:latest \
	-f ./dataimport/dataimport.Dockerfile \
	./dataimport
```
- The -f flag is used to specify the path to the docker file.
- Tags (-t) are aliases used to label images
- The last argument of the command is used to specify the directory - the scope of available files. Thus, the Docker Build command will not be able to get files that are not in this directory.


The `docker build` command downloads the base image, executes all the commands specified in the dockerfile, copies the required files, and sets the required properties.
Command output: 
```
...
Step 13/13 : ENTRYPOINT ["/usr/bin/bash"]
 ---> Running in 00dc99ab8fc8
Removing intermediate container 00dc99ab8fc8
 ---> 91fd63d47cc8
Successfully built 91fd63d47cc8
Successfully tagged aquaculture/datadownload:v1
Successfully tagged aquaculture/datadownload:latest
```

You can verify that the image was created with the following command:
```
$ docker images aquaculture/datadownload
REPOSITORY                 TAG       IMAGE ID       CREATED         SIZE
aquaculture/datadownload   latest    91fd63d47cc8   2 minutes ago   1.27GB
aquaculture/datadownload   v1        91fd63d47cc8   2 minutes ago   1.27GB
```
All images are identified by an id (image id) and a tag, so you can use both in commands.


## Let's check out the docker file 

Any dockerfile is a list of commands where the command is written in upper case (`<COMMAND NAME> <arguments>`):

``` Dockerfile
FROM python:3.10.1-slim-bullseye # The base image to be used
# Docker images support inheritance, so the output container can use settings and applications from base images. 

LABEL NM <nmaltsev@argans.eu> # metadata for the annotation  (optional)
ARG WORK_DIR="/opt" # argument declaration. Image can be parameterized at build time


WORKDIR $WORK_DIR #the WORKDIR command is used to specify the default directory used in a container. 
COPY ./src ./ # copying files from the repository (which is limited to the area specified by the last argument of the docker build command) to the container. All files from the src directory in the repository will be copied to the /opt directory (the default directory specified by the WORKDIR command)
COPY ./secret/main.cdsapirc $CDSAPIRC_DIR/.cdsapirc # example of copying one file to a specific directory in container

RUN python -m pip install --upgrade pip # executing the shell command inside the image. All created files will be available inside the container
RUN mkdir -p /media/share && \
    mkdir -p /media/share/data/IBI/{eastward_Water_current,northward_Water_current,Salinity,Phosphate,Ammonium,Nitrate,Temperature}

# RUN chmod +x ./main.sh
# CMD ["./main.sh"] # You can run any application as the main process in a container.
# or
# CMD ["bash"]
ENTRYPOINT ["/usr/bin/bash"]
```



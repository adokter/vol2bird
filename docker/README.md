### Docker container for vol2bird algorithm
Copyright 2016 Adriaan Dokter (University of Amsterdam)

If you are interested in using this software through collaboration, please contact me at a.m.dokter@uva.nl

This repository contains code to generate vertical profiles of birds (VPBs), following the methods described in:

*Bird migration flight altitudes studied by a network of operational weather radars*
Dokter A.M., Liechti F., Stark H., Delobbe L., Tabary P., Holleman I.
J. R. Soc. Interface, **8**, 30–43, 2011.
DOI: [10.1098/rsif.2010.0116](https://doi.org/10.1098/rsif.2010.0116)

This folder contains scripts for building and running vol2bird in a [Docker](https://www.docker.com/) container:
* `Dockerfile`: Dockerfile to make a Docker image for vol2bird
* `make_docker_image.sh`: Bash script to build and compact Docker image
* `start_container.sh`: starts the Docker container
* `stop_container.sh`: stops the Docker container
* `run_vol2bird.sh`: runs vol2bird on a started Docker container

To process a radar file `yourfile` in directory `your\data\durectory` using Docker:
```
# automatically pull and start the vol2bird docker image from DockerHub:
./start_container.sh your/data/directory
# goto your data directory
cd your/data/directory
# process your file:
run_vol2bird.sh yourfile
# when done, stop the container:
stop_container.sh
```

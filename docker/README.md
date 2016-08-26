### docker container for vol2bird algorithm
Copyright 2016 Adriaan Dokter (University of Amsterdam)

If you are interested in using this software through collaboration, please contact me at a.m.dokter@uva.nl

This repository contains code to generate vertical profiles of birds (VPBs), following the methods described in:

*Bird migration flight altitudes studied by a network of operational weather radars*
Dokter A.M., Liechti F., Stark H., Delobbe L., Tabary P., Holleman I.
J. R. Soc. Interface, **8**, 30–43, 2011.
DOI: [10.1098/rsif.2010.0116](https://doi.org/10.1098/rsif.2010.0116)

This folder contains:
* `Dockerfile`: Dockerfile to make a Docker image for vol2bird
* `make_docker_image.sh`: Bash script to build and compact Docker image
* `start_container.sh`: starts the Docker container
* `start_container.sh`: stops the Docker container
* `run_vol2bird.sh`: runs vol2bird on a started Docker container

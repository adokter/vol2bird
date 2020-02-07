#!/bin/bash

# build the docker image
docker build -f Dockerfile.vol2bird-mistnet -t vol2bird-mistnet/uncompacted .

echo "compacting the image..."
# get the image ID
ID=$(docker run -d vol2bird-mistnet/uncompacted /bin/bash)
# export and import
# all ENV and CMD statements are deleted in this process
# put them back with --change options
docker export $ID | docker import --change="ENV LD_LIBRARY_PATH=/opt/radar/lib:/opt/radar/rave/lib:/opt/radar/rsl/lib:/opt/radar/vol2bird/lib:/opt/radar/libtorch/lib" --change="ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/radar/rsl/bin:/opt/radar/vol2bird/bin" --change="CMD vol2bird" - vol2bird-mistnet

# stop running instances of vol2bird container
docker stop vol2bird
docker rm vol2bird

echo "cleaning up untagged images (1/2) ..."
untaggedimg=`docker images | grep "^<none>" | awk '{print $3}'`
if [ ! -z "$untaggedimg" ]; then
  docker rmi $untaggedimg
fi


echo "cleaning up ALL containers ..."
untaggedcont=`docker ps -aq`
if [ ! -z "$untaggedcont" ]; then
  docker rm $untaggedcont
fi


echo "cleaning up untagged images (2/2) ..."
untaggedimg=`docker images | grep "^<none>" | awk '{print $3}'`
if [ ! -z "$untaggedimg" ]; then
  docker rmi $untaggedimg
fi

echo "done"

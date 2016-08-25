#!/bin/bash

# build the docker image
docker build -t vol2bird/uncompacted .

echo "compacting the image..."
# get the image ID
ID=$(docker run -d vol2bird/uncompacted /bin/bash)
# export and import
# all ENV and CMD statements are deleted in this process
# put them back with --change options
docker export $ID | docker import --change="ENV LD_LIBRARY_PATH=/opt/radar/lib:/opt/radar/rave/lib:/opt/radar/rsl/lib:/opt/radar/vol2bird/lib" --change="ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/radar/vol2bird/bin" --change="CMD vol2bird" - vol2bird

echo "cleaning up untagged images ..."
docker rmi `docker images | grep "^<none>" | awk '{print $3}'`

read -p "Clean up ALL docker containers on this machine? " -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
   docker rm `docker ps -a -q`
fi

# repeat cleanup of images
docker rmi `docker images | grep "^<none>" | awk '{print $3}'`

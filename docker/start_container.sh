#!/bin/bash
if [ ! -d "$1" ]; then
    echo "usage: $0 <directory>"
    echo "   expects a data directory to mount to docker container as first argument"   
else 
    DATADIR=`realpath $1`
 
    # start docker container
    docker run -v $DATADIR:/data -d --name vol2bird adokter/vol2bird tail -f /dev/null
fi

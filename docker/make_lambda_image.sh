#!/bin/bash
docker build -f Dockerfile.lambda -t vol2bird/lambda_uncompressed .
ID=$(docker run -d vol2bird/lambda_uncompressed /bin/bash)
docker export $ID | docker import --change="ENV LD_LIBRARY_PATH=/opt/radar/lib:/opt/radar/rave/lib:/opt/radar/rsl/lib:/opt/radar/vol2bird/lib" --change="ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/radar/vol2bird/bin" - adokter/lambda_vol2bird
docker push adokter/lambda_vol2bird

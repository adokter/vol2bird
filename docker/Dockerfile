FROM ubuntu:16.04
MAINTAINER Adriaan Dokter


# install library for parsing options
RUN apt-get update && apt-get install -y libconfuse-dev
# install HDF5, Hierarchichal Data Format library
RUN apt-get install --no-install-recommends -y libhdf5-dev 
# install git, for fetching repositories from Github
RUN apt-get install -y git
# install compiler (gcc, make, etc)
RUN apt-get install --no-install-recommends -y gcc \
    && apt-get install --no-install-recommends -y make
# install zlib (gzip archiving library)
RUN apt-get install --no-install-recommends -y zlib1g-dev
# install python
RUN apt-get install --no-install-recommends -y python-dev
# install numpy 
RUN apt-get install --no-install-recommends -y python-numpy
# install proj4 library
RUN apt-get install --no-install-recommends -y libproj-dev
# we need flex, otherwise configure script of RSL library does not function properly
RUN apt-get install --no-install-recommends -y flex
# we need wget to download the RSL tarball
RUN apt-get install --no-install-recommends -y wget

# get a copy of hlhdf:
# configure and build hlhdf
# strange Docker conflict when attempting to install in /opt/radar/hlhdf, therefore in root radar instead
RUN git clone git://git.baltrad.eu/hlhdf.git \
    && cd hlhdf && ./configure --prefix=/opt/radar --with-hdf5=/usr/include/hdf5/serial,/usr/lib/x86_64-linux-gnu/hdf5/serial \
    && make && make install && cd .. && rm -rf hlhdf

# get a copy of rave:
# cd into rave source directory and configure
RUN git clone git://git.baltrad.eu/rave.git \
    && cd rave && ./configure --prefix=/opt/radar/rave --with-hlhdf=/opt/radar \
    && make && make install && cd .. && rm -rf rave



# get a copy of RSL:
RUN wget --password=guest ftp://trmm-fc.gsfc.nasa.gov/software/rsl-v1.49.tar.gz \
    && tar -xf rsl-v1.49.tar.gz && rm rsl-v1.49.tar.gz && cd rsl-v1.49 \
    && ./configure --prefix=/opt/radar/rsl && make && make install \
    && cd .. && rm -rf rsl-v1.49

# get a copy of vol2bird
# configure vol2bird
RUN git clone https://github.com/adokter/vol2bird.git \
    && cd vol2bird && ./configure --prefix=/opt/radar/vol2bird --with-rave=/opt/radar/rave --with-rsl=/opt/radar/rsl \
    && make && make install && cd .. && rm -rf vol2bird

# clean up
RUN apt-get remove -y build-essential && apt-get remove -y python-numpy && apt-get remove -y python-dev \
    && apt-get remove -y git && apt-get remove -y flex && apt-get remove -y wget\
    && apt-get clean && apt -y autoremove && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# prepare mount points
RUN mkdir data

# set the paths to installed libraries and executables
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/radar/lib:/opt/radar/rave/lib:/opt/radar/rsl/lib:/opt/radar/vol2bird/lib
ENV PATH=${PATH}:/opt/radar/vol2bird/bin

CMD vol2bird

FROM ubuntu:16.04
MAINTAINER Adriaan Dokter

# installs using apt-get:
# * libconfuse: library for parsing options
# * libhdf5: HDF5, Hierarchichal Data Format library
# * git, for fetching repositories from Github
# * compiler (gcc, make, etc)
# * zlib (gzip archiving library)
# * python
# * numpy 
# * proj4 library
# * flex, otherwise configure script of RSL library does not function properly
RUN apt-get update && apt-get install --no-install-recommends -y libconfuse-dev \
    libhdf5-dev gcc make zlib1g-dev python-dev python-numpy libproj-dev flex file \
    && apt-get install -y git && apt-get install -y libgsl-dev

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
RUN git clone https://github.com/adokter/rsl.git && cd rsl \
    && ./configure --prefix=/opt/radar/rsl && make AUTOCONF=: AUTOHEADER=: AUTOMAKE=: ACLOCAL=: \
    && make install AUTOCONF=: AUTOHEADER=: AUTOMAKE=: ACLOCAL=: \
    && cd .. && rm -rf rsl

# get a copy of vol2bird
# configure vol2bird
RUN git clone https://github.com/adokter/vol2bird.git \
    && cd vol2bird && ./configure --prefix=/opt/radar/vol2bird --with-rave=/opt/radar/rave --with-rsl=/opt/radar/rsl \
    --with-gsl=/usr/include/gsl,/usr/lib/x86_64-linux-gnu \
    && make && make install && cd .. && rm -rf vol2bird

# clean up
RUN apt-get remove -y git gcc make python-numpy -y python-dev flex \
    && apt-get clean && apt -y autoremove && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# prepare mount points
RUN mkdir data

# set the paths to installed libraries and executables
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/radar/lib:/opt/radar/rave/lib:/opt/radar/rsl/lib:/opt/radar/vol2bird/lib
ENV PATH=${PATH}:/opt/radar/vol2bird/bin

CMD vol2bird

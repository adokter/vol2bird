FROM ubuntu:18.04
MAINTAINER Adriaan Dokter

# installs using apt-get:
# * libcsv: library for parsing csv files
# * libconfuse: library for parsing options
# * libhdf5: HDF5, Hierarchichal Data Format library
# * libgsl: the GNU Scientific Library
# * git, for fetching repositories from Github
# * git-lfs, for fetching large file repositories from Github
# * wget for downloading files, specifically libtorch
# * unzip 
# * compiler (gcc, g++, make, cmake, etc)
# * zlib (gzip archiving library)
# * libbz2 (bzip2 archiving library)
# * python
# * numpy 
# * proj4 library
# * flexold, otherwise configure script of RSL library does not function properly
RUN apt-get update && apt-get install --no-install-recommends -y libcsv-dev libconfuse-dev \
    libhdf5-dev gcc g++ wget unzip make cmake zlib1g-dev python-dev python-numpy libproj-dev flex-old file \
    && apt-get install -y git git-lfs && apt-get install -y libgsl-dev && apt-get install -y libbz2-dev

# get a copy of hlhdf:
# configure and build hlhdf
# strange Docker conflict when attempting to install in /opt/radar/hlhdf, therefore in root radar instead
RUN git clone git://github.com/adokter/hlhdf.git \
    && cd hlhdf && ./configure --prefix=/opt/radar --with-hdf5=/usr/include/hdf5/serial,/usr/lib/x86_64-linux-gnu/hdf5/serial \
    && make && make install && cd .. && rm -rf hlhdf

# get a copy of rave:
# cd into rave source directory and configure
# using a clone from git://git.baltrad.eu/rave.git \
RUN git clone https://github.com/adokter/rave.git \
    && cd rave && ./configure --prefix=/opt/radar/rave --with-hlhdf=/opt/radar \
    && make && make install && cd .. && rm -rf rave

# get a copy of iris2odim:
RUN git clone https://github.com/adokter/iris2odim.git \
    && cd iris2odim && export RAVEROOT=/opt/radar \
    && make && make install && cd .. && rm -rf iris2odim

# get a copy of RSL:
RUN git clone https://github.com/adokter/rsl.git && cd rsl \
    && ./configure --prefix=/opt/radar/rsl \
    && make AUTOCONF=: AUTOHEADER=: AUTOMAKE=: ACLOCAL=: \
    && make install AUTOCONF=: AUTOHEADER=: AUTOMAKE=: ACLOCAL=: \
    && cd .. && rm -rf rsl

# get a copy of vol2bird
# configure vol2bird
RUN git clone https://github.com/adokter/vol2bird.git \
    && cd vol2bird && ./configure --prefix=/opt/radar/vol2bird --with-rave=/opt/radar/rave --with-rsl=/opt/radar/rsl \
    --with-gsl=/usr/include/gsl,/usr/lib/x86_64-linux-gnu \
    && make && make install && cd .. && rm -rf vol2bird

# clean up
RUN apt-get remove -y git git-lfs gcc g++ wget unzip make cmake python-numpy -y python-dev flex-old \
    && apt-get clean && apt -y autoremove && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# prepare mount points
RUN mkdir data

# set the paths to installed libraries and executables
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/radar/lib:/opt/radar/rave/lib:/opt/radar/rsl/lib:/opt/radar/vol2bird/lib:/usr/lib/x86_64-linux-gnu
ENV PATH=${PATH}:/opt/radar/vol2bird/bin:/opt/radar/rsl/bin

CMD vol2bird

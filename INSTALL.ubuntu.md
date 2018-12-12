# Install instructions for Ubuntu 16.04 LTS

```
# the directory in which we will install 
RADAR_ROOT_DIR=${PWD}

# prepare a directory structure:
mkdir ${RADAR_ROOT_DIR}/opt
mkdir ${RADAR_ROOT_DIR}/src

# installs using apt-get:
# * libconfuse: library for parsing options
# * libhdf5: HDF5, Hierarchichal Data Format library
# * git, for fetching repositories from Github
# * compiler (gcc, make, etc)
# * zlib (gzip archiving library)
# * python2.7
# * numpy 
# * proj4 library
# * flex, otherwise configure script of RSL library does not function properly
sudo apt-get update && apt-get install -y libconfuse-dev \
    libhdf5-dev gcc make zlib1g-dev python-dev python-numpy libproj-dev flex file git libgsl-dev

# change to the source directory...
cd ${RADAR_ROOT_DIR}/src

# get a copy of hlhdf:
# configure and build hlhdf
git clone git://git.baltrad.eu/hlhdf.git \
    && cd hlhdf && ./configure --prefix=${RADAR_ROOT_DIR}/opt/hlhdf \
    --with-hdf5=/usr/include/hdf5/serial,/usr/lib/x86_64-linux-gnu/hdf5/serial \
    && make
sudo make install

cd ${RADAR_ROOT_DIR}/src

# get a copy of rave:
# cd into rave source directory and configure
sudo git clone git://git.baltrad.eu/rave.git \
    && cd rave && ./configure --prefix=${RADAR_ROOT_DIR}/opt/rave --with-hlhdf=${RADAR_ROOT_DIR}/opt/hlhdf \
    && make 
sudo make install

cd ${RADAR_ROOT_DIR}/src 

# get a copy of RSL:
# RSL installation is optional, only required for reading US NEXRAD data
sudo git clone https://github.com/adokter/rsl.git && cd rsl \
    && cd decode_ar2v && ./configure --prefix=/usr && make && make install && cd .. \
    && ./configure --prefix=${RADAR_ROOT_DIR}/opt/rsl && make AUTOCONF=: AUTOHEADER=: AUTOMAKE=: ACLOCAL=: \
    && make install AUTOCONF=: AUTOHEADER=: AUTOMAKE=: ACLOCAL=: 
    
cd ${RADAR_ROOT_DIR}/src 

# get a copy of vol2bird
# (if not installing RSL, remove --with-rsl flag below)
git clone https://github.com/adokter/vol2bird.git \
    && cd vol2bird && ./configure --prefix=${RADAR_ROOT_DIR}/opt/vol2bird \
    --with-rave=${RADAR_ROOT_DIR}/opt/rave --with-rsl=${RADAR_ROOT_DIR}/opt/rsl \
    --with-gsl=/usr/include/gsl,/usr/lib/x86_64-linux-gnu \
    && make
sudo make install

cd ${RADAR_ROOT_DIR}

# set the paths to installed libraries and executables
# set these each time you run vol2bird, or add to .bashrc
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${RADAR_ROOT_DIR}/opt/hlhdf/lib:${RADAR_ROOT_DIR}/opt/rave/lib:${RADAR_ROOT_DIR}/opt/rsl/lib:${RADAR_ROOT_DIR}/opt/vol2bird/lib
export PATH=${PATH}:${RADAR_ROOT_DIR}/opt/vol2bird/bin
```

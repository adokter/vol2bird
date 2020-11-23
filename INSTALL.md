# Install instructions for Ubuntu 18.04 LTS and Mac OSX

```
# the directory in which we will install
# (${PWD} = the current working directory)
RADAR_ROOT_DIR=${PWD}

# prepare a directory structure:
mkdir ${RADAR_ROOT_DIR}/opt
mkdir ${RADAR_ROOT_DIR}/src

## we need to install the following dependencies:
# * libconfuse: library for parsing options
# * libhdf5: HDF5, Hierarchichal Data Format library
# * libgsl: the GNU Scientific Library
# * git, for fetching repositories from Github
# * wget for downloading files, specifically libtorch
# * unzip 
# * compiler (gcc, g++, make, cmake, etc)
# * zlib (gzip archiving library)
# * libbz2 (bzip2 archiving library)
# * python
# * numpy 
# * proj4 library
# * flexold, otherwise configure script of RSL library does not function properly

# On Ubuntu/Linux apt-get can be used to install dependencies:
apt-get update && apt-get install --no-install-recommends -y libconfuse-dev \
    libhdf5-dev gcc g++ wget unzip make cmake zlib1g-dev python-dev python-numpy libproj-dev flex-old file \
    && apt-get install -y git && apt-get install -y libgsl-dev && apt-get install -y libbz2-dev
# On Mac OSX Macports (https://www.macports.org/) can be used to install dependencies:
sudo port selfupdate
sudo port install python27 zlib hdf5 libproj4 proj47 libconfuse py27-numpy gsl bzip2 flex

# On Mac OSX, install XCode command line tools
xcode-select --install

# change to the source directory...
cd ${RADAR_ROOT_DIR}/src

# get a copy of hlhdf and configure:
# we need to point the --with-hdf5 flag of the configure script
# to the location of the hdf5 headers and binaries
# On Mac OSX find out where the headers live on your system with mdfind
mdfind -name hdf5.h
mdfind -name libhdf5.a  (static)
# for OSX 10.10.5 using Macports the locations are /opt/local/include/,/opt/local/lib/
# for Ubuntu 18.04 the locations are /usr/include/hdf5/serial,/usr/lib/x86_64-linux-gnu/hdf5/serial
# configure hlhdf:
git clone git://github.com/adokter/hlhdf.git \
    && cd hlhdf && ./configure --prefix=${RADAR_ROOT_DIR}/opt/hlhdf \
	--with-hdf5=/usr/include/hdf5/serial,/usr/lib/x86_64-linux-gnu/hdf5/serial

# On Mac OSX only, edit the file def.mk in the root directory of the hlhdf package,
# and change the LDSHARED variable from -bundle into -dynamiclib
	
# build and install hlhdf:	
make
sudo make install

cd ${RADAR_ROOT_DIR}/src

# get a copy of rave:
# cd into rave source directory and configure
git clone https://github.com/adokter/rave.git \
    && cd rave 
# on Mac OSX you have to specify the location of your proj.4 library
mdfind -name projects.h
# using Macports the location is /opt/local/lib/proj47
export PROJ4ROOT=/opt/local/lib/proj47
# now we're ready to configure the install
./configure --prefix=${RADAR_ROOT_DIR}/opt/rave  --with-proj=${PROJ4ROOT} --with-hlhdf=${RADAR_ROOT_DIR}/opt/hlhdf	
# build and install:
make
sudo make install


cd ${RADAR_ROOT_DIR}/src 

# (optional) get a copy of iris2odim:
# adds support for Vaisala IRIS format
git clone https://github.com/adokter/iris2odim.git \
    && cd iris2odim && export RAVEROOT=${RADAR_ROOT_DIR}/opt/ \
    && make 
sudo make install RAVEROOT=${RADAR_ROOT_DIR}/opt/

cd ${RADAR_ROOT_DIR}/src 

# (optional) get a copy of RSL:
# RSL installation is optional, only required for reading US NEXRAD data
git clone https://github.com/adokter/rsl.git && cd rsl \
    && ./configure --prefix=${RADAR_ROOT_DIR}/opt/rsl && make AUTOCONF=: AUTOHEADER=: AUTOMAKE=: ACLOCAL=:
sudo make install AUTOCONF=: AUTOHEADER=: AUTOMAKE=: ACLOCAL=: 
    
cd ${RADAR_ROOT_DIR}/src 

# (optional) install functionality to run MistNet rain segmentation model
# get a copy of libtorch:
# on Ubuntu/Linux:
wget https://download.pytorch.org/libtorch/cpu/libtorch-shared-with-deps-1.3.0%2Bcpu.zip \
    && unzip libtorch-shared-with-deps-1.3.0+cpu.zip -d ${RADAR_ROOT_DIR}/opt/ \
# on Mac OSX:
wget https://download.pytorch.org/libtorch/cpu/libtorch-macos-1.3.0.zip \
    && unzip libtorch-macos-1.3.0.zip -d ${RADAR_ROOT_DIR}/opt/
# change to radar root directory
cd ${RADAR_ROOT_DIR}
# get a copy of the MistNet model
mkdir MistNet && cd MistNet && wget http://mistnet.s3.amazonaws.com/mistnet_nexrad.pt

cd ${RADAR_ROOT_DIR}/src 

# get a copy of vol2bird
git clone https://github.com/adokter/vol2bird.git \
    && cd vol2bird
# (if not installing Vaisala IRIS support, remove --with-iris flag below)
# (if not installing RSL, remove --with-rsl flag below)
# (if not installing MistNet, remove --with-libtorch flag below)

# set locations of gsl and confuse libraries:
# on Mac OSX:
export GSLROOT=/opt/local
export CONFUSEROOT=/opt/local
# on Ubuntu/Linux:
export GSLROOT=/usr/include/gsl,/usr/lib/x86_64-linux-gnu
# --with-confuse can be omitted on Ubuntu/Linux, detected automatically
# configure:
./configure --prefix=${RADAR_ROOT_DIR}/opt/vol2bird \
    --with-iris \
    --with-rave=${RADAR_ROOT_DIR}/opt/rave \
	--with-rsl=${RADAR_ROOT_DIR}/opt/rsl \
    --with-libtorch=${RADAR_ROOT_DIR}/opt/libtorch \
    --with-gsl=${GSLROOT} \
	--with-confuse=${CONFUSEROOT}
# build and install:
make
sudo make install

cd ${RADAR_ROOT_DIR}

# set the paths to installed libraries and executables
# set these each time you run vol2bird, or add to .bashrc
# to load automatically in each new shell session
export PATH=${PATH}:${RADAR_ROOT_DIR}/opt/rsl/bin:${RADAR_ROOT_DIR}/opt/vol2bird/bin
# On Ubuntu/Linux:
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${RADAR_ROOT_DIR}/opt/hlhdf/lib:${RADAR_ROOT_DIR}/opt/rave/lib:${RADAR_ROOT_DIR}/opt/rsl/lib:${RADAR_ROOT_DIR}/opt/vol2bird/lib:${RADAR_ROOT_DIR}/opt/libtorch/lib
# On Mac OSX:
export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${RADAR_ROOT_DIR}/opt/hlhdf/lib:${RADAR_ROOT_DIR}/opt/rave/lib:${RADAR_ROOT_DIR}/opt/rsl/lib:${RADAR_ROOT_DIR}/opt/vol2bird/lib:${RADAR_ROOT_DIR}/opt/libtorch/lib

# (optional) to run MistNet, add these two lines to your options.conf file and place it in the directory from which you run vol2bird:
# MISTNET_PATH=${RADAR_ROOT_DIR}/MistNet/mistnet_nexrad.pt
# USE_MISTNET=TRUE
```

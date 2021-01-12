# Install instructions for Ubuntu Linux and Mac OSX

```
# specify the directory in which we will install by setting the RADAR_ROOT_DIR variable
# make sure it is an absolute (not relative) path
RADAR_ROOT_DIR=/your/directory/goes/here
# note that whenever you exit your shell the value of RADAR_ROOT_DIR (and other shell
# variables defined further down) will be erased, 
# and you will need to define them again.

# make sure the installation directory exists and
# change directories to it:
cd ${RADAR_ROOT_DIR}
# print the current installation directory:
pwd

# prepare a directory structure:
# the opt subfolder will contain the final library files and binaries after
# compilation and installation:
mkdir ${RADAR_ROOT_DIR}/opt
# the src subfolder will contain the source code for the libraries and binaries
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
# this add the gcc compiler for compiling c and c++ code
xcode-select --install
# Note on Jan 2021:
# some Mac users have noticed compatability issues between the latest XCode version 12
# and the libtorch library required for the MistNet segmentation model, and had
# to downgrade to Xcode version 11

# change to the source directory:
cd ${RADAR_ROOT_DIR}/src

# get a copy of hlhdf and configure:
# we need to point the --with-hdf5 flag of the configure script
# to the location of the hdf5 headers and binaries
# On Mac OSX find out where the headers and library live on your system
# for OSX 10.10.5 using Macports the locations are typically /opt/local/include/,/opt/local/lib/
# for Ubuntu 18.04 the locations are /usr/include/hdf5/serial,/usr/lib/x86_64-linux-gnu/hdf5/serial
# clone the repo from github
git clone git://github.com/adokter/hlhdf.git
# enter the hlhdf directory
cd hlhdf
# configure hlhdf, using the correct --with-hdf5 flag (here for Ubuntu 18.04)
./configure --prefix=${RADAR_ROOT_DIR}/opt/hlhdf \
	--with-hdf5=/usr/include/hdf5/serial,/usr/lib/x86_64-linux-gnu/hdf5/serial
	
# On Mac OSX only, edit the file def.mk in the root directory of the hlhdf package,
# and change the LDSHARED variable from -bundle into -dynamiclib
	
# build and install hlhdf:	
make
sudo make install

# go back the main src directory
cd ${RADAR_ROOT_DIR}/src

# download a copy of RAVE:
git clone https://github.com/adokter/rave.git
# cd into rave source directory
cd rave 
# on Mac OSX you have to specify the location of your proj.4 library
# using Macports the location is typically /opt/local/lib/proj49
# note we need the latest version of library proj.4, not library proj (the reversioned name for version 5 and up)
# store the location of the proj.4 library in variable PROJ4ROOT (for future reference):
export PROJ4ROOT=/opt/local/lib/proj49
# store the install directory of RAVE in variable RAVEROOT (for future reference):
export RAVEROOT=${RADAR_ROOT_DIR}/opt/
# now we're ready to configure the install:
./configure --prefix=${RAVEROOT}/rave  --with-proj=${PROJ4ROOT} --with-hlhdf=${RADAR_ROOT_DIR}/opt/hlhdf	
# build and install:
make
sudo make install

# go back to the src directory of the main installation directory
cd ${RADAR_ROOT_DIR}/src 

# (optional) get a copy of iris2odim:
# adds support for Vaisala IRIS RAW format
# download iris2odim:
git clone https://github.com/adokter/iris2odim.git
# cd into the iris2odim directory:
cd iris2odim
# build and install:
make 
sudo make install RAVEROOT=${RAVEROOT}

# go back to the main src directory:
cd ${RADAR_ROOT_DIR}/src 

# (optional) get a copy of RSL:
# RSL installation is optional, only required for reading US NEXRAD data
# download RSL:
git clone https://github.com/adokter/rsl.git
# enter the RSL directory
cd rsl
# configure RSL
./configure --prefix=${RADAR_ROOT_DIR}/opt/rsl
# build RSL:
make AUTOCONF=: AUTOHEADER=: AUTOMAKE=: ACLOCAL=:
# install RSL:
sudo make install AUTOCONF=: AUTOHEADER=: AUTOMAKE=: ACLOCAL=: 

# go back to the main src directory:
cd ${RADAR_ROOT_DIR}/src 

# (optional) install functionality to run MistNet rain segmentation model
# get a copy of libtorch:
# on Ubuntu/Linux:
wget https://download.pytorch.org/libtorch/cpu/libtorch-shared-with-deps-1.7.1%2Bcpu.zip \
    && unzip libtorch-shared-with-deps-1.7.1+cpu.zip -d ${RADAR_ROOT_DIR}/opt/ \
# on Mac OSX:
wget https://download.pytorch.org/libtorch/cpu/libtorch-macos-1.7.1.zip \
    && unzip libtorch-macos-1.7.1.zip -d ${RADAR_ROOT_DIR}/opt/
# change to radar root directory
cd ${RADAR_ROOT_DIR}
# get a copy of the MistNet model
mkdir MistNet && cd MistNet && wget http://mistnet.s3.amazonaws.com/mistnet_nexrad.pt

# go back to the main src directory:
cd ${RADAR_ROOT_DIR}/src 

# almost there, get a copy of vol2bird:
git clone https://github.com/adokter/vol2bird.git
# enter the vol2bird directory:
cd vol2bird
# set locations of gsl and confuse libraries:
# on Mac OSX (when installed with Macports):
export GSLROOT=/opt/local
export CONFUSEROOT=/opt/local
# on Ubuntu/Linux:
export GSLROOT=/usr/include/gsl,/usr/lib/x86_64-linux-gnu

# (if not installing Vaisala IRIS support, remove --with-iris flag below)
# (if not installing RSL, remove --with-rsl flag below)
# (if not installing MistNet, remove --with-libtorch flag below)
# --with-confuse can be omitted on Ubuntu/Linux, detected automatically
# configure:
./configure --prefix=${RADAR_ROOT_DIR}/opt/vol2bird \
    --with-iris \
    --with-rave=${RADAR_ROOT_DIR}/opt/rave \
    --with-rsl=${RADAR_ROOT_DIR}/opt/rsl \
    --with-libtorch=${RADAR_ROOT_DIR}/opt/libtorch \
    --with-gsl=${GSLROOT} \
    --with-confuse=${CONFUSEROOT}
# build vol2bird:
make
# install vol2bird
sudo make install
# run vol2bird.sh to check if it works (should give vol2bird version X.X message and no errors):
./vol2bird.sh

# note that the vol2bird.sh sets the DYLD_LIBRARY_PATH (for Mac) and LD_LIBRARY_PATH (for Linux) variables,
# as you can see when printing the contents of the file:
cat ./vol2bird.sh

# Preferably these path definitions are also added to your ~/.bashrc
# or ~/.bash_profile file for bash shell, or ~/.zshenv for zsh shell
# such that the paths are loaded automatically in each new shell session
# (simply copy the export statements into these shell configuration files)

# In addition, it is useful to add the paths of the vol2bird and rsl binaries
# such that you can call the vol2bird command from any directory.
# Simply update the PATH variable in your shell configuration file:
# (make sure you replace ${RADAR_ROOT_DIR} with its value specified in the beginning)
export PATH=${PATH}:${RADAR_ROOT_DIR}/opt/rsl/bin:${RADAR_ROOT_DIR}/opt/vol2bird/bin

# (optional) to run MistNet from the command line, add these two lines to your options.conf
# file and place it in the directory from which you run vol2bird:
# (not needed when running vol2bird from bioRad)
# (make sure you replace ${RADAR_ROOT_DIR} with its value specified in the beginning)
MISTNET_PATH=${RADAR_ROOT_DIR}/MistNet/mistnet_nexrad.pt
USE_MISTNET=TRUE
```

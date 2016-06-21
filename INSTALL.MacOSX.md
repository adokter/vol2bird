```bash
# instructions for Mac OSX systems
#
# we will use Macports (https://www.macports.org/) to install dependencies

# cd to project directory

# define the root directory of the project
RADAR_ROOT_DIR=${PWD}

# prepare a directory structure:
mkdir ${RADAR_ROOT_DIR}/opt
mkdir ${RADAR_ROOT_DIR}/src

# install python 2 (not python 3)
sudo port install python27

# install XCode command line tools
xcode-select --install

# install zlib (gzip archiving library)
sudo port install zlib

# install HDF5, Hierarchichal Data Format library
sudo port install hdf5

# install library that can handle different projections
sudo port install libproj4
sudo port install proj47

# install library for parsing options
sudo port install libconfuse

# install Numpy
sudo port install py27-numpy 

# change to the source directory...
cd ${RADAR_ROOT_DIR}/src

# you'll need a copy of hlhdf:
git clone git://git.baltrad.eu/hlhdf.git

# change directory into the clone of hlhdf
cd hlhdf

# configure hlhdf
# we need to point to the location of the hdf5 headers and binaries
# find out where the headers live on your system with 
mdfind -name hdf5.h
# for OSX 10.10.5 using Macports the location is /opt/local/include/ 
# 
# and the same for the binaries:
mdfind -name libhdf5.a  (static)
# for OSX 10.10.5 using Macports the location is /opt/local/lib/


./configure --prefix=${RADAR_ROOT_DIR}/opt/hlhdf --with-hdf5=/opt/local/include/,/opt/local/lib

# Edit the file def.mk in the root directory of the hlhdf package, and change the LDSHARED variable from -bundle into -dynamiclib

# compile hlhdf5 in the local directory
make 

# (optional) compile hlhdf5 tests and run them
make test

# copy the binaries to where we want them (${RADAR_ROOT_DIR}/opt/hlhdf):
make install

# cd to the source directory
cd ${RADAR_ROOT_DIR}/src

# get a copy of rave:
git clone git://git.baltrad.eu/rave.git

# cd into rave source directory
cd rave

# sanity checks
python -c 'import numpy'
python -c 'import _pyhl'

# set some environment variables
export RAVEROOT=${RADAR_ROOT_DIR}/opt/rave
export HLDIR=${RADAR_ROOT_DIR}/opt/hlhdf

# find the location of the proj4 headers and binaries
# for OSX 10.10.5 using Macports the location is /opt/local/lib/proj47
mdfind -name projects.h
export PROJ4ROOT=/opt/local/lib/proj47

# now we're ready to configure the install
./configure --prefix=${RADAR_ROOT_DIR}/opt/rave  --with-proj=${PROJ4ROOT} --with-hlhdf=${HLDIR}

make 

#
make test

#
make install

#
echo "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:${RADAR_ROOT_DIR}/opt/rave/lib" >>${RADAR_ROOT_DIR}/setup-env
echo "export PATH=\${PATH}:${RADAR_ROOT_DIR}/opt/rave/bin" >>${RADAR_ROOT_DIR}/setup-env

# cd back to the project root slash source
cd ${RADAR_ROOT_DIR}/src

# get a copy of vol2bird
git clone https://github.com/adokter/vol2bird.git

# change directory into it
cd vol2bird

# configure vol2bird 
./configure --prefix=${RADAR_ROOT_DIR}/opt/vol2bird --with-rave=${RADAR_ROOT_DIR}/opt/rave --with-confuse=/opt/local

#
make

#
make install

# add vol2bird stuff to PATH
echo "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:${RADAR_ROOT_DIR}/opt/vol2bird/lib" >>${RADAR_ROOT_DIR}/setup-env
echo "export PATH=\${PATH}:${RADAR_ROOT_DIR}/opt/vol2bird/bin" >>${RADAR_ROOT_DIR}/setup-env

# activate the environment (do this first every time you use vol2bird)
cd ${RADAR_ROOT_DIR}
. setup-env

```

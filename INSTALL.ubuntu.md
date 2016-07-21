```bash
# instructions for Ubuntu-based systems

# cd to project directory

# define the root diretory of the project
RADAR_ROOT_DIR=${PWD}

# prepare a directory structure:
mkdir ${RADAR_ROOT_DIR}/opt
mkdir ${RADAR_ROOT_DIR}/src

# install python 2 (not python 3)
sudo apt-get install python

# we'll use Pythopn virtualenv to manage the python version, and this project's dependencies
sudo apt-get install python-virtualenv

# install dependencies for rave (g++, gcc, make, etc)
sudo apt-get install build-essential

# install zlib (gzip archiving library)
sudo apt-get install zlib1g-dev

# install Tcl and Tk for user interface
sudo apt-get install tcl-dev
sudo apt-get install tk-dev

# install HDF5, Hierarchichal Data Format library
sudo apt-get install libhdf5-dev

# install dependencies for pycurl
sudo apt-get install libcurl4-gnutls-dev
sudo apt-get install libgnutls-dev

# install library that can handle different projections
sudo apt-get install libproj-dev

# (optional) install package for xmlparsing of projection definitions
sudo apt-get install libexpat1-dev

# (optional) install package for auto-generating documentation
sudo apt-get install doxygen
sudo apt-get install texlive-font-utils

# install library for parsing options
sudo apt-get install libconfuse-dev

# install python inotify
sudo apt-get install python-pyinotify

# update the system database to be able to locate the files later:
sudo updatedb

# create a python2 virtualenv
virtualenv -p /usr/bin/python2.7 ${RADAR_ROOT_DIR}/.venv

# activate the python virtual environment if necessary (don't forget the leading dot):
. ${RADAR_ROOT_DIR}/.venv/bin/activate

# install libicu for unicode support
sudo apt-get install libicu55 libicu-dev

# update pip to the latest version
pip install --upgrade pip

# install Numpy into the virtual environment
pip install numpy

# install Python Imaging Library (Pillow) into the virtual environment
pip install Pillow

# install python inotify for watching changes on files
pip install pyinotify

# change to the source directory...
cd ${RADAR_ROOT_DIR}/src

# you'll need a copy of hlhdf:
git clone git://git.baltrad.eu/hlhdf.git

# change directory into the clone of hlhdf
cd hlhdf

# configure hlhdf
# we need to point to the location of the hdf5 headers and binaries
# find out where the headers live on your system with
locate hdf5.h
# for Ubuntu 14.04, the location is /usr/include/
# for Ubuntu 16.04, the location is /usr/include/hdf5/serial (might be different on your system)
#
# and the same for the binaries:
locate libhdf5.a  (static)
locate libhdf5.so (dynamic)
# for Ubuntu 14.04, the location is /usr/lib/x86_64-linux-gnu/
# for Ubuntu 16.04, the location is /usr/lib/x86_64-linux-gnu/hdf5/serial/ (might be different on your system)

./configure --prefix=${RADAR_ROOT_DIR}/opt/hlhdf --with-hdf5=/usr/include/,/usr/lib/x86_64-linux-gnu

# compile hlhdf5 in the local directory
make

# (optional) compile hlhdf5 tests and run them
make test

# copy the binaries to where we want them (${RADAR_ROOT_DIR}/opt/hlhdf):
make install

# create a script that sets up the environment so that Python can find _pyhlmodule.so and libhlhdf.so
echo ". ${RADAR_ROOT_DIR}/.venv/bin/activate" >${RADAR_ROOT_DIR}/setup-env
echo "export PYTHONPATH=\${PYTHONPATH}:${RADAR_ROOT_DIR}/opt/hlhdf/lib" >>${RADAR_ROOT_DIR}/setup-env
echo "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:${RADAR_ROOT_DIR}/opt/hlhdf/lib" >>${RADAR_ROOT_DIR}/setup-env

# activate the environment (do this first every time you use vol2bird)
cd ${RADAR_ROOT_DIR}
. setup-env

# cd to the source directory
cd ${RADAR_ROOT_DIR}/src

# get a copy of rave:
git clone git://git.baltrad.eu/rave.git

# cd into rave source directory
cd rave

# rave needs keyczar to manage keys
pip install python-keyczar

# rave needs jprops to run PGF plugins
# you can skip this if you don't need the PGF plugin
pip install jprops

# rave needs sqlalchemy to run PGF plugins
# you can skip this if you don't need the PGF plugin
pip install sqlalchemy
pip install sqlalchemy_migrate

# rave needs psycopg2 to run PGF plugins
# you can skip this if you don't need the PGF plugin
pip install psycopg2

# install package for downloading
pip install pycurl


# sanity checks
python -c 'import numpy'
python -c 'import _pyhl'

# set some environment variables
export RAVEROOT=${RADAR_ROOT_DIR}/opt/rave
export NUMPYDIR=${RADAR_ROOT_DIR}/.venv/lib/python2.7/site-packages/numpy/core/include/numpy
export HLDIR=${RADAR_ROOT_DIR}/opt/hlhdf
export PROJ4ROOT=/usr

# for building extensions, python needs to know about config; virtualenv does not set that
# up, so we'll symlink to the system python ourselves
# (64-bit)
ln -s /usr/lib/python2.7/config-x86_64-linux-gnu ${RADAR_ROOT_DIR}/.venv/lib/python2.7/config
# (32-bit)
# ln -s /usr/lib/python2.7/config ${RADAR_ROOT_DIR}/.venv/lib/python2.7/config

# now we're ready to configure the install
./configure --prefix=${RADAR_ROOT_DIR}/opt/rave --with-expat

#
make

#
make test

#
make install

#
echo "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:${RADAR_ROOT_DIR}/opt/rave/lib" >>${RADAR_ROOT_DIR}/setup-env
echo "export PATH=\${PATH}:${RADAR_ROOT_DIR}/opt/rave/bin" >>${RADAR_ROOT_DIR}/setup-env

#
echo "${RADAR_ROOT_DIR}/opt/rave/Lib" > ${RADAR_ROOT_DIR}/.venv/lib/python2.7/site-packages/rave.pth


# cd back to the project root slash source
cd ${RADAR_ROOT_DIR}/src

# get a copy of vol2bird
git clone https://github.com/adokter/vol2bird.git

# change directory into it
cd vol2bird

# configure vol2bird
./configure --prefix=${RADAR_ROOT_DIR}/opt/vol2bird --with-rave=${RADAR_ROOT_DIR}/opt/rave

#
make

# optional: run tests
# note that we need to set some directories explicitly here, as there is no good
# way to autodetect these, and the vol2bird tests default to assuming you have RAVE
# installed as part of a full standard Baltrad installation
make test PREFIX=${RADAR_ROOT_DIR}/opt RAVEPYTHON=${RADAR_ROOT_DIR}/.venv2/bin/python

#
make install

# add vol2bird stuff to PATH
echo "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:${RADAR_ROOT_DIR}/opt/vol2bird/lib" >>${RADAR_ROOT_DIR}/setup-env
echo "export PATH=\${PATH}:${RADAR_ROOT_DIR}/opt/vol2bird/bin" >>${RADAR_ROOT_DIR}/setup-env


# activate the environment (do this first every time you use vol2bird)
cd ${RADAR_ROOT_DIR}
. setup-env

```


# Instructions for installing node-installer

```
# install instructions for Lubuntu 14.04 64-bit

sudo apt-get install python
sudo apt-get install zlib1g-dev
sudo apt-get install postgresql-client-9.3
sudo apt-get install postgresql-client-common
sudo apt-get install libpq-dev
sudo apt-get install libssl-dev
sudo apt-get install libsasl2-dev
sudo apt-get install libncurses5-dev
sudo apt-get install libdb-dev
sudo apt-get install libicu52
sudo apt-get install libicu-dev
sudo apt-get install libbz2-dev
sudo apt-get install openssl libssl-dev
sudo apt-get install doxygen
sudo apt-get install libpng-dev
sudo apt-get install libfreetype6 libfreetype6-dev
sudo apt-get install libssl-dev openssl dpkg-dev


# find if you have Java installed
java -version

# edit /etc/postgresql/9.3/main/pg_hba.conf:
change the line that says:
local   all             all                                     peer
to
local   all             all                                     md5

# restart the postgresql service
sudo service postgresql restart

# add a system user, set its password to baltrad
sudo adduser baltrad
sudo passwd baltrad

# double 'sudo':
sudo sudo -u postres psql -c "CREATE USER baltrad with PASSWORD 'baltrad';CREATE DATABASE baltrad with OWNER baltrad;"

sudo mkdir /opt/software
sudo mkdir /opt/baltrad
sudo chown -R baltrad:baltrad /opt/software
sudo chown -R baltrad:baltrad /opt/baltrad

# become the baltrad user:
sudo su - baltrad

#
cd /opt/software

# use git to download the baltrad node installer
git clone git://git.baltrad.eu/node-installer.git

# cd into the newly created directory
cd node-installer

# check if you already have an instance of tomcat running:
ps -ef | grep tomcat
# if you have one, stop it with (if your tomcat lives in /etc/init.d/):
sudo /etc/init.d/tomcat7 stop



# And now, the 'one-line installer'. Very easy...not.
./setup --tomcatpwd=baltrad --nodename=nl.esciencecenter.balthazar --prefix=/opt/baltrad \
--jdkhome=/usr/lib/jvm/java-7-openjdk-amd64 \
--with-psql=/usr/include/postgresql/,/usr/lib/ \
--dbpwd=baltrad \
--with-rave \
--experimental \
install


# if something goes wrong during installation, just rm -rf the /opt/baltrad
# directory (or the directory that you used as prefix), then repeat steps.

# importing your own local data from file goes with the program /opt/baltrad/rave/bin/odim_injector. First add these directories to the LD_LIBRARY_PATH:
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/baltrad/rave/lib:/opt/baltrad/hlhdf/lib/:/opt/baltrad/third_party/lib

# there's still some weirdness with permissions... we haven't figured out exactly how it should be. Currently
# running everything as root (bad idea):
(.venv)daisycutter@daisycutter-NLeSC:/opt/baltrad$ sudo bash -c '. /opt/baltrad/etc/bltnode.rc; odim_injector -i /home/daisycutter/projects/birdradar2016/odim-injector-watched'
# then copy an ODIM-HDF5 file into the directory, it should get picked up and show up in the messages tab of the Baltrad web interface.


```


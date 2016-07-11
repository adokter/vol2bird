include /opt/baltrad/rave/mkf/def.mk

RAVE_ROOT_DIR=      /opt/baltrad/rave

# Install prefix
prefix=             /opt/baltrad/vol2bird

# Most configuration variables are inherited from
# raves and hlhdfs configuration files so you can
# basically uncomment any of the ones below if you want.

CC=                 		/usr/bin/clang
CCOPTS=             -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes
LDFLAGS=            -L/opt/local/lib -Wl,-headerpad_max_install_names -L/opt/local/lib/db48 

# (from hlhdf)
LDSHARED=           /usr/bin/clang -dynamiclib -undefined dynamic_lookup -L/opt/local/lib -Wl,-headerpad_max_install_names -L/opt/local/lib/db48 

SHARED_FLAG=        

RAVE_INCLUDE_FLAG=  -I/opt/baltrad/rave/include
RAVE_LIB_FLAG=      -L/opt/baltrad/rave/lib

PYTHON_INCLUDE_FLAG=-I/opt/local/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7
NUMPY_INCLUDE_FLAG= -I/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/include/numpy

PROJ_INCLUDE_FLAG=  -I/opt/local/lib/proj47/include
PROJ_LIBRARY_FLAG=  -L/opt/local/lib/proj47/lib

EXPAT_INCLUDE_FLAG= 
EXPAT_LIBRARY_FLAG= 
EXPAT_SUPPRESSED=   yes

HLHDF_LIBRARY_FLAG= -L/opt/baltrad/HLHDF//lib
HLHDF_INCLUDE_FLAG= -I/opt/baltrad/HLHDF//include

# confuse

CONFUSE_LIBRARY_FLAG=-L/opt/local/lib
CONFUSE_INCLUDE_FLAG=-I/opt/local/include

# Special flag to be used for printouts of the necessary LD_LIBRARY_PATH
#
LD_PRINTOUT=        /opt/baltrad/vol2bird/lib:/opt/baltrad/rave/lib:/opt/baltrad/HLHDF//lib:/opt/local/lib

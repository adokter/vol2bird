#!/bin/sh
############################################################
# Description: Script that executes a Python script with proper
# settings in this build.
#
# Author(s):   Anders Henja
#
# Copyright:   Swedish Meteorological and Hydrological Institute, 2011
#
# History:  2011-08-30 Created by Anders Henja
############################################################
SCRFILE=`python -c "import os;print os.path.abspath(\"$0\")"`
SCRIPTPATH=`dirname "$SCRFILE"`

DEF_MK_FILE="${SCRIPTPATH}/../def.mk"

if [ ! -f "${DEF_MK_FILE}" ]; then
  echo "configure has not been run"
  exit 255
fi

RESULT=0

# RAVES MKF FILE
RAVE_ROOT_DIR=`fgrep RAVE_ROOT_DIR "${DEF_MK_FILE}" | sed -e"s/\(RAVE_ROOT_DIR=[ \t]*\)//"`
RAVE_ROOT_MKFILE="$RAVE_ROOT_DIR/mkf/def.mk"

# HLHDFS MKF FILE
HLHDF_MKFFILE=`fgrep HLHDF_HLDEF_MK_FILE "${RAVE_ROOT_MKFILE}" | sed -e"s/\(HLHDF_HLDEF_MK_FILE=[ \t]*\)//"`

# Get HDF5s ld path from hlhdfs mkf file
HDF5_LDPATH=`fgrep HDF5_LIBDIR "${HLHDF_MKFFILE}" | sed -e"s/\(HDF5_LIBDIR=[ \t]*-L\)//"`

# Get HLHDFs libpath from raves mkf file
HLHDF_LDPATH=`fgrep HLHDF_LIB_DIR "${DEF_MK_FILE}" | sed -e"s/\(HLHDF_LIB_DIR=[ \t]*\)//"`

BNAME=`python -c 'from distutils import util; import sys; print "lib.%s-%s" % (util.get_platform(), sys.version[0:3])'`

RBPATH="${SCRIPTPATH}/../pyvol2bird"
RAVE_LDPATH="${RAVE_ROOT_DIR}/lib"
RACK_LDPATH="${SCRIPTPATH}/../lib"
XRUNNERPATH="${SCRIPTPATH}/../test/lib"

# Special hack for mac osx.
ISMACOS=no
case `uname -s` in
 Darwin*)
   ISMACOS=yes
   ;;
 darwin*)
   ISMACOS=yes
   ;;
esac

if [ "x$ISMACOS" = "xyes" ]; then
  if [ "$DYLD_LIBRARY_PATH" != "" ]; then
    export DYLD_LIBRARY_PATH="${RACK_LDPATH}:${RAVE_LDPATH}:${HLHDF_LDPATH}:${HDF5_LDPATH}:${LD_LIBRARY_PATH}"
  else
    export DYLD_LIBRARY_PATH="${RACK_LDPATH}:${RAVE_LDPATH}:${HLHDF_LDPATH}:${HDF5_LDPATH}"
  fi
else
  if [ "$LD_LIBRARY_PATH" != "" ]; then
    export LD_LIBRARY_PATH="${RACK_LDPATH}:${RAVE_LDPATH}:${HLHDF_LDPATH}:${HDF5_LDPATH}:${LD_LIBRARY_PATH}"
  else
    export LD_LIBRARY_PATH="${RACK_LDPATH}:${RAVE_LDPATH}:${HLHDF_LDPATH}:${HDF5_LDPATH}"
  fi
fi

export BEAMBPATH="${RAVE_ROOT_DIR}/Lib:${HLHDF_LDPATH}:${RBPATH}:${XRUNNERPATH}"

if test "${PYTHONPATH}" != ""; then
  export PYTHONPATH="${BEAMBPATH}:${PYTHONPATH}"
else
  export PYTHONPATH="${BEAMBPATH}"
fi

# Syntax: run_python_script <pyscript> [<dir> - if script should be executed in a particular directory]

NARGS=$#
PYSCRIPT=
DIRNAME=
if [ $NARGS -eq 1 ]; then
  PYSCRIPT=`python -c "import os;print os.path.abspath(\"$1\")"`
elif [ $NARGS -eq 2 ]; then
  PYSCRIPT=`python -c "import os;print os.path.abspath(\"$1\")"`
  DIRNAME="$2"
elif [ $NARGS -eq 0 ]; then
  # Do nothing
  PYSCRIPT=
  DIRNAME=
else
  echo "Unknown command"
  exit 255
fi

if [ "$DIRNAME" != "" ]; then
  cd "$DIRNAME"
fi

if [ "$PYSCRIPT" != "" ]; then
  #valgrind --leak-check=full --show-reachable=yes 
  python "$PYSCRIPT"
else
  python
fi

VAL=$?
if [ $VAL != 0 ]; then
  RESULT=$VAL
fi

# EXIT WITH A STATUS CODE, 0 == OK, ANY OTHER VALUE = FAIL
exit $RESULT


#!/bin/bash
if [ $# = 0 ]; then
    echo "usage: $0 <polar volume> [<ODIM hdf5 profile output> [<ODIM hdf5 volume output>]]"
    echo "   runs vol2bird from docker container"
    echo "   expects OPERA ODIM hdf5 input format, see http://www.eumetnet.eu/opera-software"
    echo "   or input formats compatible with RSL, see http://trmm-fc.gsfc.nasa.gov/trmm_gv/software/rsl"
else
    docker exec vol2bird bash -c "cd data && vol2bird $1 $2 $3"
fi


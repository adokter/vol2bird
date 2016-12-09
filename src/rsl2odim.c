/*
 * Copyright 2016 - Adriaan Dokter
 *
 * This program converts RSL library compatible radar volumes to ODIM hdf5
 *
 */


#include <stdio.h>
#include "rave_io.h"
#include "polarvolume.h"
#include "libvol2bird.h"

int main(int argc, char** argv) {
//    cfg_t* cfg;

    // print default message when no input arguments
    if (argc == 1) {
        fprintf(stderr,"usage: %s <RSL polar volume input> <ODIM hdf5 volume output>\n",argv[0]);
        fprintf(stderr,"   expects input formats compatible with RSL, see <http://trmm-fc.gsfc.nasa.gov/trmm_gv/software/rsl>\n");
        fprintf(stderr,"   outputs OPERA ODIM hdf5 input format, see <http://www.eumetnet.eu/opera-software>\n");
        return -1;
    }
    // check to see if we have the right number of input arguments
    if (argc != 3) {
        fprintf(stderr, "Error: two arguments required: RSL volume in, and ODIM volume out\n");
        return -1;
    }
    
    // the polar volume file that the user provided as input
    char* fileVolIn = argv[1];
    // the volume file that the user specified as output
    const char* fileVolOut = argv[2];
        
    // read in data up to a distance of alldata.misc.rCellMax
    // we do not read in the full volume for speed/memory
    PolarVolume_t* volume = NULL;
    volume = vol2birdGetVolume(fileVolIn, NAN);
    
    if (volume == NULL) {
        fprintf(stderr,"Error: failed to read radar volume\n");
        return -1;
    }
    
    saveToODIM((RaveCoreObject*) volume, fileVolOut);
    
    RAVE_OBJECT_RELEASE(volume);
     
    return 0;

}
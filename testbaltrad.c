/*
 * Copyright 2013 Netherlands eScience Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */



#include <stdio.h>
#include <time.h>
#include "rave_io.h"
#include "polarvolume.h"
#include "libvol2bird.h"


int main(int argc, char** argv) {


    // check to see if we have the right number of input arguments
    if (argc != 2) {
        fprintf(stderr, "Only one argument is allowed\n");
        return -1;
    }
    
    // ------------------------------------------------------------- //
    //                initialization of variables                    //
    // ------------------------------------------------------------- //

    // the filename that the user provided as input
    char* filename = argv[1];
    
    // read the input file and assign it to a generic rave object
    RaveIO_t* raveio = RaveIO_open(filename);
    
    if (RaveIO_getObjectType(raveio) == Rave_ObjectType_PVOL) {

        // initialize array used for performance analysis
        struct timespec ts = { 0 };

        // the if statement above tests whether we are dealing with a 
        // PVOL object, so we can safely cast the generic object to
        // the PolarVolume_t type:
        PolarVolume_t* volume = NULL;
        volume = (PolarVolume_t*) RaveIO_getObject(raveio);

        // initialize volbird library
        int initSuccessful = vol2birdSetUp(volume) == 0;
        
        if (initSuccessful == FALSE) {
            return -1;
        }

        // call vol2bird's main routine
        vol2birdCalcProfiles();
        
        
        // ------------------------------------------------------------------- //
        //  example of how the getters can be used to get at the profile data  //
        // ------------------------------------------------------------------- //
        
        {  // getter example scope begin

            int iProfileType;
            
            for (iProfileType = 1; iProfileType <= 3;iProfileType++) {

                int nRowsProfile = vol2birdGetNRowsProfile();
                int nColsProfile = vol2birdGetNColsProfile();

                fprintf(stderr, "\n--------------------------\n\n");
                
                float *profileCopy;
                profileCopy = (float*) malloc(sizeof(float) * nRowsProfile * nColsProfile);

                profileCopy = vol2birdGetProfile(iProfileType);
                
                int iRowProfile;
                int iColProfile;
                int iCopied = 0;
                
                for (iRowProfile = 0; iRowProfile < nRowsProfile; iRowProfile++) {
                    for (iColProfile = 0; iColProfile < nColsProfile; iColProfile++) {
                        fprintf(stderr," %10.2f",profileCopy[iCopied]);
                        iCopied += 1;
                    }
                    fprintf(stderr,"\n");
                }
                
                profileCopy = NULL;
                free((void*) profileCopy);

            }
        } // getter example scope end



        // ------------------------------------------------------------------- //
        //                 end of the getter example section                   //
        // ------------------------------------------------------------------- //            


        // tear down vol2bird, give memory back
        vol2birdTearDown();

        // output some performance data
        clock_gettime(CLOCK_REALTIME, &ts);
        double nSeconds = ((double) ts.tv_nsec)/1e9;
        fprintf(stderr, "Processing done in %.2f seconds\n",nSeconds);

    }


    RAVE_OBJECT_RELEASE(raveio);
    
    return 0;

}







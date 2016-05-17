/*
 * Copyright 2015 
 *
 * This program calculates Vertical Profiles of Birds (VPBs) as described in
 *
 * Bird migration flight altitudes studied by a network of operational weather radars
 * Dokter A.M., Liechti F., Stark H., Delobbe L., Tabary P., Holleman I.
 * J. R. Soc. Interface, 8, 30–43, 2011
 * DOI: 10.1098/rsif.2010.0116
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
//    cfg_t* cfg;
    vol2bird_t alldata;

    // print default message when no input arguments
    if (argc == 1) {
        fprintf(stderr,"usage: %s <ODIM hdf5 volume> [<ODIM hdf5 profile output>] \n",argv[0]);
        fprintf(stderr,"   Version 0.2.2 (17-May-2016)\n");
        fprintf(stderr,"   expects OPERA ODIM hdf5 input format, see http://www.eumetnet.eu/opera-software\n\n");
        fprintf(stderr,"   Output fields to stdout:\n");
        fprintf(stderr,"   Date    - date in UTC\n");
        fprintf(stderr,"   Time    - time in UTC\n");
        fprintf(stderr,"   U       - speed component west to east [m/s]\n");
        fprintf(stderr,"   V       - speed component north to south [m/s]\n");
        fprintf(stderr,"   W       - vertical speed (unreliable!) [m/s]\n");
        fprintf(stderr,"   Speed   - horizontal speed [m/s]\n");
        fprintf(stderr,"   Direc   - direction [degrees, clockwise from north]\n");
        fprintf(stderr,"   Stdev   - VVP radial velocity standard deviation direction [m/s]\n");
        fprintf(stderr,"   Gap     - Angular data gap detected [T/F]\n");
        fprintf(stderr,"   dBZ     - Bird reflectivity factor [dBZ]\n");
        fprintf(stderr,"   eta     - Bird reflectivity [cm^2/km^3]\n");
        fprintf(stderr,"   DensBird- Bird density [birds/km^3]\n");
        fprintf(stderr,"   dBZAll  - Total reflectivity factor (bio+meteo scattering) [dBZ]\n");
        fprintf(stderr,"   nPts    - number of points VVP bird velocity analysis\n");
        fprintf(stderr,"   nPtsZ   - number of points bird density estimate \n");
        fprintf(stderr,"   nPtsAll - number of points VVP velocity Stdev analysis\n");
        return -1;
    }
    // check to see if we have the right number of input arguments
    if (argc > 3) {
        fprintf(stderr, "Only one or two arguments are allowed\n");
        return -1;
    }
    
    // ------------------------------------------------------------- //
    //                initialization of variables                    //
    // ------------------------------------------------------------- //

    // the filename that the user provided as input
    const char* filename = argv[1];
    const char* fileout;
    
    if (argc == 3){
        fileout = argv[2];
    }
    else{
        fileout = NULL;
    }
  
    
    // read the input file and assign it to a generic rave object
    RaveIO_t* raveio = RaveIO_open(filename);

    if (raveio == (RaveIO_t*) NULL){
        fprintf(stderr, "critical error, cannot open file %s\n", filename);
        return -1;
    }

    if (RaveIO_getObjectType(raveio) == Rave_ObjectType_PVOL) {

        // initialize array used for performance analysis
        //struct timespec ts = { 0 };

        // the if statement above tests whether we are dealing with a 
        // PVOL object, so we can safely cast the generic object to
        // the PolarVolume_t type:
        PolarVolume_t* volume = NULL;
        volume = (PolarVolume_t*) RaveIO_getObject(raveio);

        int configSuccessful = vol2birdLoadConfig(&alldata) == 0;

        if (configSuccessful == FALSE) {
            return -1;
        }

        // initialize volbird library
        int initSuccessful = vol2birdSetUp(volume, &alldata) == 0;
        
        if (initSuccessful == FALSE) {
            return -1;
        }

        // call vol2bird's main routine
        vol2birdCalcProfiles(&alldata);
        
        
        // ------------------------------------------------------------------- //
        //  example of how the getters can be used to get at the profile data  //
        // ------------------------------------------------------------------- //
        const char* date;
        const char* time;
        const char* source;

        date = PolarVolume_getDate(volume);
        time = PolarVolume_getTime(volume);
        source = PolarVolume_getSource(volume);

        {  // getter example scope begin

                int nRowsProfile = vol2birdGetNRowsProfile(&alldata);
                int nColsProfile = vol2birdGetNColsProfile(&alldata);

                fprintf(stdout, "# vol2bird vertical profile\n");
                fprintf(stdout, "# source: %s\n",source);
                fprintf(stdout, "# ODIM HDF5 input: %s\n",filename);
                printf("# Date   Time Heig    U      V       W   Speed Direc StdDev Gap dBZ     eta DensBird dBZAll   n   ndBZ  nAll nAlldBZ\n");
               
                float *profileBio;
                float *profileAll;

                profileBio = vol2birdGetProfile(1, &alldata);
                profileAll = vol2birdGetProfile(3, &alldata);
                
                int iRowProfile;
                int iCopied = 0;
                
                for (iRowProfile = 0; iRowProfile < nRowsProfile; iRowProfile++) {
                    iCopied=iRowProfile*nColsProfile;
                    printf("%8s %.4s ",date,time);
                    printf("%4.f %6.2f %6.2f %7.2f %5.2f %5.1f %6.2f %1c %6.2f %6.1f %6.2f %6.2f %5.f %5.f %5.f %5.f\n",
                    profileBio[0+iCopied],
                    profileBio[2+iCopied],profileBio[3+iCopied],
                    profileBio[4+iCopied],profileBio[5+iCopied],
                    profileBio[6+iCopied],profileAll[7+iCopied],
                    profileBio[8+iCopied] == TRUE ? 'T' : 'F',
                    profileBio[9+iCopied],profileBio[11+iCopied],
                    profileBio[12+iCopied],profileAll[9+iCopied],
                    profileBio[10+iCopied],profileBio[13+iCopied],
                    profileAll[10+iCopied],profileAll[13+iCopied]);
                }
                
                profileAll = NULL;
                profileBio = NULL;
                free((void*) profileAll);
                free((void*) profileBio);

            //}
        } // getter example scope end



        // ------------------------------------------------------------------- //
        //                 end of the getter example section                   //
        // ------------------------------------------------------------------- //            
            
        //map vol2bird profile data to Rave profile object
        mapDataToRave(volume, &alldata);
        
        //save rave profile to ODIM hdf5 file
        if (fileout != NULL){
            int result;
            result = saveToODIM(alldata.vp, fileout); 
            if (result == FALSE){
                fprintf(stderr, "critical error, cannot write file %s\n", fileout);
                return -1;
            }
        }
        
        // tear down vol2bird, give memory back
        vol2birdTearDown(&alldata);

        // output some performance data
        //clock_gettime(CLOCK_REALTIME, &ts);
        //double nSeconds = ((double) ts.tv_nsec)/1e9;
        //fprintf(stderr, "Processing done in %.2f seconds\n",nSeconds);

    }

    RAVE_OBJECT_RELEASE(raveio);
    
    return 0;

}







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
#include "constants.h"


int main(int argc, char** argv) {
//    cfg_t* cfg;
    vol2bird_t alldata;

    // print default message when no input arguments
    if (argc == 1) {
        fprintf(stderr,"usage: %s <polar volume> [<ODIM hdf5 profile output> [<ODIM hdf5 volume output>]]\n",argv[0]);
        fprintf(stderr,"   Version %s (%s)\n", VERSION, VERSIONDATE);
        fprintf(stderr,"   expects OPERA ODIM hdf5 input format, see <http://www.eumetnet.eu/opera-software>\n");
        fprintf(stderr,"   or input formats compatible with RSL, see <http://trmm-fc.gsfc.nasa.gov/trmm_gv/software/rsl>\n\n");
        fprintf(stderr,"   Output fields to stdout:\n");
        fprintf(stderr,"   date      - date [UTC]\n");
        fprintf(stderr,"   time      - time [UTC]\n");
        fprintf(stderr,"   HGHT      - height above mean sea level [m]. Alt. bin from HGHT to HGHT+interval)\n");
        fprintf(stderr,"   u         - speed component west to east [m/s]\n");
        fprintf(stderr,"   v         - speed component north to south [m/s]\n");
        fprintf(stderr,"   w         - vertical speed (unreliable!) [m/s]\n");
        fprintf(stderr,"   ff        - horizontal speed [m/s]\n");
        fprintf(stderr,"   dd        - direction [degrees, clockwise from north]\n");
        fprintf(stderr,"   sd_vvp    - VVP radial velocity standard deviation [m/s]\n");
        fprintf(stderr,"   gap       - Angular data gap detected [T/F]\n");
        fprintf(stderr,"   dbz       - Bird reflectivity factor [dBZ]\n");
        fprintf(stderr,"   eta       - Bird reflectivity [cm^2/km^3]\n");
        fprintf(stderr,"   dens      - Bird density [birds/km^3]\n");
        fprintf(stderr,"   DBZH      - Total reflectivity factor (bio+meteo scattering) [dBZ]\n");
        fprintf(stderr,"   n         - number of points VVP bird velocity analysis (u,v,w,ff,dd)\n");
        fprintf(stderr,"   n_dbz     - number of points bird density estimate (dbz,eta,dens)\n");
        fprintf(stderr,"   n_all     - number of points VVP st.dev. estimate (sd_vvp)\n");
        fprintf(stderr,"   n_dbz_all - number of points total reflectivity estimate (DBZH)\n\n");
        fprintf(stderr,"   Report bugs to: a.m.dokter@uva.nl\n");
        fprintf(stderr,"   vol2bird home page: <http://github.com/adokter/vol2bird>\n");
        return -1;
    }
    // check to see if we have the right number of input arguments
    if (argc > 4) {
        fprintf(stderr, "Up to three arguments are allowed\n");
        return -1;
    }
    
    // ------------------------------------------------------------- //
    //                initialization of variables                    //
    // ------------------------------------------------------------- //

    // the polar volume file that the user provided as input
    char* fileVolIn = argv[1];
    // the (optional) vertical profile file that the user specified as output
    const char* fileVpOut;
    // the (optional) vertical profile file that the user specified as output
    const char* fileVolOut;

    if (argc == 3){
        fileVpOut = argv[2];
        fileVolOut = NULL;
    }
    else if (argc == 4){
        fileVpOut = argv[2];
        fileVolOut = argv[3];
    }
    else{
        fileVpOut = NULL;
        fileVolOut = NULL;
    }
    
    // read configuration options
    int configSuccessful = vol2birdLoadConfig(&alldata) == 0;

    if (configSuccessful == FALSE) {
        fprintf(stderr,"Error: failed to load configuration\n");
        return -1;
    }
    
    // read in data up to a distance of alldata.misc.rCellMax
    // we do not read in the full volume for speed/memory
    PolarVolume_t* volume = NULL;
    
    volume = vol2birdGetVolume(fileVolIn, alldata.misc.rCellMax,1);
         
    if (volume == NULL) {
        fprintf(stderr,"Error: failed to read radar volume\n");
        return -1;
    }
    
    // loading static clutter map upon request
    if (alldata.options.useClutterMap){
        int clutterSuccessful = vol2birdLoadClutterMap(volume, alldata.options.clutterMap,alldata.misc.rCellMax) == 0;
        
        if (clutterSuccessful == FALSE) {
            fprintf(stderr,"Error: failed to load static clutter map '%s', aborting\n",alldata.options.clutterMap);
            return -1;
        }

    }
    
    // resample the volume upon request
    if (alldata.options.resample) {
        PolarVolume_t* volume_orig = volume;
        volume = PolarVolume_resample(volume,alldata.options.resampleRscale,
                 alldata.options.resampleNbins,alldata.options.resampleNrays);
        if (volume == NULL) {
            fprintf(stderr,"Error: volume resampling failed\n");
            return -1;
        }
        RAVE_OBJECT_RELEASE(volume_orig);
    }

    // initialize volbird library
    int initSuccessful = vol2birdSetUp(volume, &alldata) == 0;

    if (initSuccessful == FALSE) {
        fprintf(stderr,"Error: failed to initialize vol2bird\n");
        return -1;
    }

    // output (optionally de-aliased) volume
    if (fileVolOut != NULL){
        saveToODIM((RaveCoreObject*) volume, fileVolOut);
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
            
            fprintf(stdout, "# vol2bird Vertical Profile of Birds (VPB)\n");
            fprintf(stdout, "# source: %s\n",source);
            fprintf(stdout, "# ODIM HDF5 input: %s\n",fileVolIn);
            printf("# date   time HGHT    u      v       w     ff    dd  sd_vvp gap dbz     eta   dens   DBZH   n   n_dbz n_all n_dbz_all\n");
           
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
                nanify(profileBio[2+iCopied]),nanify(profileBio[3+iCopied]),
                nanify(profileBio[4+iCopied]),nanify(profileBio[5+iCopied]),
                nanify(profileBio[6+iCopied]),nanify(profileAll[7+iCopied]),
                profileBio[8+iCopied] == TRUE ? 'T' : 'F',
                nanify(profileBio[9+iCopied]),nanify(profileBio[11+iCopied]),
                nanify(profileBio[12+iCopied]),nanify(profileAll[9+iCopied]),
                nanify(profileBio[10+iCopied]),nanify(profileBio[13+iCopied]),
                nanify(profileAll[10+iCopied]),nanify(profileAll[13+iCopied]));
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
    if (fileVpOut != NULL){
        int result;
        result = saveToODIM((RaveCoreObject*) alldata.vp, fileVpOut); 
        if (result == FALSE){
            fprintf(stderr, "critical error, cannot write file %s\n", fileVpOut);
            return -1;
        }
    }
    
    // tear down vol2bird, give memory back
    vol2birdTearDown(&alldata);
    RAVE_OBJECT_RELEASE(volume);


    // output some performance data
    //clock_gettime(CLOCK_REALTIME, &ts);
    //double nSeconds = ((double) ts.tv_nsec)/1e9;
    //fprintf(stderr, "Processing done in %.2f seconds\n",nSeconds);


    
    return 0;

}







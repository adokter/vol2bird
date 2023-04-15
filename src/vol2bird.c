/** vol2bird main executable
 * @file vol2bird.c
 * @author Adriaan Dokter & Netherlands eScience Centre
 * @date 2015-05-05
 */

/*
 * Copyright 2017 Adriaan Dokter
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
#include <getopt.h>
#include <string.h>
#include "rave_io.h"
#include "polarvolume.h"
#include "libvol2bird.h"
#include "constants.h"
#include "hlhdf.h"
#include "hlhdf_debug.h"
#include "rave_debug.h"
#include <vertical_profile.h>

void usage(char* programName, int verbose){
    fprintf(stderr,"vol2bird version %s (%s)\n", VERSION, VERSIONDATE);
    fprintf(stderr,"   usage: %s <polar volume> [<ODIM hdf5 profile output> [<ODIM hdf5 volume output>]]\n",programName);
    fprintf(stderr,"   usage: %s -i <polar volume or scan> [-i <polar scan> [-i <polar scan>] ...] [-o <ODIM hdf5 profile output>] [-p <ODIM hdf5 volume output>] [-c <vol2bird configuration file>]\n",programName);
    fprintf(stderr,"   usage: %s --help\n", programName);

    if(verbose){

        fprintf(stderr,"\n   Supported radar data formats:\n");
        fprintf(stderr,"   * OPERA ODIM hdf5 input format, see <https://www.eumetnet.eu/wp-content/uploads/2019/01/ODIM_H5_v23.pdf> [enabled]\n");
        fprintf(stderr,"   * input formats compatible with RSL, see <http://trmm-fc.gsfc.nasa.gov/trmm_gv/software/rsl>");
        #ifdef RSL 
        fprintf(stderr, " [enabled]\n");
        #endif
        #ifndef RSL 
        fprintf(stderr, " [disabled]\n");
        #endif
        fprintf(stderr,"   * Vaisala Sigmet IRIS format, see <ftp://ftp.sigmet.com/outgoing/manuals/IRIS_Programmers_Manual.pdf>");
        #ifdef  IRIS
        fprintf(stderr, " [enabled]\n\n");
        #endif
        #ifndef IRIS
        fprintf(stderr, " [disabled]\n\n");
        #endif

        fprintf(stderr, "   Support for MistNet:");
        #ifdef MISTNET
        fprintf(stderr, " [enabled]\n\n");
        #endif
        #ifndef MISTNET
        fprintf(stderr, " [disabled]\n\n");
        #endif
 
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
        fprintf(stderr,"   Report bugs at: http://github.com/adokter/vol2bird/issues \n");
        fprintf(stderr,"   vol2bird home page: <http://github.com/adokter/vol2bird>\n");
    }
}


int main(int argc, char** argv) {
//    cfg_t* cfg;
    vol2bird_t alldata;
    
    // make sure executable and library version match
    if (strcmp(VERSION,libvol2bird_version()) != 0){
        fprintf(stderr,"Error: incorrect vol2bird installation: executable version (%s) has to match shared library version (%s).\n",VERSION,libvol2bird_version());
        return -1;
    }

    // print default message when no input arguments
    if (argc == 1) {
        usage(argv[0], 0);
        return -1;
    }
    
    // number of input files specified on command line
    int nInputFiles = 0;
    // the polar volume file that the user provided as input
    char* fileIn[INPUTFILESMAX];
    // the (optional) vertical profile file that the user specified as output
    const char* fileVpOut = NULL;
    // the (optional) vertical profile file that the user specified as output
    const char* fileVolOut = NULL;
    // the (optional) options.conf file path that the user specified as input
    const char* optionsFile = NULL;

    // determine whether we deal with legacy command line format (0) or getopt command line format (1)
    int commandLineFormat = 0;
    for (int i=0; i<argc; i++){
        if(strcmp("-i",argv[i])==0 || strcmp("--input",argv[i])==0  || \
           strcmp("-o",argv[i])==0 || strcmp("--output",argv[i])==0 || \
           strcmp("-p",argv[i])==0 || strcmp("--pvol",argv[i])==0   || \
           strcmp("-c",argv[i])==0 || strcmp("--config",argv[i])==0 || \
           strcmp("-h",argv[i])==0 || strcmp("--help",argv[i])==0   || \
           strcmp("-v",argv[i])==0 || strcmp("--version",argv[i])==0)
           {
            commandLineFormat = 1;
        }
    }

    // interpret legacy command line input
    if (commandLineFormat == 0){
        // check to see if we have the right number of input arguments
        if (argc > 4) {
            fprintf(stderr, "Error: Invalid command line arguments\n");
            usage(argv[0], 0);
            return -1;
        }
        
        // ------------------------------------------------------------- //
        //                initialization of variables                    //
        // ------------------------------------------------------------- //

        // the polar volume file that the user provided as input
        fileIn[0] = argv[1];
        nInputFiles = 1;

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
    }
    
    else{ // interpret command line input using getopt library
    
        int c;

        while (1) {
            static struct option long_options[] =
            {
                  /* These options don’t set a flag.
                     We distinguish them by their indices. */
                {"help",    no_argument,       0, 'h'},
                {"version", no_argument,       0, 'v'},
                {"input",   required_argument, 0, 'i'},
                {"output",  required_argument, 0, 'o'},
                {"pvol",    required_argument, 0, 'p'},
                {"config",  required_argument, 0, 'c'},
                {0, 0, 0, 0}
            };
            
            /* getopt_long stores the option index here. */
            int option_index = 0;
    
            c = getopt_long (argc, argv, "hvi:o:p:c:",
                           long_options, &option_index);

            /* Detect the end of the options. */
            if (c == -1) break;

            switch (c){
                case 0:
                /* If this option set a flag, do nothing else now. */
                    if (long_options[option_index].flag != 0)
                        break;
                    printf ("option %s", long_options[option_index].name);
                    if (optarg)
                        printf (" with arg %s", optarg);
                    printf ("\n");
                    break;

                case 'h':
                    usage(argv[0],1);
                    return -1;
                    break;

                case 'v':
                    fprintf(stdout,"%s version %s (%s)\n", argv[0], VERSION, VERSIONDATE);
                    return -1;
                    break;

                case 'i':
                    if (nInputFiles < INPUTFILESMAX){
                        fileIn[nInputFiles] = optarg;
                        nInputFiles++;
                    }
                    else{
                        fprintf(stderr, "Warning: too many input files, ignoring file %s ...\n", optarg);
                    }
                    break;

                case 'o':
                    fileVpOut = optarg;
                    break;

                case 'p':
                    fileVolOut = optarg;
                    break;

                case 'c':
                    optionsFile = optarg;
                    break;

                case '?':
                    /* getopt_long already printed an error message. */
                    break;

                default:
                    abort ();
            }
        }

        /* Print any remaining command line arguments (not options). */
        if (optind < argc) {
            printf ("unknown function argument(s): ");
            while (optind < argc)
                printf ("%s ", argv[optind++]);
            putchar ('\n');
        }        
    }

    // check that input files exist
    for (int i=0; i<nInputFiles; i++){
        if(!isRegularFile(fileIn[i])){
            fprintf(stderr, "Error: input file '%s' does not exist.\n", fileIn[i]);
            return -1;
        }
    }
    
    // check that options files exist
    if(optionsFile != NULL){
        if(!isRegularFile(optionsFile)){
            fprintf(stderr, "Error: configuration file '%s' does not exist.\n", optionsFile);
            return -1;
        }        
    }

    // Initialize hlhdf library
    HL_init();
    Rave_initializeDebugger();
    Rave_setDebugLevel(RAVE_WARNING);

    // Make rave and hlhdf library print debugging error messages
    //HL_setDebugLevel(HLHDF_SPEWDEBUG);
    //Rave_initializeDebugger();
    //Rave_setDebugLevel(RAVE_WARNING);
    //Rave_setDebugLevel(RAVE_INFO);
        
    // store the input filename TODO: add other input files
    strcpy(alldata.misc.filename_pvol, fileIn[0]);
    if (fileVpOut != NULL){
        strcpy(alldata.misc.filename_vp,fileVpOut);
    }
    else{
        strcpy(alldata.misc.filename_vp,"");
    }
    
    // read configuration options
    int configSuccessful = vol2birdLoadConfig(&alldata, optionsFile) == 0;

    if (configSuccessful == FALSE) {
        fprintf(stderr,"Error: failed to load configuration\n");
        return -1;
    }
    
    // read in data up to a distance of alldata.misc.rCellMax
    // we do not read in the full volume for speed/memory
    PolarVolume_t* volume = NULL;

    //FIXME maximum range specification not implemented for RSL / NEXRAD
    volume = vol2birdGetVolume(fileIn, nInputFiles, 1000000,1);
    //volume = vol2birdGetVolume(fileIn, nInputFiles, alldata.misc.rCellMax,1);

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
    //  using getter functions to access at the profile data               //
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
            fprintf(stdout, "# polar volume input: %s\n",fileIn[0]);
            if (alldata.misc.vcp > 0) fprintf(stdout, "# volume coverage pattern (VCP): %i\n", alldata.misc.vcp);
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
    } // getter scope end



    // ------------------------------------------------------------------- //
    //                 end of the getter example section                   //
    // ------------------------------------------------------------------- //            
        
    //map vol2bird profile data to Rave profile object
    mapDataToRave(volume, &alldata);

    int num_rows = VerticalProfile_getNumberOfRows(alldata.vp);
    int num_cols = VerticalProfile_getNumberOfColumns(alldata.vp);
    printf("Number of rows: %d\n", num_rows);
    printf("Number of columns: %d\n", num_cols);
    
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







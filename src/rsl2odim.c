/** convert NEXRAD radar volumes to ODIM hdf5
 * @file rsl2odim.c
 * @author Adriaan Dokter
 * @date 2016
 */

#include <stdio.h>
#include <float.h>
#include <getopt.h>
#include "rave_io.h"
#include "polarvolume.h"
#include "libvol2bird.h"
#include "constants.h"

void usage(char* programName, int verbose){
    fprintf(stderr,"rsl2odim version %s (%s)\n", VERSION, VERSIONDATE);
    fprintf(stderr,"usage: %s <RSL polar volume input> <ODIM hdf5 volume output>\n",programName);
    fprintf(stderr,"usage: %s -i <polar volume or scan> [-i <polar scan> [-i <polar scan>] ...] -o <ODIM hdf5 volume output>\n",programName);
    fprintf(stderr,"usage: %s --help\n", programName);

    if (verbose){
        fprintf(stderr,"rsl2odim version %s (%s)\n", VERSION, VERSIONDATE);
        fprintf(stderr,"\n   Supported radar data formats:\n");
        fprintf(stderr,"   * OPERA ODIM hdf5 input format, see <http://www.eumetnet.eu/opera-software> [enabled]\n");
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
    }
}

int main(int argc, char** argv) {
//    cfg_t* cfg;
    vol2bird_t alldata;

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
    const char* fileVolOut = NULL;

    // determine whether we deal with legacy command line format (0) or getopt command line format (1)
    int commandLineFormat = 0;
    for (int i=0; i<argc; i++){
        if(strcmp("-i",argv[i])==0 || strcmp("--input",argv[i])==0  || \
           strcmp("-o",argv[i])==0 || strcmp("--output",argv[i])==0 || \
           strcmp("-h",argv[i])==0 || strcmp("--help",argv[i])==0   || \
           strcmp("-v",argv[i])==0 || strcmp("--version",argv[i])==0)
           {
            commandLineFormat = 1;
        }
    }

    // interpret legacy command line input
    if (commandLineFormat == 0){
        // check to see if we have the right number of input arguments
        if (argc != 3) {
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
   
        // the volume file that the user specified as output
        fileVolOut = argv[2];

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
                {0, 0, 0, 0}
            };
            
            /* getopt_long stores the option index here. */
            int option_index = 0;
    
            c = getopt_long (argc, argv, "hvi:o:",
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
                    fileVolOut = optarg;
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

    // check that we have an output file specified
    if(fileVolOut == NULL){
        fprintf(stderr, "Error: no output file specified\n");
        return -1;
    }

    // check that we have at least one input file specified
    if(nInputFiles < 1){
        fprintf(stderr, "Error: no input file(s) specified\n");
        return -1;
    }

    // check that input files exist
    for (int i=0; i<nInputFiles; i++){
        if(!isRegularFile(fileIn[i])){
            fprintf(stderr, "Error: input file '%s' does not exist.\n", fileIn[i]);
            return -1;
        }
    }

    // read configuration options
    int configSuccessful = vol2birdLoadConfig(&alldata) == 0;

    if (configSuccessful == FALSE) {
        fprintf(stderr,"Error: failed to load configuration\n");
        return -1;
    }

    // read in data for the full range of distances.
    PolarVolume_t* volume = NULL;
    volume = vol2birdGetVolume(fileIn, nInputFiles, 1000000, 0);
    fprintf(stderr, "hoi1\n");

    if (volume == NULL) {
        fprintf(stderr,"Error: failed to read radar volume\n");
        return -1;
    }
    
	if(alldata.options.useMistNet){
         // initialize volbird library to run MistNet
        int initSuccessful = vol2birdSetUp(volume, &alldata) == 0;

    fprintf(stderr, "hoi2\n");
        if (initSuccessful == FALSE) {
            fprintf(stderr,"Error: failed to initialize vol2bird\n");
            return -1;
        }
	}

    saveToODIM((RaveCoreObject*) volume, fileVolOut);
    
    fprintf(stderr, "hoi3\n");
    // tear down vol2bird, give memory back
    if(alldata.options.useMistNet) vol2birdTearDown(&alldata);
    fprintf(stderr, "hoi4\n");
    RAVE_OBJECT_RELEASE(volume);
     
    fprintf(stderr, "hoi5\n");
    return 0;

}

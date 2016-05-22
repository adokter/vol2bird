//-------------------------------------------------------//
//              hard-coded constants                     //
//-------------------------------------------------------//
#define AREACELL                  4
#define CHISQMIN                  1e-5
#define CLUTPERCCELL              0.5
#define DBZCLUTTER                -10.0
#define DBZMIN                    0.0
#define FRINGEDIST                5000.0
#define NBINSGAP                  8
#define NOBSGAPMIN                5
#define NDBZMIN                   25
#define NEIGHBORS                 5
#define NTEXBINAZIM               3
#define NTEXBINRANG               3
#define NTEXMIN                   4
#define REFRACTIVE_INDEX_OF_WATER 0.964
#define STDEVCELL                 5.0
#define VDIFMAX                   10.0
#define VRADMIN                   1.0
#define NODETECT                  NAN
#define NODATA                    NAN
#define PROGRAM                   "vol2bird"
#define VERSION                   "0.2.2"
#define VERSIONDATE               "17-May-2016"

//-------------------------------------------------------//
//              user options defaults                    //
//-------------------------------------------------------//
// the number of layers in an altitude profile
#define NLAYER                       30
// the width/thickness of a layer [m]
#define HLAYER                       200.0f
//the minimum range [m] used for constructing the bird density profile
#define RANGEMIN                     5000.0f
//the maximum range [m] used for constructing the bird density profile
#define RANGEMAX                     25000.0f
//the minimum azimuth [degrees] used for constructing the bird density profile
#define AZIMMIN                      0.0f
//the maximum range [degrees] used for constructing the bird density profile
#define AZIMMAX                      360.0f
//the minimum scan elevation [degrees] used for constructing the bird density profile
#define ELEVMIN                      0.0f
//the maximum scan elevation used for constructing the bird density profile
#define ELEVMAX                      90.0f
//the wavelength [cm] of the radar
#define RADAR_WAVELENGTH_CM          5.3f
//whether a static clutter map is used
#define USE_STATIC_CLUTTER_DATA      0
//print options to stderr
#define PRINT_OPTIONS                0
//FIXME: add description
#define VERBOSE_OUTPUT_REQUIRED       0
//print dbz to stderr
#define PRINT_DBZ                    0
//print vrad to stderr
#define PRINT_VRAD                   0 
//print cell to stderr
#define PRINT_CELL                   0 
//print cell properties to stderr
#define PRINT_CELL_PROP              0 
//print texture to stderr
#define PRINT_TEXTURE                0 
//print clutter to stderr
#define PRINT_CLUT                   0 
//print profile data to stderr
#define PRINT_PROFILE                0 
//whether or not to print the 'points' array
#define PRINT_POINTS_ARRAY           0 
//Whether or not to fit a model to the observed vrad
#define FIT_VRAD                     1
//Whether to export bird profile as JSON
#define EXPORT_BIRD_PROFILE_AS_JSON  0
//Minimum Nyquist velocity
#define MIN_NYQUIST_VELOCITY         20.0f
//VVP Radial velocity standard deviation threshold
#define STDEV_BIRD                   2.0f
//Bird radar cross section [cm^2]
#define SIGMA_BIRD                   11.0f
//Maximum mean reflectivity factor for cells containing birds
#define DBZCELL                      15.0f
//Maximum reflectivity factor for gates containing birds
#define DBZMAX                       20.0f
//reflectivity quantity to use, one of DBZH, DBZV, TH, TV
#define DBZTYPE                      "DBZH" 
//for a range gate to contribute it should have a valid radial velocity
#define REQUIRE_VRAD                 0
//FIXME add description
#define EXPORT_BIRD_PROFILE_AS_JSON  0
//-------------------------------------------------------------//
//  hard-coded constants assigned to vol2birdConstants struct  //
//-------------------------------------------------------------//

#include <float.h>

// when analyzing cells, AREACELL determines the minimum size of a
// cell to be considered in the rest of the analysis [km^2]
#define AREACELL 0.5
// initialization value of rain segmentation field (CELL)
#define CELLINIT -1
// minimum standard deviation of the fit
#define CHISQMIN 1e-5
// cells with clutter fractions above this value are likely not birds
#define CLUTPERCCELL 0.5
// threshold value (on the external static clutter map!) above which gates are excluded as clutter
#define CLUTTERVALUEMIN 0.1
// each weather cell identified by findWeatherCells() is grown by a distance
// equal to 'fringeDist' using a region-growing approach
#define FRINGEDIST 5000.0
// when determining whether there are enough vrad observations in
// each direction, use NBINSGAP sectors
#define NBINSGAP 8
// there should be at least NOBSGAPMIN vrad observations in each sector
#define NOBSGAPMIN 5
// when calculating the altitude-layer averaged dbz, there should
// be at least NDBZMIN valid data points
#define NDBZMIN 25
// the minimum number of direct neighbors with dbz value above
// dbzThresMin as used in findWeatherCells()
#define NEIGHBORS 5
// vrad's texture is calculated based on the local neighborhood. The
// neighborhood size in the azimuth direction is equal to NTEXBINAZIM
// static int nAzimNeighborhood;
#define NTEXBINAZIM 3
// vrad's texture is calculated based on the local neighborhood. The
// neighborhood size in the range direction is equal to NTEXBINRANG
// static int nRangNeighborhood;
#define NTEXBINRANG 3
// the minimum number of neighbors for the texture value to be
// considered valid, as used in calcTexture()
#define NTEXMIN 4
// the refractive index of water
#define REFRACTIVE_INDEX_OF_WATER 0.964
// standard refraction coefficient (4/3)
#define REFRACTION_COEFFICIENT 1.333333f
// earth's radius (taken from NARR GRIB file)
#define EARTH_RADIUS 6371200 
// range gates up to a distance of RANGE_MAX+RCELLMAX_OFFSET are read into memory
// the extra offset allows for the raincell search to extend somewhat further
// than the maximum range used in the profile generation (RANGE_MAX).
#define RCELLMAX_OFFSET 5000.0f
// smallest range bin size to accept in metres
#define RSCALEMIN 10
// after fitting the vrad data, throw out any vrad observations that are more that VDIFMAX away
// from the fitted value, since these are likely outliers
#define VDIFMAX 10.0
// When analyzing cells, radial velocities lower than VRADMIN are treated as clutter
#define VRADMIN 1.0

//-------------------------------------------------------//
//       hard-coded options for use RSL library          //
//-------------------------------------------------------//

#ifdef RSL
// By how much the elevation angle of sweeps of different quantities
// is allowed to differ for them to be included into the same scan object.
// Applies to data read with the RSL library only
#define ELEVTOL 0.3
// Offsets and gains for encoding reflectivity, velocity
// and correlation coefficient
// Once encoded values should be positive for storage in RAVE objects
// Used with data read by the RSL library only
#define RSL_OFFSET_DBZ 0
#define RSL_GAIN_DBZ 1
#define RSL_OFFSET_VRAD 0
#define RSL_GAIN_VRAD 1
#define RSL_OFFSET_RHOHV 0
#define RSL_GAIN_RHOHV 1
// Encoded values reserved for nodata and undetects.
// Should be positive integers larger than zero
#define RSL_NODATA -1000
#define RSL_UNDETECT -999
#endif

//-------------------------------------------------------//
//            MistNet hard-coded options                 //
//-------------------------------------------------------//
// resolution of the Cartesian grid in meter for Mistnet
#define MISTNET_RESOLUTION 500
// X and Y dimension of the Cartesian grid for Mistnet,
// including a 4 pixel padding around the image, i.e. 8
// additional pixels.
#define MISTNET_DIMENSION 608
// number of pixels in MISTNET_DIMENSION used as a bleed 
#define MISTNET_BLEED 8
// number of MistNet elevation scans expected
#define MISTNET_N_ELEV 5
// predict a pixel as rain if the class probability for rain exceeds this threshold
#define MISTNET_WEATHER_THRESHOLD 0.45
// predict a pixel as rain if the average class probability for rain across
// all five elevations at that spatial location exceeds this threshold
#define MISTNET_SCAN_AVERAGE_WEATHER_THRESHOLD 0.45
// MistNet tensor indices for the background, biology and weather class probabilities
#define MISTNET_BACKGROUND_INDEX 0
#define MISTNET_BIOLOGY_INDEX 1
#define MISTNET_WEATHER_INDEX 2
// value reserved in the weather cell map for pixels identified by MistNet as weather
#define MISTNET_WEATHER_CELL_VALUE 2

//-------------------------------------------------------//
//             other hard-coded options                  //
//-------------------------------------------------------//

// maximum number of input files
#define INPUTFILESMAX 50
// Raw value used for gates or layers void of data (never ra-diated)
#define UNDETECT -999
// Raw value used for gates or layers when below the measurement detection threshold
// or when information could not be retrieved (radiated but nothing detected or calculated)
#define NODATA -1000
// name under which the calculated texture quantity will be stored
#define TEXNAME "VTEX"
// name under which the calculated raincell masking quantity will be stored
#define CELLNAME "CELL"
// name of the parameter containing the static cluttermap
#define CLUTNAME "OCCULT"
// Name of the program, to be stored as task attribute in ODIM
#define PROGRAM "vol2bird"
// Version of the program, to be stored as task_version attribute in ODIM
#define VERSION "0.5.0.9199"
// Date of latest version of the program
#define VERSIONDATE "31-Mar-2023"


//-------------------------------------------------------//
//  user options defaults (to be set in options.conf)    //
//-------------------------------------------------------//

// the number of layers in an altitude profile
#define NLAYER 25
// the width/thickness of a layer [m]
#define HLAYER 200.0f
// the minimum range [m] used for constructing the bird density profile
#define RANGEMIN 5000.0f
// the maximum range [m] used for constructing the bird density profile
#define RANGEMAX 35000.0f
// the minimum azimuth [degrees] used for constructing the bird density profile
#define AZIMMIN 0.0f
// the maximum range [degrees] used for constructing the bird density profile
#define AZIMMAX 360.0f
// the minimum scan elevation [degrees] used for constructing the bird density profile
#define ELEVMIN 0.0f
// the maximum scan elevation used for constructing the bird density profile
#define ELEVMAX 90.0f
// the wavelength [cm] of the radar
#define RADAR_WAVELENGTH_CM 5.3f
// whether a static clutter map is used
#define USE_CLUTTERMAP 0
// clutter map path and filename
#define CLUTTERMAP ""
// print options to stderr
#define PRINT_OPTIONS 0
// name of (optional) environmental variable containing path to user configuration file
#define OPTIONS_CONF "OPTIONS_CONF"
// default user configuration file name to search for in working directory
#define OPTIONS_FILE "options.conf"
// FIXME: add description
#define VERBOSE_OUTPUT_REQUIRED 0
// print dbz to stderr
#define PRINT_DBZ 0
// print aliased and dealiased vrad pairs to stderr
#define PRINT_DEALIAS 0
// print vrad to stderr
#define PRINT_VRAD 0
// print rhohv to stderr
#define PRINT_RHOHV 0
// print cell to stderr
#define PRINT_CELL 0
// print cell properties to stderr
#define PRINT_CELL_PROP 0
// print texture to stderr
#define PRINT_TEXTURE 0
// print clutter to stderr
#define PRINT_CLUT 0
// print profile data to stderr
#define PRINT_PROFILE 0
// whether or not to print the 'points' array
#define PRINT_POINTS_ARRAY 0
// Whether or not to fit a model to the observed vrad
#define FIT_VRAD 1
// Whether to export bird profile as JSON
#define EXPORT_BIRD_PROFILE_AS_JSON 0
// Scans with Nyquist velocity lower than this value are excluded
#define MIN_NYQUIST_VELOCITY 5.0f
// When all scans have nyquist velocity higher than this value, dealiasing is suppressed
#define MAX_NYQUIST_DEALIAS 25.0f
// when analyzing cells, only cells for which the stddev of vrad
// (aka the texture) is less than cellStdDevMax are considered in the
// rest of the analysis
#define STDEV_CELL 5.0f
// default VVP Radial velocity standard deviation threshold C-band (< 7.5 cm)
#define STDEV_BIRD 2.0f
// default VVP Radial velocity standard deviation threshold S-band (>= 7.5 cm)
#define STDEV_BIRD_S 1.0f
// Bird radar cross section [cm^2]
#define SIGMA_BIRD 11.0f
// Maximum mean reflectivity [cm^2/km^3] for cells containing birds
#define ETACELL 11500.0f
// Maximum reflectivity [cm^2/km^3] for single gates containing birds
#define ETAMAX 36000.0f
// minimum dbz of a gate to be considered for inclusion in a weather cell
#define DBZMIN 0.0
// reflectivity quantity to use, one of "DBZH", "DBZV", "TH", "TV"
#define DBZTYPE "DBZH"
// for a range gate to contribute it should have a valid radial velocity
#define REQUIRE_VRAD 0
// whether we should dealias the radial velocities
#define DEALIAS_VRAD 1
// whether we should dealias all data once (default), or dealias for each profile individually
#define DEALIAS_RECYCLE 1
// Test dealiasing field velocities up to VMAX m/s 
#define DEALIAS_VMAX 50.0
// Test field velocities increase in steps VMAX/VAF
#define DEALIAS_VAF 15.0
// Test field directions increase by 360/NF degrees
#define DEALIAS_NF 12.0
// whether you want to export the vertical bird profile as JSON
#define EXPORT_BIRD_PROFILE_AS_JSON 0
// whether to use dual-pol moments for filtering meteorological echoes
#define DUALPOL 1
// whether to use single-pol moments for filtering meteorological echoes
#define SINGLEPOL 1
// correlation coefficients higher than this threshold will be classified as precipitation
#define RHOHVMIN 0.95f
// whether to resample the input polar volume
#define RESAMPLE 0
// resampled range gate length in m
#define RESAMPLE_RSCALE 500.0f
// resampled number of range bins
#define RESAMPLE_NBINS 100
// resampled number of azimuth bins
#define RESAMPLE_NRAYS 360
// whether to use mistnet segmentation model
#define USE_MISTNET 0
// elevations to use in Cartesian projection for Mistnet
#define MISTNET_ELEVS "{0.5, 1.5, 2.5, 3.5, 4.5}"
// use only the specified elevation scans for mistnet to calculate profile if TRUE
// otherwise, use all available elevation scans
#define MISTNET_ELEVS_ONLY 1
// location of mistnet model in pytorch format
#define MISTNET_PATH "/MistNet/mistnet_nexrad.pt"
// initializing value of mistnet tensor
#define MISTNET_INIT 0
// require that radial velocity and spectrum width pixels rendered as mistnet input
// have a valid corresponding reflectivity value
#define MISTNET_REQUIRE_DBZ 0

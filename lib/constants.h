//-------------------------------------------------------------//
//  hard-coded constants assigned to vol2birdConstants struct  //
//-------------------------------------------------------------//

#include <float.h>

// when analyzing cells, AREAMIN determines the minimum size of a
// cell to be considered in the rest of the analysis
#define AREACELL 4
// minimum standard deviation of the fit
#define CHISQMIN 1e-5
// cells with clutter fractions above this value are likely not birds
#define CLUTPERCCELL 0.5
// threshold dbz value for excluding gates as clutter (static clutter only)
#define DBZCLUTTER -10.0
// minimum dbz for inclusion in a weather cell
#define DBZMIN 0.0
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
// minimum standard deviation of the VVP fit
#define STDEVCELL 5.0
// after fitting the vrad data, throw out any vrad observations that are more that VDIFMAX away
// from the fitted value, since these are likely outliers
#define VDIFMAX 10.0
// When analyzing cells, radial velocities lower than VRADMIN are treated as clutter
#define VRADMIN 1.0

//-------------------------------------------------------//
//             other hard-coded options                  //
//-------------------------------------------------------//

// Raw value used for gates or layers void of data (never ra-diated)
#define UNDETECT FLT_MAX
// Raw value used for gates or layers when below the measurement detection threshold
// or when information could not be retrieved (radiated but nothing detected or calculated)
#define NODATA -FLT_MAX
// Name of the program, to be stored as task attribute in ODIM
#define PROGRAM "vol2bird"
// Version of the program, to be stored as task_version attribute in ODIM
#define VERSION "0.2.2"
// Date of latest version of the program
#define VERSIONDATE "17-May-2016"

//-------------------------------------------------------//
//  user options defaults (to be set in options.conf)    //
//-------------------------------------------------------//

// the number of layers in an altitude profile
#define NLAYER 30
// the width/thickness of a layer [m]
#define HLAYER 200.0f
// the minimum range [m] used for constructing the bird density profile
#define RANGEMIN 5000.0f
// the maximum range [m] used for constructing the bird density profile
#define RANGEMAX 25000.0f
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
#define USE_STATIC_CLUTTER_DATA 0
// print options to stderr
#define PRINT_OPTIONS 0
// FIXME: add description
#define VERBOSE_OUTPUT_REQUIRED 0
// print dbz to stderr
#define PRINT_DBZ 0
// print vrad to stderr
#define PRINT_VRAD 0
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
// Minimum Nyquist velocity
#define MIN_NYQUIST_VELOCITY 20.0f
// VVP Radial velocity standard deviation threshold
#define STDEV_BIRD 2.0f
// Bird radar cross section [cm^2]
#define SIGMA_BIRD 11.0f
// Maximum mean reflectivity factor for cells containing birds
#define DBZCELL 15.0f
// Maximum reflectivity factor for gates containing birds
#define DBZMAX 20.0f
// reflectivity quantity to use, one of DBZH, DBZV, TH, TV
#define DBZTYPE "DBZH"
// for a range gate to contribute it should have a valid radial velocity
#define REQUIRE_VRAD 0
// whether you want to export the vertical bird profile as JSON
#define EXPORT_BIRD_PROFILE_AS_JSON 0

/*
 * Copyright 2015 Adriaan Dokter & Netherlands eScience Centre
 * If you want to use this software, please contact me at a.m.dokter@uva.nl
 *
 * This program calculates Vertical Profiles of Birds (VPBs) as described in
 *
 * Bird migration flight altitudes studied by a network of operational weather radars
 * Dokter A.M., Liechti F., Stark H., Delobbe L., Tabary P., Holleman I.
 * J. R. Soc. Interface, 8, 30–43, 2011
 * DOI: 10.1098/rsif.2010.0116
 *
 */

#include <confuse.h>
#include <polarvolume.h>
#include <vertical_profile.h>

// ****************************************************************************
// Definition of standard parameters.
// ****************************************************************************

#define DEG2RAD (0.017453293) // Degrees to radians.
#define RAD2DEG (57.29578)    // Radians to degrees.

// ****************************************************************************
// Definition of general macros:
// ****************************************************************************

#define XABS(x) (((x) < 0) ? (-(x)) : (x))
#define ROUND(x) (((x) > 0) ? (int)((x) + 0.5) : (int)((x)-0.5))
#define SQUARE(x) ((x) * (x))

// ****************************************************************************
// Defined (default) parameters.
// ****************************************************************************

#ifndef PI
#define PI (3.14159265358979323846)
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

// ****************************************************************************
//  Structure for containing SCAN metadata:
// ****************************************************************************

struct cellprop {
    int iRangOfMax;
    int iAzimOfMax;
    float dbzAvg;
    float texAvg;
    float cv;
    int nGates;
    int nGatesClutter;
    float dbzMax;
    int index;
    int drop;
};

struct scanmeta {
    float heig;            // Height of radar antenna in km.
    float elev;            // Elevation of scan in deg.
    int nRang;             // Number of range bins in scan.
    int nAzim;             // Number of azimuth rays in scan.
    float rangeScale;      // Size of range bins in scan in km.
    float azimScale;       // Size of azimuth steps in scan in deg.
    float valueOffset;     // Offset value of quantity contained by scan.
    float valueScale;      // Scale of value of quantity contained by scan.
    unsigned char missing; // Missing value of quantity contained by scan.
};

typedef struct cellprop CELLPROP;
typedef struct scanmeta SCANMETA;

// *****************************************************************************
// Structures for internal use
// *****************************************************************************

// ------------------------------------------------------------- //
//              vol2bird options from options.conf               //
// ------------------------------------------------------------- //

struct vol2birdOptions {
    int nLayers;              /* the number of layers in an altitude profile */
    float layerThickness;     /* the width/thickness of a layer [m] */
    float rangeMin;           /* the minimum range [m] used for constructing the bird density profile */
    float rangeMax;           /* the maximum range [m] used for constructing the bird density profile */
    float azimMin;            /* the minimum azimuth [degrees] used for constructing the bird density profile */
    float azimMax;            /* the maximum azimuth [degrees] used for constructing the bird density profile */
    float elevMin;            /* the minimum scan elevation [degrees] used for constructing the bird density profile */
    float elevMax;            /* the maximum scan elevation [degrees] used for constructing the bird density profile */
    float radarWavelength;    /* the default wavelength [cm] of the radar if it is not included in the metadata */
    int useStaticClutterData; /* whether a static clutter map is used */
    int printOptions;         /* print options to stderr */
    int printDbz;             /* print dbz to stderr */
    int printVrad;            /* print vrad to stderr */
    int printCell;            /* print cell to stderr */
    int printCellProp;        /* print cell properties to stderr */
    int printTex;             /* print texture to stderr */
    int printClut;            /* print clutter to stderr */
    int printProfileVar;      /* print profile data to stderr */
    int printPointsArray;     /* whether or not to print the 'points' array */
    int fitVrad;              /* Whether or not to fit a model to the observed vrad */
    int exportBirdProfileAsJSONVar; /* whether you want to export the vertical bird profile as JSON */
    float minNyquist;               /* Minimum Nyquist velocity [m/s] to include a scan; */
    float birdRadarCrossSection;    /* Bird radar cross section [cm^2] */
    float dbzMax;                   /* Maximum reflectivity factor of reflectivity gates containing birds */
    float cellDbzMin;               /* Maximum mean reflectivity factor of cells of birds */
    float stdDevMinBird;            /* Minimum VVP radial velocity standard deviation for layer containing birds*/
    char dbzType[10];               /* Preferred dBZ quantity to use */
    int requireVrad;                /* require range gates to have a valid radial velocity measurement */
    int dealiasVrad;                /* dealias radial velocities using torus mapping method by Haase et al. */
};
typedef struct vol2birdOptions vol2birdOptions_t;

// ------------------------------------------------------------- //
//              vol2bird options from constants.h                //
// ------------------------------------------------------------- //

struct vol2birdConstants {
    // after fitting the vrad data, throw out any vrad observations that are more that VDIFMAX away
    // from the fitted value, since these are likely outliers
    float absVDifMax;
    // when analyzing cells, AREAMIN determines the minimum size of a
    // cell to be considered in the rest of the analysis
    int nGatesCellMin;
    // cells with clutter fractions above this value are likely not birds
    float cellClutterFractionMax;
    // when analyzing cells, only cells for which the average dbz is
    // more than DBZCELL are considered in the rest of the analysis
    float cellDbzMin;
    // when analyzing cells, only cells for which the stddev of vrad
    // (aka the texture) is less than cellStdDevMax are considered in the
    // rest of the analysis
    float cellStdDevMax;
    // minimum standard deviation of the VVP fit
    float chisqMin;
    // threshold dbz value for excluding gates as clutter (static clutter only)
    float clutterValueMin;
    // maximum dbz used in calculation of profile dbzAvg
    float dbzMax;
    // minimum dbz for inclusion in a cell
    float dbzThresMin;
    // each weather cell identified by findWeatherCells() is grown by a distance
    // equal to 'fringeDist' using a region-growing approach
    float fringeDist;
    // the refractive index of water
    float refracIndex;
    // When analyzing cells, radial velocities lower than VRADMIN are treated as clutter
    float vradMin;
    // when determining whether there are enough vrad observations in
    // each direction, use NBINSGAP sectors
    int nBinsGap;
    // when calculating the altitude-layer averaged dbz, there should
    // be at least NDBZMIN valid data points
    int nPointsIncludedMin;
    // the minimum number of direct neighbors with dbz value above
    // dbzThresMin as used in findWeatherCells()
    int nNeighborsMin;
    // there should be at least NOBSGAPMIN vrad observations in each
    // sector
    int nObsGapMin;
    // vrad's texture is calculated based on the local neighborhood. The
    // neighborhood size in the azimuth direction is equal to NTEXBINAZIM
    // static int nAzimNeighborhood;
    int nAzimNeighborhood;
    // vrad's texture is calculated based on the local neighborhood. The
    // neighborhood size in the range direction is equal to NTEXBINRANG
    // static int nRangNeighborhood;
    int nRangNeighborhood;
    // the minimum number of neighbors for the texture value to be
    // considered valid, as used in calcTexture()
    int nCountMin;
};

typedef struct vol2birdConstants vol2birdConstants_t;

// ------------------------------------------------------------- //
//               information about the 'points' array            //
// ------------------------------------------------------------- //

// The data needed for calculating bird densities are collected
// in one big array, 'points'. Although this array is one
// variable, it is partitioned into 'nLayers' parts. The parts
// are not equal in size, therefore we need to keep track of
// where the data pertaining to a certain altitude bin can be
// written. The valid range of indexes into 'points' are stored
// in arrays 'indexFrom' and 'indexTo'.

struct vol2birdPoints {
    // the 'points' array has this many pseudo-columns
    int nColsPoints;
    // the 'points' array has this many rows
    int nRowsPoints;
    // the psuedo-column in 'points' that holds the azimuth angle
    int azimAngleCol;
    // the psuedo-column in 'points' that holds the elevation angle
    int elevAngleCol;
    // the psuedo-column in 'points' that holds the dbz value
    int dbzValueCol;
    // the psuedo-column in 'points' that holds the vrad value
    int vradValueCol;
    // the psuedo-column in 'points' that holds the cell value
    int cellValueCol;
    // the psuedo-column in 'points' that holds the gate classification code
    int gateCodeCol;
    // the 'points' array itself
    float* points; // Is allocated in vol2birdSetUp() and freed in vol2birdTearDown()
    // for a given altitude layer in the profile, only part of the 'points'
    // array is relevant. The 'indexFrom' and 'indexTo' arrays keep track
    // which rows in 'points' pertains to a given layer
    int* indexFrom; // Is allocated in vol2birdSetUp() and freed in vol2birdTearDown()
    int* indexTo;   // Is allocated in vol2birdSetUp() and freed in vol2birdTearDown()
    // nPointsWritten stores the number of points that was copied from one
    // of the scan elevations to the 'points' array; it should therefore
    // never exceed indexTo[i]-indexFrom[i]
    int* nPointsWritten; // Is allocated in vol2birdSetUp() and freed in vol2birdTearDown()
};
typedef struct vol2birdPoints vol2birdPoints_t;

// ------------------------------------------------------------- //
//          information about the flagfields of 'gateCode'       //
// ------------------------------------------------------------- //

struct vol2birdFlags {
    // the bit in 'gateCode' that says whether this gate is true in the static
    // clutter map (which we don't have yet TODO)
    int flagPositionStaticClutter;
    // the bit in 'gateCode' that says whether this gate is part of the
    // calculated cluttermap (without fringe)
    int flagPositionDynamicClutter;
    // the bit in 'gateCode' that says whether this gate is part of the
    // fringe of the calculated cluttermap
    int flagPositionDynamicClutterFringe;
    // the bit in 'gateCode' that says whether this gate has reflectivity data
    // but no corresponding radial velocity data
    int flagPositionVradMissing;
    // the bit in 'gateCode' the psuedo-columnsays whether this gate's dbz value is too
    // high to be due to birds, it must be caused by something else
    int flagPositionDbzTooHighForBirds;
    // the bit in 'gateCode' the psuedo-columnsays whether this gate's radial velocity is
    // close to zero. These gates are all discarded to exclude ground
    // clutter, which often has a radial velocity near zero.
    int flagPositionVradTooLow;
    // the bit in 'gateCode' that says whether this gate passed the VDIFMAX test
    int flagPositionVDifMax;
    // the bit in 'gateCode' that says whether the gate's azimuth angle was too low
    int flagPositionAzimTooLow;
    // the bit in 'gateCode' that says whether the gate's azimuth angle was too high
    int flagPositionAzimTooHigh;
};
typedef struct vol2birdFlags vol2birdFlags_t;

// ------------------------------------------------------------- //
//              information about the 'profile' array            //
// ------------------------------------------------------------- //

struct vol2birdProfiles {
    // the number of different types of profile we're making
    int nProfileTypes;
    // how many rows there are in a profile
    int nRowsProfile;
    // columns in profile contain
    // [altmin,altmax,u,v,w,hSpeed,hDir,chi,hasGap,dbzAvg,nPointsCopied,reflectivity,birdDensity]
    int nColsProfile;
    // the profile array itself
    float* profile; // Is allocated in vol2birdSetUp() and freed in vol2birdTearDown()
    // these next 3 profile arrays are an ugly way to make sure
    // vol2birdGetProfile() can deliver its data
    float* profile1;
    float* profile2;
    float* profile3;
    // the type of profile that was last calculated
    int iProfileTypeLast;
};
typedef struct vol2birdProfiles vol2birdProfiles_t;

// ------------------------------------------------------------- //
//                       some other variables                    //
// ------------------------------------------------------------- //

struct vol2birdMisc {
    // rCellMax is defined as rangeMax + 5000.0f  in order to avoid
    // edge effects when calculating the fringe
    float rCellMax;
    // the number of dimensions to describe the location of an echo
    // (azimuth and elevAngle) as used in the 'pointsSelection' array
    // that is passed to svdfit
    int nDims;
    // how many parameters are fitted by the svdfit procedure
    int nParsFitted;
    // the factor that is used when converting from Z to eta, calculated from radar wavelength
    float dbzFactor;
    // whether the vol2bird module has been initialized
    int initializationSuccessful;
    // whether vol2bird calculated a valid bird profile
    int vol2birdSuccessful;
    // number of scans used to calculate the profile
    int nScansUsed;
    // whether configuration was loaded successfully
    int loadConfigSuccessful;
    // during calculation of iProfileType == 3, use this array to store the
    // result of (chi < stdDevMinBird), such that it can be used later during
    // calculation of iProfileType == 1
    int* scatterersAreNotBirds; // Is allocated in vol2birdSetUp() and freed in vol2birdTearDown()
    // this string contains all the user options and constants, for storage in ODIM task_args attribute
    char task_args[3000];
};
typedef struct vol2birdMisc vol2birdMisc_t;

// structure for storing scan properties
struct vol2birdScanUse {
    // whether to use this scan in calculation of the profiles
    int useScan;
    // the reflectivity quantity used for this scan
    char dbzName[10];
    // the radial velocity quantity used for this scan
    char vradName[10];
};
typedef struct vol2birdScanUse vol2birdScanUse_t;

// root structure, containing all data
struct vol2bird {
    vol2birdOptions_t options;
    vol2birdConstants_t constants;
    vol2birdPoints_t points;
    vol2birdFlags_t flags;
    vol2birdProfiles_t profiles;
    vol2birdMisc_t misc;
    VerticalProfile_t* vp;
    cfg_t* cfg;
};
typedef struct vol2bird vol2bird_t;

// *****************************************************************************
// Public function prototypes
// *****************************************************************************

void vol2birdCalcProfiles(vol2bird_t* alldata);

float* vol2birdGetProfile(int iProfileType, vol2bird_t* alldata);

PolarVolume_t* vol2birdGetVolume(char* filename);

int vol2birdGetNColsProfile(vol2bird_t* alldata);

int vol2birdGetNRowsProfile(vol2bird_t* alldata);

void vol2birdPrintIndexArrays(vol2bird_t* alldata);

void vol2birdPrintOptions(vol2bird_t* alldata);

void vol2birdPrintPointsArray(vol2bird_t* alldata);

int vol2birdLoadConfig(vol2bird_t* alldata);

int vol2birdSetUp(PolarVolume_t* volume, vol2bird_t* alldata);

void vol2birdTearDown(vol2bird_t* alldata);

int mapDataToRave(PolarVolume_t* volume, vol2bird_t* alldata);

float nanify(float value);

int saveToODIM(VerticalProfile_t* volume, const char* filename);

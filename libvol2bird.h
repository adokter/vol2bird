//
// Copyright 2013 Netherlands eScience Center
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

#include <confuse.h>

// ****************************************************************************
// Definition of standard parameters.
// ****************************************************************************

#define DEG2RAD    (0.017453293)  // Degrees to radians.
#define RAD2DEG    (57.29578)     // Radians to degrees.

// ****************************************************************************
// Definition of general macros:
// ****************************************************************************

#define XABS(x)    (((x)<0)?(-(x)):(x))
#define ROUND(x)   (((x)>0)?(int)((x)+0.5):(int)((x)-0.5))
#define SQUARE(x)  ((x)*(x))

// ****************************************************************************
// Defined (default) parameters.
// ****************************************************************************


#ifndef PI
#define PI          (3.14159265358979323846)
#endif

#ifndef TRUE
#define TRUE        1
#endif

#ifndef FALSE
#define FALSE       0
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
    float heig;              // Height of radar antenna in km.
    float elev;              // Elevation of scan in deg.
    int nRang;               // Number of range bins in scan.
    int nAzim;               // Number of azimuth rays in scan.
    float rangeScale;        // Size of range bins in scan in km.
    float azimScale;         // Size of azimuth steps in scan in deg.
    float valueOffset;       // Offset value of quantity contained by scan.
    float valueScale;        // Scale of value of quantity contained by scan.
    unsigned char missing;   // Missing value of quantity contained by scan.
};


typedef struct cellprop CELLPROP;
typedef struct scanmeta SCANMETA;

// *****************************************************************************
// Structures for internal use
// *****************************************************************************

struct vol2birdOptions
{
	int nLayers;			/* the number of layers in an altitude profile */
	float layerThickness;		/* the width/thickness of a layer [m] */
	float rangeMin;			/* the minimum range [m] used for constructing the bird density profile */
	float rangeMax;			/* the maximum range [m] used for constructing the bird density profile */
	float azimMin;			/* the minimum azimuth [degrees] used for constructing the bird density profile */
	float azimMax;			/* the maximum azimuth [degrees] used for constructing the bird density profile */
	float radarWavelength;		/* the default wavelength [cm] of the radar if it is not included in the metadata */
	int useStaticClutterData;	/* whether a static clutter map is used */
	int printOptions;		/* print options to stderr */
	int printDbz;			/* print dbz to stderr */
	int printVrad;			/* print vrad to stderr */
	int printCell;			/* print cell to stderr */
	int printCellProp;		/* print cell properties to stderr */
	int printTex;			/* print texture to stderr */
	int printClut;			/* print clutter to stderr */
	int printProfileVar;		/* print profile data to stderr */
	int printPointsArray;		/* whether or not to print the 'points' array */
	int fitVrad;			/* Whether or not to fit a model to the observed vrad */
	int exportBirdProfileAsJSONVar; /* */
	float minNyquist;		/* Minimum Nyquist velocity [m/s] to include a scan; to excluded velocity scans too heavily folded */
};
typedef struct vol2birdOptions vol2birdOptions_t;



struct vol2birdConstants
{
	int nGatesCellMin;
	float cellClutterFractionMax;
	float cellDbzMin;
	float chisqMin;
	float clutterValueMin;
	float dbzMax;
	float dbzThresMin;
	float fringeDist;
	int nBinsGap;
	int nPointsIncludedMin;
	int nNeighborsMin;
	int nObsGapMin;
	int nAzimNeighborhood;
	int nRangNeighborhood;
	int nCountMin; 
	float refracIndex;
	float birdRadarCrossSection;
	float cellStdDevMax;
	float stdDevMinBird;
	float absVDifMax;
	float vradMin;
};
typedef struct vol2birdConstants vol2birdConstants_t;



struct vol2birdPoints
{
	int nColsPoints;
	int nRowsPoints;
	int azimAngleCol;
	int elevAngleCol;
	int dbzValueCol;
	int vradValueCol;
	int cellValueCol;
	int gateCodeCol;
	float* points;		// Is allocated in vol2birdSetUp() and freed in vol2birdTearDown()
	int* indexFrom;		// Is allocated in vol2birdSetUp() and freed in vol2birdTearDown()
	int* indexTo;		// Is allocated in vol2birdSetUp() and freed in vol2birdTearDown()
	int* nPointsWritten;	// Is allocated in vol2birdSetUp() and freed in vol2birdTearDown()
};
typedef struct vol2birdPoints vol2birdPoints_t;



struct vol2birdFlags
{
	int flagPositionStaticClutter;
	int flagPositionDynamicClutter;
	int flagPositionDynamicClutterFringe;
	int flagPositionVradMissing;
	int flagPositionDbzTooHighForBirds;
	int flagPositionVradTooLow;
	int flagPositionVDifMax;
	int flagPositionAzimTooLow;
	int flagPositionAzimTooHigh;
};
typedef struct vol2birdFlags vol2birdFlags_t;



struct vol2birdProfiles
{
	int nProfileTypes;
	int nRowsProfile;
	int nColsProfile; 
	float* profile;		// Is allocated in vol2birdSetUp() and freed in vol2birdTearDown()
	float* profile1;
	float* profile2;
	float* profile3;
	int iProfileTypeLast;
};
typedef struct vol2birdProfiles vol2birdProfiles_t;



struct vol2birdMisc
{
	float rCellMax;
	int nDims;
	int nParsFitted;
	float dbzFactor;
	int initializationSuccessful;
	int* scatterersAreNotBirds;	// Is allocated in vol2birdSetUp() and freed in vol2birdTearDown()
};
typedef struct vol2birdMisc vol2birdMisc_t;

struct vol2birdScanUse
{
	int useScan;
	char dbzName[10];
	char vradName[10];
};
typedef struct vol2birdScanUse vol2birdScanUse_t;

struct vol2bird
{
	vol2birdOptions_t options;
	vol2birdConstants_t constants;
	vol2birdPoints_t points;
	vol2birdFlags_t flags;
	vol2birdProfiles_t profiles;
	vol2birdMisc_t misc;
};
typedef struct vol2bird vol2bird_t;





// *****************************************************************************
// Public function prototypes
// *****************************************************************************

void vol2birdCalcProfiles(vol2bird_t* alldata);

float* vol2birdGetProfile(int iProfileType, vol2bird_t *alldata);

int vol2birdGetNColsProfile(vol2bird_t *alldata);
    
int vol2birdGetNRowsProfile(vol2bird_t *alldata);

void vol2birdPrintIndexArrays(vol2bird_t* alldata);

void vol2birdPrintOptions(vol2bird_t* alldata);

void vol2birdPrintPointsArray(vol2bird_t* alldata);

int vol2birdSetUp(PolarVolume_t* volume, cfg_t* cfg, vol2bird_t* alldata);

void vol2birdTearDown(cfg_t* cfg, vol2bird_t* alldata);

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
// Public function prototypes
// *****************************************************************************

void vol2birdCalcProfiles();

void vol2birdPrintIndexArrays(void);

void vol2birdPrintOptions(void);

void vol2birdPrintPointsArray(void);

int vol2birdSetUp(PolarVolume_t* volume);

void vol2birdTearDown();

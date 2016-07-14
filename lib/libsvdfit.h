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


// ****************************************************************************
// Definition of general macros:
// ****************************************************************************

#define SIGN(x)    (((x)<0)?-1:1)
#define SQUARE(x)  ((x) * (x))
#define XYMAX(x,y) (((x)<(y))?(y):(x))
#define XYMIN(x,y) (((x)<(y))?(x):(y))

// no longer needed - already defined in dealias.h
#define DEG2RAD    (0.017453293)  // Degrees to radians.

// ****************************************************************************
// Definition of parameters for fitting
// ****************************************************************************

#define NPARSFITTEDMAX    (16)      /*Maximum number of fit parameters.*/
#define SVDTOL            (1e-5)      /*Accuracy in SV decomposition.*/



// *****************************************************************************
// Function prototypes
// *****************************************************************************

int svd_vvp1func(const float points[], const int nDims, float afunc[], const int nParsFitted);
int svbksb(float *u, float w[], float *v, int m, int n, const float b[], float x[]);
int svdcmp(float *a, int m, int n, float w[], float *v);
float svdfit(const float *points, const int nDims, const float yObs[], float yFitted[], const int nPoints,
             float parameterVector[], float avar[], const int nParsFitted);




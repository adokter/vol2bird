#include <math.h>
#ifndef NAN
#include <bits/nan.h>
#endif
#include <string.h>

#include "rave_alloc.h"

/******************************************************************************/
/*Definition of standard parameters.                                          */
/******************************************************************************/

#ifndef DEG2RAD
#define DEG2RAD (0.017453293) // Degrees to radians.
#define RAD2DEG (57.29578)    // Radians to degrees.
#endif
#define VMAX       48         /* Test field velocities up to VMAX m/s             */
#define VAF        6          /* Test field velocities increase in steps VMAX/VAF */
#define NF         8          /* Test field directions increase by 360/NF degrees */

void printDealias(const float *points, const int nDims, const float nyquist[], 
	const float vradObs[], float vradDealias[], const int nPoints, const int iProfileType, const int iLayer, const int iPass);

int dealias_points(const float *points, const int nDims, const float nyquist[], 
	const double NI_MIN, const float vo[], float vradDealias[], const int nPoints);

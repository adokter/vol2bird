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
#define VAF        4          /* Test field velocities increase in steps NI/VAF   */
#define NF         40         /* Test field directions increase by 360/NF degrees */

int dealias_points(const float *points, const int nDims, const float nyquist[], 
	const double NI_MIN, const float vo[], float vradDealias[], const int nPoints);

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
#define VMAX       48              /* Maximum velocity */
#define VAF        4               /*   */
#define NF         40              /*   */
#define MVA        8               /*   */

int dealias_points(const float *points, const int nDims, const float nyquist[], 
	const double NI, const float vo[], float vradDealias[], const int nPoints);

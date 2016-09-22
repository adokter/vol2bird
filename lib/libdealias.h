#include <math.h>
#ifndef NAN
#include <bits/nan.h>
#endif
#include <string.h>

#include "rave_alloc.h"

/******************************************************************************/
/*Definition of standard parameters.                                          */
/******************************************************************************/

#define DEG2RAD    DEG_TO_RAD      /* Degrees to radians. From PROJ.4 */
#define RAD2DEG    RAD_TO_DEG      /* Radians to degrees. From PROJ.4 */
#define VMAX       48              /* Maximum velocity */
#define VAF        4               /*   */
#define NF         40              /*   */
#define MVA        8               /*   */

int dealias_points(const float *points, const int nDims, const float nyquist[], 
	const double NI, const float vo[], float vradDealias[], const int nPoints);

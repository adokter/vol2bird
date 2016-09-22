/* --------------------------------------------------------------------
Copyright (C) 2012 Swedish Meteorological and Hydrological Institute, SMHI,

This file is part of RAVE.

RAVE is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

RAVE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with RAVE.  If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------*/
/**
 * Provides functionality for dealiasing radial wind data.
 * @file
 * @author Gunther Haase, SMHI
  * @date 2013-02-06
 */
#ifndef DEALIAS_H
#define DEALIAS_H

#include <math.h>
#ifndef NAN
#include <bits/nan.h>
#endif
#include <string.h>

#include "polarvolume.h"
#include "polarscan.h"
#include "polarscanparam.h"
#include "rave_attribute.h"
#include "rave_object.h"
#include "rave_alloc.h"
#include "rave_types.h"
#include "raveutil.h"
#include "rave_debug.h"

/******************************************************************************/
/*Definition of standard parameters.                                          */
/******************************************************************************/

#define DEG2RAD    DEG_TO_RAD      /* Degrees to radians. From PROJ.4 */
#define RAD2DEG    RAD_TO_DEG      /* Radians to degrees. From PROJ.4 */
#define VMAX       48              /* Maximum velocity */
#define VAF        4               /*   */
#define NF         40              /*   */
#define MVA        8               /*   */
#define EPSILON    1E-10           /*   */
#define EMAX       2.0             /* Maximum elevation angle [deg] */


/**
 * Returns the maximum of a vector.
 * @param[in] a - input vector
 * @param[in] n - vector length
 * @returns maximum of a vector
 */
double max_vector (double *a, int n);

/**
 * Returns the minimum of a vector.
 * @param[in] a - input vector
 * @param[in] n - vector length
 * @returns minimum of a vector
 */
double min_vector (double *a, int n);

/**
 * Determines whether a scan is dealiased or not
 * @param[in] scan - input scan
 * @param[in] quantity - the quantity to be dealiased
 * @returns int 1 if scan is dealiased, otherwise 0
 */
int dealiased_by_quantity(PolarScan_t* scan, const char* quantity);

/**
 * Determines whether a scan is dealiased or not and performs dealiasing of the VRAD parameter. Same as calling
 * dealiased_by_quantity(scan, "VRAD").
 *
 * @param[in] scan - input scan
 * @returns int 1 if scan is dealiased, otherwise 0
 */
int dealiased(PolarScan_t* scan);

/**
 * Function for dealiasing polar volume data for the specified parameter.
 * @param[in] source - input volume
 * @param[in] quantity - the quantity
 * @param[in] emax - the maximum elevation angle
 * @returns int 1 upon success, otherwise 0
 */
int dealias_pvol_by_quantity(PolarVolume_t* inobj, const char* quantity, double emax);

/**
 * Function for dealiasing polar volume data for the VRAD parameter. Same as invoking dealias_pvol_by_quantity(vol, "VRAD").
 * @param[in] source - input volume
 * @returns int 1 upon success, otherwise 0
 */
int dealias_pvol(PolarVolume_t* inobj);

/**
 * Function for dealiasing polar scan data for the specified parameter
 * @param[in] source - input scan
 * @param[in] quantity - the quantity
 * @param[in] emax - the maximum elevation angle
 * @returns int 1 upon success, otherwise 0
 */
int dealias_scan_by_quantity(PolarScan_t* inobj, const char* quantity, double emax);

/**
 * Function for dealiasing polar scan data for the VRAD parameter. Same as calling dealias_scan_by_quantity(scan, "VRAD").
 * @param[in] source - input scan
 * @returns int 1 upon success, otherwise 0
 */
int dealias_scan(PolarScan_t* inobj);

#endif /* DEALIAS_H */

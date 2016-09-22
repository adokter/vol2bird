/** Function for dealiasing weather radar radial velocities.
 * @file libdealias.c
 * @author Adriaan Dokter, adapted from dealias.c RAVE functions by Gunther Haase
 * @date 2016-09-22
 */

#include "libdealias.h"

int dealias_points(const float *points, const int nDims, const float nyquist[], 
	const double NI, const float vo[], float vradDealias[], const int nPoints){
  
	int i, j, n, m, eind;
    double vm, min1, esum, u1, v1, min2, dmy;
  
	// number of rows
	m = floor (VAF/NI*VMAX);
	// number of columns
	n = NF;

	// polarscan matrix, torus projected x coordinate, eq. 6 Haase et al. 2004 jaot
	double *x = RAVE_CALLOC ((size_t)nPoints, sizeof(double));
	// polarscan matrix, torus projected y coordinate, eq. 7 Haase et al. 2004 jaot
	double *y = RAVE_CALLOC ((size_t)nPoints, sizeof(double));
	// the observed, possibly aliased, radial velocity 
	//float *vo = RAVE_CALLOC ((size_t)nPoints, sizeof(float));
	// U-components of test velocity fields
	double *uh = RAVE_CALLOC ((size_t)(m*n), sizeof(double));
	// V-components of test velocity fields
	double *vh = RAVE_CALLOC ((size_t)(m*n), sizeof(double));
	// summed absolute differences between observed velocity field and test fields
	double *e = RAVE_CALLOC ((size_t)(m*n*nPoints), sizeof(double));
	// see Eq. 7 Haase et al. 2004 jaot
	double *xt = RAVE_CALLOC ((size_t)(m*n*nPoints), sizeof(double));
	// see Eq. 6 Haase et al. 2004 jaot
	double *yt = RAVE_CALLOC ((size_t)(m*n*nPoints), sizeof(double));
      // radial velocities of the best fitting test field
	double *vt1 = RAVE_CALLOC ((size_t)nPoints, sizeof(double));
	// ordered array of possible aliases

	// map measured data to 3F
	for (i=0; i<nPoints; i++) {
		*(x+i) = nyquist[i]/M_PI * cos(*(vo+i)*M_PI/nyquist[i]);
		*(y+i) = nyquist[i]/M_PI * sin(*(vo+i)*M_PI/nyquist[i]);
	}
	  
	// Setting up the u and v component of the test velocity fields:
	// index n=NF gives number of azimuthal directions (default n=40, i.e. steps of 360/40=9 degrees)
	// index m=VAF/NI*VMAX gives number of speeds (maximum speed is VMAX, steps of NI/VAF)
	for (i=0; i<n; i++) {
        for (j=0; j<m; j++) {
			*(uh+i*m+j) = NI/VAF*(j+1) * sin(2*M_PI/NF*i);
			*(vh+i*m+j) = NI/VAF*(j+1) * cos(2*M_PI/NF*i);
        }
	}

	for (int iPoint=0; iPoint<nPoints; iPoint++) {
        for (i=0; i<n; i++) {
			for (j=0; j<m; j++) {
				// calculate the radial velocities for the test wind fields, eq 4 in Haase et al. 2004 jaot
				vm = *(uh+i*m+j) * sin(points[nDims*iPoint]) +
					 *(vh+i*m+j) * cos(points[nDims*iPoint]);
				// Eq. 7 for the test radial wind:
				*(xt+i*m+j+iPoint*m*n) = nyquist[iPoint]/M_PI * cos(vm*M_PI/nyquist[iPoint]);
				// Eq. 6 for the test radial wind:
				*(yt+i*m+j+iPoint*m*n) = nyquist[iPoint]/M_PI * sin(vm*M_PI/nyquist[iPoint]);
			}
        }
	}  
	
	for (int iPoint=0; iPoint<nPoints; iPoint++) {
		for (i=0; i<m*n; i++) {
            // difference between the radial velocity of the test wind field, and the observed radial velocity, for a given range bin.
            *(e+i+iPoint*m*n) = fabs(*(xt+i+iPoint*m*n)-*(x+iPoint)) +
                            fabs(*(yt+i+iPoint*m*n)-*(y+iPoint));
		}
	}

	min1 = 1e32;
	eind = 0;
	u1 = 0;
	v1 = 0;
	// select the test wind with the best fit
	for (i=0; i<m*n; i++) {
		esum = 0;
		for (int iPoint=0; iPoint<nPoints; iPoint++) {
            if (!isnan(*(e+i+iPoint*m*n))) {
				esum = esum + *(e+i+iPoint*m*n);
            }
		}
		if (esum<min1) {
            min1 = esum;
            eind = i;
		}
		u1 = *(uh+eind);
		v1 = *(vh+eind);
	}

	// the radial velocity of the best fitting test velocity field:
	for (int iPoint=0; iPoint<nPoints; iPoint++) {
		*(vt1+iPoint) = u1*sin(points[nDims*iPoint]) + v1*cos(points[nDims*iPoint]);
	}

	// dealias the observed velocities using the best test velocity field
	for (int iPoint=0; iPoint<nPoints; iPoint++) {
		min2 = 1e32;
		dmy = 0;
		float diffVTest = (*(vt1+iPoint)-*(vo+iPoint));
		float dv = 0;
		for (i=0; i<MVA+1; i++) {
			// get a candidate velocity fold, i.e. integer multiple of nyquist velocity
			dv = nyquist[iPoint]*(2*i-MVA);
            // checking how many folds we have, vt1-vo is the residual between the real velocity
            // and the folded velocity; which equals a  multiple of the nyquist interval			
            dmy = fabs(dv-diffVTest);
            if ((dmy<min2) && (!isnan(dmy))) {
				// add the aliased interval to the observed velocity field, and obtain dealiased velocity
				*(vradDealias+iPoint) = *(vo+iPoint) + dv;
				min2 = dmy;
            }
		} // loop MVA
	} // loop over points
	
	return 1;
}

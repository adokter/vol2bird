/** Function for dealiasing weather radar radial velocities.
 * @file libdealias.c
 * @author Adriaan Dokter, adapted from dealias.c RAVE functions by Gunther Haase
 * @date 2016-09-22
 */

#include "libdealias.h"
#include <stdio.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

void vol2bird_err_printf(const char* fmt, ...);

void printDealias(const float *points, const int nDims, const float nyquist[], 
    const float vradObs[], float vradDealias[], const int nPoints, const int iProfileType, const int iLayer, const int iPass){
    vol2bird_err_printf("#iProfile iLayer iPass azim elev nyquist vrad vradd\n");
    for(int i=0; i<nPoints;i++){
      vol2bird_err_printf("%i %i %i %3.1f %3.1f %3.1f %3.1f %3.1f\n",
        iProfileType,iLayer,iPass,points[i*nDims],points[i*nDims+1],nyquist[i],vradObs[i],vradDealias[i]);
    }
}

double test_field(float u, float v,const float *points, const float *pointsTrigon, const int nPoints, const int nDims, double x[], double y[], const float nyquist[]){
    double xt, yt, e, vm;
    double esum = 0;
    for (int iPoint=0; iPoint<nPoints; iPoint++) {
        // calculate the radial velocities for the test wind fields, eq 4 in Haase et al. 2004 jaot
        vm = (u*pointsTrigon[3*iPoint] + v*pointsTrigon[3*iPoint+1])*pointsTrigon[3*iPoint+2];
        // Eq. 7 for the test radial wind:
        xt = nyquist[iPoint]/M_PI * cos(vm*M_PI/nyquist[iPoint]);
        // Eq. 6 for the test radial wind:
        yt = nyquist[iPoint]/M_PI * sin(vm*M_PI/nyquist[iPoint]);
        
        // summed absolute differences between observed velocity field and test fields
        e = fabs(xt-x[iPoint]) + fabs(yt-y[iPoint]);
//        fprintf(stderr,"test_field: i=%i,u=%f,v=%f,vm=%f,xt=%f,yt=%f,e=%f,sin=%f,cos=%f,cos=%f\n",iPoint,u,v,vm,xt,yt,e,pointsTrigon[3*iPoint],pointsTrigon[3*iPoint+1],pointsTrigon[3*iPoint+2]);
        if (!isnan(e)) {
            esum = esum + e;
        }
    }
    return esum;
}

double test_field_gsl(const gsl_vector *uv, void* params){
    double u,v;
    float *points = ((void **) params)[0];
    float *pointsTrigon = ((void **) params)[1];
    int nPoints = *(int *) ((void **) params)[2];
    int nDims = *(int *) ((void **) params)[3];
    double *x = ((void **) params)[4];
    double *y = ((void **) params)[5];
    float *nyquist = ((void **) params)[6];

    u = gsl_vector_get(uv, 0);
    v = gsl_vector_get(uv, 1);
 
    return test_field(u, v, points, pointsTrigon, nPoints, nDims, x, y, nyquist);
 
}

int fit_field_gsl(gsl_vector *uv, void *params){
    gsl_vector *ss;
    
    // define which minimizer to use
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    
    int iter = 0;
    int status;
    double size;
    double u1 = 0;
    double v1 = 0;
    
    // Set initial step sizes to 1
    ss = gsl_vector_alloc(2);
    gsl_vector_set_all(ss, 1);

    // Initialize method
    gsl_multimin_function minex_func;
    minex_func.n = 2;
    minex_func.f = &test_field_gsl;
    minex_func.params = params;
    s = gsl_multimin_fminimizer_alloc (T, 2);
    gsl_multimin_fminimizer_set (s, &minex_func, uv, ss);
    
    // minimize by iteration
    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        
        if (status) 
        break;

        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, 1e-2);
        
        if (status == GSL_SUCCESS)
        {
            #ifdef FPRINTFON 
            printf ("converged to minimum on iteration ");
            printf ("%d at %10.3e %10.3e f() = %7.3f size = %.3f\n", 
                iter,
                gsl_vector_get (s->x, 0), 
                gsl_vector_get (s->x, 1), 
                s->fval, size);
            #endif

            u1=gsl_vector_get(s->x, 0);
            v1=gsl_vector_get(s->x, 1);
            
            gsl_vector_set(uv, 0, u1);
            gsl_vector_set(uv, 1, v1);
        }
    }
    while (status == GSL_CONTINUE && iter < 100);

    #ifdef FPRINTFON
    fprintf(stdout,"Finished dealias at (x,y)=%f,%f at f()=%f ...\n",u1,v1,s->fval);
    #endif

    // clean up
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);
    
    if (status != GSL_SUCCESS) return 0;
    return 1;
}


int dealias_points(const float *points, const int nDims, const float nyquist[], 
    const double NI_MIN, const float vo[], float vradDealias[], const int nPoints){
  
    int i, j, n, m, eind, fitOk;
    double min1, esum, u1, v1, min2, dmy;
    
    // number of rows
    m = DEALIAS_VAF;
    // number of columns
    n = DEALIAS_NF;
    
    // max number of folds of nyquist interval to test for
    double MVA=2*ceil(DEALIAS_VMAX/(2*NI_MIN));

    // polarscan matrix, torus projected x coordinate, eq. 6 Haase et al. 2004 jaot
    double *x = RAVE_CALLOC ((size_t)nPoints, sizeof(double));
    // polarscan matrix, torus projected y coordinate, eq. 7 Haase et al. 2004 jaot
    double *y = RAVE_CALLOC ((size_t)nPoints, sizeof(double));
    // U-components of test velocity fields
    double *uh = RAVE_CALLOC ((size_t)(m*n), sizeof(double));
    // V-components of test velocity fields
    double *vh = RAVE_CALLOC ((size_t)(m*n), sizeof(double));
    // radial velocities of the best fitting test field
    double *vt1 = RAVE_CALLOC ((size_t)nPoints, sizeof(double));
    // array with trigonometric conversions of the points array
    float *pointsTrigon = RAVE_CALLOC ((size_t)(3*nPoints), sizeof(double));
    
    // map measured data to 3D
    for (i=0; i<nPoints; i++) {
        x[i] = nyquist[i]/M_PI * cos(vo[i]*M_PI/nyquist[i]);
        y[i] = nyquist[i]/M_PI * sin(vo[i]*M_PI/nyquist[i]);
    }

    // trigonometric conversion points array (compute once to speed up code)
    for (int iPoint=0; iPoint<nPoints; iPoint++) {
        pointsTrigon[3*iPoint+0] = sin(points[nDims*iPoint]*DEG2RAD);
        pointsTrigon[3*iPoint+1] = cos(points[nDims*iPoint]*DEG2RAD);
        pointsTrigon[3*iPoint+2] = cos(points[nDims*iPoint+1]*DEG2RAD);
    }

    // Setting up the u and v component of the test velocity fields:
    // index n=DEALIAS_NF gives number of azimuthal directions (default n=40, i.e. steps of 360/40=9 degrees)
    // index m=DEALIAS_VAF/NI_MIN*DEALIAS_VMAX gives number of speeds (maximum speed is DEALIAS_VMAX, steps of NI_MIN/DEALIAS_VAF)
    for (i=0; i<n; i++) {
        for (j=0; j<m; j++) {
            *(uh+i*m+j) = DEALIAS_VMAX/DEALIAS_VAF*(j+1) * sin(2*M_PI/DEALIAS_NF*i);
            *(vh+i*m+j) = DEALIAS_VMAX/DEALIAS_VAF*(j+1) * cos(2*M_PI/DEALIAS_NF*i);
        }
    }

    min1 = 1e32;
    eind = 0;
    u1 = 0;
    v1 = 0;

    gsl_vector *uv;
    uv = gsl_vector_alloc(2);
    
    void *params[7] = {(void *) points, (void *) pointsTrigon, (void *) &nPoints, (void *) &nDims, (void *) x, (void *) y, (void *) nyquist};     
   
    // try several test velocity fields for use as starting point in GSL fit

    for (i=0; i<m*n; i++) {
        
        gsl_vector_set(uv, 0, *(uh+i));
        gsl_vector_set(uv, 1, *(vh+i));
        esum = test_field_gsl(uv, &params);
                
        if (esum<min1) {
            min1 = esum;
            eind = i;
        }
        u1 = *(uh+eind);
        v1 = *(vh+eind);
    }
    gsl_vector_set(uv, 0, u1);
    gsl_vector_set(uv, 1, v1);

    
    #ifdef FPRINTFON
    fprintf(stdout,"Start dealiasing at (x,y)=%f,%f at f()=%f ...\n",u1,v1,esum);
    #endif
    
    fitOk = fit_field_gsl(uv, &params);
    if(!fitOk) goto cleanup;
        
    // the radial velocity of the best fitting test velocity field:
    for (int iPoint=0; iPoint<nPoints; iPoint++) {
        *(vt1+iPoint) = (u1*sin(points[nDims*iPoint]*DEG2RAD) + v1*cos(points[nDims*iPoint]*DEG2RAD))
                        *cos(points[nDims*iPoint+1]*DEG2RAD);
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

    cleanup:
        RAVE_FREE(x);
        RAVE_FREE(y);
        RAVE_FREE(uh);
        RAVE_FREE(vh);
        RAVE_FREE(vt1);
        RAVE_FREE(pointsTrigon);
        gsl_vector_free(uv);
        
        if(fitOk) return 1;
        else return 0;
}

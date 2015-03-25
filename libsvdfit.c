/*
 * Copyright 2013 Netherlands eScience Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */




// ****************************************************************************
// Functions for linear fitting using Singular Value Decomposition:
// ****************************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "polarvolume.h"
#include "libsvdfit.h"





int svd_vvp1func(const float points[], const int nDims, float afunc[], const int nParsFitted) {

    // ************************************************************************************
    // This function contains the basis-functions and other relevant parameters of
    // the fit model. The function returns the basis functions of the
    // multi-parameter linear fit evaluated at vector X=points[0..nDims-1] in the array
    // afunc[0..nParsFitted-1].
    // This function is intended to be supplied to the fitting routine 'svdfit'.
    // In this case, a three-parameter linear Area-VVP fit is performed on a two
    // dimensional vector space 'X' in radar coordinates: azimuth (alpha), and
    // elevation (gamma).
    //
    // FIXME this description is unclear
    //
    // *************************************************************************************

    float sinAlpha;
    float cosAlpha;
    float sinGamma;
    float cosGamma;

    if (nDims != 2) {
        fprintf(stderr, "Number of dimensions is wrong!\n");
        return -1;
    }
    if (nParsFitted != 3) {
        fprintf(stderr, "Number of parameters is wrong!\n");
        return -1;
    }

    sinAlpha = sin(points[0] * DEG2RAD);
    cosAlpha = cos(points[0] * DEG2RAD);
    sinGamma = sin(points[1] * DEG2RAD);
    cosGamma = cos(points[1] * DEG2RAD);
    
    afunc[0] = sinAlpha * cosGamma;   // u
    afunc[1] = cosAlpha * cosGamma;   // v
    afunc[2] = sinGamma;              // w

    return 0;

} //svd_vvp1func





int svdcmp(float *a,int m,int n,float w[],float *v) {

    // ************************************************************************************
    // Given a matrix a[0..m*n-1], this function computes its singular value
    // decomposition, A=U.W.V^T. The matrix U replaces A on output. The diagonal
    // matrix of singular values W is output as a vector w[0..n-1]. The matrix V
    // (not the transpose V^T) is output as v[0..n*n-1]. This function has
    // been modified from the original one in Numerical Recipes (2nd ed., paragraph
    // 2.6). Indices in the matrices now run from 0 to n-1/m-1, and allocation of
    // arrays is simplified. In addition, the calculation of 'pythagoras' has been
    // incorporated explicitly. See Numerical Recipes for more details.
    // ************************************************************************************


    //
    // FIXME it would be nice to have more meaningful variable names
    //


    float anorm;
    float c;
    float f;
    int flag;
    float g;
    float h;
    int i;
    int iIteration;
    int j;
    int jj;
    int k;
    int l;
    int nIterationsMax;
    int nm;
    float *rv1;
    float s;
    float scale;
    float x;
    float y;
    float z;

    nIterationsMax = 30;

    l = 0;
    nm = 0;

    /*Allocation of memory.*/

    rv1 = (float *)malloc(n*sizeof(float));
    if (rv1 == NULL) {
        fprintf(stderr, "Requested memory could not be allocated!\n");
        return -1;
    }

    /*Start of very stable algorithm by Forsythe et al.*/
    /*Householder reduction to bidiagonal form.*/

    g = 0.0;
    scale = 0.0;
    anorm = 0.0;
    for (i = 0; i < n; i++) {
        l = i + 1;
        rv1[i] = scale * g;
        g = 0.0;
        s = 0.0;
        scale = 0.0;
        if (i < m) {
            for (k = i; k < m; k++) {
                scale += fabs(a[i+n*k]);
            }
            if (scale != 0) {
                for (k = i; k < m; k++) {
                    a[i+n*k] /= scale;
                    s += a[i+n*k] * a[i+n*k];
                }
                f = a[i+n*i];
                g = -sqrt(s) * SIGN(f);
                h = f * g - s;
                a[i+n*i] = f - g;

                for (j = l; j < n; j++) {
                    s = 0.0;
                    for (k = i; k < m; k++) {
                        s += a[i+n*k] * a[j+n*k];
                    }
                    f = s / h;
                    for (k = i; k < m; k++) {
                        a[j+n*k] += f * a[i+n*k];
                    }
                }
                for (k = i; k < m; k++) {
                    a[i+n*k] *= scale;
                }
            }
        }
        w[i] = scale * g;
        g = 0.0;
        s = 0.0;
        scale = 0.0;
        if (i < m && i != (n-1)) {
            for (k = l; k < n; k++) {
                scale += fabs(a[k+n*i]);
            }
            if (scale != 0) {
                for (k = l; k < n; k++) {
                    a[k+n*i] /= scale;
                    s += a[k+n*i] * a[k+n*i];
                }

                f = a[l+n*i];
                g = -sqrt(s) * SIGN(f);
                h = f * g - s;
                a[l+n*i] = f - g;
                for (k = l; k < n; k++) {
                    rv1[k] = a[k+n*i] / h;
                }

                for (j = l; j < m; j++) {
                    s = 0.0;
                    for (k = l; k < n; k++) {
                        s += a[k+n*j] * a[k+n*i];
                    }

                    for (k = l; k < n; k++) {
                        a[k+n*j] += s * rv1[k];
                    }
                }

                for (k = l; k < n; k++) {
                    a[k+n*i] *= scale;
                }
            }
        }
        anorm = XYMAX(anorm,(fabs(w[i]) + fabs(rv1[i])));
    }

    /*Accumulation of right-hand transformations.*/

    for (i = n-1; i >= 0; i--) {
        if (i < n-1) {
            if (g != 0) {
                for (j = l; j < n; j++) {
                    v[i+n*j] = (a[j+n*i] / a[l+n*i]) / g;
                }
                for (j = l; j < n; j++) {
                    s = 0.0;
                    for (k = l; k < n; k++) {
                        s += a[k+n*i] * v[j+n*k];
                    }

                    for (k = l; k < n; k++) {
                        v[j+n*k] += s * v[i+n*k];
                    }
                }
            }
            for (j = l; j < n; j++) {
                v[j+n*i] = v[i+n*j] = 0.0;
            }
        }
        v[i+n*i] = 1.0;
        g = rv1[i];
        l = i;
    }

    /*Accumulation of left-hand transformations.*/

    for (i = XYMIN(m,n) - 1; i >= 0; i--) {
        l = i + 1;
        g = w[i];
        for (j = l; j < n; j++) {
            a[j+n*i] = 0.0;
        }
        if (g != 0) {
            g = 1.0 / g;
            for (j = l; j < n; j++) {
                s = 0.0;
                for (k = l; k < m; k++) {
                    s += a[i+n*k] * a[j+n*k];
                }
                f = (s / a[i+n*i]) * g;
                for (k = i; k < m; k++) {
                    a[j+n*k] += f * a[i+n*k];
                }
            }
            for (j = i; j < m; j++) {
                a[i+n*j] *= g;
            }
        }
        else {
            for (j = i; j < m; j++) {
                a[i+n*j] = 0.0;
            }
        }
        ++a[i+n*i];
    }

    /*Diagonalization of the bidiagonal form: Loop over singular values, and */
    /*over allowed iterations.*/

    for (k = n-1; k >= 0; k--) {
        for (iIteration = 1; iIteration <= nIterationsMax; iIteration++) {
            flag = 1;
            for (l = k; l >= 0; l--) {
                nm = l - 1;
                if ((float)(fabs(rv1[l]) + anorm) == anorm) {
                    flag = 0;
                    break;
                }
                if ((float)(fabs(w[nm]) + anorm) == anorm) {
                    break;
                }
            }
            if (flag != 0) {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++) {
                    f = s * rv1[i];
                    rv1[i] = c * rv1[i];
                    if ((float)(fabs(f) + anorm) == anorm) {
                        break;
                    }
                    g = w[i];
                    if (fabs(f) > fabs(g)) {
                        h = fabs(f) * sqrt(1 + SQUARE(g / f));
                    }
                    else {
                        h = (fabs(g) == 0 ? 0.0 : fabs(g) * sqrt(1 + SQUARE(f / g)));
                    }
                    w[i] = h;
                    h = 1.0 / h;
                    c = g * h;
                    s = -f * h;
                    for (j = 0; j < m; j++) {
                        y = a[nm+n*j];
                        z = a[i+n*j];
                        a[nm+n*j] = y * c + z * s;
                        a[i+n*j] = z * c - y * s;
                    }
                }
            }

            /*Convergence. Singular value is made nonnegative.*/

            z = w[k];
            if (l == k) {
                if (z < 0.0) {
                    w[k] = -z;
                    for (j = 0; j < n; j++) {
                        v[k+n*j] = -v[k+n*j];
                    }
                }
                break;
            }
            if (iIteration == nIterationsMax) {
                fprintf(stderr, "No convergence in %d svdcmp iterations!\n",nIterationsMax);
                return -1;
            }
            x = w[l];
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            if (fabs(f) > 1.0) {
                g = fabs(f) * sqrt(1 + SQUARE(1 / f));
            }
            else {
                g = sqrt(1 + SQUARE(f));
            }
            f = ((x - z) * (x + z) + h * ((y / (f + g * SIGN(f))) - h)) / x;

            /*Next QR transformation.*/

            c = 1.0;
            s = 1.0;
            for (j = l; j <= nm; j++) {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                if (fabs(f) > fabs(h)) {
                    z = fabs(f) * sqrt(1 + SQUARE(h / f));
                }
                else {
                    z = (fabs(h) == 0 ? 0.0 : fabs(h) * sqrt(1 + SQUARE(f / h)));
                }
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y *= c;
                for (jj = 0; jj < n; jj++) {
                    x = v[j+n*jj];
                    z = v[i+n*jj];
                    v[j+n*jj] = x * c + z * s;
                    v[i+n*jj] = z * c - x * s;
                }
                if (fabs(f) > fabs(h)) {
                    z = fabs(f) * sqrt(1 + SQUARE(h / f));
                }
                else {
                    z = (fabs(h) == 0 ? 0.0 : fabs(h) * sqrt(1 + SQUARE(f / h)));
                }

                /*Rotation can be arbitrary if z=0.*/

                w[j] = z;
                if (z != 0) {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = c * g + s * y;
                x = c * y - s * g;
                for (jj = 0; jj < m; jj++) {
                    y = a[j+n*jj];
                    z = a[i+n*jj];
                    a[j+n*jj] = y * c + z * s;
                    a[i+n*jj] = z * c - y * s;
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }

    /*Cleaning of memory.*/

    free(rv1);

    return 0;

} //svdcmp





float svdfit(const float *points, const int nDims, const float vradObs[], float vradFitted[], const int nPoints,
             float parameterVector[], float avar[], const int nParsFitted) {


    // ************************************************************************************************
    // This function performs a multi-dimensional linear fit using singular value
    // decomposition. See Numerical Recipes (2nd ed., paragraph 15.4) for details.
    // Given a set of data points points[0..nPoints*nDims-1], vradObs[0..nPoints-1], use Chi-square
    // minimization to determine the coefficients a[0..nParsFitted-1] of the fitting
    // function y=SUM_i[apar_i*afunc_i(points)]. The dimensionality of the input vector
    // 'points' is equal to 'nDims', generally this will be equal to 1.
    // Here we solve the fitting equations using singular value decomposition of
    // the nPoints by nParsFitted matrix. The program returns values for the array with fit
    // parameters 'parameterVector', variance matrix 'avar', and the Chi-square fitness score.
    //
    // The user supplies a function with the fit model 'funcs(points,nDims,afunc,nParsFitted)'
    // that returns the 'nParsFitted' basis functions evaluated at points[0..nDims-1] in the
    // array afunc[0..nParsFitted-1].
    //
    // ************************************************************************************************


    int iParFitted;
    int iParFittedRows;
    int iParFittedCols;
    int k;
    int iPoint;
    float afunc[nParsFitted];                // FIXME order of afunc parameters is not u,v,w?
    float singularValues[nParsFitted];
    float v[nParsFitted*nParsFitted];
    float wti[nParsFitted];
    float *u;
    float singularValueMax;
    float sum;
    float chisq;

    // Checking whether the input numbers are within bounds:
    if (nParsFitted > NPARSFITTEDMAX) {
        fprintf(stderr, "Number of fit parameters is too large!\n");
        return -1.0;
    }
    if (nPoints <= nParsFitted) {
        fprintf(stderr, "Number of data points is too small!\n");
        return -1.0;
    }

    // Allocation of memory for arrays.
    u = (float *)malloc(nPoints*nParsFitted*sizeof(float));
    if (!u) {
        fprintf(stderr, "Requested memory could not be allocated!\n");
        return -1.0;
    }

    // Filling of the design matrix of the fitting problem (u[iPoint][iParFitted]).
    for (iPoint = 0; iPoint < nPoints; iPoint++) {

        // note pointer arithmetic in this next statement:
        if (svd_vvp1func(points+nDims*iPoint,nDims,afunc,nParsFitted)) {
            return -1.0;
        }

        for (iParFitted = 0; iParFitted < nParsFitted; iParFitted++) {
            u[iParFitted+nParsFitted*iPoint] = afunc[iParFitted];
        }
    }


    // Singular value decomposition of the design matrix of the fit.
    if (svdcmp(u,nPoints,nParsFitted,singularValues,v)) {
        return -1.0;
    }

    // Removal of the singular values.
    singularValueMax = 0.0;
    for (iParFitted = 0; iParFitted < nParsFitted; iParFitted++) {
        if (singularValues[iParFitted] > singularValueMax) {
            singularValueMax = singularValues[iParFitted];
        }
    }
    for (iParFitted = 0; iParFitted < nParsFitted; iParFitted++) {
        if (singularValues[iParFitted] < SVDTOL*singularValueMax) {
            singularValues[iParFitted] = 0.0;
        }
    }



    // Calculation of fit parameters 'parameterVector' using backsubstitution with 'vradObs'.
    if (svbksb(u,singularValues,v,nPoints,nParsFitted,vradObs,parameterVector)) {
        return -1.0;
    }



    // Calculation of variances of fit parameters 'parameterVector'.
    for (iParFitted = 0; iParFitted < nParsFitted; iParFitted++) {
        wti[iParFitted] = 0.0;
        if (singularValues[iParFitted] != 0) {
            wti[iParFitted] = 1.0/(singularValues[iParFitted]*singularValues[iParFitted]);
        }
    }

    // FIXME the calculation of the variance is still a bit weird...buggy?
    for (iParFittedCols = 0; iParFittedCols < nParsFitted; iParFittedCols++) {

        avar[iParFittedCols] = 0.0;

        for (iParFittedRows = 0 ; iParFittedRows < nParsFitted; iParFittedRows++) {

            k = iParFittedCols + nParsFitted*iParFittedRows;
            avar[iParFittedCols] += v[k] * v[k] * wti[iParFittedRows];

        }
    }

    /*Calculation of vradFitted and Chi-square of the fit.*/
    chisq = 0.0;
    for (iPoint = 0; iPoint < nPoints; iPoint++) {

        // note pointer arithmetic in this next statement:
        if (svd_vvp1func(points+nDims*iPoint,nDims,afunc,nParsFitted)) {
            return -1.0;
        }
        sum = 0.0;
        for (iParFitted = 0; iParFitted < nParsFitted; iParFitted++) {
            sum += parameterVector[iParFitted] * afunc[iParFitted];
        }
        vradFitted[iPoint] = sum;

        chisq += SQUARE(vradObs[iPoint]-vradFitted[iPoint]);

    }
    chisq /= nPoints-nParsFitted;

    /*Cleaning of memory.*/

    free(u);

    return chisq;
} //svdfit






int svbksb(float *u,float w[],float *v,int m,int n,const float b[],float x[])
{
    int jj,j,i;
    float sum,*tmp;

    /*Allocation of memory.*/

    tmp=(float *)malloc(n*sizeof(float));
    if (tmp==NULL) {
        printf("Requested memory could not be allocated!\n");
        return -1;
    }

    /*First part of inversion: calculation of Tmp = (W^-1.U^T).B. Singular values */
    /*of W are discarded.*/

    for (j=0 ; j<n ; j++) {
        sum=0.0;
        if (w[j]) {
            for (i=0 ; i<m ; i++) sum+=u[j+n*i]*b[i];
            sum/=w[j];
        }
        tmp[j]=sum;
    }

    /*Second part: calculation of X = V.Tmp = (V.W^-1.U^T).B = A^-1.B.*/

    for (j=0 ; j<n ; j++) {
        sum=0.0;
        for (jj=0 ; jj<n ; jj++) sum+=v[jj+n*j]*tmp[jj];
        x[j]=sum;
    }

    /*Cleaning of memory.*/

    free(tmp);

    return 0;
} //svbksb



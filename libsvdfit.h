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
// Definition of general macros:
// ****************************************************************************

#define SIGN(x)    (((x)<0)?-1:1)
#define SQUARE(x)  ((x)*(x))
#define XYMAX(x,y) (((x)<(y))?(y):(x))
#define XYMIN(x,y) (((x)<(y))?(x):(y))

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




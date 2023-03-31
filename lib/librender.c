#include "polarvolume.h"
#include "polarscan.h"
#include "cartesian.h"
#include "rave_debug.h"
#include "rave_list.h"
#include "rave_utilities.h"
#include "rave_alloc.h"
#include "constants.h"
#include "libvol2bird.h"
#include "librender.h"
#include <string.h>
#include <math.h>

#ifdef MISTNET
#include "../libmistnet/libmistnet.h"
#endif


/**
 * FUNCTION PROTOTYPES
 **/
double distance2range(double distance,double elev);

double distance2height(double distance,double elev);

double range2distance(double range,double elev);

double range2height(double range,double elev);

Cartesian_t* polarVolumeToCartesian(PolarVolume_t* pvol, long dim, long res, double init);

RaveObjectList_t* polarVolumeToCartesianList(PolarVolume_t* pvol, long dim, long res, double init, int *nParam);

Cartesian_t* polarScanToCartesian(PolarScan_t* scan, long dim, long res, double init);

void free4DTensor(float ****tensor, int dim1, int dim2, int dim3);

float**** create4DTensor(float *array, int dim1, int dim2, int dim3, int dim4);

int addTensorToPolarVolume(PolarVolume_t* pvol, float ****tensor, int dim1, int dim2, int dim3, int dim4, long res);

int addClassificationToPolarVolume(PolarVolume_t* pvol, float ****tensor, int dim1, int dim2, int dim3, int dim4, long res);

double*** init3DTensor(int dim1, int dim2, int dim3, double init);

void free3DTensor(double ***tensor, int dim1, int dim2);

int fill3DTensor(double ***tensor, RaveObjectList_t* list, int dim1, int dim2, int dim3);

float* flatten3DTensor(double ***tensor, int dim1, int dim2, int dim3);

int polarVolumeTo3DTensor(PolarVolume_t* pvol, double ****tensor, int dim, long res, int nParam);

PolarVolume_t* PolarVolume_selectScansByElevation(PolarVolume_t* volume, float elevs[], int nElevs);

PolarVolume_t* PolarVolume_selectScansByScanUse(PolarVolume_t* volume, vol2birdScanUse_t *scanUse, int nScansUsed);

PolarScan_t* PolarVolume_getScanClosestToElevation_vol2bird(PolarVolume_t* volume, double elev);

#ifdef MISTNET
int run_mistnet(float* tensor_in, float** tensor_out, const char* model_path, int tensor_size);
#endif

#ifndef MIN
#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#endif

#ifndef MAX
#define MAX(x,y) (((x) < (y)) ? (y) : (x))
#endif

#ifndef ABS
#define ABS(x) (((x) < 0) ? (-(x)) : (x))
#endif


/**
 * FUNCTION BODIES
 **/

/**
 * Convert from ground distance and elevation to slant range.
 *
 * Uses spherical earth
 *
 * Based on RSL_get_slantr_and_h
 * See alo Doviak and Zrnic 1993 Eqs. (2.28b) and (2.28c)
 * 
 * @param distance - range along ground (great circle distance)
 * @param elev - beam elevation in radians
 * @return range along radar path in meter
 */
double distance2range(double distance,double elev){
    
    double effectiveEarthRadius = EARTH_RADIUS * REFRACTION_COEFFICIENT;
    double alpha, beta, gamma;
    double range;
    
    /*
    Law of sines of triangle ABC

       A = center of earth
       B = radar station
       C = pulse volume
      
    gamma := angle(AB,AC) = distance/effectiveEarthRadius             
    alpha := angle(BA,BC) = 90 + elev          
    beta  := angle(CB,CA) = pi - alpha - gamma 

    Law of sines says:
    effectiveEarthRadius/sin(beta) = (effectiveEarthRadius + h)/sin(alpha) = r/sin(gamma)

    We know effectiveEarthRadius, so we can solve for (effectiveEarthRadius + h) and r

    */
    gamma = distance/effectiveEarthRadius;
    alpha = PI/2+elev;
    beta = PI - alpha - gamma;

    range = effectiveEarthRadius * (sin(gamma)/sin(beta));
    
    return range;
}

/**
 * Convert from ground distance and elevation to height.
 *
 * Uses spherical earth
 *
 * Based on RSL_get_slantr_and_h
 * See alo Doviak and Zrnic 1993 Eqs. (2.28b) and (2.28c)
 * 
 * @param distance - range along ground (great circle distance)
 * @param elev - beam elevation in degrees
 * @return range along radar path in meter
 */
double distance2height(double distance,double elev){
    
    double effectiveEarthRadius = EARTH_RADIUS * REFRACTION_COEFFICIENT;
    double alpha, beta, gamma;
    double height;
    
    gamma = distance/effectiveEarthRadius;
    alpha = PI/2+elev;
    beta = PI - alpha - gamma;

    height = effectiveEarthRadius * (sin(alpha)/sin(beta)) - effectiveEarthRadius;
    
    return height;
}

/**
 * Convert from slant range and elevation to ground distance.
 *
 * Uses spherical earth
 *
 * See also:
 * Doviak and Zrnic 1993 Eqs. (2.28b) and (2.28c)
 * https://bitbucket.org/deeplycloudy/lmatools/src/3ad332f9171e/coordinateSystems.py?at=default
 * 
 * @param range - slant range along radar line of sight in meter
 * @param elev - beam elevation in radians
 * @return distance range along ground (great circle distance) in meter
 */
double range2distance(double range,double elev){
    
    double effectiveEarthRadius = EARTH_RADIUS * REFRACTION_COEFFICIENT;
    double distance;
    double height;

    height = range2height(range, elev);
    distance = effectiveEarthRadius * asin(range * cos(elev) / ( effectiveEarthRadius + height ) );
    
    return(distance);
}

/**
 * Convert from slant range and elevation to ground distance.
 *
 * Uses spherical earth
 *
 * See also:
 * Doviak and Zrnic 1993 Eqs. (2.28b) and (2.28c)
 * https://bitbucket.org/deeplycloudy/lmatools/src/3ad332f9171e/coordinateSystems.py?at=default
 * 
 * @param range - slant range along radar line of sight in meter
 * @param elev - beam elevation in radians
 * @return height above ground in meter
 */
double range2height(double range,double elev){
    
    double effectiveEarthRadius = EARTH_RADIUS * REFRACTION_COEFFICIENT;
    double height;

    height = sqrt(SQUARE(range) + SQUARE(effectiveEarthRadius) + (2 * effectiveEarthRadius * range * sin(elev))) - effectiveEarthRadius;
    
    return(height);
}


Cartesian_t* polarVolumeToCartesian(PolarVolume_t* pvol, long dim, long res, double init){
    
    RAVE_ASSERT((pvol != NULL), "pvol == NULL");
        
    // initialize scan, param, and cartesian RAVE objects
    PolarScan_t *scan = NULL;
    Cartesian_t *cartesian = NULL;
    CartesianParam_t *cartesianParam = NULL;
    RaveList_t* scanParameterNames;
    char* scanParameterName;
    char* parameterName;
    
    // create a new Cartesian grid object
    cartesian = RAVE_OBJECT_NEW(&Cartesian_TYPE);

    // copy metadata from volume
    Cartesian_setTime(cartesian, PolarVolume_getTime(pvol));
    Cartesian_setDate(cartesian, PolarVolume_getDate(pvol));
    Cartesian_setSource(cartesian, PolarVolume_getSource(pvol));

    //set cartesian product and object type
    Cartesian_setObjectType(cartesian, Rave_ObjectType_IMAGE);
    Cartesian_setProduct(cartesian, Rave_ProductType_PPI);
    
    if (cartesian == NULL){
        vol2bird_err_printf("failed to allocate memory for new cartesian object\n");
        return NULL;
    }

    //set dimensions and resolution of the grid
    Cartesian_setXSize(cartesian, dim);
    Cartesian_setYSize(cartesian, dim);
    Cartesian_setXScale(cartesian, res);
    Cartesian_setYScale(cartesian, res);
    //Cartesian_setAreaExtent(cartesian, -res*dim/2, -res*dim/2, res*dim/2, res*dim/2);
    
    int nScans;
    // determine how many scan elevations the volume object contains
    nScans = PolarVolume_getNumberOfScans(pvol);
    
    if(nScans<=0){
        vol2bird_err_printf("Error: polar volume contains no scans\n");
        return NULL;
    }

    // iterate over the selected scans in 'volume'
    for (int iScan = 0; iScan < nScans; iScan++) {
        
        // extract the scan object from the volume object
        scan = PolarVolume_getScan(pvol,iScan);
        
        double elev = PolarScan_getElangle(scan);
        
        scanParameterNames = PolarScan_getParameterNames(scan);
        
        if(RaveList_size(scanParameterNames)<=0){
            vol2bird_err_printf("Warning: ignoring scan without scan parameters\n");
            continue;            
        }
                
        for(int iParam = 0; iParam<RaveList_size(scanParameterNames); iParam++){
            // retrieve name of the scan parameter
            scanParameterName = RaveList_get(scanParameterNames, iParam);
            
            char iElevString[11];
            // copy iElev to iElevString
            sprintf(iElevString, "%d", iScan);

            // create a new scan parameter name with index for the sweep
            char *parameterNameFull = malloc(strlen(scanParameterName)+strlen(iElevString)+1);
            strcpy(parameterNameFull,scanParameterName);
            strcat(parameterNameFull,iElevString);
            parameterName = RaveUtilities_trimText(parameterNameFull, strlen(parameterNameFull));

            // create a cartesian scan parameter with the same name
            cartesianParam = Cartesian_createParameter(cartesian,parameterName,RaveDataType_DOUBLE, init);
            CartesianParam_setNodata(cartesianParam, PolarScanParam_getNodata(PolarScan_getParameter(scan, scanParameterName)));
            CartesianParam_setUndetect(cartesianParam, PolarScanParam_getUndetect(PolarScan_getParameter(scan, scanParameterName)));
            
            double range,azim,distance;
            double value;
            
            RaveValueType a;
            // loop over the grid, and fill it
            for(long x = 0; x<dim; x++){
                for(long y = 0; y<dim; y++){
                    double xx=((double)res)*((double)(x-dim/2));
                    double yy=((double)res)*((double)(y-dim/2));
                    azim=atan2(yy,xx);
                    distance=sqrt(SQUARE(xx)+SQUARE(yy));
                    range=distance2range(distance,elev);
                    a=PolarScan_getConvertedParameterValueAtAzimuthAndRange(scan,scanParameterName,azim,range,&value);
                    if(a!=RaveValueType_DATA){
                        PolarScan_getParameterValueAtAzimuthAndRange(scan,scanParameterName,azim,range,&value);
                    }
                    CartesianParam_setValue(cartesianParam, x, y, value);
                }
            }
            
            // add the cartesian scan parameter to the cartesian object
            Cartesian_addParameter(cartesian,cartesianParam);
            
            free(parameterNameFull);
            RAVE_FREE(parameterName);
            RAVE_OBJECT_RELEASE(cartesianParam);
        } // iParam
        
    } // iElev
    
    return cartesian;
}




RaveObjectList_t* polarVolumeToCartesianList(PolarVolume_t* pvol, long dim, long res, double init, int *nParam){
    
    PolarScan_t *scan = NULL;
    RaveObjectList_t* list;
    Cartesian_t *cartesian = NULL;
    
    // create a new list
    list = (RaveObjectList_t *) RAVE_OBJECT_NEW(&RaveObjectList_TYPE);
    
    int nScans;
    // determine how many scan elevations the volume object contains
    nScans = PolarVolume_getNumberOfScans(pvol);
    
    if(nScans<=0){
        vol2bird_err_printf("Error: polar volume contains no scans\n");
        return NULL;
    }

    // iterate over the selected scans in 'volume'
    for (int iScan = 0; iScan < nScans; iScan++) {
                
        // extract the scan object from the volume object
        scan = PolarVolume_getScan(pvol, iScan);
                
        cartesian = polarScanToCartesian(scan, dim, res, init);
        
        *nParam += Cartesian_getParameterCount(cartesian);
        
        RaveObjectList_add(list, (RaveCoreObject*) cartesian); 
        
        RAVE_OBJECT_RELEASE(cartesian);

        RAVE_OBJECT_RELEASE(scan);

    }
    
    return list;

}



Cartesian_t* polarScanToCartesian(PolarScan_t* scan, long dim, long res, double init){
    
    RAVE_ASSERT((scan != NULL), "scan == NULL");
        
    // initialize scan, param, and cartesian RAVE objects
    Cartesian_t *cartesian = NULL;
    CartesianParam_t *cartesianParam = NULL;
    PolarScanParam_t* polarScanParam = NULL;

    RaveList_t* scanParameterNames;
    char* scanParameterName;
    
    // create a new Cartesian grid object
    cartesian = RAVE_OBJECT_NEW(&Cartesian_TYPE);
    if (cartesian == NULL){
        vol2bird_err_printf( "failed to allocate memory for new cartesian object\n");
        return NULL;
    }

    // copy metadata from volume
    Cartesian_setTime(cartesian, PolarScan_getTime(scan));
    Cartesian_setDate(cartesian, PolarScan_getDate(scan));
    Cartesian_setSource(cartesian, PolarScan_getSource(scan));

    //set cartesian product and object type
    Cartesian_setObjectType(cartesian, Rave_ObjectType_IMAGE);
    Cartesian_setProduct(cartesian, Rave_ProductType_PPI);
    
    //set dimensions and resolution of the grid
    Cartesian_setXSize(cartesian, dim);
    Cartesian_setYSize(cartesian, dim);
    Cartesian_setXScale(cartesian, res);
    Cartesian_setYScale(cartesian, res);
    //Cartesian_setAreaExtent(cartesian, -res*dim/2, -res*dim/2, res*dim/2, res*dim/2);

    
    double elev = PolarScan_getElangle(scan);
    
    scanParameterNames = PolarScan_getParameterNames(scan);
    
    if(RaveList_size(scanParameterNames)<=0){
        vol2bird_err_printf("Warning: scan without scan parameters\n");
        RaveList_freeAndDestroy(&scanParameterNames);
        RAVE_OBJECT_RELEASE(cartesian);
        return NULL;
    }
            
    for(int iParam = 0; iParam<RaveList_size(scanParameterNames); iParam++){
        // retrieve name of the scan parameter
        scanParameterName = (char*)RaveList_get(scanParameterNames, iParam);
        polarScanParam = PolarScan_getParameter(scan, scanParameterName);
        
        // create a cartesian scan parameter with the same name
        cartesianParam = Cartesian_createParameter(cartesian, scanParameterName, RaveDataType_DOUBLE, init);
        CartesianParam_setNodata(cartesianParam, PolarScanParam_getNodata(polarScanParam));
        CartesianParam_setUndetect(cartesianParam, PolarScanParam_getUndetect(polarScanParam));

        double range,azim,distance;
        double value;
        
        RaveValueType a;
        // loop over the grid, and fill it
        for(long x = 0; x<dim; x++){
            for(long y = 0; y<dim; y++){
                double xx=((double)res)*((double)(x-dim/2));
                double yy=((double)res)*((double)(y-dim/2));
                azim=atan2(yy,xx);
                distance=sqrt(SQUARE(xx)+SQUARE(yy));
                range=distance2range(distance,elev);
                a = PolarScan_getConvertedParameterValueAtAzimuthAndRange(scan, scanParameterName, azim, range, &value);
                if(a != RaveValueType_DATA){
                    PolarScan_getParameterValueAtAzimuthAndRange(scan, scanParameterName, azim, range, &value);
                }
                CartesianParam_setValue(cartesianParam, x, y, value);
            }
        }
        
        // add the cartesian scan parameter to the cartesian object
        Cartesian_addParameter(cartesian, cartesianParam);
        
        RAVE_OBJECT_RELEASE(polarScanParam);
        RAVE_OBJECT_RELEASE(cartesianParam);
    } // iParam
    
    RaveList_freeAndDestroy(&scanParameterNames);
    return cartesian;
}


float**** create4DTensor(float *array, int dim1, int dim2, int dim3, int dim4) {
    float ****tensor = (float ****)malloc(dim1 * sizeof(float***));
    for(int i=0 ; i < dim1 ; i++) {
        tensor[i] = (float ***) malloc(dim2 * sizeof(float**));
        for(int j=0 ; j < dim2 ; j++) {
            tensor[i][j] = (float **)malloc(dim3 * sizeof(float*));
            for (int k=0 ; k < dim3 ; k++){
                tensor[i][j][k] = (float *)malloc(dim4 * sizeof(float));
                for (int l=0; l<dim4; l++){
                    tensor[i][j][k][l] = array[i * dim4 * dim3 * dim2 + j * dim4 * dim3 + k * dim4 + l];
                }                    
            }
        }
    }

    return tensor;
}



void free4DTensor(float ****tensor, int dim1, int dim2, int dim3){
	// deallocate memory
	for (int i = 0; i < dim1; i++) 
	{
		for (int j = 0; j < dim2; j++){
            for (int k = 0; k < dim3; k++){
                free(tensor[i][j][k]);
            }
            free(tensor[i][j]);
        }
		free(tensor[i]);
	}
	free(tensor);
}

double*** init3DTensor(int dim1, int dim2, int dim3, double init)
{
    double ***tensor = (double ***)malloc(dim1*sizeof(double**));
    
    if(tensor == NULL){
        vol2bird_err_printf("failed to initialize 3D tensor (1)");
#ifdef VOL2BIRD_R
        return NULL;
#else        
        exit(0);
#endif        
    }

    for (int i = 0; i < dim1; i++) {
        tensor[i] = NULL;
    }

    for (int i = 0; i< dim1; i++) {
        tensor[i] = (double **) malloc(dim2*sizeof(double *));
        
        if(tensor[i] == NULL){
            vol2bird_err_printf("failed to initialize 3D tensor (2)");
            free3DTensor(tensor, dim1, dim2);
#ifdef VOL2BIRD_R
            return NULL;
#else
            exit(0);
#endif            
        } else {
            for (int j = 0; j < dim2; j++) {
                tensor[i][j] = NULL;
            }
        }
        
        for (int j = 0; j < dim2; j++) {
            tensor[i][j] = (double *)malloc(dim3*sizeof(double));
            if(tensor[i][j] == NULL){
                vol2bird_err_printf("failed to initialize 3D tensor (3)");
	              free3DTensor(tensor, dim1, dim2);
#ifdef VOL2BIRD_R	              
                return NULL;
#else                
                exit(0);
#endif                
            }
        } // j
    } // i
    
    // assign values to allocated memory
    for (int i = 0; i < dim1; i++) {
        for (int j = 0; j < dim2; j++) {
            for (int k = 0; k < dim3; k++) {
                tensor[i][j][k] = init;
            }
        }
    }

    return tensor;
}

void free3DTensor(double ***tensor, int dim1, int dim2)
{
  if (tensor != NULL) {
    for (int i = 0; i < dim1; i++) {
      if (tensor[i] != NULL) {
        for (int j = 0; j < dim2; j++) {
          if (tensor[i][j] != NULL) {
            free(tensor[i][j]);
          }
        }
        free(tensor[i]);
      }
    }
    free(tensor);
  }
}

int fill3DTensor(double ***tensor, RaveObjectList_t* list, int dim1, int dim2, int dim3){

    int nScan = RaveObjectList_size(list);
    
    for(int iScan = 0; iScan<nScan; iScan++){
        
        Cartesian_t* cartesian = (Cartesian_t*) RaveObjectList_get(list, iScan);
        
        int nCartesianParam = Cartesian_getParameterCount(cartesian);
        long xSize	= Cartesian_getXSize(cartesian);
        long ySize	= Cartesian_getYSize(cartesian);
        
        if(dim2 != xSize){
            vol2bird_err_printf( "Error: expecting %i bins in X dimension, but found only %li\n", dim2, xSize);
            RAVE_OBJECT_RELEASE(cartesian);
            return -1;
        }

        if(dim3 != ySize){
            vol2bird_err_printf( "Error: expecting %i bins in Y dimension, but found only %li\n", dim3, ySize);
            RAVE_OBJECT_RELEASE(cartesian);
            return -1;
        }
        
        RaveList_t* cartesianParameterNames = Cartesian_getParameterNames(cartesian);
        double value;
        RaveValueType valueType;
        int dbz_count = 0;
        int vrad_count = 0;
        int wrad_count = 0;
        
        for(int iOrder = 0; iOrder < 3; iOrder++){
            for(int iCartesianParam = 0; iCartesianParam < nCartesianParam; iCartesianParam++){
                char* parameterName = (char *) RaveList_get(cartesianParameterNames, iCartesianParam);
                
                // make sure parameters are stored in order DBZ, VRAD, WRAD, RHOHV
                switch(iOrder){
                    case 0:
                       if(strncmp("DBZ",parameterName,3)!=0){
                           continue;
                       }
                       break;
                    case 1:
                       if(strncmp("VRAD",parameterName,4)!=0){
                           continue;
                       }
                       break;
                    case 2:
                       if(strncmp("WRAD",parameterName,4)!=0){
                           continue;
                       }
                       break;
                    case 3:
                       if(strncmp("RHOHV",parameterName,5)!=0){
                           // note: this case is never selected because iOrder<3
                           continue;                       
                       }
                       break;
                }
                
                CartesianParam_t* cartesianParam = Cartesian_getParameter(cartesian, parameterName);
                
                #ifdef FPRINTFON
                vol2bird_err_printf("Writing Cartesian parameter %s at index %i (scan=%i, param=%i) \n",
                    parameterName,iScan+nScan*iOrder, iScan, iOrder);
                #endif

                if(iScan+nScan*iOrder>=dim1){
                   vol2bird_err_printf( "Error: exceeding 3D tensor dimension\n");
                   RaveList_freeAndDestroy(&cartesianParameterNames);
                   RAVE_OBJECT_RELEASE(cartesian);
                   RAVE_OBJECT_RELEASE(cartesianParam);
                   return(-1);
                }
                
                if(iOrder == 0) dbz_count+=1;
                if(iOrder == 1) vrad_count+=1;
                if(iOrder == 2) wrad_count+=1;
                
                // fill tensor
                for(int x = 0; x < xSize; x++){
                    for(int y = 0; y < ySize; y++){
                        valueType = CartesianParam_getValue(cartesianParam, x, y, &value);
                        
                        if (valueType == RaveValueType_DATA){
                            // only copy radial velocity and spectrum width values that have a corresponding reflectivity value
                            // this is to account for occasional sweeps where radial velocity extends to shorter ranges than reflectivity
                            if(MISTNET_REQUIRE_DBZ && (iOrder > 0) && isnan(tensor[iScan][x][y])){
                                tensor[iScan+nScan*iOrder][x][y] = NAN;
                            }
                            else{
                                tensor[iScan+nScan*iOrder][x][y] = value;
                            }
                        }
                        else{
                            tensor[iScan+nScan*iOrder][x][y] = NAN;
                        }
                    } //y
                } //x
                
                RAVE_OBJECT_RELEASE(cartesianParam);
                
            } // iParam
        } // iOrder
        
        if(dbz_count == 0) vol2bird_err_printf( "Warning: no reflectivity data found for MistNet input scan %i, initializing with values %i instead.\n", iScan, MISTNET_INIT);
        if(vrad_count == 0) vol2bird_err_printf( "Warning: no radial velocity data found for MistNet input scan %i, initializing with values %i instead.\n", iScan, MISTNET_INIT);
        if(wrad_count == 0) vol2bird_err_printf( "Warning: no spectrum width data found for MistNet input scan %i, initializing with values %i instead.\n", iScan, MISTNET_INIT);
        RaveList_freeAndDestroy(&cartesianParameterNames);
        RAVE_OBJECT_RELEASE(cartesian);
    }  // iScan
    
    return 0;
}


float* flatten3DTensor(double ***tensor, int dim1, int dim2, int dim3){
    float* output = (float *)malloc(dim1 * dim2 * dim3 *sizeof(float));
    float* temp = output;
    for (int i=0 ; i < dim1 ; i++){
        for (int j=0 ; j < dim2 ; j++){
            for (int k=0 ; k < dim3 ; k++){
                *temp = tensor[i][j][k];
                temp++;
            }
        }
    }
    return output;
}


int polarVolumeTo3DTensor(PolarVolume_t* pvol, double ****tensor, int dim, long res, int nParam){
    //Un-comment these two lines to save a rendering to file
    //Cartesian_t *cartesian = NULL;
    //cartesian = polarVolumeToCartesian(pvol, elevs, nElevs, dim, res, 0);            
    //saveToODIM((RaveCoreObject*) cartesian, "rendering.h5");
    
    // convert polar volume to a list of Cartesian objects, one for each scan
    // store the total number of scan parameters for all scans in nCartesianParam
    int nCartesianParam = 0;

    RaveObjectList_t* list = polarVolumeToCartesianList(pvol, dim, res, 0, &nCartesianParam);

    if(list == NULL){
        vol2bird_err_printf( "Error: failed to load Cartesian objects from polar volume\n");
        return -1;
    }
    // if nParam is specified, restrict the number of output parameters to its value.
    // nParam typically equals 3x5=15, selecting DBZ, VRAD and WRAD for Misnet segmentation model input.
    if(nParam > 0){
        if(nParam < nCartesianParam){
            nCartesianParam = nParam;
        }
    }

    // initialize a 3D tensor, and fill it
    *tensor = init3DTensor(nCartesianParam,dim,dim,MISTNET_INIT);
    fill3DTensor(*tensor, list, nCartesianParam, dim, dim);

    // clean up
    RAVE_OBJECT_RELEASE(list);

    return(nCartesianParam);
}


/**
 * Return a polar volume containing a selection of scans by elevation
 * 
 * @param volume - a polar volume
 * @param elevs - array with elevation angles of scans to be selected
 * @param nElevs - length of the elevation angle array
 * @return a polar volume containing only the selected scans. Note:
 * the returned volume is NOT a copy of the input volume, both objects
 * reference the same scan objects.
 */
PolarVolume_t* PolarVolume_selectScansByElevation(PolarVolume_t* volume, float elevs[], int nElevs){
    int iScan;
    int nScans;

    PolarScan_t* scan = NULL;
    PolarVolume_t* volume_select = NULL;

    // copy the volume 
    volume_select = RAVE_OBJECT_CLONE(volume); 

    nScans = PolarVolume_getNumberOfScans(volume_select);

    if(nScans<=0){
        vol2bird_err_printf("Error: polar volume contains no scans\n");
        return volume_select;
    }
    
    // get the number of elevations.
    if(nElevs>nScans){
        vol2bird_err_printf("Warning: requesting %i elevations scans, but only %i available\n", nElevs, nScans);
    }
         
    // empty the scans in the copied volume 
    for (iScan = nScans-1; iScan>=0 ; iScan--) {
        PolarVolume_removeScan(volume_select,iScan); 
    }

    // iterate over the selected scans in 'volume' and add them to 'volume_select'
    for (int iElev = 0; iElev < nElevs; iElev++) {
        // extract the scan object from the volume object
        scan = PolarVolume_getScanClosestToElevation_vol2bird(volume,DEG2RAD*elevs[iElev]);
        if (ABS(RAD2DEG*PolarScan_getElangle(scan)-elevs[iElev]) > 0.1){
            vol2bird_err_printf("Warning: Requested elevation scan at %f degrees but selected scan at %f degrees\n",
                elevs[iElev],RAD2DEG*PolarScan_getElangle(scan));
        }
        
        // add it to the selected volume
        PolarVolume_addScan(volume_select, scan);
        RAVE_OBJECT_RELEASE(scan);
    }
    
    // sort polar volume by ascending elevation
    PolarVolume_sortByElevations(volume_select, 1);
    
    return(volume_select);
}

/**
 * Return a polar volume containing a scans selected by scanUse object. 
 * 
 * @param volume - a polar volume
 * @param scanUse - a scanUse object
 * @return a polar volume containing only the selected scans. Note:
 * the returned volume is NOT a copy of the input volume, both objects
 * reference the same scan objects.
 */
PolarVolume_t* PolarVolume_selectScansByScanUse(PolarVolume_t* volume, vol2birdScanUse_t *scanUse, int nScansUsed){
    int iScan;
    int nScans;

    PolarScan_t* scan = NULL;
    PolarVolume_t* volume_select = NULL;

    // copy the volume 
    volume_select = RAVE_OBJECT_CLONE(volume); 

    nScans = PolarVolume_getNumberOfScans(volume_select);

    if(nScans<=0){
        vol2bird_err_printf("Error: polar volume contains no scans\n");
        return volume;
    }
             
    // empty the scans in the cloned volume 
    for (iScan = nScans-1; iScan>=0 ; iScan--) {
            PolarVolume_removeScan(volume_select,iScan);
    }

    // iterate over the selected scans in 'volume' and add them to 'volume_select'
    for (int iScan = 0; iScan < nScans; iScan++) {
        // extract the scan object from the volume object
        scan = PolarVolume_getScan(volume,iScan);
        // add it to the selected volume
        if (scanUse[iScan].useScan){
            PolarVolume_addScan(volume_select, scan);
        }
        RAVE_OBJECT_RELEASE(scan);
    }

    // sort polar volume by ascending elevation
    PolarVolume_sortByElevations(volume_select, 1);
    
    return(volume_select);
}


/**
 * Return the polar scan of a volume closest to a given elevation
 * This is a replacement function for PolarVolume_getScanClosestToElevation
 * available in RAVE, which fails when there are multiple scans at the same elevation. 
 *
 * @param volume - a polar volume
 * @param elev - an elevation angle in radians
 * @return a polar scan
 */
PolarScan_t* PolarVolume_getScanClosestToElevation_vol2bird(PolarVolume_t* volume, double elev){
    int nScans;
    double elevDifference = 1000;
    double elevDifferenceCandidate = 1000;

    nScans = PolarVolume_getNumberOfScans(volume);

    PolarScan_t* scan = NULL;
    PolarScan_t* scanCandidate = NULL;

    if(nScans<=0){
        vol2bird_err_printf("Error: polar volume contains no scans\n");
        return scan;
    }

    for (int iScan = 0; iScan < nScans; iScan++) {
        // extract the scan object from the volume object
        scanCandidate = PolarVolume_getScan(volume,iScan);
        elevDifferenceCandidate = ABS(elev - PolarScan_getElangle(scanCandidate));
        
        // this happens when there are two elevation scans at the same elevation
        if(elevDifferenceCandidate == elevDifference){
            // pick the higest resolution scan
            if (PolarScan_getRscale(scanCandidate) < PolarScan_getRscale(scan)){
                RAVE_OBJECT_RELEASE(scan);
                scan = RAVE_OBJECT_COPY(scanCandidate);
            }
        }
    
        if(elevDifferenceCandidate < elevDifference){
            elevDifference = elevDifferenceCandidate;
            RAVE_OBJECT_RELEASE(scan);
            scan = RAVE_OBJECT_COPY(scanCandidate);
        }
        RAVE_OBJECT_RELEASE(scanCandidate);
    }

    return(scan);
}

int addTensorToPolarVolume(PolarVolume_t* pvol, float ****tensor, int dim1, int dim2, int dim3, int dim4, long res){
    
    RAVE_ASSERT((pvol != NULL), "pvol == NULL");

    PolarScan_t* scan = NULL;

    int nScans;
    // determine how many scan elevations the volume object contains
    nScans = PolarVolume_getNumberOfScans(pvol);
    
    if(nScans != dim2){
        vol2bird_err_printf( "Error: polar volume has %i scans, while tensor has data for %i scans.\n", nScans, dim2);
    }

   // iterate over the selected scans in 'volume'
    for (int iScan = 0; iScan < nScans; iScan++) {
        // extract the scan object from the volume object
        scan = PolarVolume_getScan(pvol,iScan);
        
        if(PolarScan_hasParameter(scan, "WEATHER")){
            vol2bird_err_printf( "Warning: scan used multiple times as MistNet input, ignoring segmentation %i/%i\n", iScan+1, MISTNET_N_ELEV);
            RAVE_OBJECT_RELEASE(scan);
            continue;
        }
        
        PolarScanParam_t *mistnetParamWeather = PolarScan_newParam(scan, "WEATHER", RaveDataType_DOUBLE);
        PolarScanParam_t *mistnetParamBiology = PolarScan_newParam(scan, "BIOLOGY", RaveDataType_DOUBLE);
        PolarScanParam_t *mistnetParamBackground= PolarScan_newParam(scan, "BACKGROUND", RaveDataType_DOUBLE);
        PolarScanParam_t *mistnetParamClassification= PolarScan_newParam(scan, CELLNAME, RaveDataType_INT);
        
        long nRang = PolarScan_getNbins(scan);
        long nAzim = PolarScan_getNrays(scan);
        double elev = PolarScan_getElangle(scan);
        double rangeScale = PolarScan_getRscale(scan);
        
        for(int iRang=0; iRang<nRang; iRang++){
            for(int iAzim=0; iAzim<nAzim; iAzim++){
                //range in meter
                double range = iRang*rangeScale;
                //azimuth in radials
                double azim = iAzim*2*PI/nAzim;
                // ground distance in meter
                double distance=range2distance(range,elev);
                // Cartesian x coordinate, with radar at center
                double xx=distance*cos(azim);
                // Cartesian y coordinate, with radar at center
                double yy=distance*sin(azim);
                // do not assign values outside the mistnet grid
                if(ABS(xx) > MISTNET_RESOLUTION * (MISTNET_DIMENSION-MISTNET_BLEED)/2) continue;
                if(ABS(yy) > MISTNET_RESOLUTION * (MISTNET_DIMENSION-MISTNET_BLEED)/2) continue;
                // Cartesian grid index x
                int x=MIN(dim3-1,MAX(0,ROUND(xx/res+dim3/2)));
                // Cartesian grid index y
                int y=MIN(dim4-1,MAX(0,ROUND(yy/res+dim4/2)));
                //
                float valueBackground=tensor[MISTNET_BACKGROUND_INDEX][iScan][x][y];
                float valueBiology=tensor[MISTNET_BIOLOGY_INDEX][iScan][x][y];
                float valueWeather=tensor[MISTNET_WEATHER_INDEX][iScan][x][y];
                float valueWeatherAvg=0;
                for(int i=0; i<nScans; i++){
                    valueWeatherAvg+=(tensor[MISTNET_WEATHER_INDEX][i][x][y]/nScans);
                }
                int valueClassification = CELLINIT;
                // post-processing prediction rules for weather, as defined in Lin et al. 2019, doi 10.1111/2041-210X.13280
                if(valueWeather > MISTNET_WEATHER_THRESHOLD || valueWeatherAvg > MISTNET_SCAN_AVERAGE_WEATHER_THRESHOLD){
                    valueClassification=MISTNET_WEATHER_CELL_VALUE;
                }
                PolarScanParam_setValue(mistnetParamBackground, iRang, iAzim, valueBackground);
                PolarScanParam_setValue(mistnetParamBiology, iRang, iAzim, valueBiology);
                PolarScanParam_setValue(mistnetParamWeather, iRang, iAzim, valueWeather);
                PolarScanParam_setValue(mistnetParamClassification, iRang, iAzim, valueClassification);                
            }            
        }
        RAVE_OBJECT_RELEASE(mistnetParamWeather);
        RAVE_OBJECT_RELEASE(mistnetParamBiology);
        RAVE_OBJECT_RELEASE(mistnetParamBackground);
        RAVE_OBJECT_RELEASE(mistnetParamClassification);
        RAVE_OBJECT_RELEASE(scan);
    }
    
    return(0);

}

int addClassificationToPolarVolume(PolarVolume_t* pvol, float ****tensor, int dim1, int dim2, int dim3, int dim4, long res){
    
    RAVE_ASSERT((pvol != NULL), "pvol == NULL");

    PolarScan_t* scan = NULL;

    int nScans;
    // determine how many scan elevations the volume object contains
    nScans = PolarVolume_getNumberOfScans(pvol);
    
   // iterate over the selected scans in 'volume'
    for (int iScan = 0; iScan < nScans; iScan++) {
        // extract the scan object from the volume object
        scan = PolarVolume_getScan(pvol,iScan);

        if(PolarScan_hasParameter(scan, CELLNAME)){
            RAVE_OBJECT_RELEASE(scan);
            continue;
        }

        PolarScanParam_t *mistnetParamClassification= PolarScan_newParam(scan, CELLNAME, RaveDataType_INT);
        
        long nRang = PolarScan_getNbins(scan);
        long nAzim = PolarScan_getNrays(scan);
        double elev = PolarScan_getElangle(scan);
        double rangeScale = PolarScan_getRscale(scan);
        
        for(int iRang=0; iRang<nRang; iRang++){
            for(int iAzim=0; iAzim<nAzim; iAzim++){
                //range in meter
                double range = iRang*rangeScale;
                //azimuth in radials
                double azim = iAzim*2*PI/nAzim;
                // ground distance in meter
                double distance=range2distance(range,elev);
                // Cartesian x coordinate, with radar at center
                double xx=distance*cos(azim);
                // Cartesian y coordinate, with radar at center
                double yy=distance*sin(azim);
                // Cartesian grid index x
                int x=MIN(dim3-1,MAX(0,ROUND(xx/res+dim3/2)));
                // Cartesian grid index y
                int y=MIN(dim4-1,MAX(0,ROUND(yy/res+dim4/2)));
                //
                float valueWeatherAvg=0;
                for(int i=0; i<dim2; i++){
                    valueWeatherAvg+=(tensor[MISTNET_WEATHER_INDEX][i][x][y]/dim2);
                }
                int valueClassification = CELLINIT;
                // post-processing prediction rules for weather, modified for scans not
                // part of the segmentation model, after Lin et al. 2019, doi 10.1111/2041-210X.13280
                if(valueWeatherAvg > MISTNET_SCAN_AVERAGE_WEATHER_THRESHOLD){
                    valueClassification=MISTNET_WEATHER_CELL_VALUE;
                }
                PolarScanParam_setValue(mistnetParamClassification, iRang, iAzim, valueClassification);                
            }            
        }
        
        PolarScan_addParameter(scan, mistnetParamClassification);
        RAVE_OBJECT_RELEASE(mistnetParamClassification);
        RAVE_OBJECT_RELEASE(scan);
    }
    
    return(0);

}

#ifdef MISTNET
// segments biology from precipitation using mistnet deep convolution net.
int segmentScansUsingMistnet(PolarVolume_t* volume, vol2birdScanUse_t *scanUse, vol2bird_t* alldata){    
    // volume with only the 5 selected elevations
    PolarVolume_t* volume_mistnet = NULL;
    PolarVolume_t* volume_select = NULL;
    int result = 0;

    volume_select = PolarVolume_selectScansByScanUse(volume, scanUse, alldata->misc.nScansUsed);
    volume_mistnet = PolarVolume_selectScansByElevation(volume_select, alldata->options.mistNetElevs, alldata->options.mistNetNElevs);
        
    if (PolarVolume_getNumberOfScans(volume_mistnet) != alldata->options.mistNetNElevs){
        vol2bird_err_printf("Error: found only %i/%i scans required by mistnet segmentation model\n",
            PolarVolume_getNumberOfScans(volume_mistnet),alldata->options.mistNetNElevs);
            RAVE_OBJECT_RELEASE(volume_select);
            RAVE_OBJECT_RELEASE(volume_mistnet);
            return -1;
    }

    // set scanUse to false for scans not entering the MistNet segmentation model
    if(alldata->options.mistNetElevsOnly){
        int printWarning = TRUE;
        for(int iScan = 0; iScan < PolarVolume_getNumberOfScans(volume); iScan++){
            PolarScan_t* scan = PolarVolume_getScan(volume, iScan);
            if(PolarVolume_indexOf(volume_mistnet, scan) == -1){
                if(printWarning) vol2bird_err_printf("Warning: Ignoring scan(s) not used as MistNet input: ");
                vol2bird_err_printf( "%i ", iScan + 1);
                printWarning = FALSE;
                scanUse[iScan].useScan = FALSE;
            }
            RAVE_OBJECT_RELEASE(scan);
        }
        if(!printWarning) vol2bird_err_printf( "...\n");
    }

    // convert polar volume into 3D tensor array
    double ***mistnetTensorInput3D = NULL;
    int nCartesianParam = polarVolumeTo3DTensor(volume_mistnet,&mistnetTensorInput3D,MISTNET_DIMENSION,MISTNET_RESOLUTION,3*alldata->options.mistNetNElevs);

    // flatten 3D tensor into a 1D array
    float *mistnetTensorInput;
    mistnetTensorInput = flatten3DTensor(mistnetTensorInput3D,3*alldata->options.mistNetNElevs,MISTNET_DIMENSION,MISTNET_DIMENSION);
    // run mistnet, which outputs a 1D array
    int mistnetTensorSize=3*alldata->options.mistNetNElevs*MISTNET_DIMENSION*MISTNET_DIMENSION;
    float *mistnetTensorOutput = (float *) malloc(mistnetTensorSize*sizeof(float));

    vol2bird_err_printf( "Running MistNet...");

    result = run_mistnet(mistnetTensorInput, &mistnetTensorOutput, alldata->options.mistNetPath, mistnetTensorSize);

    // if mistnet run failed, clean up and exit
    if(result < 0){
        if(nCartesianParam > 0){
            free(mistnetTensorInput);
            free3DTensor(mistnetTensorInput3D,nCartesianParam,MISTNET_RESOLUTION);
        }
        RAVE_OBJECT_RELEASE(volume_select);
        RAVE_OBJECT_RELEASE(volume_mistnet);
        vol2bird_err_printf( "failed\n");
        return -1;
    }

    vol2bird_err_printf( "done\n");
    // convert mistnet 1D array into a 4D tensor
    float ****mistnetTensorOutput4D = create4DTensor(mistnetTensorOutput,3,alldata->options.mistNetNElevs,MISTNET_DIMENSION,MISTNET_DIMENSION);
    // add segmentation to polar volume
    addTensorToPolarVolume(volume_mistnet, mistnetTensorOutput4D,3,alldata->options.mistNetNElevs,MISTNET_DIMENSION,MISTNET_DIMENSION,MISTNET_RESOLUTION);

    // add segmentation for scans that weren't input to the segmentation model to polar volume
    // note: all scans in 'volume_mistnet' are also contained in 'volume', i.e. its scan pointers point to the same objects
    addClassificationToPolarVolume(volume, mistnetTensorOutput4D,3,alldata->options.mistNetNElevs,MISTNET_DIMENSION,MISTNET_DIMENSION,MISTNET_RESOLUTION);

    //clean up 3D array
    if(nCartesianParam > 0){
        free(mistnetTensorInput);
        free(mistnetTensorOutput);
        free3DTensor(mistnetTensorInput3D,nCartesianParam,MISTNET_RESOLUTION);
        free4DTensor(mistnetTensorOutput4D, 3, alldata->options.mistNetNElevs, MISTNET_RESOLUTION);
    }

    RAVE_OBJECT_RELEASE(volume_select);
    RAVE_OBJECT_RELEASE(volume_mistnet);
    
    return result;
}   // segmentScansUsingMistnet
#endif

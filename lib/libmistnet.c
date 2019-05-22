#include "polarvolume.h"
#include "polarscan.h"
#include "cartesian.h"
#include "rave_debug.h"
#include "rave_list.h"
#include "rave_utilities.h"
#include "rave_alloc.h"
#include "constants.h"
#include "libvol2bird.h"


/**
 * Convert from ground distance and elevation to slant range.
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


Cartesian_t* polarVolumeToCartesian(PolarVolume_t* pvol, double elevs[], int nElevs, long dim, long res, double init){
    
    // use only these elevations
    // FIXME make this adjustable
    
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
        fprintf(stderr, "failed to allocate memory for new cartesian object\n");
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
        fprintf(stderr,"Error: polar volume contains no scans\n");
        return NULL;
    }

    // get the number of elevations.
    if(nElevs>nScans){
        fprintf(stderr,"Warning: requesting %i elevations scans, but only %i available\n", nElevs, nScans);
    }

    // iterate over the selected scans in 'volume'
    for (int iElev = 0; iElev < nElevs; iElev++) {
        
        // extract the scan object from the volume object
        scan = PolarVolume_getScanClosestToElevation(pvol,DEG2RAD*elevs[iElev],0);
        
        double elev = PolarScan_getElangle(scan);
        
        scanParameterNames = PolarScan_getParameterNames(scan);
        
        if(RaveList_size(scanParameterNames)<=0){
            fprintf(stderr,"Warning: ignoring scan without scan parameters\n");
            continue;            
        }
        
        for(int iParam = 0; iParam<RaveList_size(scanParameterNames); iParam++){
            // retrieve name of the scan parameter
            scanParameterName = RaveList_get(scanParameterNames, iParam);
            
            char iElevString[4];
            sprintf(iElevString, "%d", iElev);

            // create a new scan parameter name with index for the sweep
            char *parameterNameFull = malloc(strlen(scanParameterName)+strlen(iElevString)+1);
            strcpy(parameterNameFull,scanParameterName);
            strcat(parameterNameFull,iElevString);
            parameterName = RaveUtilities_trimText(parameterNameFull, strlen(parameterNameFull));

            // create a cartesian scan parameter with the same name
            cartesianParam = Cartesian_createParameter(cartesian,parameterName,RaveDataType_DOUBLE, init);
            
            double range,azim,distance;
            double value;
            
            RaveValueType a;
            // loop over the grid, and fill it
            for(long x = 0,xx; x<dim; x++){
                for(long y = 0,yy; y<dim; y++){
                    xx=res*(x-dim/2);
                    yy=res*(y-dim/2);
                    azim=atan2(yy,xx);
                    distance=sqrt(SQUARE(xx)+SQUARE(yy));
                    range=distance2range(distance,elev);
                    a=PolarScan_getConvertedParameterValueAtAzimuthAndRange(scan,scanParameterName,azim,range,&value);
                    if(a==RaveValueType_DATA){
                        CartesianParam_setValue(cartesianParam, x, y, value);
                        //fprintf(stderr,"no data at scan=%i,param=%s,azim=%f,range=%f,distance=%f,x=%li,y=%li,elev=%f\n",iElev,scanParameterName,azim,range,distance,xx,yy,elev*RAD2DEG);
                    }
                    else{
                        CartesianParam_setValue(cartesianParam, x, y, NAN);
                    }
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
};

double*** init3DArray(int dim1, int dim2, int dim3, double init){
    
    double ***array = (double ***)malloc(dim1*sizeof(double**));
    
    if(array == NULL){
        fprintf(stderr,"failed to initialize 3D array (1)");
        exit(0);
    }
    
    for (int i = 0; i< dim1; i++) {

        array[i] = (double **) malloc(dim2*sizeof(double *));
        
        if(array[i] == NULL){
            fprintf(stderr,"failed to initialize 3D array (2)");
            exit(0);
        }
        
        for (int j = 0; j < dim2; j++) {
            
            array[i][j] = (double *)malloc(dim3*sizeof(double));
            
            if(array[i][j] == NULL){
                fprintf(stderr,"failed to initialize 3D array (3)");
                exit(0);
            }
        } // j
    } // i
    
    	// assign values to allocated memory
	for (int i = 0; i < dim1; i++)
		for (int j = 0; j < dim2; j++)
			for (int k = 0; k < dim3; k++)
				array[i][j][k] = init;

    return array;
}

 
void free3DArray(double ***array, int dim1, int dim2){
	// deallocate memory
	for (int i = 0; i < dim1; i++) 
	{
		for (int j = 0; j < dim2; j++)
			free(array[i][j]);

		free(array[i]);
	}
	free(array);
}

int fill3DArray(double ***array, Cartesian_t* cartesian, int dim1, int dim2, int dim3, int nScanParam){
    int nCartesianParam;
    
    nCartesianParam = Cartesian_getParameterCount(cartesian);
    long xSize	= Cartesian_getXSize(cartesian);
    long ySize	= Cartesian_getYSize(cartesian);
    
    if(dim1*nScanParam != nCartesianParam){
        fprintf(stderr, "Error: expecting %i cartesian scan parameters, but found %i\n", dim1*nScanParam, nCartesianParam);
        return -1;
    }

    if(dim2 != xSize){
        fprintf(stderr, "Error: expecting %i bins in X dimension, but found only %li\n", dim2, xSize);
        return -1;
    }

    if(dim3 != ySize){
        fprintf(stderr, "Error: expecting %i bins in Y dimension, but found only %li\n", dim3, ySize);
        return -1;
    }
    
    RaveList_t* cartesianParameterNames = Cartesian_getParameterNames(cartesian);
    double value;
    RaveValueType valueType;
    
    for(int iParam = 0; iParam < nCartesianParam; iParam++){
        char* parameterName = RaveList_get(cartesianParameterNames, iParam);
        fprintf(stderr, "writing %s to 3D array\n", parameterName);
        CartesianParam_t* cartesianParam = Cartesian_getParameter(cartesian, parameterName);
        
        for(int x = 0; x < xSize; x++){
            for(int y = 0; y < ySize; y++){
                valueType = CartesianParam_getValue(cartesianParam, x, y, &value);
                if (valueType == RaveValueType_DATA){
                    array[iParam][x][y] = value;
                }
                else{
                    array[iParam][x][y] = NAN;
                }
            } //y
        } //x
        
        RAVE_OBJECT_RELEASE(cartesianParam);
    } //iScan
    
    return 0;
}

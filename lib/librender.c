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
#ifdef MISTNET
#include "../libmistnet/libmistnet.h"
#endif

/**
 * FUNCTION PROTOTYPES
 **/
double distance2range(double distance,double elev);

double distance2height(double distance,double elev);

Cartesian_t* polarVolumeToCartesian(PolarVolume_t* pvol, float elevs[], int nElevs, long dim, long res, double init);

RaveObjectList_t* polarVolumeToCartesianList(PolarVolume_t* pvol, float elevs[], int nElevs, long dim, long res, double init, int *nParam);

Cartesian_t* polarScanToCartesian(PolarScan_t* scan, long dim, long res, double init);

void free4DTensor(float ****tensor, int dim1, int dim2, int dim3);

float**** create4DTensor(float *array, int dim1, int dim2, int dim3, int dim4);

double*** init3DTensor(int dim1, int dim2, int dim3, double init);

void free3DTensor(double ***tensor, int dim1, int dim2);

int fill3DTensor(double ***tensor, RaveObjectList_t* list, int dim1, int dim2, int dim3);

float* flatten3DTensor(double ***tensor, int dim1, int dim2, int dim3);

int polarVolumeTo3DTensor(PolarVolume_t* pvol, double ****tensor, float elevs[], int nElevs, int dim, double res, int nParam);

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


Cartesian_t* polarVolumeToCartesian(PolarVolume_t* pvol, float elevs[], int nElevs, long dim, long res, double init){
    
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
            
            char iElevString[11];
            // copy iElev to iElevString
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




RaveObjectList_t* polarVolumeToCartesianList(PolarVolume_t* pvol, float elevs[], int nElevs, long dim, long res, double init, int *nParam){
    
    PolarScan_t *scan = NULL;
    RaveObjectList_t* list;
    Cartesian_t *cartesian = NULL;
    
    // create a new list
    list = (RaveObjectList_t *) RAVE_OBJECT_NEW(&RaveObjectList_TYPE);
    
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
                
        cartesian = polarScanToCartesian(scan, dim, res, init);
        
        *nParam += Cartesian_getParameterCount(cartesian);
        
        RaveObjectList_add(list, (RaveCoreObject*) cartesian); 
        
        RAVE_OBJECT_RELEASE(cartesian);

    }
    
    return list;

}



Cartesian_t* polarScanToCartesian(PolarScan_t* scan, long dim, long res, double init){
    
    RAVE_ASSERT((scan != NULL), "scan == NULL");
        
    // initialize scan, param, and cartesian RAVE objects
    Cartesian_t *cartesian = NULL;
    CartesianParam_t *cartesianParam = NULL;
    RaveList_t* scanParameterNames;
    char* scanParameterName;
    
    // create a new Cartesian grid object
    cartesian = RAVE_OBJECT_NEW(&Cartesian_TYPE);

    // copy metadata from volume
    Cartesian_setTime(cartesian, PolarScan_getTime(scan));
    Cartesian_setDate(cartesian, PolarScan_getDate(scan));
    Cartesian_setSource(cartesian, PolarScan_getSource(scan));

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

    
    double elev = PolarScan_getElangle(scan);
    
    scanParameterNames = PolarScan_getParameterNames(scan);
    
    if(RaveList_size(scanParameterNames)<=0){
        fprintf(stderr,"Warning: scan without scan parameters\n");
        return NULL;
    }
            
    for(int iParam = 0; iParam<RaveList_size(scanParameterNames); iParam++){
        // retrieve name of the scan parameter
        scanParameterName = RaveList_get(scanParameterNames, iParam);
        
        // create a cartesian scan parameter with the same name
        cartesianParam = Cartesian_createParameter(cartesian,scanParameterName,RaveDataType_DOUBLE, init);
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
                }
                else{
                    CartesianParam_setValue(cartesianParam, x, y, NAN);
                }
            }
        }
        
        // add the cartesian scan parameter to the cartesian object
        Cartesian_addParameter(cartesian,cartesianParam);
        
        RAVE_OBJECT_RELEASE(cartesianParam);
    } // iParam
    
    return cartesian;
};


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

    

double*** init3DTensor(int dim1, int dim2, int dim3, double init){
    
    double ***tensor = (double ***)malloc(dim1*sizeof(double**));
    
    if(tensor == NULL){
        fprintf(stderr,"failed to initialize 3D tensor (1)");
        exit(0);
    }
    
    for (int i = 0; i< dim1; i++) {

        tensor[i] = (double **) malloc(dim2*sizeof(double *));
        
        if(tensor[i] == NULL){
            fprintf(stderr,"failed to initialize 3D tensor (2)");
            exit(0);
        }
        
        for (int j = 0; j < dim2; j++) {
            
            tensor[i][j] = (double *)malloc(dim3*sizeof(double));
            
            if(tensor[i][j] == NULL){
                fprintf(stderr,"failed to initialize 3D tensor (3)");
                exit(0);
            }
        } // j
    } // i
    
    	// assign values to allocated memory
	for (int i = 0; i < dim1; i++)
		for (int j = 0; j < dim2; j++)
			for (int k = 0; k < dim3; k++)
				tensor[i][j][k] = init;

    return tensor;
}

void free3DTensor(double ***tensor, int dim1, int dim2){
	// deallocate memory
	for (int i = 0; i < dim1; i++) 
	{
		for (int j = 0; j < dim2; j++){
            free(tensor[i][j]);
        }
		free(tensor[i]);
	}
	free(tensor);
}



int fill3DTensor(double ***tensor, RaveObjectList_t* list, int dim1, int dim2, int dim3){

    int nScan = RaveObjectList_size(list);
    int iParam = 0;
    
    for(int iScan = 0; iScan<nScan; iScan++){
        
        Cartesian_t* cartesian = (Cartesian_t*) RaveObjectList_get(list, iScan);
        
        int nCartesianParam = Cartesian_getParameterCount(cartesian);
        long xSize	= Cartesian_getXSize(cartesian);
        long ySize	= Cartesian_getYSize(cartesian);
        
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

                fprintf(stderr,"writing %s at index %i\n",parameterName,iParam);
                
                if(iParam>=dim1){
                   fprintf(stderr, "Error: exceeding 3D tensor dimension\n"); 
                   RAVE_OBJECT_RELEASE(cartesianParam);
                   return(-1);
                }
                
                // fill tensor
                for(int x = 0; x < xSize; x++){
                    for(int y = 0; y < ySize; y++){
                        valueType = CartesianParam_getValue(cartesianParam, x, y, &value);
                        if (valueType == RaveValueType_DATA){
                            tensor[iParam][x][y] = value;
                        }
                        else{
                            tensor[iParam][x][y] = NAN;
                        }
                    } //y
                } //x
                
                RAVE_OBJECT_RELEASE(cartesianParam);
                
                //increase counter
                iParam++;
                                
            } //iScan
        }   
    }
    
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


int polarVolumeTo3DTensor(PolarVolume_t* pvol, double ****tensor, float elevs[], int nElevs, int dim, double res, int nParam){
    //Un-comment these two lines to save a rendering to file
    //Cartesian_t *cartesian = NULL;
    //cartesian = polarVolumeToCartesian(pvol, elevs, nElevs, dim, res, 0);            
    //saveToODIM((RaveCoreObject*) cartesian, "rendering.h5");
    
    // convert polar volume to a list of Cartesian objects, one for each scan
    // store the total number of scan parameters for all scans in nCartesianParam
    int nCartesianParam = 0;
    RaveObjectList_t* list = polarVolumeToCartesianList(pvol, elevs, nElevs, dim, res, 0, &nCartesianParam);
    
    if(list == NULL){
        fprintf(stderr, "Error: failed to load Cartesian objects from polar volume\n");
        return -1;
    }
    
    // if nParam is specified, restrict the number of output parameters to its value.
    // nParam typically equals 3, selecting DBZ, VRAD and WRAD for Misnet segmentation model input.
    if(nParam > 0){
        if(nParam < nCartesianParam){
            nCartesianParam = nParam;
        }
    }

    // initialize a 3D tensor, and fill it
    *tensor = init3DTensor(nCartesianParam,dim,dim,0);
    fill3DTensor(*tensor, list, nCartesianParam, dim, dim);

    // clean up
    RAVE_OBJECT_RELEASE(list);
    
    fprintf(stderr,"DONE\n");
    
    return(nCartesianParam);
}

#ifdef MISTNET

// segments biology from precipitation using mistnet deep convolution net.
int segmentScansUsingMistnet(PolarVolume_t* volume, vol2bird_t* alldata){
    
    if (PolarVolume_getNumberOfScans(volume) != alldata->options.cartesianNElevs){
        fprintf(stderr,"Error: polar volume has %i scans but segmentation model expects %i scans",
            PolarVolume_getNumberOfScans(volume),alldata->options.cartesianNElevs);
        return -1;
    }

    // convert polar volume into 3D tensor array
    double ***mistnetTensorInput3D = NULL;
    fprintf(stderr, "convert pvol to 3D tensor...\n");
    int nCartesianParam = polarVolumeTo3DTensor(volume,&mistnetTensorInput3D,alldata->options.cartesianElevs,alldata->options.cartesianNElevs,MISTNET_DIMENSION,MISTNET_RESOLUTION,3*alldata->options.cartesianNElevs);
    // flatten 3D tensor into a 1D array
    float *mistnetTensorInput;
    fprintf(stderr, "flatten 3D tensor...\n");
    mistnetTensorInput = flatten3DTensor(mistnetTensorInput3D,3*alldata->options.cartesianNElevs,MISTNET_DIMENSION,MISTNET_DIMENSION);
    // run mistnet, which outputs a 1D array
    int mistnetTensorSize=3*alldata->options.cartesianNElevs*MISTNET_DIMENSION*MISTNET_DIMENSION;
    float *mistnetTensorOutput = (float *) malloc(mistnetTensorSize*sizeof(float));
    //float**** mistnetTensorOutput4D = create4DTensor(3,alldata->options.cartesianNElevs,MISTNET_DIMENSION,MISTNET_DIMENSION);
    fprintf(stderr, "START MISTNET...");
    run_mistnet(mistnetTensorInput, &mistnetTensorOutput, MISTNET_PATH, mistnetTensorSize);
    fprintf(stderr, "done\n");
    // convert mistnet 1D array into a 4D tensor
    float ****mistnetTensorOutput4D = create4DTensor(mistnetTensorOutput,3,alldata->options.cartesianNElevs,MISTNET_DIMENSION,MISTNET_DIMENSION);
    // add segmentation to polar volume
    int result = 0;
    //XXX TODO 
    // 1) map tensor to cartesian object (not essential), or a list of cartesian objects.
    // 2) map tensor to polar volume
    
    //clean up 3D array
    if(nCartesianParam > 0){
        fprintf(stderr,"DONE WITH MISTNET, cleaning up %i params\n",nCartesianParam);
        free(mistnetTensorInput);
        free(mistnetTensorOutput);
        free3DTensor(mistnetTensorInput3D,nCartesianParam,MISTNET_RESOLUTION);
        free4DTensor(mistnetTensorOutput4D, 3, alldata->options.cartesianNElevs, MISTNET_RESOLUTION);
    }
    
    return result;
}   // segmentScansUsingMistnet

#endif

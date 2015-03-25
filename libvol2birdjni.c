



#include <jni.h>
#include "libvol2bird.c"
#include "libsvdfit.h"


JNIEXPORT jint JNICALL
Java_nl_esciencecenter_ncradar_JNIMethodsVol2Bird_analyzeCells(
JNIEnv *env,
jobject obj,
jintArray dbzImageInt,
jintArray vradImageInt,
jintArray texImageInt,
jintArray clutterImageInt,
jintArray cellImageInt,
jint dbznRang,
jint dbznAzim,
jfloat dbzElev,
jfloat dbzValueScale,
jfloat dbzValueOffset,
jfloat vradValueScale,
jfloat vradValueOffset,
jfloat clutterValueScale,
jfloat clutterValueOffset,
jfloat texValueScale,
jfloat texValueOffset,
jint nCells,
jint areaMin,
jfloat cellDbzMin,
jfloat cellStdDevMax,
jfloat cellClutterFraction,
jfloat absVradMin,
jfloat clutterValueMax,
jint cmFlagInt,
jint verboseInt
)
{

    // do some Java Native interface tricks:
    jint *dbzImageIntBody = (*env)->GetIntArrayElements(env, dbzImageInt, NULL);
    jint *vradImageIntBody = (*env)->GetIntArrayElements(env, vradImageInt, NULL);
    jint *texImageIntBody = (*env)->GetIntArrayElements(env, texImageInt, NULL);
    jint *clutterImageIntBody = (*env)->GetIntArrayElements(env, clutterImageInt, NULL);
    jint *cellImageIntBody = (*env)->GetIntArrayElements(env, cellImageInt, NULL);
    jsize nGlobal = (*env)->GetArrayLength(env, dbzImageInt);
    // end of Java Native Interface tricks



    unsigned char dbzImageBody[nGlobal];
    unsigned char vradImageBody[nGlobal];
    unsigned char texImageBody[nGlobal];
    unsigned char clutterImageBody[nGlobal];
    unsigned char cmFlag;
    unsigned char verbose;


    int nAzim = dbznAzim;
    int nRang = dbznRang;
    int iAzim;
    int iRang;
    int iGlobal;
    int nCellsValid;

    for (iAzim = 0;iAzim<nAzim;iAzim++){
        for (iRang= 0 ; iRang<nRang;iRang++){

            iGlobal = iAzim*nRang + iRang;

            if (0<=dbzImageIntBody[iGlobal] && dbzImageIntBody[iGlobal]<=255) {
                dbzImageBody[iGlobal] = (unsigned char) dbzImageIntBody[iGlobal];
            }
            else {
                fprintf(stderr,"Error converting type (dbzImageIntBody[iGlobal]).\n");
                return -1;
            }

            if (0<=vradImageIntBody[iGlobal] && vradImageIntBody[iGlobal]<=255) {
                vradImageBody[iGlobal] = (unsigned char) vradImageIntBody[iGlobal];
            }
            else {
                fprintf(stderr,"Error converting type (vradImageIntBody[iGlobal]).\n");
                return -1;
            }

            if (0<=texImageIntBody[iGlobal] && texImageIntBody[iGlobal]<=255) {
                texImageBody[iGlobal] = (unsigned char) texImageIntBody[iGlobal];
            }
            else {
                fprintf(stderr,"Error converting type (texImageIntBody[iGlobal]).\n");
                return -1;
            }

            if (0<=clutterImageIntBody[iGlobal] && clutterImageIntBody[iGlobal]<=255) {
                clutterImageBody[iGlobal] = (unsigned char) clutterImageIntBody[iGlobal];
            }
            else {
                fprintf(stderr,"Error converting type (clutterImageIntBody[iGlobal]).\n");
                return -1;
            }

        }
    }


    SCANMETA dbzMeta;
    SCANMETA vradMeta;
    SCANMETA texMeta;
    SCANMETA clutterMeta;

    dbzMeta.nRang = dbznRang;
    dbzMeta.nAzim = dbznAzim;
    dbzMeta.elev = dbzElev;

    dbzMeta.valueScale = dbzValueScale;
    dbzMeta.valueOffset = dbzValueOffset;

    vradMeta.valueScale = vradValueScale;
    vradMeta.valueOffset = vradValueOffset;

    clutterMeta.valueScale = clutterValueScale;
    clutterMeta.valueOffset = clutterValueOffset;

    texMeta.valueScale = texValueScale;
    texMeta.valueOffset = texValueOffset;


    // cast to unsigned char
    if (0<=cmFlagInt && cmFlagInt<=255) {
        cmFlag = (unsigned char) cmFlagInt;
    }
    else {
        fprintf(stderr,"Error converting type.");
        return -1;
    }


    // cast to unsigned char
    if (0<=verboseInt && verboseInt<=255) {
        verbose = (unsigned char) verboseInt;
    }
    else {
        fprintf(stderr,"Error converting type.");
        return -1;
    }



    nCellsValid = analyzeCells(&dbzImageBody[0], &vradImageBody[0], &texImageBody[0],
                               &clutterImageBody[0], &cellImageIntBody[0],
                               &dbzMeta, &vradMeta, &texMeta, &clutterMeta,
                               nCells, areaMin, cellDbzMin, cellStdDevMax, cellClutterFraction,
                               absVradMin, clutterValueMax, cmFlag,
                               verbose);


    #ifdef FPRINTFON
    int minValue = cellImageIntBody[0];
    int maxValue = cellImageIntBody[0];
    for (iGlobal = 0;iGlobal < nGlobal;iGlobal++) {
        if (cellImageIntBody[iGlobal] < minValue) {
            minValue = cellImageIntBody[iGlobal];
        }
        if (cellImageIntBody[iGlobal] > maxValue) {
            maxValue = cellImageIntBody[iGlobal];
        }
    }
    fprintf(stderr,"minimum value in cellImageIntBody array = %d.\n", minValue);
    fprintf(stderr,"maximum value in cellImageIntBody array = %d.\n", maxValue);
    #endif



    #ifdef FPRINTFON
    fprintf(stderr,"DEBUG libvol2birdjni.c before casting to jint\n");
    #endif


    // cast to the right type:
    for (iAzim = 0; iAzim < nAzim; iAzim++){
        for (iRang = 0; iRang < nRang; iRang++){

            iGlobal = iAzim*nRang + iRang;

            dbzImageIntBody[iGlobal] = (jint) dbzImageBody[iGlobal];
            vradImageIntBody[iGlobal] = (jint) vradImageBody[iGlobal];
            texImageIntBody[iGlobal] = (jint) texImageBody[iGlobal];
            clutterImageIntBody[iGlobal] = (jint) clutterImageBody[iGlobal];

        }
    }


    // do some Java Native interface tricks:
    (*env)->ReleaseIntArrayElements(env, dbzImageInt, dbzImageIntBody, JNI_ABORT);  // FIXME maybe don't use ABORT?
    (*env)->ReleaseIntArrayElements(env, vradImageInt, vradImageIntBody, JNI_ABORT);  // FIXME maybe don't use ABORT?
    (*env)->ReleaseIntArrayElements(env, texImageInt, texImageIntBody, JNI_ABORT);  // FIXME maybe don't use ABORT?
    (*env)->ReleaseIntArrayElements(env, clutterImageInt, clutterImageIntBody, JNI_ABORT);  // FIXME maybe don't use ABORT?
    (*env)->ReleaseIntArrayElements(env, cellImageInt, cellImageIntBody, 0);
    // end of Java Native Interface tricks


    #ifdef FPRINTFON
    fprintf(stderr,"DEBUG libvol2birdjni.c end of analyzeCells\n");
    #endif


    return nCellsValid;



}






JNIEXPORT jfloat JNICALL
Java_nl_esciencecenter_ncradar_JNIMethodsVol2Bird_calcDist(
JNIEnv *env,
jobject obj,
const jint range1,
const jint azim1,
const jint range2,
const jint azim2,
const jfloat rscale,
const jfloat ascale)
{

    return calcDist(range1, azim1, range2, azim2, rscale, ascale);

}




JNIEXPORT void JNICALL
Java_nl_esciencecenter_ncradar_JNIMethodsVol2Bird_calcTexture(
        JNIEnv *env,
        jobject obj,
        jintArray texImageInt,
        const jintArray dbzImageInt,
        const jintArray vradImageInt,
        const jint nRangNeighborhood,
        const jint nAzimNeighborhood,
        const jint nCountMin,
        const jfloat texOffset,
        const jfloat texScale,
        const jfloat dbzOffset,
        const jfloat dbzScale,
        const jfloat vradOffset,
        const jfloat vradScale,
        const jint vradMissing,
        const jint nRang,
        const jint nAzim) {

    // do some Java Native interface tricks:
    jint *texImageIntBody = (*env)->GetIntArrayElements(env, texImageInt, NULL);
    jint *dbzImageIntBody = (*env)->GetIntArrayElements(env, dbzImageInt, NULL);
    jint *vradImageIntBody = (*env)->GetIntArrayElements(env, vradImageInt, NULL);
    jsize nGlobal = (*env)->GetArrayLength(env, texImageInt);
    // end of Java Native Interface tricks

    unsigned char texImageUCharBody[nGlobal];
    unsigned char dbzImageUCharBody[nGlobal];
    unsigned char vradImageUCharBody[nGlobal];

    int iGlobal;

    for (iGlobal = 0; iGlobal < nGlobal; iGlobal++){

        if (0<=texImageIntBody[iGlobal] && texImageIntBody[iGlobal]<=255) {
            texImageUCharBody[iGlobal] = (unsigned char) texImageIntBody[iGlobal];
        }
        else {
            fprintf(stderr,"Error converting type (texImageIntBody[iGlobal]).\n");
            return;
        }

        if (0<=dbzImageIntBody[iGlobal] && dbzImageIntBody[iGlobal]<=255) {
            dbzImageUCharBody[iGlobal] = (unsigned char) dbzImageIntBody[iGlobal];
        }
        else {
            fprintf(stderr,"Error converting type (dbzImageIntBody[iGlobal]).\n");
            return;
        }

        if (0<=vradImageIntBody[iGlobal] && vradImageIntBody[iGlobal]<=255) {
            vradImageUCharBody[iGlobal] = (unsigned char) vradImageIntBody[iGlobal];
        }
        else {
            fprintf(stderr,"Error converting type (vradImageIntBody[iGlobal]).\n");
            return;
        }
    }

    SCANMETA texMeta;
    SCANMETA dbzMeta;
    SCANMETA vradMeta;

    texMeta.valueOffset = texOffset;
    texMeta.valueScale = texScale;

    dbzMeta.valueOffset = dbzOffset;
    dbzMeta.valueScale = dbzScale;

    vradMeta.valueOffset = vradOffset;
    vradMeta.valueScale = vradScale;
    vradMeta.nRang = nRang;
    vradMeta.nAzim = nAzim;
    vradMeta.missing = vradMissing;



    calcTexture(texImageUCharBody, vradImageUCharBody, dbzImageUCharBody,
            &texMeta, &vradMeta, &dbzMeta,
            nRangNeighborhood, nAzimNeighborhood, nCountMin);



    // cast back to integer type:
    for (iGlobal = 0; iGlobal < nGlobal; iGlobal++) {
        texImageIntBody[iGlobal] = (jint) texImageUCharBody[iGlobal];
        vradImageIntBody[iGlobal] = (jint) vradImageUCharBody[iGlobal];
        dbzImageIntBody[iGlobal] = (jint) dbzImageUCharBody[iGlobal];
    }



    // do some Java Native interface tricks:
    (*env)->ReleaseIntArrayElements(env, texImageInt, texImageIntBody, 0);
    (*env)->ReleaseIntArrayElements(env, dbzImageInt, dbzImageIntBody, JNI_ABORT);
    (*env)->ReleaseIntArrayElements(env, vradImageInt, vradImageIntBody, JNI_ABORT);
    // end of Java Native Interface tricks

}




JNIEXPORT void JNICALL
Java_nl_esciencecenter_ncradar_JNIMethodsVol2Bird_classifyGates(
JNIEnv *env,
jobject obj,
const jint dbznRang,
const jint dbznAzim,
const jfloat dbzRangeScale,
const jfloat dbzElev,
const jfloat dbzHeig,
const jfloat dbzValueScale,
const jfloat dbzValueOffset,
const jfloat dbzAzimScale,
const jint dbzMissing,
const jfloat vradValueScale,
const jfloat vradValueOffset,
const jint vradMissing,
const jint rawReflMissing,
const jfloat clutterValueScale,
const jfloat clutterValueOffset,
jintArray cellImage,
const jintArray dbzImageInt,
const jintArray vradImageInt,
const jintArray rawReflImageInt,
const jintArray clutterImageInt,
jfloatArray zdata,
jintArray nzdata,
const jfloat rangeMin,
const jfloat rangeMax,
const jfloat HLAYER,
const jfloat XOFFSET,
const jfloat XSCALE,
const jfloat XMEAN,
const jfloat heightOfInterest,
const jfloat azimMin,
const jfloat azimMax,
const jfloat absVradMin,
const jfloat dbzClutter,
const jfloat dbzMin,
const jfloat dBZx,
const jfloat DBZNOISE,
const jint iLayer,
const jint clutterFlagInt,
const jint rawReflFlagInt,
const jint xflagInt
)
{


//    fprintf(stderr, "dbznRang = %d\n", dbznRang);
//    fprintf(stderr, "dbznAzim = %d\n", dbznAzim);
//    fprintf(stderr, "dbzRangeScale = %f\n", dbzRangeScale);
//    fprintf(stderr, "dbzElev = %f\n", dbzElev);
//    fprintf(stderr, "dbzHeig = %f\n", dbzHeig);
//    fprintf(stderr, "dbzValueScale = %f\n", dbzValueScale);
//    fprintf(stderr, "dbzValueOffset = %f\n", dbzValueOffset);
//    fprintf(stderr, "dbzAzimScale = %f\n", dbzAzimScale);
//    fprintf(stderr, "dbzMissing = %d\n", dbzMissing);
//    fprintf(stderr, "vradValueScale = %f\n", vradValueScale);
//    fprintf(stderr, "vradValueOffset = %f\n", vradValueOffset);
//    fprintf(stderr, "vradMissing = %d\n", vradMissing);
//    fprintf(stderr, "rawReflMissing = %d\n", rawReflMissing);
//    fprintf(stderr, "clutterValueScale = %f\n", clutterValueScale);
//    fprintf(stderr, "clutterValueOffset = %f\n", clutterValueOffset);
//    fprintf(stderr, "cellImage = %p\n", cellImage);


    // do some Java Native interface tricks:
    jint *cellImageBody = (*env)->GetIntArrayElements(env, cellImage, NULL);
    jint *dbzImageIntBody = (*env)->GetIntArrayElements(env, dbzImageInt, NULL);
    jint *vradImageIntBody = (*env)->GetIntArrayElements(env, vradImageInt, NULL);
    jint *rawReflImageIntBody = (*env)->GetIntArrayElements(env, rawReflImageInt, NULL);
    jint *clutterImageIntBody = (*env)->GetIntArrayElements(env, clutterImageInt, NULL);
    jfloat *zdataBody = (*env)->GetFloatArrayElements(env, zdata, NULL);
    jint *nzdataBody = (*env)->GetIntArrayElements(env, nzdata, NULL);
    jsize nGlobal = (*env)->GetArrayLength(env, cellImage);
    // end of Java Native Interface tricks


    unsigned char dbzImageBody[nGlobal];
    unsigned char vradImageBody[nGlobal];
    unsigned char rawReflImageBody[nGlobal];
    unsigned char clutterImageBody[nGlobal];
    unsigned char clutterFlag;
    unsigned char rawReflFlag;
    unsigned char xflag;

    int iGlobal;


    for (iGlobal = 0; iGlobal < nGlobal; iGlobal++){

            if (0<=dbzImageIntBody[iGlobal] && dbzImageIntBody[iGlobal]<=255) {
                dbzImageBody[iGlobal] = (unsigned char) dbzImageIntBody[iGlobal];
            }
            else {
                fprintf(stderr,"Error converting type (dbzImageIntBody[iGlobal]).\n");
                return;
            }

            if (0<=vradImageIntBody[iGlobal] && vradImageIntBody[iGlobal]<=255) {
                vradImageBody[iGlobal] = (unsigned char) vradImageIntBody[iGlobal];
            }
            else {
                fprintf(stderr,"Error converting type (vradImageIntBody[iGlobal]).\n");
                return;
            }

            if (0<=rawReflImageIntBody[iGlobal] && rawReflImageIntBody[iGlobal]<=255) {
                rawReflImageBody[iGlobal] = (unsigned char) rawReflImageIntBody[iGlobal];
            }
            else {
                fprintf(stderr,"Error converting type (rawReflImageIntBody[iGlobal]).\n");
                return;
            }

            if (0<=clutterImageIntBody[iGlobal] && clutterImageIntBody[iGlobal]<=255) {
                clutterImageBody[iGlobal] = (unsigned char) clutterImageIntBody[iGlobal];
            }
            else {
                fprintf(stderr,"Error converting type (clutterImageIntBody[iGlobal]).\n");
                return;
            }

    }

    SCANMETA dbzMeta;
    SCANMETA vradMeta;
    SCANMETA rawReflMeta;
    SCANMETA clutterMeta;

    dbzMeta.nRang = dbznRang;
    dbzMeta.nAzim = dbznAzim;
    dbzMeta.rangeScale = dbzRangeScale;
    dbzMeta.elev = dbzElev;
    dbzMeta.heig = dbzHeig;
    dbzMeta.valueScale = dbzValueScale;
    dbzMeta.valueOffset = dbzValueOffset;
    dbzMeta.azimScale = dbzAzimScale;
    dbzMeta.missing = dbzMissing;

    vradMeta.valueScale = vradValueScale;
    vradMeta.valueOffset = vradValueOffset;
    vradMeta.missing = vradMissing;

    rawReflMeta.missing = rawReflMissing;

    clutterMeta.valueScale = clutterValueScale;
    clutterMeta.valueOffset = clutterValueOffset;



    // cast to unsigned char
    if (0<=clutterFlagInt && clutterFlagInt<=255) {
        clutterFlag = (unsigned char) clutterFlagInt;
    }
    else {
        fprintf(stderr,"Error converting type.\n");
        return;
    }

    // cast to unsigned char
    if (0<=rawReflFlagInt && rawReflFlagInt<=255) {
        rawReflFlag = (unsigned char) rawReflFlagInt;
    }
    else {
        fprintf(stderr,"Error converting type.\n");
        return;
    }

    // cast to unsigned char
    if (0<=xflagInt && xflagInt<=255) {
        xflag = (unsigned char) xflagInt;
    }
    else {
        fprintf(stderr,"Error converting type.\n");
        return;
    }


    classifyGates(dbzMeta, vradMeta, rawReflMeta,
            clutterMeta, cellImageBody, &dbzImageBody[0], &vradImageBody[0],
            &rawReflImageBody[0], &clutterImageBody[0],
            &zdataBody[0], &nzdataBody[0],
            rangeMin, rangeMax, HLAYER, XOFFSET,
            XSCALE, XMEAN, heightOfInterest,
            azimMin, azimMax, absVradMin, dbzClutter, dbzMin,
            dBZx, DBZNOISE,
            iLayer, clutterFlag, rawReflFlag, xflag);

    // cast back to integer type:
    for (iGlobal = 0; iGlobal < nGlobal; iGlobal++) {
        dbzImageIntBody[iGlobal] = (jint) dbzImageBody[iGlobal];
        vradImageIntBody[iGlobal] = (jint) vradImageBody[iGlobal];
        rawReflImageIntBody[iGlobal] = (jint) rawReflImageBody[iGlobal];
        clutterImageIntBody[iGlobal] = (jint) clutterImageBody[iGlobal];
    }


    // do some Java Native interface tricks:
    (*env)->ReleaseIntArrayElements(env, cellImage, cellImageBody, JNI_ABORT);
    (*env)->ReleaseIntArrayElements(env, dbzImageInt, dbzImageIntBody, JNI_ABORT);
    (*env)->ReleaseIntArrayElements(env, vradImageInt, vradImageIntBody, JNI_ABORT);
    (*env)->ReleaseIntArrayElements(env, rawReflImageInt, rawReflImageIntBody, JNI_ABORT);
    (*env)->ReleaseIntArrayElements(env, clutterImageInt, clutterImageIntBody, JNI_ABORT);
    (*env)->ReleaseFloatArrayElements(env, zdata, zdataBody, JNI_ABORT);
    (*env)->ReleaseIntArrayElements(env, nzdata, nzdataBody, JNI_ABORT);
    // end of Java Native Interface tricks

    return;

}



JNIEXPORT jint JNICALL
Java_nl_esciencecenter_ncradar_JNIMethodsVol2Bird_detNumberOfGates(
        JNIEnv *env,
        jobject obj,
        const jint iLayer,
        const jfloat layerThickness,
        const jfloat rangeMin,
        const jfloat rangeMax,
        const jfloat rangeScale,
        const jfloat elevAngle,
        const jint nRang,
        const jint nAzim,
        const jfloat radarHeight)
{

    int nRecordsMax;

    nRecordsMax = detNumberOfGates(iLayer, layerThickness,
                              rangeMin, rangeMax,
                              rangeScale, elevAngle,
                              nRang, nAzim,
                              radarHeight);


    return nRecordsMax;

}




JNIEXPORT jint JNICALL
Java_nl_esciencecenter_ncradar_JNIMethodsVol2Bird_findCells(
        JNIEnv *env,
        jobject obj,
        const jintArray dbzImageInt,
        jintArray cellImageInt,
        const jint dbzMissing,
        const jint dbznAzim,
        const jint dbznRang,
        const jfloat dbzValueOffset,
        const jfloat dbzValueScale,
        const jfloat dbzRangeScale,
        const jfloat dbzThresMin,
        const jfloat rCellMax)
{

    int nAzim = dbznAzim;
    int nRang = dbznRang;
    int iAzim;
    int iRang;
    int iGlobal;
    int nGlobal;


    nGlobal = nAzim * nRang;

    // do some Java Native interface tricks:

    jint *dbzImageIntBody = (*env)->GetIntArrayElements(env, dbzImageInt, NULL);
    unsigned char dbzImageBody[nGlobal];
    for (iAzim = 0; iAzim < nAzim; iAzim++){
         for (iRang= 0 ; iRang<nRang;iRang++){
             iGlobal = iAzim*nRang + iRang;
             if (0<=dbzImageIntBody[iGlobal] && dbzImageIntBody[iGlobal]<=255) {
                 dbzImageBody[iGlobal] = (unsigned char) dbzImageIntBody[iGlobal];
             }
             else {
                 fprintf(stderr,"Error converting type (dbzImageIntBody[iGlobal]).\n");
                 return -1;
             }
         }
     }


    jint *cellImageIntBody = (*env)->GetIntArrayElements(env, cellImageInt, NULL);
    // end of Java Native Interface tricks


    // Allocating and initializing memory for cell properties.
    SCANMETA dbzMeta;

    dbzMeta.missing = dbzMissing;
    dbzMeta.nAzim = dbznAzim;
    dbzMeta.nRang = dbznRang;
    dbzMeta.valueOffset = dbzValueOffset;
    dbzMeta.valueScale = dbzValueScale;
    dbzMeta.rangeScale = dbzRangeScale;

    int nCells;

    nCells = findCells(&dbzImageBody[0], &cellImageIntBody[0],
                       &dbzMeta,
                       dbzThresMin,
                       rCellMax);


    for (iAzim = 0; iAzim < nAzim; iAzim++){
        for (iRang= 0 ; iRang<nRang;iRang++){

            iGlobal = iAzim*nRang + iRang;
            dbzImageIntBody[iGlobal] = (jint) dbzImageBody[iGlobal];

        }
    }



    // do some Java Native interface tricks:
    (*env)->ReleaseIntArrayElements(env, cellImageInt, cellImageIntBody, 0);
    (*env)->ReleaseIntArrayElements(env, dbzImageInt, dbzImageIntBody, JNI_ABORT);
    // end of Java Native Interface tricks


    return nCells;


}




JNIEXPORT jint JNICALL
Java_nl_esciencecenter_ncradar_JNIMethodsVol2Bird_findNearbyGateIndex(
JNIEnv *env,
jobject obj,
const jint nAzimParent,
const jint nRangParent,
const jint iParent,
const jint nAzimChild,
const jint nRangChild,
const jint iChild)
{

    return findNearbyGateIndex(nAzimParent, nRangParent, iParent,
                               nAzimChild,  nRangChild,  iChild);

}








JNIEXPORT void JNICALL
Java_nl_esciencecenter_ncradar_JNIMethodsVol2Bird_fringeCells(
JNIEnv *env,
jobject obj,
jintArray cellImage,
const jint nRang,
const jint nAzim,
const jfloat aScale,
const jfloat rScale,
const jfloat fringeDist)
{

    // do some Java Native interface tricks:
    jint *cellImageBody = (*env)->GetIntArrayElements(env, cellImage, NULL);
    // end of Java Native Interface tricks

    fringeCells(cellImageBody,nRang,nAzim,aScale,rScale,fringeDist);

    // do some Java Native interface tricks:
    (*env)->ReleaseIntArrayElements(env, cellImage, cellImageBody, 0);
    // end of Java Native Interface tricks


}






JNIEXPORT int JNICALL
Java_nl_esciencecenter_ncradar_JNIMethodsVol2Bird_getListOfSelectedGates(
JNIEnv *env,
jobject obj,
jint nRang,
jint nAzim,
jfloat rangeScale,
jfloat azimuthScale,
jfloat elevAngle,
jint missing,
jfloat radarHeight,
jfloat vradValueOffset,
jfloat vradValueScale,
jintArray vradImageInt,
jfloat dbzValueOffset,
jfloat dbzValueScale,
jintArray dbzImageInt,
jfloatArray listOfAzimuths,
jfloatArray listOfElevAngles,
jfloatArray listOfVradObs,
jfloatArray listOfDbzObs,
jintArray listOfCellIds,
jintArray cellImage,
jfloat rangeMin,
jfloat rangeMax,
jfloat altitudeMin,
jfloat altitudeMax,
jfloat absVradMin,
jint iData,
jint nPoints)
{

    int iAzim;
    int iRang;
    int iGlobal;
    int iPoint;

    // do some Java Native interface tricks:
    jint *vradImageIntBody = (*env)->GetIntArrayElements(env, vradImageInt, NULL);
    jint *dbzImageIntBody = (*env)->GetIntArrayElements(env, dbzImageInt, NULL);
    jint *cellImageBody = (*env)->GetIntArrayElements(env, cellImage, NULL);
    jsize nGlobal = (*env)->GetArrayLength(env, vradImageInt);

    jfloat *listOfAzimuthsBody = (*env)->GetFloatArrayElements(env, listOfAzimuths, NULL);
    jfloat *listOfElevAnglesBody = (*env)->GetFloatArrayElements(env, listOfElevAngles, NULL);
    jfloat *listOfVradObsBody = (*env)->GetFloatArrayElements(env, listOfVradObs, NULL);
    jfloat *listOfDbzObsBody = (*env)->GetFloatArrayElements(env, listOfDbzObs, NULL);
    jint *listOfCellIdsBody = (*env)->GetIntArrayElements(env, listOfCellIds, NULL);
    // end of Java Native Interface tricks

    unsigned char vradImageBody[nGlobal];
    unsigned char dbzImageBody[nGlobal];

    for (iAzim = 0; iAzim < nAzim; iAzim++){
         for (iRang = 0; iRang < nRang; iRang++){
             iGlobal = iAzim*nRang + iRang;
             if (0<=vradImageIntBody[iGlobal] && vradImageIntBody[iGlobal]<=255) {
                 vradImageBody[iGlobal] = (unsigned char) vradImageIntBody[iGlobal];
             }
             else {
                 fprintf(stderr,"Error converting type (vradImageIntBody[iGlobal]).\n");
                 return -1;
             }
             if (0<=dbzImageIntBody[iGlobal] && dbzImageIntBody[iGlobal]<=255) {
                 dbzImageBody[iGlobal] = (unsigned char) dbzImageIntBody[iGlobal];
             }
             else {
                 fprintf(stderr,"Error converting type (dbzImageIntBody[iGlobal]).\n");
                 return -1;
             }
         }
     }


    SCANMETA vradMeta;

    vradMeta.nRang = nRang;
    vradMeta.nAzim = nAzim;
    vradMeta.rangeScale = rangeScale;
    vradMeta.azimScale = azimuthScale;
    vradMeta.elev = elevAngle;
    if (0<=missing && missing<=255) {
        vradMeta.missing = (unsigned char) missing;
    } else {
        fprintf(stderr,"Error converting type (vradImageIntBody[iGlobal]).\n");
        return -1;
    }
    vradMeta.heig = radarHeight;
    vradMeta.valueOffset = vradValueOffset;
    vradMeta.valueScale = vradValueScale;


    SCANMETA dbzMeta;
    dbzMeta.valueOffset = dbzValueOffset;
    dbzMeta.valueScale = dbzValueScale;

    fprintf(stderr, "nPoints = %d\n", nPoints);


    getListOfSelectedGates(vradMeta, &vradImageBody[0],
                             dbzMeta, &dbzImageBody[0],
                             &cellImageBody[0],
                             rangeMin, rangeMax,
                             altitudeMin, altitudeMax,
                             absVradMin, iData,
                             &nPoints, &listOfAzimuthsBody[0], &listOfElevAnglesBody[0], &listOfVradObsBody[0],
                             &listOfDbzObsBody[0], &listOfCellIdsBody[0]);






    for (iPoint = 0; iPoint < nPoints ; iPoint++) {
        if (listOfAzimuthsBody[iPoint] != -1.0f) {
            fprintf(stderr, "%5d %10.3f %10.3f %10.3f %10.3f %4d\n", iPoint, listOfAzimuthsBody[iPoint], listOfElevAnglesBody[iPoint],
                    listOfVradObsBody[iPoint], listOfDbzObsBody[iPoint], listOfCellIdsBody[iPoint]);
        }
    }
    fprintf(stderr, "nPoints = %d\n", nPoints);

    // cast to the right type:
    for (iAzim = 0; iAzim < nAzim; iAzim++){
        for (iRang = 0; iRang < nRang; iRang++){

            iGlobal = iAzim*nRang + iRang;
            vradImageIntBody[iGlobal] = (jint) vradImageBody[iGlobal];

        }
    }


    // do some Java Native interface tricks:

    (*env)->ReleaseFloatArrayElements(env, listOfAzimuths, listOfAzimuthsBody, 0);
    (*env)->ReleaseFloatArrayElements(env, listOfElevAngles, listOfElevAnglesBody, 0);
    (*env)->ReleaseFloatArrayElements(env, listOfVradObs, listOfVradObsBody, 0);
    (*env)->ReleaseFloatArrayElements(env, listOfDbzObs, listOfDbzObsBody, 0);
    (*env)->ReleaseIntArrayElements(env, listOfCellIds, listOfCellIdsBody, 0);

    (*env)->ReleaseIntArrayElements(env, cellImage, cellImageBody, JNI_ABORT);
    (*env)->ReleaseIntArrayElements(env, vradImageInt, vradImageIntBody, JNI_ABORT);

    // end of Java Native Interface tricks

    return nPoints;


}









JNIEXPORT void JNICALL
Java_nl_esciencecenter_ncradar_JNIMethodsVol2Bird_sortCells(
JNIEnv *env,
jobject obj,
jintArray cellPropIRangOfMax,
jintArray cellPropIAzimOfMax,
jfloatArray cellPropDbzAvg,
jfloatArray cellPropTexAvg,
jfloatArray cellPropCv,
jfloatArray cellPropArea,
jfloatArray cellPropClutterArea,
jfloatArray cellPropDbzMax,
jintArray cellPropIndex,
jintArray cellPropDrop,
const jint nCells
)
{

    // do some Java Native interface tricks:

    // these arrays have nCells elements:
    jint *cellPropIRangOfMaxBody = (*env)->GetIntArrayElements(env, cellPropIRangOfMax, NULL);
    jint *cellPropIAzimOfMaxBody = (*env)->GetIntArrayElements(env, cellPropIAzimOfMax, NULL);
    jfloat *cellPropDbzAvgBody = (*env)->GetFloatArrayElements(env, cellPropDbzAvg, NULL);
    jfloat *cellPropTexAvgBody = (*env)->GetFloatArrayElements(env, cellPropTexAvg, NULL);
    jfloat *cellPropCvBody = (*env)->GetFloatArrayElements(env, cellPropCv, NULL);
    jfloat *cellPropAreaBody = (*env)->GetFloatArrayElements(env, cellPropArea, NULL);
    jfloat *cellPropClutterAreaBody = (*env)->GetFloatArrayElements(env, cellPropClutterArea, NULL);
    jfloat *cellPropDbzMaxBody = (*env)->GetFloatArrayElements(env, cellPropDbzMax, NULL);
    jint *cellPropIndexBody = (*env)->GetIntArrayElements(env, cellPropIndex, NULL);
    jint *cellPropDropBody = (*env)->GetIntArrayElements(env, cellPropDrop, NULL);
    // end of Java Native Interface tricks

    int iCell;

    CELLPROP *cellProp;
    /*Allocating and initializing memory for cell properties.*/
    cellProp = (CELLPROP *)malloc(nCells*sizeof(CELLPROP));

    // construct the CELLPROP struct
    for (iCell = 0; iCell<nCells;iCell++) {
        cellProp[iCell].iRangOfMax = cellPropIRangOfMaxBody[iCell];
        cellProp[iCell].iAzimOfMax = cellPropIAzimOfMaxBody[iCell];
        cellProp[iCell].dbzAvg = cellPropDbzAvgBody[iCell];
        cellProp[iCell].texAvg = cellPropTexAvgBody[iCell];
        cellProp[iCell].cv = cellPropCvBody[iCell];
        cellProp[iCell].area = cellPropAreaBody[iCell];
        cellProp[iCell].clutterArea = cellPropClutterAreaBody[iCell];
        cellProp[iCell].dbzMax = cellPropDbzMaxBody[iCell];
        cellProp[iCell].index = cellPropIndexBody[iCell];
        cellProp[iCell].drop = cellPropDropBody[iCell];
    }

    for (iCell=0;iCell<nCells;iCell++) {
        fprintf(stderr,"B: cellProp[%d].area = %d\n",iCell,cellProp[iCell].area);
    }

    sortCells(cellProp,nCells);

    for (iCell=0;iCell<nCells;iCell++) {
        fprintf(stderr,"A: cellProp[%d].area = %d\n",iCell,cellProp[iCell].area);
    }

    // deconstruct the CELLPROP struct
    for (iCell = 0; iCell<nCells;iCell++) {
        cellPropIRangOfMaxBody[iCell] = cellProp[iCell].iRangOfMax;
        cellPropIAzimOfMaxBody[iCell] = cellProp[iCell].iAzimOfMax;
        cellPropDbzAvgBody[iCell] = cellProp[iCell].dbzAvg;
        cellPropTexAvgBody[iCell] = cellProp[iCell].texAvg;
        cellPropCvBody[iCell] = cellProp[iCell].cv;
        cellPropAreaBody[iCell] = cellProp[iCell].area;
        cellPropClutterAreaBody[iCell] = cellProp[iCell].clutterArea;
        cellPropDbzMaxBody[iCell] = cellProp[iCell].dbzMax;
        cellPropIndexBody[iCell] = cellProp[iCell].index;
        cellPropDropBody[iCell] = cellProp[iCell].drop;

    }



    // do some more Java Native Interface tricks:
    (*env)->ReleaseIntArrayElements(env, cellPropIRangOfMax, cellPropIRangOfMaxBody, 0);
    (*env)->ReleaseIntArrayElements(env, cellPropIAzimOfMax, cellPropIAzimOfMaxBody, 0);
    (*env)->ReleaseFloatArrayElements(env, cellPropDbzAvg, cellPropDbzAvgBody, 0);
    (*env)->ReleaseFloatArrayElements(env, cellPropTexAvg, cellPropTexAvgBody, 0);
    (*env)->ReleaseFloatArrayElements(env, cellPropCv, cellPropCvBody, 0);
    (*env)->ReleaseFloatArrayElements(env, cellPropArea, cellPropAreaBody, 0);
    (*env)->ReleaseFloatArrayElements(env, cellPropClutterArea, cellPropClutterAreaBody, 0);
    (*env)->ReleaseFloatArrayElements(env, cellPropDbzMax, cellPropDbzMaxBody, 0);
    (*env)->ReleaseIntArrayElements(env, cellPropIndex, cellPropIndexBody, 0);
    (*env)->ReleaseIntArrayElements(env, cellPropDrop, cellPropDropBody, 0);
    // end of Java Native Interface tricks


    // free up the memory that we allocated earlier
    free(cellProp);
    cellProp = NULL;



}








JNIEXPORT jint JNICALL
Java_nl_esciencecenter_ncradar_JNIMethodsVol2Bird_updateMap(
JNIEnv *env,
jobject obj,
jintArray cellImage,
jintArray cellPropIRangOfMax,
jintArray cellPropIAzimOfMax,
jfloatArray cellPropDbzAvg,
jfloatArray cellPropTexAvg,
jfloatArray cellPropCv,
jfloatArray cellPropArea,
jfloatArray cellPropClutterArea,
jfloatArray cellPropDbzMax,
jintArray cellPropIndex,
jintArray cellPropDrop,
const jint nCells,
const jint nGlobal,
const jint minCellArea)
{

    // do some Java Native interface tricks:
    // arrays with nGlobal elements:
    jint *cellImageBody = (*env)->GetIntArrayElements(env, cellImage, NULL);
    // arrays with nCells elements:
    jint *cellPropIRangOfMaxBody = (*env)->GetIntArrayElements(env, cellPropIRangOfMax, NULL);
    jint *cellPropIAzimOfMaxBody = (*env)->GetIntArrayElements(env, cellPropIAzimOfMax, NULL);
    jfloat *cellPropDbzAvgBody = (*env)->GetFloatArrayElements(env, cellPropDbzAvg, NULL);
    jfloat *cellPropTexAvgBody = (*env)->GetFloatArrayElements(env, cellPropTexAvg, NULL);
    jfloat *cellPropCvBody = (*env)->GetFloatArrayElements(env, cellPropCv, NULL);
    jfloat *cellPropAreaBody = (*env)->GetFloatArrayElements(env, cellPropArea, NULL);
    jfloat *cellPropClutterAreaBody = (*env)->GetFloatArrayElements(env, cellPropClutterArea, NULL);
    jfloat *cellPropDbzMaxBody = (*env)->GetFloatArrayElements(env, cellPropDbzMax, NULL);
    jint *cellPropIndexBody = (*env)->GetIntArrayElements(env, cellPropIndex, NULL);
    jint *cellPropDropBody = (*env)->GetIntArrayElements(env, cellPropDrop, NULL);
    // end of Java Native Interface tricks

    int iCell;
    int nCellsValid;

    CELLPROP *cellProp;
    /*Allocating and initializing memory for cell properties.*/
    cellProp = (CELLPROP *)malloc(nCells*sizeof(CELLPROP));

    // construct the CELLPROP struct
    for (iCell = 0; iCell<nCells;iCell++) {
        cellProp[iCell].iRangOfMax = cellPropIRangOfMaxBody[iCell];
        cellProp[iCell].iAzimOfMax = cellPropIAzimOfMaxBody[iCell];
        cellProp[iCell].dbzAvg = cellPropDbzAvgBody[iCell];
        cellProp[iCell].texAvg = cellPropTexAvgBody[iCell];
        cellProp[iCell].cv = cellPropCvBody[iCell];
        cellProp[iCell].area = cellPropAreaBody[iCell];
        cellProp[iCell].clutterArea = cellPropClutterAreaBody[iCell];
        cellProp[iCell].dbzMax = cellPropDbzMaxBody[iCell];
        cellProp[iCell].index = cellPropIndexBody[iCell];
        cellProp[iCell].drop = cellPropDropBody[iCell];
    }

    nCellsValid = updateMap(cellImageBody,nGlobal,cellProp,nCells,minCellArea);


    // deconstruct the CELLPROP struct
    for (iCell = 0; iCell<nCells;iCell++) {
        cellPropIRangOfMaxBody[iCell] = cellProp[iCell].iRangOfMax;
        cellPropIAzimOfMaxBody[iCell] = cellProp[iCell].iAzimOfMax;
        cellPropDbzAvgBody[iCell] = cellProp[iCell].dbzAvg;
        cellPropTexAvgBody[iCell] = cellProp[iCell].texAvg;
        cellPropCvBody[iCell] = cellProp[iCell].cv;
        cellPropAreaBody[iCell] = cellProp[iCell].area;
        cellPropClutterAreaBody[iCell] = cellProp[iCell].clutterArea;
        cellPropDbzMaxBody[iCell] = cellProp[iCell].dbzMax;
        cellPropIndexBody[iCell] = cellProp[iCell].index;
        cellPropDropBody[iCell] = cellProp[iCell].drop;
    }


    // do some Java Native interface tricks:
    // arrays with nGlobal elements:
    (*env)->ReleaseIntArrayElements(env, cellImage, cellImageBody, 0);
    // arrays with nCells elements:
    (*env)->ReleaseIntArrayElements(env, cellPropIRangOfMax, cellPropIRangOfMaxBody, 0);
    (*env)->ReleaseIntArrayElements(env, cellPropIAzimOfMax, cellPropIAzimOfMaxBody, 0);
    (*env)->ReleaseFloatArrayElements(env, cellPropDbzAvg, cellPropDbzAvgBody, 0);
    (*env)->ReleaseFloatArrayElements(env, cellPropTexAvg, cellPropTexAvgBody, 0);
    (*env)->ReleaseFloatArrayElements(env, cellPropCv, cellPropCvBody, 0);
    (*env)->ReleaseFloatArrayElements(env, cellPropArea, cellPropAreaBody, 0);
    (*env)->ReleaseFloatArrayElements(env, cellPropClutterArea, cellPropClutterAreaBody, 0);
    (*env)->ReleaseFloatArrayElements(env, cellPropDbzMax, cellPropDbzMaxBody, 0);
    (*env)->ReleaseIntArrayElements(env, cellPropIndex, cellPropIndexBody, 0);
    (*env)->ReleaseIntArrayElements(env, cellPropDrop, cellPropDropBody, 0);
    // end of Java Native Interface tricks

    free(cellProp);
    cellProp = NULL;

    return nCellsValid;


}







JNIEXPORT jfloat JNICALL
Java_nl_esciencecenter_ncradar_JNIMethodsVol2Bird_svdfit(
JNIEnv *env,
jobject obj,
jfloatArray points,
jint nDims,
jfloatArray yObs,
jfloatArray yFitted,
jint nPoints,
jfloatArray parameterVector,
jfloatArray avar,
jint nParsFitted
)
{

    // do some Java Native interface tricks:
    jfloat *pointsBody = (*env)->GetFloatArrayElements(env, points, NULL);
    jfloat *yObsBody = (*env)->GetFloatArrayElements(env, yObs, NULL);
    jfloat *yFittedBody = (*env)->GetFloatArrayElements(env, yFitted, NULL);
    jfloat *parameterVectorBody = (*env)->GetFloatArrayElements(env, parameterVector, NULL);
    jfloat *avarBody = (*env)->GetFloatArrayElements(env, avar, NULL);
    // end of Java Native Interface tricks

    float chisq;

    chisq = svdfit(&pointsBody[0], nDims, &yObsBody[0], &yFittedBody[0], nPoints,
            &parameterVectorBody[0], &avarBody[0], nParsFitted);

    // do some Java Native interface tricks:
    (*env)->ReleaseFloatArrayElements(env, points, pointsBody, JNI_ABORT);
    (*env)->ReleaseFloatArrayElements(env, yObs, yObsBody, JNI_ABORT);
    (*env)->ReleaseFloatArrayElements(env, yFitted, yFittedBody, 0);
    (*env)->ReleaseFloatArrayElements(env, parameterVector, parameterVectorBody, 0);
    (*env)->ReleaseFloatArrayElements(env, avar, avarBody, 0);
    // end of Java Native Interface tricks


    return chisq;



}





JNIEXPORT jint JNICALL
Java_nl_esciencecenter_ncradar_JNIMethodsVol2Bird_svdcmp(
JNIEnv *env,
jobject obj,
jfloatArray arrayA,
jint nRows,
jint nCols,
jfloatArray arrayW,
jfloatArray arrayV
)
{

    if (nCols>nRows) {
        fprintf(stderr,"Behavior undefined for matrices that have more columns than rows.\n");
    };

    // do some Java Native interface tricks:
    jfloat *arrayABody = (*env)->GetFloatArrayElements(env, arrayA, NULL);
    jfloat *arrayWBody = (*env)->GetFloatArrayElements(env, arrayW, NULL);
    jfloat *arrayVBody = (*env)->GetFloatArrayElements(env, arrayV, NULL);
    // end of Java Native Interface tricks

    int status;

//    int iRow;
//    int iCol;
//    int iGlobal;
//
//    for (iRow = 0; iRow < nRows; iRow++) {
//        for (iCol = 0; iCol < nCols; iCol++) {
//            iGlobal = iRow*nCols + iCol;
//            fprintf(stderr, "A[%d,%d] = %f\n",iRow,iCol,arrayABody[iGlobal]);
//        }
//    }
//
//    fprintf(stderr,"\n");
//
//    for (iCol = 0; iCol < nCols; iCol++) {
//        fprintf(stderr, "W[%d] = %f\n",iCol,arrayWBody[iCol]);
//    }
//
//    fprintf(stderr,"\n");
//
//    int iCol1;
//    int iCol2;
//
//    for (iCol1 = 0; iCol1 < nCols; iCol1++) {
//        for (iCol2 = 0; iCol2 < nCols; iCol2++) {
//            iGlobal = iCol1*nCols + iCol2;
//            fprintf(stderr, "V[%d,%d] = %f\n",iCol1,iCol2,arrayVBody[iGlobal]);
//        }
//    }
//
//    fprintf(stderr,"\n");

    status = svdcmp(&arrayABody[0], nRows, nCols, &arrayWBody[0], &arrayVBody[0]);

//    fprintf(stderr,"= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = \n");
//
//    for (iRow = 0; iRow < nRows; iRow++) {
//        for (iCol = 0; iCol < nCols; iCol++) {
//            iGlobal = iRow*nCols + iCol;
//            fprintf(stderr, "U[%d,%d] = %f\n",iRow,iCol,arrayABody[iGlobal]);
//        }
//    }
//
//    fprintf(stderr,"\n");
//
//    for (iCol = 0; iCol < nCols; iCol++) {
//        fprintf(stderr, "W[%d] = %f\n",iCol,arrayWBody[iCol]);
//    }
//
//    fprintf(stderr,"\n");
//
//
//    for (iCol1 = 0; iCol1 < nCols; iCol1++) {
//        for (iCol2 = 0; iCol2 < nCols; iCol2++) {
//            iGlobal = iCol1*nCols + iCol2;
//            fprintf(stderr, "V[%d,%d] = %f\n",iCol1,iCol2,arrayVBody[iGlobal]);
//        }
//    }


    // do some Java Native interface tricks:
    (*env)->ReleaseFloatArrayElements(env, arrayA, arrayABody, 0);
    (*env)->ReleaseFloatArrayElements(env, arrayW, arrayWBody, 0);
    (*env)->ReleaseFloatArrayElements(env, arrayV, arrayVBody, 0);
    // end of Java Native Interface tricks

    return status;

}







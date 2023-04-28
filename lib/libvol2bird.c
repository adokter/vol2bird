/** vol2bird main functions
 * @file libvol2bird.c
 * @author Adriaan Dokter & Netherlands eScience Centre
 * @date 2015-05-05
 */

/*
 * Copyright 2015 Adriaan Dokter & Netherlands eScience Centre
 * If you want to use this software, please contact me at a.m.dokter@uva.nl
 *
 * This program calculates Vertical Profiles of Birds (VPBs) as described in
 *
 * Bird migration flight altitudes studied by a network of operational weather radars
 * Dokter A.M., Liechti F., Stark H., Delobbe L., Tabary P., Holleman I.
 * J. R. Soc. Interface, 8, 30–43, 2011
 * DOI: 10.1098/rsif.2010.0116
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#ifndef NOCONFUSE
#include <confuse.h>
#endif
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <unistd.h>
#include <vertical_profile.h>
#include "rave_io.h"
#include "rave_debug.h"
#include "polarvolume.h"
#include "polarscan.h"
#include "libvol2bird.h"
#include "libsvdfit.h"
#include "constants.h"
#undef RAD2DEG // to suppress redefine warning, also defined in dealias.h
#undef DEG2RAD // to suppress redefine warning, also defined in dealias.h
#include "libdealias.h"
#include <assert.h>
#include "librender.h"


#ifdef RSL
#include "rsl.h"
#include "librsl.h"
#endif

#ifdef IRIS
#include "iris2odim.h"
#endif


// non-public function prototypes (local to this file/translation unit)

static int analyzeCells(PolarScan_t *scan, vol2birdScanUse_t scanUse, const int nCells, int dualpol, vol2bird_t *alldata);

static float calcDist(const int range1, const int azim1, const int range2, const int azim2, const float rscale, const float ascale);

static void calcTexture(PolarScan_t *scan, vol2birdScanUse_t scanUse, vol2bird_t* alldata);

static void classifyGatesSimple(vol2bird_t* alldata);

static void constructPointsArray(PolarVolume_t* volume, vol2birdScanUse_t *scanUse, vol2bird_t* alldata);

static int detNumberOfGates(const int iLayer, const float rangeScale, const float elevAngle,
                            const int nRang, const int nAzim, const float radarHeight, vol2bird_t* alldata);

static int detSvdfitArraySize(PolarVolume_t* volume, vol2birdScanUse_t *scanUse, vol2bird_t* alldata);

static vol2birdScanUse_t *determineScanUse(PolarVolume_t* volume, vol2bird_t* alldata);

static void exportBirdProfileAsJSON(vol2bird_t* alldata);

static int findWeatherCells(PolarScan_t *scan, const char* quantity, float quantityThreshold,
        int selectAboveThreshold, int iCellStart, int initialize, vol2bird_t* alldata);

static int findNearbyGateIndex(const int nAzimParent, const int nRangParent, const int iParent,
                        const int nAzimChild,  const int nRangChild,  const int iChild, int *iAzimReturn, int *iRangReturn);

static void fringeCells(PolarScan_t* scan, vol2bird_t* alldata);

CELLPROP* getCellProperties(PolarScan_t* scan, vol2birdScanUse_t scanUse, const int nCells, vol2bird_t* alldata);

static int getListOfSelectedGates(PolarScan_t* scan, vol2birdScanUse_t scanUse,
                                      const float altitudeMin, const float altitudeMax,
                                      float* points_local, int iRowPoints, int nColsPoints_local, vol2bird_t* alldata);

static int hasAzimuthGap(const float *points_local, const int nPoints, vol2bird_t* alldata);

static int includeGate(const int iProfileType, const int iQuantityType, const unsigned int gateCode, vol2bird_t* alldata);

const char* libvol2bird_version(void);

static int verticalProfile_AddCustomField(VerticalProfile_t* self, RaveField_t* field, const char* quantity);

static int profileArray2RaveField(vol2bird_t* alldata, int idx_profile, int idx_quantity, const char* quantity, RaveDataType raveType);

static int mapVolumeToProfile(VerticalProfile_t* vp, PolarVolume_t* volume);

PolarScanParam_t* PolarScan_newParam(PolarScan_t *scan, const char *quantity, RaveDataType type);

int PolarVolume_dealias(PolarVolume_t* pvol);

int PolarVolume_getEndDateTime(PolarVolume_t* pvol, char** EndDate, char** EndTime);

int PolarVolume_getStartDateTime(PolarVolume_t* pvol, char** StartDate, char** StartTime);

const char* PolarVolume_getEndTime(PolarVolume_t* pvol);

const char* PolarVolume_getEndDate(PolarVolume_t* pvol);

const char* PolarVolume_getStartDate(PolarVolume_t* pvol);

const char* PolarVolume_getStartTime(PolarVolume_t* pvol);

double PolarVolume_getWavelength(PolarVolume_t* pvol);

PolarScan_t* PolarScan_resample(PolarScan_t* scan, double rscale_proj, long nbins_proj, long nrays_proj);

PolarScanParam_t* PolarScanParam_resample(PolarScanParam_t* param, double rscale, double rscale_proj, long nbins_proj, long nrays_proj);

static void printCellProp(CELLPROP* cellProp, float elev, int nCells, int nCellsValid, vol2bird_t *alldata);

static void printGateCode(char* flags, const unsigned int gateCode);

static void printImage(PolarScan_t* scan, const char* quantity);

static int printMeta(PolarScan_t* scan, const char* quantity);

static void printProfile(vol2bird_t* alldata);

static int removeDroppedCells(CELLPROP *cellProp, const int nCells);

static int selectCellsToDrop(CELLPROP *cellProp, int nCells, int dualpol, vol2bird_t* alldata);

static int selectCellsToDrop_singlePol(CELLPROP *cellProp, int nCells, vol2bird_t* alldata);

static int selectCellsToDrop_dualPol(CELLPROP *cellProp, int nCells, vol2bird_t* alldata);

static void sortCellsByArea(CELLPROP *cellProp, const int nCells);

static void updateFlagFieldsInPointsArray(const float* yObs, const float* yFitted, const int* includedIndex, 
                                          const int nPointsIncluded, float* points_local, vol2bird_t* alldata);

static int updateMap(PolarScan_t* scan, CELLPROP *cellProp, const int nCells, vol2bird_t* alldata);

#ifdef IRIS
PolarVolume_t* vol2birdGetIRISVolume(char* filenames[], int nInputFiles);
#endif

PolarVolume_t* vol2birdGetODIMVolume(char* filenames[], int nInputFiles);

#ifdef VOL2BIRD_R
int check_mistnet_loaded_c(void);
#endif

static vol2bird_printfun vol2bird_internal_printf_fun = vol2bird_default_print;

static vol2bird_printfun vol2bird_internal_err_printf_fun = vol2bird_default_err_print;

// non-public function declarations (local to this file/translation unit)

void vol2bird_printf(const char* fmt, ...)
{
  va_list ap;
  int n;
  char msg[65536];
  va_start(ap, fmt);
  n = vsnprintf(msg, 1024, fmt, ap);
  va_end(ap);
  if (n >= 0 && n <= 65536) {
    vol2bird_internal_printf_fun(msg);
  } else {
    vol2bird_internal_printf_fun("vol2bird_printf failed when printing message");
  }
}

void vol2bird_err_printf(const char* fmt, ...)
{
  va_list ap;
  int n;
  char msg[65536];
  va_start(ap, fmt);
  n = vsnprintf(msg, 1024, fmt, ap);
  va_end(ap);
  if (n >= 0 && n <= 65536) {
    vol2bird_internal_err_printf_fun(msg);
  } else {
    vol2bird_internal_err_printf_fun("vol2bird_err_printf failed when printing message");
  }
}

void vol2bird_default_print(const char* msg)
{
#ifndef NO_VOL2BIRD_PRINTF
  fprintf(stdout, "%s", msg);
#endif
}

void vol2bird_default_err_print(const char* msg)
{
#ifndef NO_VOL2BIRD_PRINTF
  fprintf(stderr, "%s", msg);
#endif
}

void vol2bird_set_printf(vol2bird_printfun fun)
{
  if (fun != NULL) {
    vol2bird_internal_printf_fun = fun;
  }
}

void vol2bird_set_err_printf(vol2bird_printfun fun)
{
  if (fun != NULL) {
    vol2bird_internal_err_printf_fun = fun;
  }
}

static int analyzeCells(PolarScan_t *scan, vol2birdScanUse_t scanUse, const int nCells, int dualpol, vol2bird_t *alldata) {

    // ----------------------------------------------------------------------------------- // 
    //  This function analyzes the cellImage array found by the 'findWeatherCells'         //
    //  procedure. Small cells are rejected and the cells are re-numbered according        //
    //  to size. The final number of cells in cellImage is returned as an integer.         // 
    // ----------------------------------------------------------------------------------- //

    CELLPROP *cellProp;
    int nCellsValid;
    long nAzim;
    long nRang;

    nCellsValid = nCells;
    nRang = PolarScan_getNbins(scan);
    nAzim = PolarScan_getNrays(scan);
    nCellsValid = 0;
    
    if(!PolarScan_hasParameter(scan, scanUse.cellName)){
        vol2bird_err_printf("no CELL quantity in polar scan, aborting analyzeCells()\n");
        return 0;
    }
    
    // first deal with the case that no weather cells were detected by findWeatherCells
    if (nCells == 0) {
        for (int iAzim = 0; iAzim < nAzim; iAzim++) {
            for (int iRang = 0; iRang < nRang; iRang++) {
                PolarScan_setParameterValue(scan, scanUse.cellName, iRang, iAzim, -1);
            }
        }
        return nCellsValid;
    }
    
    cellProp = getCellProperties(scan, scanUse, nCells, alldata);

    selectCellsToDrop(cellProp, nCells, dualpol, alldata);    
    
    // sorting cell properties according to cell area. Drop small cells from map
    nCellsValid = updateMap(scan, cellProp, nCells, alldata);
        
    // printing of cell properties to stderr
    if (alldata->options.printCellProp == TRUE) {
        printCellProp(cellProp, (float) PolarScan_getElangle(scan), nCells, nCellsValid, alldata);
    } // endif (printCellProp == TRUE)

    free(cellProp);

    return nCellsValid;
} // analyzeCells



static float calcDist(const int iRang1, const int iAzim1, const int iRang2, 
    const int iAzim2, const float rangScale, const float azimScaleDeg) {

    // -------------------------------------------------------- //
    // This function calculates the distance between two gates  //
    // -------------------------------------------------------- //

    float range1;
    float range2;
    float azimuth1;
    float azimuth2;

    range1 = (float) iRang1 * rangScale;
    range2 = (float) iRang2 * rangScale;

    azimuth1 = (float) iAzim1 * azimScaleDeg * (float) DEG2RAD;
    azimuth2 = (float) iAzim2 * azimScaleDeg * (float) DEG2RAD;

    return (float) sqrt((range1 * range1) +
                (range2 * range2) -
                2 * (range1 * range2) * cos(azimuth1-azimuth2));


} // calcDist



static void calcTexture(PolarScan_t *scan, vol2birdScanUse_t scanUse, vol2bird_t* alldata) {


    // --------------------------------------------------------------------------------------- //
    // This function computes a texture parameter based on a block of (nRangNeighborhood x     //
    // nAzimNeighborhood) pixels. The texture parameter equals the local standard deviation    //
    // in the radial velocity field                                                            //
    // --------------------------------------------------------------------------------------- //


    int iRang, iRangLocal;
    int iAzim, iAzimLocal;
    long nRang;
    long nAzim;
    int iNeighborhood;
    int nNeighborhood;
    int count;
    double vradValGlobal;
    double vradValLocal;
    double dbzValGlobal;
    double dbzValLocal;
    double vradMissingValue;
    double vradUndetectValue;
    double dbzMissingValue;
    double dbzUndetectValue;
    double texMissingValue;
    double vmoment1;
    double vmoment2;
    double dbz;
    double tex;
    int iGlobal;
    int iLocal;
    double texOffset;
    double texScale;
    double dbzOffset;
    double dbzScale;
    double vradOffset;
    double vradScale;
    double vRadDiff;

    nRang = PolarScan_getNbins(scan);
    nAzim = PolarScan_getNrays(scan);

    PolarScanParam_t* texImage = PolarScan_getParameter(scan, scanUse.texName);
    PolarScanParam_t* vradImage = PolarScan_getParameter(scan, scanUse.vradName);
    PolarScanParam_t* dbzImage = PolarScan_getParameter(scan, scanUse.dbzName);

    if(scanUse.useScan != 1){
      vol2bird_err_printf("Error: scanUse unequal to 1 (%i), this scan should not be used\n",scanUse.useScan);
    }
    if (texImage == NULL) {
      vol2bird_err_printf("Error: Couldn't fetch texture parameter for texture calculation\n");
    }
    if (vradImage == NULL) {
      vol2bird_err_printf("Error: Couldn't fetch radial velocity parameter for texture calculation\n");
    }
    if (dbzImage == NULL) {
      vol2bird_err_printf("Error: Couldn't fetch reflectivity parameter for texture calculation\n");
    }

    dbzOffset = PolarScanParam_getOffset(dbzImage);
    dbzScale = PolarScanParam_getGain(dbzImage);
    dbzMissingValue = PolarScanParam_getNodata(dbzImage);
    dbzUndetectValue = PolarScanParam_getUndetect(dbzImage);

    vradOffset = PolarScanParam_getOffset(vradImage);
    vradScale = PolarScanParam_getGain(vradImage);
    vradMissingValue = PolarScanParam_getNodata(vradImage);
    vradUndetectValue = PolarScanParam_getUndetect(vradImage);

    texOffset = PolarScanParam_getOffset(texImage);
    texScale = PolarScanParam_getGain(texImage);
    texMissingValue = PolarScanParam_getNodata(texImage);

    nNeighborhood = (int)(alldata->constants.nRangNeighborhood * alldata->constants.nAzimNeighborhood);

    for (iAzim = 0; iAzim < nAzim; iAzim++) {
        for (iRang = 0; iRang < nRang; iRang++) {

            iGlobal = iRang + iAzim * nRang;

            // count number of direct neighbors above threshold
            count = 0;
            vmoment1 = 0;
            vmoment2 = 0;

            dbz = 0;

            for (iNeighborhood = 0; iNeighborhood < nNeighborhood; iNeighborhood++) {

                iLocal = findNearbyGateIndex(nAzim,nRang,iGlobal,alldata->constants.nAzimNeighborhood,
                         alldata->constants.nRangNeighborhood,iNeighborhood,&iAzimLocal,&iRangLocal);

                #ifdef FPRINTFON
                vol2bird_err_printf("iLocal = %d; ",iLocal);
                #endif

                if (iLocal < 0) {
                    // iLocal less than zero are error codes
                    continue;
                }
                
                PolarScanParam_getValue(vradImage, iRang, iAzim, &vradValGlobal);
                PolarScanParam_getValue(vradImage, iRangLocal, iAzimLocal, &vradValLocal);
                PolarScanParam_getValue(dbzImage, iRang, iAzim, &dbzValGlobal);
                PolarScanParam_getValue(dbzImage, iRangLocal, iAzimLocal, &dbzValLocal);

                
                if (vradValLocal == vradMissingValue || dbzValLocal == dbzMissingValue ||
                    vradValLocal == vradUndetectValue || dbzValLocal == dbzUndetectValue) {
                    continue;
                }

                vRadDiff = vradOffset + vradScale * (vradValGlobal - vradValLocal);
                vmoment1 += vRadDiff;
                vmoment2 += SQUARE(vRadDiff);

                dbz += dbzOffset + dbzScale * dbzValLocal;

                count++;

            }

            vmoment1 /= count;
            vmoment2 /= count;
            dbz /= count;

            // when not enough neighbors, continue
            if (count < alldata->constants.nCountMin) {
                PolarScanParam_setValue(texImage,iRang,iAzim,texMissingValue);
            }
            else {

                tex = sqrt(XABS(vmoment2-SQUARE(vmoment1)));

                float tmpTex = (tex - texOffset) / texScale;
                if (-FLT_MAX <= tmpTex && tmpTex <= FLT_MAX) {
                    PolarScanParam_setValue(texImage,iRang,iAzim,(double) tmpTex);
                }
                else {
                    vol2bird_err_printf("Error casting texture value of %f to float type at texImage[%d]. Aborting.\n",tmpTex,iGlobal);
                    goto done;
                }
                

                #ifdef FPRINTFON
                vol2bird_err_printf(
                        "\n(C) count = %d; nCountMin = %d; vmoment1 = %f; vmoment2 = %f; tex = %f; texBody[%d] = %f\n",
                        count, alldata->constants.nCountMin, vmoment1, vmoment2, tex,
                        iGlobal, tmpTex);
                #endif

            } //else
        } //for
    } //for
done:
    RAVE_OBJECT_RELEASE(texImage);
    RAVE_OBJECT_RELEASE(vradImage);
    RAVE_OBJECT_RELEASE(dbzImage);
} // calcTexture



static void classifyGatesSimple(vol2bird_t* alldata) {
    
    int iPoint;
    
    for (iPoint = 0; iPoint < alldata->points.nRowsPoints; iPoint++) {
    
        const float azimValue = alldata->points.points[iPoint * alldata->points.nColsPoints + alldata->points.azimAngleCol];
        const float dbzValue = alldata->points.points[iPoint * alldata->points.nColsPoints + alldata->points.dbzValueCol];        
        const float vradValue = alldata->points.points[iPoint * alldata->points.nColsPoints + alldata->points.vradValueCol];
        const int cellValue = (int) alldata->points.points[iPoint * alldata->points.nColsPoints + alldata->points.cellValueCol];
        const float clutValue = alldata->points.points[iPoint * alldata->points.nColsPoints + alldata->points.clutValueCol];

        unsigned int gateCode = 0;
        
        if (alldata->options.useClutterMap && clutValue > alldata->options.clutterValueMin) {
            // this gate is true in the static clutter map
            gateCode |= 1<<(alldata->flags.flagPositionStaticClutter);
        }

        if (cellValue > 1) {
            // this gate is part of the cluttermap (without fringe)
            gateCode |= 1<<(alldata->flags.flagPositionDynamicClutter);
        }

        if (cellValue == 1) {
            // this gate is part of the fringe of the cluttermap
            gateCode |= 1<<(alldata->flags.flagPositionDynamicClutterFringe);
        }

        if ((isnan(vradValue) == TRUE) || (isnan(dbzValue) == TRUE)) {
            // this gate has no valid radial velocity data
            gateCode |= 1<<(alldata->flags.flagPositionVradMissing);
        }

        if (dbzValue > alldata->misc.dbzMax) {
            // this gate's dbz value is too high to be due to birds, it must be 
            // caused by something else
            gateCode |= 1<<(alldata->flags.flagPositionDbzTooHighForBirds);
        }

        if (fabs(vradValue) < alldata->constants.vradMin) {
            // this gate's radial velocity is too low; Excluded because possibly clutter
            gateCode |= 1<<(alldata->flags.flagPositionVradTooLow);
        }

        if (FALSE) {
            // this flag is set later on by updateFlagFieldsInPointsArray()
            // flagPositionVDifMax
        }

        if (alldata->options.azimMin < alldata->options.azimMax){
            if ((azimValue < alldata->options.azimMin) || (azimValue > alldata->options.azimMax)) {
                // the user can specify to exclude gates based on their azimuth;
                // this clause is for gates that have too low azimuth
                gateCode |= 1<<(alldata->flags.flagPositionAzimOutOfRange);
            }
        }
        else{
            if ((azimValue < alldata->options.azimMin) && (azimValue > alldata->options.azimMax)) {
                // the user can specify to exclude gates based on their azimuth;
                // this clause is for gates that have too low azimuth
                gateCode |= 1<<(alldata->flags.flagPositionAzimOutOfRange);
            }
        }

        alldata->points.points[iPoint * alldata->points.nColsPoints + alldata->points.gateCodeCol] = (float) gateCode;
        
    }

    return;
}



static void constructPointsArray(PolarVolume_t* volume, vol2birdScanUse_t* scanUse, vol2bird_t* alldata) {
    
        // iterate over the scans in 'volume'
        int iScan;
        int nScans;
        
        // determine how many scan elevations the volume object contains
        nScans = PolarVolume_getNumberOfScans(volume);


        for (iScan = 0; iScan < nScans; iScan++) {
            if (scanUse[iScan].useScan == 1)
            {
                // extract the scan object from the volume object
                PolarScan_t* scan = PolarVolume_getScan(volume, iScan);

                PolarScanParam_t *cellScanParam = NULL;
                PolarScanParam_t *texScanParam = NULL;
                
                // check that CELL parameter is not present, which might be after running MistNet
                if (!PolarScan_hasParameter(scan, CELLNAME)){
                    cellScanParam = PolarScan_newParam(scan, scanUse[iScan].cellName, RaveDataType_INT);
                }
                // only when dealing with normal (non-dual pol) data, generate a vrad texture field
                if (alldata->options.singlePol){
                    // ------------------------------------------------------------- //
                    //                      calculate vrad texture                   //
                    // ------------------------------------------------------------- //

                    texScanParam = PolarScan_newParam(scan, scanUse[iScan].texName, RaveDataType_DOUBLE);

                    calcTexture(scan, scanUse[iScan], alldata);					
                }

                int nCells = -1;

                // ------------------------------------------------------------- //
                //        find (weather) cells in the reflectivity image         //
                // ------------------------------------------------------------- //
				
                if (alldata->options.dualPol && !alldata->options.useMistNet){
                    
                    if (alldata->options.singlePol){
						
                        // first pass: single pol rain filtering
                        nCells = findWeatherCells(scan,scanUse[iScan].dbzName,alldata->options.dbzThresMin,TRUE,2,TRUE,alldata);
                        // first pass: single pol analysis of precipitation cells
                        analyzeCells(scan, scanUse[iScan], nCells, FALSE, alldata);
                        // second pass: dual pol precipitation filtering
                        nCells = findWeatherCells(scan,scanUse[iScan].rhohvName,
                                    alldata->options.rhohvThresMin,TRUE,nCells+1,FALSE,alldata);
                    }
                    else{
                        nCells = findWeatherCells(scan,scanUse[iScan].rhohvName,
                                    alldata->options.rhohvThresMin,TRUE,2,TRUE,alldata);						
                    }

                }

                if (!alldata->options.dualPol && !alldata->options.useMistNet){
                    
                    nCells = findWeatherCells(scan,scanUse[iScan].dbzName,alldata->options.dbzThresMin,TRUE,2,TRUE,alldata);

                }
                
                if (alldata->options.useMistNet){
                    nCells = 2;
                }
                
                if (nCells<0){
                    vol2bird_err_printf("Error: findWeatherCells exited with errors\n");
                    RAVE_OBJECT_RELEASE(scan);
                    RAVE_OBJECT_RELEASE(cellScanParam);
                    RAVE_OBJECT_RELEASE(texScanParam);
                    return;
                }
                
                if (alldata->options.printCellProp == TRUE) {
                    vol2bird_err_printf("(%d/%d): found %d cells.\n",iScan+1, nScans, nCells);
                }
                
                // ------------------------------------------------------------- //
                //                      analyze cells                            //
                // ------------------------------------------------------------- //
                if (!alldata->options.useMistNet){
                    nCells=analyzeCells(scan, scanUse[iScan], nCells, alldata->options.dualPol, alldata);
                }
                // ------------------------------------------------------------- //
                //                     calculate fringe                          //
                // ------------------------------------------------------------- //
    
                fringeCells(scan, alldata); 
                // ------------------------------------------------------------- //
                //            print selected outputs to stderr                   //
                // ------------------------------------------------------------- //
    
                if (alldata->options.printDbz == TRUE) {
                    vol2bird_err_printf("product = dbz\n");
                    printMeta(scan,scanUse[iScan].dbzName);
                    printImage(scan,scanUse[iScan].dbzName);
                }
                if (alldata->options.printVrad == TRUE) {
                    vol2bird_err_printf("product = vrad\n");
                    printMeta(scan,scanUse[iScan].vradName);
                    printImage(scan,scanUse[iScan].vradName);
                }
                if (alldata->options.printRhohv == TRUE) {
                    vol2bird_err_printf("product = rhohv\n");
                    printMeta(scan,scanUse[iScan].rhohvName);
                    printImage(scan,scanUse[iScan].rhohvName);
                }
                if (alldata->options.printTex == TRUE) {
                    vol2bird_err_printf("product = tex\n");
                    printMeta(scan,scanUse[iScan].texName);
                    printImage(scan,scanUse[iScan].texName);
                }
                if (alldata->options.printCell == TRUE) {
                    vol2bird_err_printf("product = cell\n");
                    printMeta(scan,scanUse[iScan].cellName);
                    printImage(scan,scanUse[iScan].cellName);
                }
                if (alldata->options.printClut == TRUE) { 
                    vol2bird_err_printf("product = clut\n");
                    printMeta(scan,scanUse[iScan].clutName);
                    printImage(scan,scanUse[iScan].clutName);
                }
                            
                // ------------------------------------------------------------- //
                //    fill in the appropriate elements in the points array       //
                // ------------------------------------------------------------- //
    
                int iLayer;
                
                for (iLayer = 0; iLayer < alldata->options.nLayers; iLayer++) {
                    
                    float altitudeMin = iLayer * alldata->options.layerThickness;
                    float altitudeMax = (iLayer + 1) * alldata->options.layerThickness;
                    int iRowPoints = alldata->points.indexFrom[iLayer] + alldata->points.nPointsWritten[iLayer];
                        
                    int n = getListOfSelectedGates(scan, scanUse[iScan], altitudeMin, altitudeMax, 
                        &(alldata->points.points[0]), iRowPoints, alldata->points.nColsPoints, alldata);
                    
                    alldata->points.nPointsWritten[iLayer] += n;

                    if (alldata->points.indexFrom[iLayer] + alldata->points.nPointsWritten[iLayer] > alldata->points.indexTo[iLayer]) {
                        vol2bird_err_printf("Problem occurred: writing over existing data\n");
                        RAVE_OBJECT_RELEASE(scan);
                        RAVE_OBJECT_RELEASE(cellScanParam);
                        RAVE_OBJECT_RELEASE(texScanParam);
                        return;
                    }
    
                } // endfor (iLayer = 0; iLayer < nLayers; iLayer++)

                // ------------------------------------------------------------- //
                //                         clean up                              //
                // ------------------------------------------------------------- //
    
                // free previously malloc'ed arrays                
                RAVE_OBJECT_RELEASE(scan);
                RAVE_OBJECT_RELEASE(texScanParam);
                RAVE_OBJECT_RELEASE(cellScanParam);
            }
        } // endfor (iScan = 0; iScan < nScans; iScan++)
}



static int detNumberOfGates(const int iLayer, 
                     const float rangeScale, const float elevAngle,
                     const int nRang, const int nAzim,
                     const float radarHeight, vol2bird_t* alldata) {

    // Determine the number of gates that are within the limits set
    // by (rangeMin,rangeMax) as well as by (iLayer*layerThickness,
    // (iLayer+1)*layerThickness).

    int nGates;
    int iRang;

    float layerHeight;
    float range;
    float beamHeight;


    nGates = 0;
    layerHeight = (iLayer + 0.5) * alldata->options.layerThickness;
    for (iRang = 0; iRang < nRang; iRang++) {
        range = (iRang + 0.5) * rangeScale;
        if (range < alldata->options.rangeMin || range > alldata->options.rangeMax) {
            // the gate is too close to the radar, or too far away
            continue;
        }
        beamHeight = range2height(range, elevAngle) + radarHeight;
        if (fabs(layerHeight - beamHeight) > 0.5*alldata->options.layerThickness) {
            // the gate is not close enough to the altitude layer of interest
            continue;
        }

        #ifdef FPRINTFON
        vol2bird_err_printf("iRang = %d; range = %f; beamHeight = %f\n",iRang,range,beamHeight);
        #endif

        nGates += nAzim;

    } // for iRang

    return nGates;

} // detNumberOfGates()





static int detSvdfitArraySize(PolarVolume_t* volume, vol2birdScanUse_t* scanUse, vol2bird_t* alldata) {
    
    int iScan;
    int nScans = PolarVolume_getNumberOfScans(volume);

    int iLayer;
    int nRowsPoints_local = 0;
    
    int* nGates = malloc(sizeof(int) * alldata->options.nLayers);
    int* nGatesAcc = malloc(sizeof(int) * alldata->options.nLayers);
    
    for (iLayer = 0; iLayer < alldata->options.nLayers; iLayer++) {
        nGates[iLayer] = 0;
        nGatesAcc[iLayer] = 0;
    }

    for (iScan = 0; iScan < nScans; iScan++) {
        if (scanUse[iScan].useScan == 1)
        {
            for (iLayer = 0; iLayer < alldata->options.nLayers; iLayer++) {
    
                PolarScan_t* scan = PolarVolume_getScan(volume, iScan);
                
                int nRang = (int) PolarScan_getNbins(scan);
                int nAzim = (int) PolarScan_getNrays(scan);
                float elevAngle = (float) PolarScan_getElangle(scan);
                float rangeScale = (float) PolarScan_getRscale(scan);
                float radarHeight = (float) PolarScan_getHeight(scan);

                nGates[iLayer] += detNumberOfGates(iLayer, 
                    rangeScale, elevAngle, nRang, 
                    nAzim, radarHeight, alldata);
                    
                RAVE_OBJECT_RELEASE(scan);

            }
        }
    }

    for (iLayer = 0; iLayer < alldata->options.nLayers; iLayer++) {
        if (iLayer == 0) {
            nGatesAcc[iLayer] = nGates[iLayer];
        }
        else {
            nGatesAcc[iLayer] = nGatesAcc[iLayer-1] + nGates[iLayer];
        }
    }
    nRowsPoints_local = nGatesAcc[alldata->options.nLayers-1];

    for (iLayer = 0; iLayer < alldata->options.nLayers; iLayer++) {
        if (iLayer == 0) {
            alldata->points.indexFrom[iLayer] = 0;
        }
        else {
            alldata->points.indexFrom[iLayer] = nGatesAcc[iLayer-1];
        }
    }
    for (iLayer = 0; iLayer < alldata->options.nLayers; iLayer++) {
        alldata->points.indexTo[iLayer] = nGatesAcc[iLayer];
    }

    free((void*) nGates);
    free((void*) nGatesAcc);

    return nRowsPoints_local;
    
}  // detSvdfitArraySize()




// Determine whether to use scan
static vol2birdScanUse_t* determineScanUse(PolarVolume_t* volume, vol2bird_t* alldata)
{
    RAVE_ASSERT((volume != NULL), "pvol == NULL");
    
	RaveAttribute_t *attr;
	PolarScan_t *scan;
	PolarScanParam_t *param = NULL;
	int result, nScans, iScan, nScansUsed;
	int noNyquist=0;
	vol2birdScanUse_t *scanUse;
	double nyquist, nyquistMin = DBL_MAX, nyquistMinUsed = DBL_MAX, nyquistMax = 0;
 
	
	// Read number of scans
	nScans = PolarVolume_getNumberOfScans(volume);
    
    // variable that will contain the number of scans to be used as determined by this function
    nScansUsed = 0;
	
	// Allocate memory for useScan variable
	scanUse = (vol2birdScanUse_t *) malloc(nScans * sizeof(vol2birdScanUse_t));
	
    // check that correlation coefficient is present
    // if not, revert to single pol mode
    if (alldata->options.dualPol){
        int dualPolPresent = FALSE;
        
        for (iScan = 0; iScan < nScans; iScan++)
        {
            scan = PolarVolume_getScan(volume, iScan);
            if (PolarScan_hasParameter(scan, "RHOHV")){
                dualPolPresent = TRUE;
            }
            RAVE_OBJECT_RELEASE(scan);
        }
        if (!dualPolPresent){
            vol2bird_err_printf("Warning: no dual-pol moments found, switching to SINGLE POL mode\n");
            alldata->options.dualPol = FALSE;
        }
    }
    
	for (iScan = 0; iScan < nScans; iScan++)
	{
		// Initialize useScan and result
		scanUse[iScan].useScan = FALSE;
		result = 0;
		
		scan = PolarVolume_getScan(volume, iScan);

        // check that radial velocity parameter is present
		if (PolarScan_hasParameter(scan, "VRAD")){
			sprintf(scanUse[iScan].vradName,"VRAD");	
			scanUse[iScan].useScan = TRUE;
		}
		else{
			if (PolarScan_hasParameter(scan, "VRADH")){
				sprintf(scanUse[iScan].vradName,"VRADH");	
				scanUse[iScan].useScan = TRUE;
			}
			else{
				if (PolarScan_hasParameter(scan, "VRADV")){
					sprintf(scanUse[iScan].vradName,"VRADV");	
					scanUse[iScan].useScan = TRUE;
				}
			}
		}
		if (scanUse[iScan].useScan == FALSE){
		  vol2bird_err_printf("Warning: radial velocity missing, dropping scan %i ...\n",iScan+1);
		}

        // check that reflectivity parameter is present
        if (scanUse[iScan].useScan){
            if(PolarScan_hasParameter(scan,alldata->options.dbzType)){
                strcpy(scanUse[iScan].dbzName,alldata->options.dbzType);	
            }
            else{
                vol2bird_err_printf("Warning: requested reflectivity factor '%s' missing, searching for alternatives ...\n",alldata->options.dbzType);
                if(PolarScan_hasParameter(scan, "DBZH")){
                    sprintf(scanUse[iScan].dbzName,"DBZH");	
                }
                else{
                    if (PolarScan_hasParameter(scan, "DBZV")){
                        sprintf(scanUse[iScan].dbzName,"DBZV");	
                    }
                    else{	
                        scanUse[iScan].useScan = FALSE;
                    }
	        }
            }
            if (scanUse[iScan].useScan == FALSE){
                vol2bird_err_printf("Warning: reflectivity factor missing, dropping scan %i ...\n",iScan+1);
            }
        }
        
        // check that correlation coefficient is present in addition to radial velocity and reflectivity
        if (scanUse[iScan].useScan & alldata->options.dualPol){
            if (PolarScan_hasParameter(scan, "RHOHV")){
                sprintf(scanUse[iScan].rhohvName,"RHOHV");	
                scanUse[iScan].useScan = TRUE;
            }
            else{
                vol2bird_err_printf("Warning: correlation coefficient missing, dropping scan %i ...\n",iScan+1);
                scanUse[iScan].useScan = FALSE;
            }
        }
        
        // check whether spectrum width parameter is present, and store its name
        strcpy(scanUse[iScan].wradName,"");
		if (PolarScan_hasParameter(scan, "WRAD")){
			sprintf(scanUse[iScan].wradName,"WRAD");	
		}
		else{
			if (PolarScan_hasParameter(scan, "WRADH")){
				sprintf(scanUse[iScan].wradName,"WRADH");	
			}
			else{
				if (PolarScan_hasParameter(scan, "WRADV")){
					sprintf(scanUse[iScan].vradName,"WRADV");	
				}
			}
		}
        
        // check that elevation is not too high or too low
        if (scanUse[iScan].useScan)
        {
            // drop scans with elevations outside the range set by user
            double elev = 360*PolarScan_getElangle(scan) / 2 / PI;
            if (elev < alldata->options.elevMin || elev > alldata->options.elevMax)
            {
                scanUse[iScan].useScan = FALSE;
                  vol2bird_err_printf("Warning: elevation (%.1f deg) outside valid elevation range (%.1f-%.1f deg), dropping scan %i ...\n",\
                    elev,alldata->options.elevMin,alldata->options.elevMax,iScan+1);
            }
        }
        
        // check that range bin size has correct units
        if (scanUse[iScan].useScan)
        {
            // drop scans with range bin sizes below 1 meter
            double rscale = PolarScan_getRscale(scan);
            if (rscale < RSCALEMIN || rscale == 0)
            {
                scanUse[iScan].useScan = FALSE;
                  vol2bird_err_printf("Warning: range bin size (%.2f metre) too small, dropping scan %i ...\n",\
                    rscale,iScan+1);
            }
        }
        
        // check that Nyquist velocity is not too low
        // retrieve Nyquist velocity if not present at scan level
        if (scanUse[iScan].useScan)
        {
            // Read Nyquist interval from scan how group
            attr = PolarScan_getAttribute(scan, "how/NI");
            result = 0;
            if (attr != (RaveAttribute_t *) NULL) result = RaveAttribute_getDouble(attr, &nyquist);
            RAVE_OBJECT_RELEASE(attr);

            // Read Nyquist interval from top level how group
            if (result == 0)
            {	
                // set flag that we found no nyquist interval at scan level
                noNyquist = 1;
                // proceed to top level how group
                attr = PolarVolume_getAttribute(volume, "how/NI");
                if (attr != (RaveAttribute_t *) NULL) result = RaveAttribute_getDouble(attr, &nyquist);
                RAVE_OBJECT_RELEASE(attr);
            }
            
            // Derive Nyquist interval from the offset attribute of the dataset
            // when no NI attribute was found
            if (result == 0)
            {
                // in case that we'll be using dealiased velocities, but also the original is available
                // we want to get the offset attribute of the original
                if(alldata->options.dealiasVrad && PolarScan_hasParameter(scan, "VRADDH") && PolarScan_hasParameter(scan, "VRADH")){
                    param = PolarScan_getParameter(scan, "VRADH");
                    nyquist = fabs(PolarScanParam_getOffset(param));
                }
                else{
                    param = PolarScan_getParameter(scan, scanUse[iScan].vradName);
                    nyquist = fabs(PolarScanParam_getOffset(param));
                }
                vol2bird_err_printf("Warning: Nyquist interval attribute not found for scan %i, using radial velocity offset (%.1f m/s) instead \n",iScan+1,nyquist);
                RAVE_OBJECT_RELEASE(param);
            }
            
            // Set useScan to 0 if no Nyquist interval is available or if it is too low
            // only check for nyquist interval when we are NOT dealiasing the velocities
            if (nyquist < alldata->options.minNyquist){
                scanUse[iScan].useScan = 0;
                vol2bird_err_printf("Warning: Nyquist velocity (%.1f m/s) too low, dropping scan %i ...\n",nyquist,iScan+1);
            }
            
            // if Nyquist interval (NI) attribute was missing at scan level, add it now
            if (noNyquist){
                RaveAttribute_t* attr_NI = RaveAttributeHelp_createDouble("how/NI", (double) nyquist);
                if (attr_NI == NULL && nyquist > 0){
                    vol2bird_err_printf("warning: no valid Nyquist attribute could be added to scan\n");
                }
                else{    
                    PolarScan_addAttribute(scan, attr_NI);
                }
                RAVE_OBJECT_RELEASE(attr_NI);
            }
            
            if (nyquist < nyquistMin){
                nyquistMin=nyquist;
            }
            if (nyquist < nyquistMinUsed && nyquist > alldata->options.minNyquist){
                nyquistMinUsed=nyquist;
            }
            if (nyquist > nyquistMax){
                nyquistMax=nyquist;
            }
            
        }
        
        if (scanUse[iScan].useScan){

            // copy names of scans that will be generated by the program
            sprintf(scanUse[iScan].texName,TEXNAME);	
            sprintf(scanUse[iScan].clutName,CLUTNAME);	
            sprintf(scanUse[iScan].cellName,CELLNAME);	

            nScansUsed+=1;
        }
        RAVE_OBJECT_RELEASE(scan);
    }
        
        
    alldata->misc.nyquistMin = nyquistMin;
    alldata->misc.nyquistMinUsed = nyquistMinUsed;
    alldata->misc.nyquistMax = nyquistMax;
    
    // FIXME: better to make scanUse a struct that contains both the array of vol2birdScanUse_t objects
    // and the number nScansUSed, which now is stored ad hoc under alldata->misc
    alldata->misc.nScansUsed = nScansUsed;
    if(nScansUsed == 0){
        alldata->misc.vol2birdSuccessful = FALSE;
        return (vol2birdScanUse_t*) NULL;
    }
    
    // Return the array scanUse
    return scanUse;
}



static void exportBirdProfileAsJSON(vol2bird_t *alldata) {
    
    // produces valid JSON according to http://jsonlint.com/
    
    if (alldata->misc.initializationSuccessful==FALSE) {
        vol2bird_err_printf("You need to initialize vol2bird before you can use it. Aborting.\n");
        return;
    }
    
    if (alldata->profiles.iProfileTypeLast != 1) {
        vol2bird_err_printf("Export method expects profile 1, but found %d. Aborting.",alldata->profiles.iProfileTypeLast);
        return; 
    }
    
    int iLayer;


    FILE *f = fopen("vol2bird-profile1.json", "w");
    if (f == NULL)
    {
        vol2bird_printf("Error opening file 'vol2bird-profile1.json'!\n");
#ifdef VOL2BIRD_R        
        return;
#else
	exit(1);
#endif        
    }
    
    fprintf(f,"[\n");
    for (iLayer = 0;iLayer < alldata->options.nLayers; iLayer += 1) {
        
        fprintf(f,"   {\n");
        
        {
            char varName[] = "HGHT";
            float val = alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile +  0];
            if (isnan(val) == TRUE) {
                fprintf(f,"    \"%s\":null,\n",varName);
            }
            else {
                fprintf(f,"    \"%s\":%.2f,\n",varName,val);
            }
        }
        
        {
            char varName[] = "HGHT_max";
            float val = alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile +  1];
            if (isnan(val) == TRUE) {
                fprintf(f,"    \"%s\":null,\n",varName);
            }
            else {
                fprintf(f,"    \"%s\":%.2f,\n",varName,val);
            }
        }
            
        {
            char varName[] = "u";
            float val = alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile +  2];
            if (isnan(val) == TRUE) {
                fprintf(f,"    \"%s\":null,\n",varName);
            }
            else {
                fprintf(f,"    \"%s\":%.2f,\n",varName,val);
            }
        }
        
        {    
            char varName[] = "v";
            float val = alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile +  3];
            if (isnan(val) == TRUE) {
                fprintf(f,"    \"%s\":null,\n",varName);
            }
            else {
                fprintf(f,"    \"%s\":%.2f,\n",varName,val);
            }
        }
        
        {    
            char varName[] = "w";
            float val = alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile +  4];
            if (isnan(val) == TRUE) {
                fprintf(f,"    \"%s\":null,\n",varName);
            }
            else {
                fprintf(f,"    \"%s\":%.2f,\n",varName,val);
            }
        }
        
        {
            char varName[] = "ff";
            float val = alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile +  5];
            if (isnan(val) == TRUE) {
                fprintf(f,"    \"%s\":null,\n",varName);
            }
            else {
                fprintf(f,"    \"%s\":%.2f,\n",varName,val);
            }
        }
            
        {
            char varName[] = "dd";
            float val = alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile +  6];
            if (isnan(val) == TRUE) {
                fprintf(f,"    \"%s\":null,\n",varName);
            }
            else {
                fprintf(f,"    \"%s\":%.2f,\n",varName,val);
            }
        }
        
        {    
            char varName[] = "sd_vvp";
            float val = alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile +  7];
            if (isnan(val) == TRUE) {
                fprintf(f,"    \"%s\":null,\n",varName);
            }
            else {
                fprintf(f,"    \"%s\":%.2f,\n",varName,val);
            }
        }
            
        {
            char varName[] = "gap";
            float val = alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile +  8];
            if (isnan(val) == TRUE) {
                fprintf(f,"    \"%s\":null,\n",varName);
            }
            else {
                fprintf(f,"    \"%s\":%s,\n",varName,val == TRUE ? "true" : "false");
            }
        }

        {            
            char varName[] = "dbz";
            float val = alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile +  9];
            if (isnan(val) == TRUE) {
                fprintf(f,"    \"%s\":null,\n",varName);
            }
            else {
                fprintf(f,"    \"%s\":%.2f,\n",varName,val);
            }
        }
            
        {
            char varName[] = "n";
            float val = alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile +  10];
            if (isnan(val) == TRUE) {
                fprintf(f,"    \"%s\":null,\n",varName);
            }
            else {
                fprintf(f,"    \"%s\":%d,\n",varName,(int) val);
            }
        }
            
        {
            char varName[] = "eta";
            float val = alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile +  11];
            if (isnan(val) == TRUE) {
                fprintf(f,"    \"%s\":null,\n",varName);
            }
            else {
                fprintf(f,"    \"%s\":%.2f,\n",varName,val);
            }
        }
        
        {    
            char varName[] = "dens";
            float val = alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile +  12];
            if (isnan(val) == TRUE) {
                fprintf(f,"    \"%s\":null\n",varName);
            }
            else {
                fprintf(f,"    \"%s\":%.2f\n",varName,val);
            }
        }

        {
            char varName[] = "n_dbz";
            float val = alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile +  13];
            if (isnan(val) == TRUE) {
                fprintf(f,"    \"%s\":null,\n",varName);
            }
            else {
                fprintf(f,"    \"%s\":%d,\n",varName,(int) val);
            }
        }

        fprintf(f,"   }");
        if (iLayer < alldata->options.nLayers - 1) {
            fprintf(f,",");
        }
        fprintf(f,"\n");
    }
    fprintf(f,"]\n");

    fclose(f);
}



static int findWeatherCells(PolarScan_t *scan, const char* quantity, float quantityThreshold,
        int selectAboveThreshold, int iCellStart, int initialize, vol2bird_t* alldata) {

    //  ----------------------------------------------------------------------------- //
    //  This function detects the cells in 'quantityImage' using a threshold value of      //
    //  'dbzThresMin' and a non-recursive algorithm which looks for neighboring       //
    //  pixels above that threshold. On return the marked cells are contained by      //
    //  'cellImage'. The number of detected cells/highest index value is returned.    //
    //  ----------------------------------------------------------------------------- //


    int iCellIdentifier;
    int nCells;
    int iRang,iRangLocal;
    int nRang;
    int iAzim, iAzimLocal;
    int nAzim;
    int iNeighborhood;
    int nNeighborhood;
    int count;
    int cellImageInitialValue;

    float quantityThres;

    int iGlobal;
    int iGlobalOther;
    int nGlobal;
    int iLocal;

    double quantityMissing;
    double quantityUndetect;
    int quantitynAzim;
    int quantitynRang;
    double quantityValueOffset;
    double quantityValueScale;
    double quantityValueGlobal, quantityValueLocal;
    double cellValueGlobal, cellValueLocal, cellValueOther;


    float quantityRangeScale;


    // These next two local variables 'nAzimNeighborhood' and 
    // 'nRangNeighborhood' shadow two static global variables of the 
    // same name, but that is in fact what you want. This is because the
    // value of 'nAzimNeighborHood' (function scope) is not equal to 
    // 'nAzimNeighborHood' (file scope) and because 'nRangNeighborHood' 
    // (function scope) is not equal to 'nRangNeighborHood' (file 
    // scope). The variable names are the same though, because they play
    // the same role in 'findNearbyGateIndex()'.
    int nAzimNeighborhood_local;
    int nRangNeighborhood_local;
    
    int nHalfNeighborhood;

    int cellIdentifierGlobal;
    int cellIdentifierGlobalOther;
    int iGlobalInner;


    #ifdef FPRINTFON
    int dbg = 0;
    #endif

    if (scan==NULL) {
        vol2bird_err_printf("Input scan is NULL.\n");
        return -1;
    }
    
    PolarScanParam_t *scanParam = PolarScan_getParameter(scan, quantity);
    PolarScanParam_t *cellParam = PolarScan_getParameter(scan, CELLNAME);

    if (scanParam == NULL || cellParam == NULL) {
        RAVE_OBJECT_RELEASE(scanParam);
        RAVE_OBJECT_RELEASE(cellParam);
        vol2bird_err_printf("%s and/or CELL quantities not found in polar scan\n", quantity);
        return -1;
    }

    // get direct pointer to the data block
    int* cellParamData = (int*) PolarScanParam_getData(cellParam);
    
    quantityMissing = PolarScanParam_getNodata(scanParam);
    quantityUndetect = PolarScanParam_getUndetect(scanParam);
    quantitynAzim = PolarScan_getNrays(scan);
    quantitynRang = PolarScan_getNbins(scan);
    quantityValueOffset = PolarScanParam_getOffset(scanParam);
    quantityValueScale = PolarScanParam_getGain(scanParam);
    quantityRangeScale = PolarScan_getRscale(scan);

    nAzim = quantitynAzim;
    nRang = quantitynRang;

    nGlobal = nAzim * nRang;
    
    // We use a neighborhood of 3x3, because in this function we are
    // interested in a cell's direct neighbors. See also comment at
    // variable definitions of 'nAzimNeighborhood' and 'nRangNeighborhood'
    // above. 
    nAzimNeighborhood_local = 3;
    nRangNeighborhood_local = 3;
    
    nNeighborhood = nAzimNeighborhood_local * nRangNeighborhood_local;
    nHalfNeighborhood = (nNeighborhood - 1)/2;


    if (scanParam != NULL) {
        quantityThres = (float) ((quantityThreshold - quantityValueOffset) / quantityValueScale);
    }

    cellImageInitialValue = CELLINIT;
	if(initialize){
		for (int iAzim = 0; iAzim < nAzim; iAzim++) {
			for (int iRang = 0; iRang < nRang; iRang++) {
				PolarScanParam_setValue(cellParam, iRang, iAzim, cellImageInitialValue);
			}
		}
	}

    // If threshold value is equal to missing value, produce a warning
    if (quantityThres == quantityMissing) {
        vol2bird_err_printf("Warning: in function findWeatherCells, quantityThres equals quantityMissing\n");
    }

    // ----------------------------------------------------------------------- //
    // Labeling of groups of connected pixels using horizontal, vertical, and  //
    // diagonal connections. The algorithm is described in 'Digital Image      //
    // Processing' by Gonzales and Woods published by Addison-Wesley.          //
    // ----------------------------------------------------------------------- //

    // Typically the first cell will have iCellIdentifier = 2, because we reserve 1 for fringe
    // to be added by function fringeCells
    iCellIdentifier = iCellStart;

    for (iAzim = 0; iAzim < nAzim; iAzim++) {
        for (iRang = 0; iRang < nRang; iRang++) {

            iGlobal = iRang + iAzim * nRang;

            if ((float)(iRang + 1) * quantityRangeScale > alldata->misc.rCellMax) {
                continue;
            }
            else {
                #ifdef FPRINTFON
                vol2bird_err_printf("iGlobal = %d\niRang + 1 = %d\n"
                "quantityRangeScale = %f\n"
                "rCellMax = %f\n"
                "(iRang + 1) * quantityRangeScale = %f\n"
                "((iRang + 1) * quantityRangeScale > rCellMax) = %d\n"
                "dbg=%d\n",iGlobal,iRang + 1,quantityRangeScale,alldata->misc.rCellMax,
                (iRang + 1) * quantityRangeScale,
                ((iRang + 1) * quantityRangeScale > alldata->misc.rCellMax),dbg);
                
                dbg++;
                #endif
            }

            #ifdef FPRINTFON
            vol2bird_err_printf("iGlobal = %d\n",iGlobal);
            #endif
            
            PolarScanParam_getValue(scanParam, iRang, iAzim, &quantityValueGlobal);
            PolarScanParam_getValue(cellParam, iRang, iAzim, &cellValueGlobal);


            if (quantityValueGlobal == quantityMissing || quantityValueGlobal == quantityUndetect) {

                #ifdef FPRINTFON
                vol2bird_err_printf("quantityImage[%d] == quantityMissing\n",iGlobal);
                #endif

                continue;
            }
            
            if (selectAboveThreshold && quantityValueGlobal < (float) quantityThres){
                continue;
            }
            if (!selectAboveThreshold && quantityValueGlobal > (float) quantityThres){
                continue;
            }

            // count number of direct neighbors above threshold
            count = 0;
            for (iNeighborhood = 0; iNeighborhood < nNeighborhood; iNeighborhood++) {

                iLocal = findNearbyGateIndex(nAzim,nRang,iGlobal,nAzimNeighborhood_local,nRangNeighborhood_local,iNeighborhood,&iAzimLocal,&iRangLocal);

                if (iLocal < 0) {
                    // iLocal below zero are error codes
                    continue;
                }
                
                PolarScanParam_getValue(scanParam, iRangLocal, iAzimLocal, &quantityValueLocal);

                if (quantityValueLocal > quantityThres) {
                    count++;
                }

            }
            // when not enough qualified neighbors, continue
            if (count - 1 < alldata->constants.nNeighborsMin) {
                continue;
            }

            #ifdef FPRINTFON
            vol2bird_err_printf("iGlobal = %d, count = %d\n",iGlobal,count);
            #endif


            // Looking for horizontal, vertical, forward diagonal, and backward diagonal connections.
            for (iNeighborhood = 0; iNeighborhood < nHalfNeighborhood; iNeighborhood++) {

                iLocal = findNearbyGateIndex(nAzim,nRang,iGlobal,nAzimNeighborhood_local,nRangNeighborhood_local,iNeighborhood,&iAzimLocal,&iRangLocal);
//                fprintf(stderr, "iLocal %i - args %i %i %i %i %i %i %i %i\n",iLocal,nAzim,nRang,iGlobal,nAzimNeighborhood_local,nRangNeighborhood_local,iNeighborhood,iAzimLocal,iRangLocal);
                if (iLocal < 0) {
                    // iLocal less than zero are error codes
                    continue;
                }
                
                PolarScanParam_getValue(cellParam, iRangLocal, iAzimLocal, &cellValueLocal);

                // no connection found, go to next pixel within neighborhood
                if ((int) cellValueLocal == cellImageInitialValue) {
                    continue;
                }

                // if pixel still unassigned, assign same iCellIdentifier as connection
                if ((int) cellValueGlobal == cellImageInitialValue) {
                    PolarScanParam_setValue(cellParam, iRang, iAzim, cellValueLocal);
                    cellValueGlobal = cellValueLocal;
                }
                else {
                    // if connection found but pixel is already assigned a different iCellIdentifier:
                    if (cellValueGlobal != cellValueLocal) {
                        // merging cells detected: replace all other occurences by value of connection:
                        for (iGlobalOther = 0; iGlobalOther < nGlobal; iGlobalOther++) {
                            if (cellParamData[iGlobalOther] == cellValueGlobal) {
                                cellParamData[iGlobalOther] = cellValueLocal;
                            }
                            // note: not all iCellIdentifier need to be used eventually
                        }
                    }
                }
            }

            // When no connections are found, give a new number.
            if ((int) cellValueGlobal == cellImageInitialValue) {

                #ifdef FPRINTFON
                vol2bird_err_printf("new cell found...assigning number %d\n",iCellIdentifier);
                #endif
                
                PolarScanParam_setValue(cellParam, iRang, iAzim, iCellIdentifier);
                iCellIdentifier++;
            }

        } // (iRang=0; iRang<nRang; iRang++)
    } // for (iAzim=0; iAzim<nAzim; iAzim++)


    // check whether a cell crosses the border of the array (remember that iAzim=0 is
    // adjacent to iAzim=nAzim-1):
    iAzim = 0;

    for (iRang = 0; iRang < nRang; iRang++) {

        iGlobal = iRang + iAzim * nRang;
        // index 1 in a 3x3 child array refers to the cell that is a direct neighbor of
        // iGlobal, but on the other side of the array (because the polar plot is wrapped
        // in the azimuth dimension):
        iGlobalOther = findNearbyGateIndex(nAzim,nRang,iGlobal,3,3,1,&iAzimLocal,&iRangLocal);
        PolarScanParam_getValue(cellParam, iRangLocal, iAzimLocal, &cellValueOther);

        #ifdef FPRINTFON
        vol2bird_err_printf("iGlobal = %d, iGlobalOther = %d\n",iGlobal,iGlobalOther);
        #endif

        cellIdentifierGlobal = cellValueGlobal;
        cellIdentifierGlobalOther = cellValueOther;
        if (cellIdentifierGlobal != cellImageInitialValue && cellIdentifierGlobalOther != cellImageInitialValue ) {
            // adjacent gates, both part of a cell -> assign them the same identifier, i.e. assign
            // all elements of cellImage that are equal to cellImage[iGlobalOther] the value of
            // cellImage[iGlobal]

            for (iGlobalInner = 0; iGlobalInner < nGlobal; iGlobalInner++) {
                if (cellParamData[iGlobalInner] == cellIdentifierGlobalOther) {
                    cellParamData[iGlobalInner] = cellIdentifierGlobal;
                }
            }
        }
    }

    // Returning number of detected cells (including fringe/clutter)
    nCells = iCellIdentifier;

    RAVE_OBJECT_RELEASE(scanParam);
    RAVE_OBJECT_RELEASE(cellParam);

    return nCells;
} // findWeatherCells



static int findNearbyGateIndex(const int nAzimParent, const int nRangParent, const int iParent,
                        const int nAzimChild,  const int nRangChild,  const int iChild, int *iAzimReturn, int *iRangReturn) {



    if (nRangChild%2 != 1) {

        #ifdef FPRINTFON
        vol2bird_err_printf("nRangChild must be an odd integer number.\n");
        #endif

        return -1;
    }

    if (nAzimChild%2 != 1) {

        #ifdef FPRINTFON
        vol2bird_err_printf("nAzimChild must be an odd integer number.\n");
        #endif

        return -2;
    }

    if (iChild > nAzimChild * nRangChild - 1) {

        #ifdef FPRINTFON
        vol2bird_err_printf("iChild is outside the child window.\n");
        #endif

        return -3;
    }


    const int iAzimParent = iParent / nRangParent;
    const int iRangParent = iParent % nRangParent;

    const int iAzimChild = iChild / nRangChild;
    const int iRangChild = iChild % nRangChild;

    // the azimuth dimension is wrapped (polar plot); iAzim = 0 is adjacent to iAzim = nAzim-1:
    *iAzimReturn = (iAzimParent - nAzimChild/2 + iAzimChild + nAzimParent) % nAzimParent;
    *iRangReturn = iRangParent - nRangChild/2 + iRangChild;


    if (*iRangReturn > nRangParent - 1) {

        #ifdef FPRINTFON
        vol2bird_err_printf("iChild is outside the parent array on the right-hand side.\n");
        #endif

        return -4;
    }
    if (*iRangReturn < 0) {

        #ifdef FPRINTFON
        vol2bird_err_printf("iChild is outside the parent array on the left-hand side.\n");
        #endif

        return -5;
    }
    
    return *iAzimReturn * nRangParent + *iRangReturn;
    
} // findNearbyGateIndex


static void fringeCells(PolarScan_t* scan, vol2bird_t* alldata) {

    // -------------------------------------------------------------------------- //
    // This function enlarges cells in cellImage by an additional fringe.         //
    // First a block around each pixel is searched for pixels within a distance   //
    // equal to 'fringeDist'.                                                     //
    // -------------------------------------------------------------------------- //

    if(!PolarScan_hasParameter(scan, CELLNAME)){
        vol2bird_err_printf("no CELL quantity in polar scan, aborting fringeCells()\n");
        return;
    }

    int nRang = (int) PolarScan_getNbins(scan);
    int nAzim = (int) PolarScan_getNrays(scan);
    float aScale = 360.0f/ nAzim;
    float rScale = (float) PolarScan_getRscale(scan);
    PolarScanParam_t *cellParam = PolarScan_getParameter(scan,CELLNAME);
    int *cellImage = (int *) PolarScanParam_getData(cellParam);

    int iRang;
    int iAzim;
    int iNeighborhood;
    int nNeighborhood;
    int iRangLocal;
    int iAzimLocal;
    int iLocal;
    int rBlock;
    int aBlock;
    int isEdge;
    int iGlobal;
    float theDist;
    int nAzimChild;
    int nRangChild;

    float actualRange;
    float circumferenceAtActualRange;


    rBlock = ROUND(alldata->constants.fringeDist / rScale);
    for (iAzim = 0; iAzim < nAzim; iAzim++) {
        for (iRang = 0; iRang < nRang; iRang++) {

            iGlobal = iRang + iAzim * nRang;

            if (cellImage[iGlobal] <= 1) {
                continue; // with the next iGlobal; already fringe or not in cellImage
            }

            // determine whether current pixel is a pixel on the edge of a cell
            isEdge = FALSE;
            for (iNeighborhood = 0; iNeighborhood < 9 && isEdge == FALSE; iNeighborhood++) {

                iLocal = findNearbyGateIndex(nAzim,nRang,iGlobal,3,3,iNeighborhood,&iAzimLocal,&iRangLocal);

                if (iLocal < 0) {
                    // iLocal less than zero are error codes
                    continue; // with the next iGlobal
                }

                if (cellImage[iLocal] < 1) { //FIXME: Strictly this should be <=1, but then much slower
                    isEdge = TRUE;           //Now only pixels without any bordering fringe are 'fringed', giving very similar result (but improvement welcome)
                }

            }
            
            if (isEdge == FALSE) {
                continue; // with the next iGlobal
            }

            actualRange = (iRang+0.5) * rScale;
            circumferenceAtActualRange = 2 * PI * actualRange;
            aBlock = (alldata->constants.fringeDist / circumferenceAtActualRange) * nAzim;


            #ifdef FPRINTFON
            vol2bird_err_printf("actualRange = %f\n",actualRange);
            vol2bird_err_printf("circumferenceAtActualRange = %f\n",circumferenceAtActualRange);
            vol2bird_err_printf("fringeDist / circumferenceAtActualRange = %f\n",alldata->constants.fringeDist / circumferenceAtActualRange);
            vol2bird_err_printf("aBlock = %d\n", aBlock);
            vol2bird_err_printf("rBlock = %d\n", rBlock);
            #endif

            nAzimChild = 2 * aBlock + 1;
            nRangChild = 2 * rBlock + 1;
            nNeighborhood = nAzimChild * nRangChild;

            for (iNeighborhood = 0; iNeighborhood < nNeighborhood; iNeighborhood++) {

                iLocal = findNearbyGateIndex(nAzim,nRang,iGlobal,nAzimChild,nRangChild,iNeighborhood,&iAzimLocal,&iRangLocal);

                if (iLocal < 0) {
                    // iLocal less than zero are error codes
                    continue;
                }

//                iAzimLocal = iLocal / nRang;
//                iRangLocal = iLocal % nRang;

                // if not within range or already in cellImage or already a fringe, do nothing
                theDist = calcDist(iRang, iAzim, iRangLocal, iAzimLocal, rScale, aScale);
                if (theDist > alldata->constants.fringeDist || cellImage[iLocal] >= 1) {
                    continue; // with the next iGlobal
                }
                // include pixel (iRangLocal,iAzimLocal) in fringe
                cellImage[iLocal] = 1;
                
            } // (iNeighborhood = 0; iNeighborhood < nNeighborhood; iNeighborhood++)
        } // (iRang = 0; iRang < nRang; iRang++)
    } // (iAzim = 0; iAzim < nAzim; iAzim++)

    RAVE_OBJECT_RELEASE(cellParam);
    return;

} // fringeCells


CELLPROP* getCellProperties(PolarScan_t* scan, vol2birdScanUse_t scanUse, const int nCells, vol2bird_t* alldata){    
    int iCell;
    int iGlobal;
    int iRang;
    int iAzim;
    long nAzim;
    long nRang;
    double rScale;
    double aScale;
    float nGatesValid;
    double dbzValue;
    double texValue = 0;
    double vradValue;
    double clutterValue = alldata->options.clutterValueMin;
    double cellValue;
    
    PolarScanParam_t *dbzParam = PolarScan_getParameter(scan,scanUse.dbzName);
    PolarScanParam_t *vradParam = PolarScan_getParameter(scan,scanUse.vradName);
    PolarScanParam_t *texParam = PolarScan_getParameter(scan,scanUse.texName);
    PolarScanParam_t *cellParam = PolarScan_getParameter(scan,scanUse.cellName);
    PolarScanParam_t *clutParam = PolarScan_getParameter(scan,scanUse.clutName);

    nRang = PolarScan_getNbins(scan);
    nAzim = PolarScan_getNrays(scan);
    rScale = PolarScan_getRscale(scan);
    aScale = (360.0/nAzim)*PI/180; // in radials
    
    // Allocating and initializing memory for cell properties.
    CELLPROP* cellProp;
    cellProp = (CELLPROP *)malloc(nCells*sizeof(CELLPROP));
    
    if (!cellProp) {
        vol2bird_err_printf("Requested memory could not be allocated in getCellProperties!\n");
        return cellProp;
    }
    
    for (iCell = 0; iCell < nCells; iCell++) {
        cellProp[iCell].iRangOfMax = -1;
        cellProp[iCell].iAzimOfMax = -1;
        cellProp[iCell].nGates = 0;
        cellProp[iCell].nGatesClutter = 0;
        cellProp[iCell].area = 0;
        cellProp[iCell].dbzAvg = NAN;
        cellProp[iCell].texAvg = NAN;
        cellProp[iCell].dbzMax = NAN;
        cellProp[iCell].index = iCell;
        cellProp[iCell].drop = TRUE;
        cellProp[iCell].cv = NAN;
    }

    // Calculation of cell properties.
    RaveValueType typeDbz, typeVrad, typeTex, typeCell;
    typeTex = RaveValueType_DATA;
    for (iAzim = 0; iAzim < nAzim; iAzim++) {
        for (iRang = 0; iRang < nRang; iRang++) {

            iGlobal = iRang + iAzim * nRang;

            typeDbz=PolarScanParam_getConvertedValue(dbzParam, iRang, iAzim, &dbzValue);
            typeVrad=PolarScanParam_getConvertedValue(vradParam, iRang, iAzim, &vradValue);
            if (clutParam != NULL) PolarScanParam_getConvertedValue(clutParam, iRang, iAzim, &clutterValue);
            if (texParam != NULL) typeTex=PolarScanParam_getConvertedValue(texParam, iRang, iAzim, &texValue);
            typeCell = PolarScanParam_getConvertedValue(cellParam, iRang, iAzim, &cellValue);
	    
            iCell = (int) cellValue;

            // Note: this also throws out all nodata/undetect values for dbzValue
            if (iCell<0 || typeCell != RaveValueType_DATA) {
                continue;
            }

            #ifdef FPRINTFON
            vol2bird_err_printf("dbzValue = %f; vradValue = %f; clutterValue = %f; texValue = %f\n",dbzValue,vradValue,clutterValue,texValue);
            vol2bird_err_printf("iGlobal = %d, iCell = %d\n",iGlobal,iCell);
            #endif

            cellProp[iCell].nGates += 1;
            cellProp[iCell].area += rScale * rScale * iRang * sin(aScale)/(1000*1000);
            cellProp[iCell].drop = FALSE;

            // low radial velocities are treated as clutter, not included in calculation cell properties
            if ((fabs(vradValue) < alldata->constants.vradMin) && (typeVrad == RaveValueType_DATA)){

                cellProp[iCell].nGatesClutter += 1;

                #ifdef FPRINTFON
                vol2bird_err_printf("iGlobal = %d: vrad too low...treating as clutter\n",iGlobal);
                #endif

                continue;
            }

            if (typeVrad != RaveValueType_DATA || typeDbz != RaveValueType_DATA || typeTex != RaveValueType_DATA){

                cellProp[iCell].nGatesClutter += 1;

                continue;
            }


            // pixels in clutter map not included in calculation cell properties
            if (alldata->options.useClutterMap == TRUE){
                if (clutterValue > alldata->options.clutterValueMin){
                    cellProp[iCell].nGatesClutter += 1;
                    continue;
                }
            }

            if (isnan(cellProp[iCell].dbzMax) || (dbzValue > cellProp[iCell].dbzMax)) {

                #ifdef FPRINTFON
                vol2bird_err_printf("%d: new dbzMax value of %f found for this cell (%d).\n",iGlobal,dbzValue,iCell);
                #endif

                cellProp[iCell].dbzMax = dbzValue;
                cellProp[iCell].iRangOfMax = iGlobal%nRang;
                cellProp[iCell].iAzimOfMax = iGlobal/nRang;
            }

            if (isnan(cellProp[iCell].dbzAvg)) {
                cellProp[iCell].dbzAvg = dbzValue;
            } 
            else {
                // note that we average logarithmic values - done like this
                // to limit the contribution of high reflectivity outliers to the average
                cellProp[iCell].dbzAvg += dbzValue;
            }

            if (isnan(cellProp[iCell].texAvg)) {
                cellProp[iCell].texAvg = texValue;
            } 
            else {
                cellProp[iCell].texAvg += texValue;
            }
            
        } // for (iRang = 0; iRang < nRang; iRang++)
    } // for (iAzim = 0; iAzim < nAzim; iAzim++)


    for (iCell = 0; iCell < nCells; iCell++) {
        nGatesValid = cellProp[iCell].nGates - cellProp[iCell].nGatesClutter;
        if (nGatesValid > 0){
            cellProp[iCell].dbzAvg /= nGatesValid;
            cellProp[iCell].texAvg /= nGatesValid;
            cellProp[iCell].cv = cellProp[iCell].texAvg / cellProp[iCell].dbzAvg;
        }
    }
    
    RAVE_OBJECT_RELEASE(dbzParam);
    RAVE_OBJECT_RELEASE(vradParam);
    RAVE_OBJECT_RELEASE(texParam);
    RAVE_OBJECT_RELEASE(cellParam);
    RAVE_OBJECT_RELEASE(clutParam);

    return cellProp;
} // getCellProperties



static int getListOfSelectedGates(PolarScan_t* scan, vol2birdScanUse_t scanUse, const float altitudeMin, const float altitudeMax,
                           float* points_local, int iRowPoints, int nColsPoints_local, vol2bird_t* alldata) {

    // ------------------------------------------------------------------- //
    // Write combinations of an azimuth angle, an elevation angle, an      // 
    // observed vrad value, an observed dbz value, and a cell identifier   //
    // value into an external larger list.                                 //
    // ------------------------------------------------------------------- //

    int iAzim;
    int iRang;
    int nRang;
    int nAzim;
    int nPointsWritten_local;

    float gateHeight;
    float gateRange;
    float gateAzim;
    float rangeScale;
    float azimuthScale;
    float elevAngle;
    float radarHeight;
    double vradValue;
    double dbzValue;
    double cellValue;
    double clutValue = NAN;
    double nyquist = 0;
    
    nPointsWritten_local = 0;
    
    RaveValueType vradValueType, dbzValueType;
    
    nRang = (int) PolarScan_getNbins(scan);
    nAzim = (int) PolarScan_getNrays(scan);
    rangeScale = (float) PolarScan_getRscale(scan);
    azimuthScale = 360.0f/nAzim;
    elevAngle = (float) PolarScan_getElangle(scan);
    radarHeight = (float) PolarScan_getHeight(scan);
    RaveAttribute_t* attr = PolarScan_getAttribute(scan, "how/NI");
    if (attr != (RaveAttribute_t *) NULL){ 
        RaveAttribute_getDouble(attr, &nyquist);
    }
    RAVE_OBJECT_RELEASE(attr);
    
    PolarScanParam_t* vradParam = PolarScan_getParameter(scan,scanUse.vradName);
    PolarScanParam_t* dbzParam = PolarScan_getParameter(scan,scanUse.dbzName);
    PolarScanParam_t* cellParam = PolarScan_getParameter(scan,scanUse.cellName);
    PolarScanParam_t* clutParam = NULL;

    if (alldata->options.useClutterMap){
        clutParam = PolarScan_getParameter(scan,scanUse.clutName);
    }

    for (iRang = 0; iRang < nRang; iRang++) {

        // so gateRange represents a distance along the view direction (not necessarily horizontal)
        gateRange = ((float) iRang + 0.5f) * rangeScale;

        // note that "sin(elevAngle*DEG2RAD)" is equivalent to = "cos((90 - elevAngle)*DEG2RAD)":
        gateHeight = range2height(gateRange, elevAngle) + radarHeight;

        if (gateRange < alldata->options.rangeMin || gateRange > alldata->options.rangeMax) {
            // the current gate is either 
            // (1) too close to the radar; or
            // (2) too far away.
            continue;
        }
        if (gateHeight < altitudeMin || gateHeight > altitudeMax) {
            // if the height of the middle of the current gate is too far away from
            // the requested height, continue with the next gate
            continue;
        }

        // the gates at this range and elevation angle are within bounds,
        // include their data in the 'points' array:

        for (iAzim = 0; iAzim < nAzim; iAzim++) {

            gateAzim = ((float) iAzim + 0.5f) * azimuthScale;
            vradValueType = PolarScanParam_getConvertedValue(vradParam, iRang, iAzim, &vradValue);
            dbzValueType = PolarScanParam_getConvertedValue(dbzParam, iRang, iAzim, &dbzValue);
            PolarScanParam_getValue(cellParam, iRang, iAzim, &cellValue);
            if (alldata->options.useClutterMap){
                PolarScanParam_getValue(clutParam, iRang, iAzim, &clutValue);
            }

            // in the points array, store missing reflectivity values as the lowest possible reflectivity
            // this is to treat undetects as absence of scatterers
            if (dbzValueType != RaveValueType_DATA){
                dbzValue = NAN;
            }

            // in the points array, store missing vrad values as NAN
            // this is necessary because different scans may have different missing values
            if (vradValueType != RaveValueType_DATA){
                vradValue = NAN;
            }

            // store the location as a range, azimuth angle, elevation angle combination
            points_local[iRowPoints * nColsPoints_local + alldata->points.rangeCol] = gateRange;
            points_local[iRowPoints * nColsPoints_local + alldata->points.azimAngleCol] = gateAzim;
            points_local[iRowPoints * nColsPoints_local + alldata->points.elevAngleCol] = elevAngle * RAD2DEG;

            // also store the dbz value --useful when estimating the bird density
            points_local[iRowPoints * nColsPoints_local + alldata->points.dbzValueCol] = (float) dbzValue;
            
            // store the corresponding observed vrad value
            points_local[iRowPoints * nColsPoints_local + alldata->points.vradValueCol] = (float) vradValue;

            // store the corresponding cellImage value
            points_local[iRowPoints * nColsPoints_local + alldata->points.cellValueCol] = (float) cellValue;

            // set the gateCode to zero for now
            points_local[iRowPoints * nColsPoints_local + alldata->points.gateCodeCol] = (float) 0;

            // store the corresponding observed nyquist velocity
            points_local[iRowPoints * nColsPoints_local + alldata->points.nyquistCol] = (float) nyquist;

            // store the corresponding observed vrad value for now (to be dealiased later)
            points_local[iRowPoints * nColsPoints_local + alldata->points.vraddValueCol] = (float) vradValue;

            // store the corresponding observed vrad value for now (to be dealiased later)
            points_local[iRowPoints * nColsPoints_local + alldata->points.clutValueCol] = (float) clutValue;

            // raise the row counter by 1
            iRowPoints += 1;
            
            // raise number of points written by 1
            nPointsWritten_local += 1;

        }  //for iAzim
    } //for iRang

    RAVE_OBJECT_RELEASE(vradParam);
    RAVE_OBJECT_RELEASE(dbzParam);
    RAVE_OBJECT_RELEASE(cellParam);
    RAVE_OBJECT_RELEASE(clutParam);

    return nPointsWritten_local;
} // getListOfSelectedGates


int vol2birdLoadClutterMap(PolarVolume_t* volume, char* file, float rangeMax){
    PolarVolume_t* clutVol = NULL;

    clutVol = vol2birdGetVolume(&file, 1, rangeMax,1);
            
    if(clutVol == NULL){
        vol2bird_err_printf( "Error: function loadClutterMap: failed to load file '%s'\n",file);
        return -1;
    }
    
    int nClutScans = PolarVolume_getNumberOfScans(clutVol);

    if(nClutScans < 1){
        vol2bird_err_printf( "Error: function loadClutterMap: no clutter map data found in file '%s'\n",file);
        RAVE_OBJECT_RELEASE(clutVol);
        return -1;
    }

    // iterate over the scans in 'volume'
    int iScan;
    int nScans;
    
    // determine how many scan elevations the volume object contains
    nScans = PolarVolume_getNumberOfScans(volume);

    for (iScan = 0; iScan < nScans; iScan++) {
        // initialize the scan object
        PolarScan_t* scan = NULL;
        PolarScan_t* clutScan = NULL;
        PolarScanParam_t* param = NULL;
        PolarScanParam_t* param_proj = NULL;
        
        int result;
        double elev, rscale;
        
        // extract the scan object from the volume object
        scan = PolarVolume_getScan(volume, iScan);
        
        // extract the cluttermap scan parameter closest in elevation
        elev = PolarScan_getElangle(scan);

        //FIXME: here PolarVolume_getScanClosestToElevation_vol2bird leads to a segmentation fault. It finds the correct scan, but
        //any operation here on the pointer leads to segfault...
        //clutScan = PolarVolume_getScanClosestToElevation_vol2bird(clutVol,elev);
        clutScan = PolarVolume_getScanClosestToElevation(clutVol,elev,0);

        param = PolarScan_getParameter(clutScan,CLUTNAME);
        
        if(param == NULL){
            vol2bird_err_printf( "Error in loadClutterMap: no scan parameter %s found in file %s\n", CLUTNAME,file);
            RAVE_OBJECT_RELEASE(scan);
            RAVE_OBJECT_RELEASE(clutScan);
            RAVE_OBJECT_RELEASE(clutVol);
            return -1;
        }
        
        // project the clutter map scan parameter to the correct dimensions
        rscale = PolarScan_getRscale(clutScan);
        param_proj = PolarScanParam_project_on_scan(param, scan, rscale);
        
        // add the clutter map scan parameter to the polar volume
        result = PolarScan_addParameter(scan, param_proj);

        if(result == 0){
            vol2bird_err_printf( "Warning in loadClutterMap: failed to add cluttermap for scan %i\n",iScan+1);
        }
        
        RAVE_OBJECT_RELEASE(scan);
        RAVE_OBJECT_RELEASE(clutScan);
        RAVE_OBJECT_RELEASE(param);
        RAVE_OBJECT_RELEASE(param_proj);
    }
    
    RAVE_OBJECT_RELEASE(clutVol);
    
    return 0;
}


// adds a scan parameter to the scan
PolarScanParam_t* PolarScan_newParam(PolarScan_t *scan, const char *quantity, RaveDataType type){
    if (scan == NULL){
        vol2bird_err_printf( "error in PolarScan_newParam(): cannat create a new polar scan parameter for scan NULL pointer\n");
        return NULL;        
    }
    
    if (PolarScan_hasParameter(scan, quantity)){
        vol2bird_err_printf( "Parameter %s already exists in polar scan\n", quantity);
        return NULL;
    }

    PolarScanParam_t *scanParam = RAVE_OBJECT_NEW(&PolarScanParam_TYPE);

    if (scanParam == NULL){
        vol2bird_err_printf( "failed to allocate memory for new polar scan parameter\n");
        RAVE_OBJECT_RELEASE(scanParam);
        return NULL;
    }
    
    PolarScanParam_createData(scanParam, PolarScan_getNbins(scan), PolarScan_getNrays(scan), type);
    PolarScanParam_setQuantity(scanParam,quantity);
    PolarScanParam_setNodata(scanParam,NODATA);
    PolarScanParam_setUndetect(scanParam,UNDETECT);
    PolarScanParam_setOffset(scanParam,0);
    PolarScanParam_setGain(scanParam,1);
    
    // initialize all values to NODATA
    // (NOTE: this fails for negative values when type == RaveDataType_FLOAT, results in 0 being set. RAVE bug?)
    double nodata = PolarScanParam_getNodata(scanParam);
    for(int iRang = 0; iRang < PolarScan_getNbins(scan); iRang++){
        for(int iAzim = 0; iAzim < PolarScan_getNrays(scan); iAzim++){
            
            PolarScanParam_setValue(scanParam, iRang, iAzim, nodata);
        }
    }
        
    PolarScan_addParameter(scan, scanParam);
    
 //   RAVE_OBJECT_RELEASE(scanParam);
    
 //   scanParam = PolarScan_getParameter(scan, quantity);
    
    return scanParam;
}


long datetime2long(char* date, char* time){

    //concatenate date and time into a datetime string
    char *datetime = malloc(strlen(date)+strlen(time)+1);
    long result;
    strcpy(datetime,date);
    strcat(datetime,time);

    //convert datetime string to a decimal long
    char *eptr;
    long ldatetime;
    ldatetime = strtol(datetime, &eptr, 10);

    // check for conversion errors
    if (ldatetime == 0){
        #ifdef FPRINTFON
            vol2bird_err_printf("Conversion error occurred\n");
        #endif
        result = 0;
    }
    else{
        result = ldatetime;
    }
    
    free(datetime);
    
    return(result);
}

const char* PolarVolume_getStartDate(PolarVolume_t* pvol){
    RAVE_ASSERT((pvol != NULL), "pvol == NULL");
    char* date = (char *) PolarVolume_getDate(pvol);
    char* time = (char *) PolarVolume_getTime(pvol);
    char* result = (char*) NULL;
    if(PolarVolume_getStartDateTime(pvol, &date, &time)==0){
        result = date;
    }
    return(result);  
}

const char* PolarVolume_getStartTime(PolarVolume_t* pvol){
    RAVE_ASSERT((pvol != NULL), "pvol == NULL");
    char* date = (char *) PolarVolume_getDate(pvol);
    char* time = (char *) PolarVolume_getTime(pvol);
    char* result = (char*) NULL;
    if(PolarVolume_getStartDateTime(pvol, &date, &time)==0){
        result = time;
    }
    return(result);  
}

const char* PolarVolume_getEndDate(PolarVolume_t* pvol){
    RAVE_ASSERT((pvol != NULL), "pvol == NULL");
    char* date = (char *) PolarVolume_getDate(pvol);
    char* time = (char *) PolarVolume_getTime(pvol);
    char* result = (char*) NULL;
    if(PolarVolume_getEndDateTime(pvol, &date, &time)==0){
        result = date;
    }
    return(result);  
}

const char* PolarVolume_getEndTime(PolarVolume_t* pvol){
    RAVE_ASSERT((pvol != NULL), "pvol == NULL");
    char* date = (char *) PolarVolume_getDate(pvol);
    char* time = (char *) PolarVolume_getTime(pvol);
    char* result = (char*) NULL;
    if(PolarVolume_getEndDateTime(pvol, &date, &time)==0){
        result = time;
    }
    return(result);  
}

int PolarVolume_getStartDateTime(PolarVolume_t* pvol, char** StartDate, char** StartTime)
{
    RAVE_ASSERT((pvol != NULL), "pvol == NULL");
    
    int result = -1;
    
    // Initialize datetimes
    long StartDateTime = (long) 99999999999999;
    
    int nScans;
    char* date;
    char* time;

    // Read number of scans
    nScans = PolarVolume_getNumberOfScans(pvol);
    
    // find the start date and time
    for (int iScan = 0; iScan < nScans; iScan++)
    {
        PolarScan_t* scan = PolarVolume_getScan(pvol, iScan);
        // useScan is set to zero if the quantity is not found in the scan
        if (scan != (PolarScan_t *) NULL){
            date = (char *) PolarScan_getStartDate(scan);
            time = (char *) PolarScan_getStartTime(scan);

            long datetime = datetime2long(date, time);
            
            //continue if no valid datetime can be constructed
            if (datetime == 0){
                RAVE_OBJECT_RELEASE(scan);
                continue;
            }
            
            if (datetime < StartDateTime){
                StartDateTime = datetime;
                *StartDate = date;
                *StartTime = time;
                // success, we found a valid start date time
                result = 0;
            }
        }
        RAVE_OBJECT_RELEASE(scan);
    }
    
    return result;
}

int PolarVolume_getEndDateTime(PolarVolume_t* pvol, char** EndDate, char** EndTime)
{
    RAVE_ASSERT((pvol != NULL), "pvol == NULL");

    int result = -1;

    // Initialize datetimes
    long EndDateTime = 00000000000000;
    
    int nScans;
    char* date;
    char* time;

    // Read number of scans
    nScans = PolarVolume_getNumberOfScans(pvol);
    
     // find the end date and time
    for (int iScan = 0; iScan < nScans; iScan++)
    {
        PolarScan_t* scan = PolarVolume_getScan(pvol, iScan);
        // useScan is set to zero if the quantity is not found in the scan
        if (scan != (PolarScan_t *) NULL){
            date = (char *) PolarScan_getEndDate(scan);
            time = (char *) PolarScan_getEndTime(scan);

            long datetime = datetime2long(date, time);
            
            //continue if no valid datetime can be constructed
            if ((date == NULL || time == NULL || datetime == 0)){
                RAVE_OBJECT_RELEASE(scan);
                continue;
            }
            
            if (datetime > EndDateTime){
                EndDateTime = datetime;
                *EndDate = date;
                *EndTime = time;
                // success, we found a valid start date time
                result = 0;
            }
        }
        RAVE_OBJECT_RELEASE(scan);
    }
    return result;
}


double PolarVolume_getWavelength(PolarVolume_t* pvol)
{
    RAVE_ASSERT((pvol != NULL), "pvol == NULL");
    
    double value = 0;
    int speed_of_light = 299792458;

    RaveAttribute_t* attr = PolarVolume_getAttribute(pvol, "how/wavelength");
    if (attr != (RaveAttribute_t *) NULL){
        RaveAttribute_getDouble(attr, &value);
    }
    else{
        attr = PolarVolume_getAttribute(pvol, "how/frequency");
        if (attr != (RaveAttribute_t *) NULL){
            RaveAttribute_getDouble(attr, &value);
            // convert frequency in Hz to wavelength in cm
            value = 100*speed_of_light/value;
        }
        else{
            // wavelength attribute was not found in the root /how attribute
            // check whether we can find it under /dataset1/how 
            PolarScan_t* scan = PolarVolume_getScan(pvol, 1);
            if (scan != (PolarScan_t *) NULL){
                attr = PolarScan_getAttribute(scan, "how/wavelength");
                if (attr != (RaveAttribute_t *) NULL){
                    RaveAttribute_getDouble(attr, &value);
                    vol2bird_err_printf( "Warning: using radar wavelength stored for scan 1 (%f cm) for all scans ...\n", value);
                }
                else{
                    attr = PolarScan_getAttribute(scan, "how/frequency");
                    if (attr != (RaveAttribute_t *) NULL){
                        RaveAttribute_getDouble(attr, &value);
                        // convert frequency in Hz to wavelength in cm
                        value = 100*speed_of_light/value;
                        vol2bird_err_printf( "Warning: using radar frequency stored for scan 1 (%f Hz) for all scans ...\n", value);
                    }
                }
            }
            RAVE_OBJECT_RELEASE(scan);
        }
    }
    RAVE_OBJECT_RELEASE(attr);
    return value;
}

double PolarVolume_setWavelength(PolarVolume_t* pvol, double wavelength)
{
    RAVE_ASSERT((pvol != NULL), "pvol == NULL");
    
    double value = 0;

    RaveAttribute_t* attr = PolarVolume_getAttribute(pvol, "how/wavelength");
    if (attr != (RaveAttribute_t *) NULL){
        RaveAttribute_getDouble(attr, &value);
    }
    else{
        // wavelength attribute was not found in the root /how attribute
        // check whether we can find it under /dataset1/how 
        PolarScan_t* scan = PolarVolume_getScan(pvol, 1);
        if (scan != (PolarScan_t *) NULL){
            attr = PolarScan_getAttribute(scan, "how/wavelength");
            if (attr != (RaveAttribute_t *) NULL){
                RaveAttribute_getDouble(attr, &value);
                vol2bird_err_printf( "Warning: using radar wavelength stored for scan 1 (%f cm) for all scans ...\n", value);
            }
        }
        RAVE_OBJECT_RELEASE(scan);
    }
    RAVE_OBJECT_RELEASE(attr);
    return value;
}

PolarScanParam_t* PolarScanParam_project_on_scan(PolarScanParam_t* param, PolarScan_t* scan, double rscale){
    PolarScanParam_t* param_proj = NULL;
    
    double rscale_proj = PolarScan_getRscale(scan);
    long nbins_proj = PolarScan_getNbins(scan);
    long nrays_proj = PolarScan_getNrays(scan);
    
    param_proj = PolarScanParam_resample(param, rscale, rscale_proj, nbins_proj, nrays_proj);
    
    return param_proj;
}

PolarVolume_t* PolarVolume_resample(PolarVolume_t* volume, double rscale_proj, long nbins_proj, long nrays_proj){
    int iScan;
    int nScans;
    
    nScans = PolarVolume_getNumberOfScans(volume);
    
    PolarScan_t* scan = NULL;
    PolarScan_t* scan_proj = NULL;

    // copy the volume
    PolarVolume_t* volume_proj = RAVE_OBJECT_CLONE(volume);
        
    // empty the scans in the copied volume
    for (iScan = nScans-1; iScan>=0 ; iScan--) { 
        PolarVolume_removeScan(volume_proj, iScan);
    }
   
    // iterate over the scans in source volume
    for (iScan = 0; iScan < nScans; iScan++) {
        scan = PolarVolume_getScan(volume, iScan);
        scan_proj = PolarScan_resample(scan, rscale_proj, nbins_proj, nrays_proj);
        PolarVolume_addScan(volume_proj, scan_proj);
        RAVE_OBJECT_RELEASE(scan_proj);
        RAVE_OBJECT_RELEASE(scan);
    }
    
    return volume_proj;
}

PolarScan_t* PolarScan_resample(PolarScan_t* scan, double rscale_proj, long nbins_proj, long nrays_proj){
    int iParam;
    
    // determine how many parameters the scan contains
    RaveList_t* ParamNames = PolarScan_getParameterNames(scan);
    int nParams = RaveList_size(ParamNames);

    // initialize the scan object
    PolarScan_t* scan_proj = NULL;
    PolarScanParam_t* param = NULL;
    PolarScanParam_t* param_proj = NULL;

    scan_proj = RAVE_OBJECT_CLONE(scan);
    PolarScan_removeAllParameters(scan_proj);
        
    double rscale = PolarScan_getRscale(scan);
    long nbins = PolarScan_getNbins(scan);
    long nrays = PolarScan_getNrays(scan);
    double elev = PolarScan_getElangle(scan)*180/PI;
    
    if (rscale > rscale_proj){
        vol2bird_err_printf( "Warning: requested range gate size (rscale=%3.1f m) too small for %2.1f degree scan, using %4.1f m\n", rscale_proj, elev, rscale);
        rscale_proj = rscale;
    }

    if (nbins < nbins_proj){
        vol2bird_err_printf( "Warning: requested number of range bins (Nbins=%li) too large for %3.1f degree scan, using %li bins\n", nbins_proj, elev, nbins);
        nbins_proj = nbins;
    }

    if (nrays < nrays_proj){
        vol2bird_err_printf( "Warning: requested number of azimuth rays (Nrays=%li) too large for %3.1f degree scan, using %li rays\n", nrays_proj, elev, nrays);
        nrays_proj = nrays;
    }
    
    // update scan object with new rscale
    PolarScan_setRscale(scan_proj, rscale_proj);
    
    // iterate over the parameters in scan
    for (iParam = 0; iParam < nParams; iParam++) {
        // extract the scan object from the volume object
        param = PolarScan_getParameter(scan, RaveList_get(ParamNames,iParam));
        // resample the parameter
        param_proj = PolarScanParam_resample(param, rscale, rscale_proj, nbins_proj, nrays_proj);
        // add parameter to scan
        PolarScan_addParameter(scan_proj, param_proj);
        // release the parameter
        RAVE_OBJECT_RELEASE(param);
        RAVE_OBJECT_RELEASE(param_proj);
    }
    RAVE_OBJECT_RELEASE(param);
    RAVE_OBJECT_RELEASE(param_proj);

    RaveList_freeAndDestroy(&ParamNames);
    
    return scan_proj;
}

PolarScanParam_t* PolarScanParam_resample(PolarScanParam_t* param, double rscale, double rscale_proj, long nbins_proj, long nrays_proj){
    PolarScanParam_t* param_proj = NULL;
    
    long nrays = PolarScanParam_getNrays(param);
    
    double bin_scaling = rscale_proj/rscale;
    double ray_scaling = (double) nrays/nrays_proj;
    
    param_proj = RAVE_OBJECT_NEW(&PolarScanParam_TYPE);
    
    if (param_proj == NULL || !PolarScanParam_createData(param_proj, nbins_proj, nrays_proj, RaveDataType_DOUBLE)) {
        RAVE_ERROR0("Failed to create resampled polar scan parameter");
        goto error;
    }
    
    // copy the metadata
    PolarScanParam_setQuantity(param_proj, PolarScanParam_getQuantity(param));
    PolarScanParam_setOffset(param_proj,PolarScanParam_getOffset(param));
    PolarScanParam_setGain(param_proj,PolarScanParam_getGain(param));
    PolarScanParam_setNodata(param_proj,PolarScanParam_getNodata(param));
    PolarScanParam_setUndetect(param_proj,PolarScanParam_getUndetect(param));
    
    double value;
    RaveValueType valueType;
    
    // project onto new grid
    for(int iRay=0; iRay<nrays_proj; iRay++){
        for(int iBin=0; iBin<nbins_proj; iBin++){
            // initialize to nodata
            PolarScanParam_setValue(param_proj, iBin, iRay, PolarScanParam_getNodata(param));
            // read data from the scan parameter
            valueType = PolarScanParam_getValue(param, round(iBin*bin_scaling - 0.499999), round(iRay*ray_scaling - 0.499999), &value);
            // write data from the source scan parameter, to the newly projected scan paramter
            if (valueType != RaveValueType_UNDEFINED){
                PolarScanParam_setValue(param_proj, iBin, iRay, value);
            }
        }
    }
    error:
        return param_proj;
}



static int hasAzimuthGap(const float* points_local, const int nPoints, vol2bird_t* alldata) {

    int hasGap;
    int nObs[alldata->constants.nBinsGap];
    int iPoint;
    int iBinGap;
    int iBinGapNext;
    float azimuth;

    hasGap = FALSE;

    // Initialize histogram
    for (iBinGap = 0; iBinGap < alldata->constants.nBinsGap; iBinGap++) {
        nObs[iBinGap] = 0;
    }

    // Collect histogram data
    for (iPoint = 0; iPoint < nPoints; iPoint++) {
        azimuth = points_local[iPoint*alldata->misc.nDims];
        iBinGap = ((int) floor((azimuth / 360.0) * alldata->constants.nBinsGap)) % alldata->constants.nBinsGap;
        nObs[iBinGap]++;
    }

    // Detect adjacent bins in which the number of azimuth observations 
    // is less than the minimum required number
    for (iBinGap = 0; iBinGap < alldata->constants.nBinsGap; iBinGap++) {
        
        iBinGapNext = (iBinGap + 1) % alldata->constants.nBinsGap;
        
        if (nObs[iBinGap] < alldata->constants.nObsGapMin && nObs[iBinGapNext] < alldata->constants.nObsGapMin) {
            hasGap = TRUE;
        }
    }

    return hasGap;
    
} // hasAzimuthGap


const char* libvol2bird_version(void){
    return VERSION;
}


static int includeGate(const int iProfileType, const int iQuantityType, const unsigned int gateCode, vol2bird_t* alldata) {
    
    int doInclude = TRUE;
    
    if (gateCode & 1<<(alldata->flags.flagPositionStaticClutter)) {
        
        // i.e. flag 0 in gateCode is true
        // this gate is true in the static clutter map
        
        switch (iProfileType) {
            case 1 : 
                doInclude = FALSE;
                break;
            case 2 : 
                doInclude = FALSE;
                break;
            case 3 : 
                doInclude = FALSE;
                break;
            default :
                vol2bird_err_printf( "Something went wrong; behavior not implemented for given iProfileType.\n");
        }
    }
    
    if (gateCode & 1<<(alldata->flags.flagPositionDynamicClutter)) {
        
        // i.e. flag 1 in gateCode is true
        // this gate is part of the cluttermap (without fringe)
        
        switch (iProfileType) {
            case 1 : 
                doInclude = FALSE;
                break;
            case 2 : 
                break;
            case 3 : 
                break;
            default :
                vol2bird_err_printf( "Something went wrong; behavior not implemented for given iProfileType.\n");
        }
    }
        
    if (gateCode & 1<<(alldata->flags.flagPositionDynamicClutterFringe)) {
        
        // i.e. flag 2 in gateCode is true
        // this gate is part of the fringe of the cluttermap
        
        switch (iProfileType) {
            case 1 : 
                doInclude = FALSE;
                break;
            case 2 : 
                doInclude = FALSE;
                break;
            case 3 : 
                break;
            default :
                vol2bird_err_printf( "Something went wrong; behavior not implemented for given iProfileType.\n");
        }
    }
    
    if (iQuantityType && (gateCode & 1<<(alldata->flags.flagPositionVradMissing))) {
        
        // i.e. flag 3 in gateCode is true
        // this gate has reflectivity data but no corresponding radial velocity data
        // and iQuantityType != 0, i.e. we are not dealing with reflectivity quantities
        
        switch (iProfileType) {
            case 1 : 
		doInclude = FALSE;
                break;
            case 2 : 
                doInclude = FALSE;
                break;
            case 3 : 
                doInclude = FALSE;
                break;
            default :
                vol2bird_err_printf( "Something went wrong; behavior not implemented for given iProfileType.\n");
        }
    }

    if (!iQuantityType && (gateCode & 1<<(alldata->flags.flagPositionVradMissing)) && alldata->options.requireVrad) {
        
        // i.e. flag 3 in gateCode is true
        // this gate has reflectivity data but no corresponding radial velocity data
        // and iQuantityType == 0, i.e. we are dealing with reflectivity quantities
        // and requireVrad is true, i.e. we exclude gates without radial velocity data
        
        switch (iProfileType) {
            case 1 : 
                doInclude = FALSE;
                break;
            case 2 : 
                doInclude = FALSE;
                break;
            case 3 : 
                doInclude = FALSE;
                break;
            default :
                vol2bird_err_printf( "Something went wrong; behavior not implemented for given iProfileType.\n");
        }
    }
    
    if (gateCode & 1<<(alldata->flags.flagPositionDbzTooHighForBirds)) {
        
        // i.e. flag 4 in gateCode is true
        // this gate's dbz value is too high to be due to birds, it must be 
        // caused by something else

        
        switch (iProfileType) {
            case 1 : 
                doInclude = FALSE;
                break;
            case 2 : 
                break;
            case 3 : 
                break;
            default :
                vol2bird_err_printf( "Something went wrong; behavior not implemented for given iProfileType.\n");
        }
    }
    
    if (gateCode & 1<<(alldata->flags.flagPositionVradTooLow)) {
        
        // i.e. flag 5 in gateCode is true
        // this gate's radial velocity is very low, and therefore excluded as potential clutter.
        
        switch (iProfileType) {
            case 1 : 
                doInclude = FALSE;
                break;
            case 2 : 
                doInclude = FALSE;
                break;
            case 3 : 
                doInclude = FALSE;
                break;
            default :
                vol2bird_err_printf( "Something went wrong; behavior not implemented for given iProfileType.\n");
        }
    }
    
    if (iQuantityType && (gateCode & 1<<(alldata->flags.flagPositionVDifMax))) {

        // i.e. iQuantityType !=0, we are dealing with a selection for svdfit.
        // i.e. flag 6 in gateCode is true
        // after the first svdfit, this gate's fitted vRad was more than 
        // VDIFMAX away from the observed vRad for that gate. It is therefore
        // considered an outlier
        
        switch (iProfileType) {
            case 1 : 
                doInclude = FALSE;
                break;
            case 2 : 
                doInclude = FALSE;
                break;
            case 3 : 
                doInclude = FALSE;
                break;
            default :
                vol2bird_err_printf( "Something went wrong; behavior not implemented for given iProfileType.\n");
        }
    }



    if (!iQuantityType && (gateCode & 1<<(alldata->flags.flagPositionAzimOutOfRange))) {

        // i.e. iQuantityType == 0, we are NOT dealing with a selection for svdfit, but with a selection of reflectivities.
	// Azimuth selection does not apply to svdfit, because svdfit requires data at all azimuths
        // i.e. flag 7 in gateCode is true
        // the user can specify to exclude gates based on their azimuth;
        // this clause is for gates that have too low azimuth
        
        switch (iProfileType) {
            case 1 : 
                doInclude = FALSE;
                break;
            case 2 : 
                doInclude = FALSE;
                break;
            case 3 : 
                doInclude = FALSE;
                break;
            default :
                vol2bird_err_printf( "Something went wrong; behavior not implemented for given iProfileType.\n");
        }
    }


    return doInclude;

} // includeGate


/**
 * Function name: isRegularFile
 * Intent: determines whether the given path is to a regular file
 * Note: also returns true on existing directories
 */
int isRegularFile(const char *path) {
    return (access(path, F_OK) != -1);
} /* end function is_regular_file */

#ifndef NOCONFUSE
static int readUserConfigOptions(cfg_t** cfg, const char * optsConfFilename) {


    cfg_opt_t opts[] = {
        CFG_FLOAT("HLAYER",HLAYER, CFGF_NONE),
        CFG_INT("NLAYER",NLAYER, CFGF_NONE),
        CFG_FLOAT("RANGEMIN",RANGEMIN, CFGF_NONE),
        CFG_FLOAT("RANGEMAX",RANGEMAX, CFGF_NONE),
        CFG_FLOAT("AZIMMIN",AZIMMIN, CFGF_NONE),
        CFG_FLOAT("AZIMMAX",AZIMMAX, CFGF_NONE),
        CFG_FLOAT("ELEVMIN",ELEVMIN, CFGF_NONE),
        CFG_FLOAT("ELEVMAX",ELEVMAX, CFGF_NONE),
        CFG_FLOAT("RADAR_WAVELENGTH_CM",RADAR_WAVELENGTH_CM,CFGF_NONE),
        CFG_BOOL("USE_CLUTTERMAP",USE_CLUTTERMAP,CFGF_NONE),
        CFG_STR("CLUTTERMAP",CLUTTERMAP,CFGF_NONE),
        CFG_FLOAT("CLUTTERVALUEMIN",CLUTTERVALUEMIN, CFGF_NONE),
        CFG_BOOL("VERBOSE_OUTPUT_REQUIRED",VERBOSE_OUTPUT_REQUIRED,CFGF_NONE),
        CFG_BOOL("PRINT_DBZ",PRINT_DBZ,CFGF_NONE),
        CFG_BOOL("PRINT_DEALIAS",PRINT_DEALIAS,CFGF_NONE),
        CFG_BOOL("PRINT_VRAD",PRINT_VRAD,CFGF_NONE),
        CFG_BOOL("PRINT_RHOHV",PRINT_RHOHV,CFGF_NONE),
        CFG_BOOL("PRINT_CELL",PRINT_CELL,CFGF_NONE),
        CFG_BOOL("PRINT_CELL_PROP",PRINT_CELL_PROP,CFGF_NONE),
        CFG_BOOL("PRINT_TEXTURE",PRINT_TEXTURE,CFGF_NONE),
        CFG_BOOL("PRINT_CLUT",PRINT_CLUT,CFGF_NONE),
        CFG_BOOL("PRINT_OPTIONS",PRINT_OPTIONS,CFGF_NONE),
        CFG_BOOL("FIT_VRAD",FIT_VRAD,CFGF_NONE),
        CFG_BOOL("PRINT_PROFILE",PRINT_PROFILE,CFGF_NONE),
        CFG_BOOL("PRINT_POINTS_ARRAY",PRINT_POINTS_ARRAY,CFGF_NONE),
        CFG_FLOAT("MIN_NYQUIST_VELOCITY",MIN_NYQUIST_VELOCITY,CFGF_NONE),
        CFG_FLOAT("MAX_NYQUIST_DEALIAS",MAX_NYQUIST_DEALIAS,CFGF_NONE),
        /* initialize STDEV_BIRD to -FLT_MAX, final initialization will depend on radar wavelength */
        CFG_FLOAT("STDEV_BIRD",-FLT_MAX,CFGF_NONE),
        CFG_FLOAT("STDEV_CELL",STDEV_CELL,CFGF_NONE),
        CFG_FLOAT("SIGMA_BIRD",SIGMA_BIRD,CFGF_NONE),
        CFG_FLOAT("ETAMAX",ETAMAX,CFGF_NONE),
        CFG_FLOAT("ETACELL",ETACELL,CFGF_NONE),
        CFG_STR("DBZTYPE",DBZTYPE,CFGF_NONE),
        CFG_BOOL("REQUIRE_VRAD",REQUIRE_VRAD,CFGF_NONE),
        CFG_BOOL("DEALIAS_VRAD",DEALIAS_VRAD,CFGF_NONE),
        CFG_BOOL("DEALIAS_RECYCLE",DEALIAS_RECYCLE,CFGF_NONE),
        CFG_BOOL("EXPORT_BIRD_PROFILE_AS_JSON",FALSE,CFGF_NONE),
        CFG_BOOL("DUALPOL",DUALPOL,CFGF_NONE),
        CFG_BOOL("SINGLEPOL",SINGLEPOL,CFGF_NONE),
        CFG_FLOAT("DBZMIN",DBZMIN,CFGF_NONE),
        CFG_FLOAT("RHOHVMIN",RHOHVMIN,CFGF_NONE),
        CFG_BOOL("RESAMPLE",RESAMPLE,CFGF_NONE),
        CFG_FLOAT("RESAMPLE_RSCALE",RESAMPLE_RSCALE,CFGF_NONE),
        CFG_INT("RESAMPLE_NBINS",RESAMPLE_NBINS,CFGF_NONE),
        CFG_INT("RESAMPLE_NRAYS",RESAMPLE_NRAYS,CFGF_NONE),
        CFG_FLOAT_LIST("MISTNET_ELEVS", MISTNET_ELEVS, CFGF_NONE),
        CFG_BOOL("MISTNET_ELEVS_ONLY", MISTNET_ELEVS_ONLY, CFGF_NONE),
        CFG_BOOL("USE_MISTNET", USE_MISTNET, CFGF_NONE),
        CFG_STR("MISTNET_PATH",MISTNET_PATH,CFGF_NONE),
        CFG_END()
    };
    
    (*cfg) = cfg_init(opts, CFGF_NONE);
    int result = cfg_parse((*cfg), optsConfFilename);

    if (result == CFG_FILE_ERROR){
       vol2bird_err_printf( "Warning: no user configuration file '%s' found. Using default settings ...\n", optsConfFilename);
    }
    else{
       vol2bird_err_printf( "Loaded user configuration file '%s' ...\n", optsConfFilename);
    }

    if (result == CFG_PARSE_ERROR) {
        return 1;
    }
   
    return 0;

} // readUserConfigOptions
    
#endif


// copies shared metadata from rave polar volume to rave vertical profile
static int mapVolumeToProfile(VerticalProfile_t* vp, PolarVolume_t* volume){
    //assert that the volume and vertical profile are defined
    RAVE_ASSERT((vp != NULL), "vp == NULL");
    RAVE_ASSERT((volume != NULL), "volume == NULL");
    
    //copy the metadata
    VerticalProfile_setTime(vp,PolarVolume_getTime(volume));
    VerticalProfile_setDate(vp,PolarVolume_getDate(volume));
    VerticalProfile_setSource(vp,PolarVolume_getSource(volume));
    VerticalProfile_setLongitude(vp,PolarVolume_getLongitude(volume));
    VerticalProfile_setLatitude(vp,PolarVolume_getLatitude(volume));
    VerticalProfile_setHeight(vp,PolarVolume_getHeight(volume));
   
    return 0;
}

int mapDataToRave(PolarVolume_t* volume, vol2bird_t* alldata) {
    int result = 0;
    //assert that the vertical profile is defined
    RAVE_ASSERT((alldata->vp != NULL), "vp == NULL");
    
    // copy shared metadata from volume to profile
    mapVolumeToProfile(alldata->vp, volume);

    // copy vol2bird profile data to RAVE profile
    VerticalProfile_setLevels(alldata->vp,alldata->options.nLayers);
    VerticalProfile_setInterval(alldata->vp,alldata->options.layerThickness);
    VerticalProfile_setMinheight(alldata->vp, 0);
    VerticalProfile_setMaxheight(alldata->vp, alldata->options.nLayers * alldata->options.layerThickness);
    
    //intialize attributes for /how and /what
    RaveAttribute_t* attr_beamwidth = RaveAttributeHelp_createDouble("how/beamwidth", PolarVolume_getBeamwidth(volume)*180/PI);
    RaveAttribute_t* attr_wavelength = RaveAttributeHelp_createDouble("how/wavelength", alldata->options.radarWavelength);
    RaveAttribute_t* attr_rcs_bird = RaveAttributeHelp_createDouble("how/rcs_bird", alldata->options.birdRadarCrossSection);
    RaveAttribute_t* attr_sd_vvp_thresh = RaveAttributeHelp_createDouble("how/sd_vvp_thresh", alldata->options.stdDevMinBird);
    RaveAttribute_t* attr_dealiased = RaveAttributeHelp_createLong("how/dealiased", alldata->options.dealiasVrad);
    RaveAttribute_t* attr_task = RaveAttributeHelp_createString("how/task", PROGRAM);
    RaveAttribute_t* attr_task_version = RaveAttributeHelp_createString("how/task_version", VERSION);
    RaveAttribute_t* attr_task_args = RaveAttributeHelp_createString("how/task_args", alldata->misc.task_args);
    RaveAttribute_t* attr_comment = RaveAttributeHelp_createString("how/comment", "");
    RaveAttribute_t* attr_minrange = RaveAttributeHelp_createDouble("how/minrange", alldata->options.rangeMin/1000);
    RaveAttribute_t* attr_maxrange = RaveAttributeHelp_createDouble("how/maxrange", alldata->options.rangeMax/1000);
    RaveAttribute_t* attr_minazim = RaveAttributeHelp_createDouble("how/minazim", alldata->options.azimMin);
    RaveAttribute_t* attr_maxazim = RaveAttributeHelp_createDouble("how/maxazim", alldata->options.azimMax);
    RaveAttribute_t* attr_cluttermap = RaveAttributeHelp_createString("how/clutterMap", "");
    RaveAttribute_t* attr_filename_pvol = RaveAttributeHelp_createString("how/filename_pvol", alldata->misc.filename_pvol);
    RaveAttribute_t* attr_filename_vp = RaveAttributeHelp_createString("how/filename_vp", alldata->misc.filename_vp);
    RaveAttribute_t* attr_vcp = RaveAttributeHelp_createLong("how/vcp", alldata->misc.vcp);

    //add /how and /what attributes to the vertical profile object
    VerticalProfile_addAttribute(alldata->vp, attr_beamwidth);
    VerticalProfile_addAttribute(alldata->vp, attr_wavelength);
    VerticalProfile_addAttribute(alldata->vp, attr_rcs_bird);
    VerticalProfile_addAttribute(alldata->vp, attr_sd_vvp_thresh);
    VerticalProfile_addAttribute(alldata->vp, attr_dealiased);
    VerticalProfile_addAttribute(alldata->vp, attr_task);
    VerticalProfile_addAttribute(alldata->vp, attr_task_version);
    VerticalProfile_addAttribute(alldata->vp, attr_task_args);
    VerticalProfile_addAttribute(alldata->vp, attr_comment);
    VerticalProfile_addAttribute(alldata->vp, attr_minrange);
    VerticalProfile_addAttribute(alldata->vp, attr_maxrange);
    VerticalProfile_addAttribute(alldata->vp, attr_minazim);
    VerticalProfile_addAttribute(alldata->vp, attr_maxazim);
    VerticalProfile_addAttribute(alldata->vp, attr_cluttermap);
    VerticalProfile_addAttribute(alldata->vp, attr_filename_pvol);
    VerticalProfile_addAttribute(alldata->vp, attr_filename_vp);
    VerticalProfile_addAttribute(alldata->vp, attr_vcp);
      
    //-------------------------------------------//
    //   map the profile data to rave fields     //
    //-------------------------------------------//
    
    //layer specification:
    profileArray2RaveField(alldata, 1, 0, "HGHT", RaveDataType_DOUBLE);

    //bird-specific quantities:
    profileArray2RaveField(alldata, 1, 5, "ff", RaveDataType_DOUBLE);
    profileArray2RaveField(alldata, 1, 6, "dd", RaveDataType_DOUBLE);
    profileArray2RaveField(alldata, 1, 2, "u", RaveDataType_DOUBLE);
    profileArray2RaveField(alldata, 1, 3, "v", RaveDataType_DOUBLE);
    profileArray2RaveField(alldata, 1, 4, "w", RaveDataType_DOUBLE);
    profileArray2RaveField(alldata, 1, 8, "gap", RaveDataType_INT);
    profileArray2RaveField(alldata, 1, 9, "dbz", RaveDataType_DOUBLE);    
    profileArray2RaveField(alldata, 1, 11, "eta", RaveDataType_DOUBLE);
    profileArray2RaveField(alldata, 1, 12, "dens", RaveDataType_DOUBLE);        
    profileArray2RaveField(alldata, 1, 10, "n", RaveDataType_LONG);
    profileArray2RaveField(alldata, 1, 13, "n_dbz", RaveDataType_LONG);    

    //quantities calculated from all scatterers:
    profileArray2RaveField(alldata, 3, 7, "sd_vvp", RaveDataType_DOUBLE);    
    profileArray2RaveField(alldata, 3, 9, alldata->options.dbzType, RaveDataType_DOUBLE);    
    profileArray2RaveField(alldata, 3, 10, "n_all", RaveDataType_LONG);    
    profileArray2RaveField(alldata, 3, 13, "n_dbz_all", RaveDataType_LONG);        

    //some unused quantities for later reference:
    //profileArray2RaveField(alldata, 1, 2, "u", RaveDataType_DOUBLE);
    //profileArray2RaveField(alldata, 1, 3, "v", RaveDataType_DOUBLE);
     //initialize start and end date attributes to the vertical profile object
    RaveAttribute_t* attr_startdate = RaveAttributeHelp_createString("how/startdate", PolarVolume_getStartDate(volume));
    RaveAttribute_t* attr_starttime = RaveAttributeHelp_createString("how/starttime", PolarVolume_getStartTime(volume));
    RaveAttribute_t* attr_enddate = RaveAttributeHelp_createString("how/enddate", PolarVolume_getEndDate(volume));
    RaveAttribute_t* attr_endtime = RaveAttributeHelp_createString("how/endtime", PolarVolume_getEndTime(volume));


    //add the start and end date attributes to the vertical profile object
    VerticalProfile_addAttribute(alldata->vp, attr_startdate);
    VerticalProfile_addAttribute(alldata->vp, attr_starttime);
    VerticalProfile_addAttribute(alldata->vp, attr_enddate);
    VerticalProfile_addAttribute(alldata->vp, attr_endtime);

    RAVE_OBJECT_RELEASE(attr_beamwidth);
    RAVE_OBJECT_RELEASE(attr_wavelength);
    RAVE_OBJECT_RELEASE(attr_rcs_bird);
    RAVE_OBJECT_RELEASE(attr_sd_vvp_thresh);
    RAVE_OBJECT_RELEASE(attr_dealiased);
    RAVE_OBJECT_RELEASE(attr_task);
    RAVE_OBJECT_RELEASE(attr_task_version);
    RAVE_OBJECT_RELEASE(attr_task_args);
    RAVE_OBJECT_RELEASE(attr_comment);
    RAVE_OBJECT_RELEASE(attr_minrange);
    RAVE_OBJECT_RELEASE(attr_maxrange);
    RAVE_OBJECT_RELEASE(attr_minazim);
    RAVE_OBJECT_RELEASE(attr_maxazim);
    RAVE_OBJECT_RELEASE(attr_cluttermap);
    RAVE_OBJECT_RELEASE(attr_filename_pvol);
    RAVE_OBJECT_RELEASE(attr_filename_vp);
    RAVE_OBJECT_RELEASE(attr_vcp);

    RAVE_OBJECT_RELEASE(attr_startdate);
    RAVE_OBJECT_RELEASE(attr_starttime);
    RAVE_OBJECT_RELEASE(attr_enddate);
    RAVE_OBJECT_RELEASE(attr_endtime);
    result=1;

    return result;
    
}



// this function replaces NODATA and UNDETECT float values to NA and NAN
double nanify(double value){
    double output = value;
    if(value == NODATA) output = NAN;
    if(value == UNDETECT) output = NAN;
    return output;
} // nanify

void nanify_str(char* buff, const char* fmt, double v) {
  if (v == NODATA) {
    strcpy(buff, "na");
  } else if (v == UNDETECT) {
    strcpy(buff, "nan");
  } else {
    sprintf(buff, fmt, v);
  }
}

char* nanify_vpts(float value, const char* fmt) {
  char* output = malloc(15 + 1); // Allocate enough memory for a 15-character float string plus null terminator
  if (value == NODATA) {
    strcpy(output, "");
  } else if (value == UNDETECT) {
    strcpy(output, "NaN");
  } else {
    sprintf(output, fmt, value);
  }
  output[15] = '\0'; // Ensure output is null-terminated
  return output;
}

void create_profile_printout_str(char* printbuffer, int buflen, const char* date, const char* time,
    float HGHT, float u, float v, float w, float ff, float dd,
    float sd_vvp, char gap, float dbz, float eta, float dens, float DBZH,
    float n, float n_dbz, float n_all, float n_dbz_all)
{
  char s_HGHT[16], s_u[16], s_v[16], s_w[16], s_ff[16], s_dd[16];
  char s_sd_vvp[16], s_dbz[16], s_eta[16], s_dens[16], s_DBZH[16];
  char s_n[16], s_n_dbz[16], s_n_all[16], s_n_dbz_all[16];
  memset(printbuffer, 0, sizeof(char)*buflen);
  sprintf(s_HGHT, "%4.f", HGHT);
  nanify_str(s_u, "%6.2f", u);
  nanify_str(s_v, "%6.2f", v);
  nanify_str(s_w, "%7.2f", w);
  nanify_str(s_ff, "%5.2f", ff);
  nanify_str(s_dd, "%5.1f", dd);
  nanify_str(s_sd_vvp, "%6.2f", sd_vvp);
  nanify_str(s_dbz, "%6.2f", dbz);
  nanify_str(s_eta, "%6.1f", eta);
  nanify_str(s_dens, "%6.2f", dens);
  nanify_str(s_DBZH, "%6.2f", DBZH);
  nanify_str(s_n, "%5.f", n);
  nanify_str(s_n_dbz, "%5.f", n_dbz);
  nanify_str(s_n_all, "%5.f", n_all);
  nanify_str(s_n_dbz_all, "%5.f", n_dbz_all);
  sprintf(printbuffer, "%8s %.4s %4s %6s %6s %7s %5s %5s %6s %1c %6s %6s %6s %6s %5s %5s %5s %5s", date, time, s_HGHT,
      s_u, s_v, s_w, s_ff, s_dd, s_sd_vvp, gap, s_dbz, s_eta, s_dens, s_DBZH, s_n, s_n_dbz, s_n_all, s_n_dbz_all);
}

static int profileArray2RaveField(vol2bird_t* alldata, int idx_profile, int idx_quantity, const char* quantity, RaveDataType raveType){
    int result = 0;
    float* profile;
    
    RaveField_t* field = RAVE_OBJECT_NEW(&RaveField_TYPE);
    if (RaveField_createData(field, 1, alldata->options.nLayers, raveType) == 0){
        vol2bird_err_printf("Error pre-allocating field '%s'.\n", quantity);
        return -1;
    }
    
    switch (idx_profile) {
        case 1 :
            profile=alldata->profiles.profile1;
            break; 
        case 2 :
            profile=alldata->profiles.profile2;
            break; 
        case 3 :
            profile=alldata->profiles.profile3;
            break;
        default:
            vol2bird_err_printf( "Something is wrong this should not happen.\n");
            goto done;
    }
    
    int iRowProfile;
    int nColProfile = alldata->profiles.nColsProfile;
    for (iRowProfile = 0; iRowProfile < alldata->profiles.nRowsProfile; iRowProfile++) {
        RaveField_setValue(field, 0, iRowProfile, profile[idx_quantity +iRowProfile * nColProfile]);
    }
    
    result = verticalProfile_AddCustomField(alldata->vp, field, quantity);
    
    done:
        RAVE_OBJECT_RELEASE(field);
    
        return result;
}


static int verticalProfile_AddCustomField(VerticalProfile_t* self, RaveField_t* field, const char* quantity)
{
    int result = 0;
    RAVE_ASSERT((self != NULL), "self == NULL");
    RaveAttribute_t* attr = RaveAttributeHelp_createString("what/quantity", quantity);
    RaveAttribute_t* attr_gain = RaveAttributeHelp_createDouble("what/gain", 1.0);
    RaveAttribute_t* attr_offset = RaveAttributeHelp_createDouble("what/offset", 0.0);
    RaveAttribute_t* attr_nodata = RaveAttributeHelp_createDouble("what/nodata", NODATA);
    RaveAttribute_t* attr_undetect = RaveAttributeHelp_createDouble("what/undetect", UNDETECT);

    if (attr == NULL || !RaveField_addAttribute(field, attr)) {
        RAVE_ERROR0("Failed to add what/quantity attribute to field");
        goto done;
    }
    if (attr_gain == NULL || !RaveField_addAttribute(field, attr_gain)) {
        RAVE_ERROR0("Failed to add what/gain attribute to field");
        goto done;
    }
    if (attr_offset == NULL || !RaveField_addAttribute(field, attr_offset)) {
        RAVE_ERROR0("Failed to add what/offset attribute to field");
        goto done;
    }
    if (attr_nodata == NULL || !RaveField_addAttribute(field, attr_nodata)) {
        RAVE_ERROR0("Failed to add what/nodata attribute to field");
        goto done;
    }
    if (attr_undetect == NULL || !RaveField_addAttribute(field, attr_undetect)) {
        RAVE_ERROR0("Failed to add what/undetect attribute to field");
        goto done;
    }
    result = VerticalProfile_addField(self, field);
    
    done:
        RAVE_OBJECT_RELEASE(attr);
        RAVE_OBJECT_RELEASE(attr_gain);
        RAVE_OBJECT_RELEASE(attr_offset);
        RAVE_OBJECT_RELEASE(attr_nodata);
        RAVE_OBJECT_RELEASE(attr_undetect);
        return result;
}

const char *missing_values[] = {"", "NA", "NaN"};

//----------------------------------------------------------//
//   vpts exchange format https://aloftdata.eu/vpts-csv    //
//---------------------------------------------------------//
    

const field_t fields[] = {

    {
        .name = "radar",
        .description = "Radar identifier.",
        .type = "string",
        .example = "KBGM",
        .required = true,
        .constraints = {
            .required = true
        }
    },
    {
        .name = "datetime",
        .description = "Nominal date and time of the measurement, as an ISO 8601 formatted string in UTC.",
        .type = "string",
        .format = "%Y-%m-%dT%H:%M:%SZ",
        .example = "2016-09-01T00:02:00Z",
        .required = true,
        .constraints = {
            .required = true
        }
    },
    {
        .name = "height",
        .name_alternatives = (const char *[]){"HGHT", "bin_lower", NULL},
        .description = "Lower bound of the altitude bin in m above sea level.",
        .type = "integer",
        .example = "600",
        .required = true,
        .constraints = {
            .required = true,
            .minimum = -200,
            .maximum = 25000
        }
    },
    {
        .name = "u",
        .description = "Ground speed component west to east in m/s.",
        .type = "number",
        .example = "4.14",
        .required = false,
        .constraints = {
            .minimum = -100,
            .maximum = 100
        }
    },
    {
        .name = "v",
        .description = "Ground speed component north to south in m/s.",
        .type = "number",
        .example = "3.84",
        .required = false,
        .constraints = {
            .minimum = -100,
            .maximum = 100
        }
    },
    {
        .name = "w",
        .description = "Vertical speed in m/s.",
        .type = "number",
        .example = "12.17"
    },
    {
        .name = "ff",
        .name_alternatives = "speed",
        .description = "Horizontal ground speed in m/s.",
        .type = "number",
        .example = "5.65",
        .constraints = {
            .minimum = 0,
            .maximum = 100
        }
    },
    {
        .name = "dd",
        .name_alternatives = "direction",
        .description = "Ground speed direction in degrees clockwise from north.",
        .type = "number",
        .example = "47.2",
        .constraints = {
            .minimum = 0,
            .maximum = 360
        }
    },
    {
        .name = "sd_vvp",
        .name_alternatives = "rmse",
        .description = "VVP radial velocity standard deviation in m/s.",
        .type = "number",
        .example = "2.8",
        .constraints = {
            .minimum = 0,
            .maximum = 100
        }
    },
    {
        .name = "gap",
        .description = "Angular data gap detected.",
        .type = "boolean",
        .example = "FALSE"
    },
    {
        .name = "eta",
        .name_alternatives = "linear_eta",
        .description = "Animal reflectivity in cm^2/km^3.",
        .example = "46.9",
        .type = "number",
        .constraints = {
            .minimum = 0,
            .maximum = INFINITY
        }
    },
    {
        .name = "dens",
        .description = "Animal density in animals/km^3.",
        .type = "number",
        .example = "4.263636363636364",
        .constraints = {
            .minimum = 0,
            .maximum = INFINITY
        }
    },
    {
        .name = "dbz",
        .description = "Animal reflectivity factor in dBZ.",
        .type = "number",
        .example = "3.36",
        .constraints = {
            .minimum = -INFINITY,
            .maximum = 100
        }
    },
    {
        .name = "dbz_all",
        .name_alternatives = "DBZH",
        .description = "Total reflectivity factor (bio + meteo scattering) in dBZ.",
        .type = "number",
        .example = "0.5",
        .constraints = {
            .minimum = -INFINITY,
            .maximum = 100
        }
    },
    {
        .name = "n",
        .description = "Number of data points used for the ground speed estimates (quantities `u`, `v`, `w`, `ff`, `dd`).",
        .type = "integer",
        .example = "9006",
        .constraints = {
            .minimum = 0
        }
    },    
    {
        .name = "n_dbz",
        .description = "Number of data points used for reflectivity-based estimates (quantities `dbz`, `eta`, `dens`).",
        .type = "integer",
        .example = "13442",
        .constraints = {
            .minimum = 0
        }
    },
    {
        .name = "n_all",
        .description = "Number of data points used for the radial velocity standard deviation estimate (quantity `sd_vvp`).",
        .type = "integer",
        .example = "65947",
        .constraints = {
            .minimum = 0
        }
    },
    {
        .name = "n_dbz_all",
        .name_alternatives = "nbins",
        .description = "Number of data points used for the total reflectivity estimate (quantity `dbz_all`).",
        .type ="integer",
        .example = "104455",
        .constraints = {
            .minimum = 0
        }
    },
    {
        .name = "rcs",
        .description = "Radar cross section per bird in cm^2.",
        .type = "number",
        .example = "11",
        .constraints = {
            .minimum = 1e-15,
            .maximum = INFINITY
        }
    },
    {
        .name = "sd_vvp_threshold",
        .description = "Lower threshold in radial velocity standard deviation (profile quantity `sd_vvp`) in m/s. Biological signals with `sd_vvp < sd_vvp_threshold` are set to zero. Defaults to 2 m/s for C-band radars and 1 m/s for S-band radars if not specified.",
        .type = "number",
        .example = "2",
        .constraints = {
            .minimum = 0,
            .maximum = 100
        }
    },
    {
        .name = "vcp",
        .name_alternatives = (const char *[]){"scan_strategy", NULL},
        .description = "Volume coverage pattern, unitless. Documented on [Wikipedia](https://en.wikipedia.org/wiki/NEXRAD#Scan_strategies) for NEXRAD.",
        .type = "integer",
        .example = ""
    },
    {
        .name = "radar_latitude",
        .description = "Latitude of the radar location in decimal degrees, using the WGS84 datum. Constant for all records from the same `radar`.",
        .type = "number",
        .example = "42.19972",
        .constraints = {
            .required = true,
            .minimum = -90,
            .maximum = 90
        }
    },
    {
        .name = "radar_longitude",
        .description = "Longitude of the radar location in decimal degrees, using the WGS84 datum. Constant for all records from the same `radar`.",
        .type = "number",
        .example = "-75.98472",
        .constraints = {
            .required = true,
            .minimum = -180,
            .maximum = 180
        }
    },
    {
        .name = "radar_height",
        .description = "Height of the center of the radar antenna in m above sea level. Constant for all records from the same `radar`.",
        .type = "integer",
        .example = "519",
        .constraints = {
            .required = true,
            .minimum = -200,
            .maximum = 9000
        }
    },
        {
        .name = "radar_wavelength",
        .description = "Wavelength of the radar in cm. Constant for all records from the same `radar`. Most C-band radars operate at approximately 5.3 cm wavelength, most S-band radars at 10.6 cm.",
        .type = "number",
        .example = "10.6",
        .constraints = {
            .required = true,
            .minimum = 0.1,
            .maximum = 100
        }
    },
    {
        .name = "source_file",
        .description = "URL or path to the source file from which the data were derived.",
        .type = "string",
        .example = "s3://noaa-nexrad-level2/2016/09/01/KBGM/KBGM20160901_000212_V06",
        .constraints = {
            .pattern = "^(?=^[^./~])(^((?!\\.{2}).)*$).*$"
        }
    }
};

//-----------------------------------------//
//   functions to validate vpts fields    //
//---------------------------------------//
    

int is_datetime(const char *value, const char *format) {
    // Check if the value is in the correct format
    struct tm tm;
    return (strptime(value, format, &tm) != NULL);
}

int is_boolean(const char *value) {
    // Check if the value is either TRUE or FALSE
    return (strcmp(value, "TRUE") == 0 || strcmp(value, "FALSE") == 0);
}

bool is_float(const char *value) {
    // check if the string can be converted to a floating-point number
    char *endptr;
    strtod(value, &endptr);
    return endptr != value && *endptr == '\0';
}

bool is_integer(int value) {
    // an integer is always a valid value
    return true;
}
bool is_number(double value) {
    // check if the given value is a number (not NaN or infinity)
    return isfinite(value) || isnan(value);
}

bool is_string(const char *value) {
    return (value != NULL);
}

typedef union {
    int i;
    double d;
    char* c;
} VptsValue;


bool validate_value(const field_t field, const VptsValue value) {
    if (strcmp(field.type, "string") == 0) {
        // Check if the value is a string
        if (!is_string(value.c)) {
            return false;
        }
        // No additional checks needed for string type
        return true;
    }
    else if (strcmp(field.type, "number") == 0) {
        // Check if the value is a number
        if (!isfinite(value.d)) {
            return false;
        }
        double d_value = value.d;

        // Check if minimum value is set and validate
        if (!isnan(field.constraints.minimum)) {
            double min_value = field.constraints.minimum;
            if (d_value < min_value) {
                printf("Value for field '%s' is below minimum value of %f: %f\n", field.name, min_value, d_value);
                return false;
            }
        }

        // Check if maximum value is set and validate
        if (!isnan(field.constraints.maximum)) {
            double max_value = field.constraints.maximum;
            if (max_value != INFINITY && d_value > max_value) {
                printf("Value for field '%s' is above maximum value of %f: %f\n", field.name, max_value, d_value);
                return false;
            }
        }
    }
    else if (strcmp(field.type, "integer") == 0) {
        // Check if the value is an integer
        if (!is_integer(value.i)) {
            return false;
        }
        int i_value = value.i;

        // Check if the value is within the minimum and maximum bounds
        if (!isnan(field.constraints.minimum) && i_value < (int)field.constraints.minimum) {
            printf("Value for field '%s' is below minimum value of %d: %d\n", field.name, (int)field.constraints.minimum, i_value);
            return false;
        }

        if (!isnan(field.constraints.maximum)) {
            if (i_value > (int)field.constraints.maximum) {
                printf("Value for field '%s' is above maximum value of %d: %d\n", field.name, (int)field.constraints.maximum, i_value);
                return false;
            }
        }
    }
    else {
        // Invalid type
        return false;
    }
    return true;
}


void validate_fields(const field_t fields[], int num_fields, const VptsValue values[]) {
    for (int i = 0; i < num_fields; i++) {
        if (fields[i].required) {
            const VptsValue value = values[i];
            if (value.c == NULL || !is_number(value.d) || !is_integer(value.i)) {
                printf("WARNING! Missing value for required field: '%s'\n", fields[i].name);
            }
            else if (!validate_value(fields[i], value)) {
                printf("WARNING! Missing or invalid value for required field '%s' of type '%s'\n", fields[i].name, fields[i].type);
            }
        }
    }
}


int saveToODIM(RaveCoreObject* object, const char* filename){
    
    //define new Rave IO instance
    RaveIO_t* raveio = RAVE_OBJECT_NEW(&RaveIO_TYPE);
    //VpOdimIO_t* raveio = RAVE_OBJECT_NEW(&VpOdimIO_TYPE);

    //save in ODIM version 2.3, to keep HGHT in unit m and
    //keep deprecated wavelength attribute, as expected by bioRad
    RaveIO_setOdimVersion(raveio, RaveIO_ODIM_Version_2_3);

    //set the object to be saved
    RaveIO_setObject(raveio, object);
    
    //save the object
    int result;
    result = RaveIO_save(raveio, filename);
    
    RAVE_OBJECT_RELEASE(raveio);

    return result;    
}


void writeCSV(char *filename, vol2bird_t* alldata, char* source, char* fileIn, char* date, char* time, PolarVolume_t* pvol){
    
    // ----------------------------------------------------------------------------------------- //
    // this function writes the vertical profile to CSV format https://aloftdata.eu/vpts-csv     //
    // ---------------------------------------------------------------------------------------- //

    double longitude, latitude;
    int height;

    longitude = PolarVolume_getLongitude(pvol) / (M_PI/180.0);
    latitude = PolarVolume_getLatitude(pvol) / (M_PI/180.0);
    height = (int)PolarVolume_getHeight(pvol);
    source = PolarVolume_getSource(pvol);

    FILE *fp;
    fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("Failed to open file %s for writing.\n", filename);
        return;
    }

    int nRowsProfile = vol2birdGetNRowsProfile(alldata);
    int nColsProfile = vol2birdGetNColsProfile(alldata);

    float *profileBio, *profileAll;

    profileBio = vol2birdGetProfile(1, alldata);
    profileAll = vol2birdGetProfile(3, alldata);

    const int num_fields = sizeof(fields) / sizeof(fields[0]);

    union VptsValue {
        int i;
        double d;
        char* c;
    };

    union VptsValue vpts_values[num_fields];

    float *rcs, *sd_vvp_thresh, *wavelength;
    int *vcp;

    rcs = &alldata->options.birdRadarCrossSection;
    sd_vvp_thresh = &alldata->options.stdDevMinBird;
    vcp = &alldata->misc.vcp;
    wavelength = &alldata->options.radarWavelength;

    // Extract the radar name from the source variable
    char* radarName = NULL;
    char* p = strstr(source, "radar_name:");
    if (p != NULL) {
        p += strlen("radar_name:");
        radarName = strtok(p, ",");
    }

    fprintf(fp, "radar, datetime, height, u,v,w,ff,dd,sd_vvp,gap,dbz,eta,dens,DBZH,n,n_dbz,n_all,n_dbz_all,rcs,sd_vvp_threshold,vcp,radar_latitude,radar_longitude,radar_height,radar_wavelenght,source_file\n");

    int iRowProfile;
    int iCopied = 0;
    for (iRowProfile = 0; iRowProfile < nRowsProfile; iRowProfile++) {
        iCopied=iRowProfile*nColsProfile;

        char datetime[24];
        sprintf(datetime, "%.4s-%.2s-%.2sT%.2s:%.2s:00Z", date, date+4, date+6, time, time+2);

        //validate field functions
        /*
        union VptsValue vpts_values[] = {
            { .c = radarName },                                           // radar*
            { .c = datetime },                                            // datetime*
            { .i = (int)nanify(profileBio[0+iCopied])},                   // height*
            { .d = nanify(profileBio[2 + iCopied])},                      // u
            { .d = nanify(profileBio[3 + iCopied])},                      // v
            { .d = nanify(profileBio[4 + iCopied])},                      // w
            { .d = nanify(profileBio[5 + iCopied])},                      // ff
            { .d = nanify(profileBio[6 + iCopied])},                      // dd
            { .d = nanify(profileBio[7 + iCopied])},                      // sd_vvp
            { .c = profileBio[8 + iCopied] == TRUE ? "TRUE" : "FALSE"},   // gap
            { .d = nanify(profileBio[11 + iCopied])},                     // eta
            { .d = nanify(profileBio[12 + iCopied])},                     // dens
            { .d = nanify(profileBio[9 + iCopied])},                      // dbz
            { .d = nanify(profileAll[9 + iCopied])},                      // DBZH
            { .d = profileBio[10 + iCopied]},                             // n
            { .d = nanify(profileBio[13 + iCopied])},                     // n_dbz
            { .d = nanify(profileAll[10 + iCopied])},                     // n_all
            { .d = nanify(profileAll[13 + iCopied])},                     // n_dbz_all
            { .d = *rcs },                                                // rcs
            { .d = *sd_vvp_thresh },                                      // sd_vvp_threshold
            { .i = *vcp },                                                // vcp
            { .d = latitude },                                            // radar_latitude
            { .d = longitude },                                           // radar_longitude
            { .i = height },                                              // radar_height
            { .d = *wavelength },                                         // radar_wavelength
            { .c = source }                                               // source_file
        };

        validate_fields(fields, num_fields, vpts_values);
        */

        //write to CSV format
        fprintf(fp,"%s,%s,%d,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%.2f,%.2f,%d,%f,%f,%d,%f,%s\n",
        radarName,                                                              //radar*
        datetime,                                                               //datetime*    
        (int)nanify(profileBio[0+iCopied]),                                     //height*
        nanify_vpts(profileBio[2 + iCopied],"%6.2f"),                           //u
        nanify_vpts(profileBio[3 + iCopied], "%6.2f"),                          //v
        nanify_vpts(profileBio[4 + iCopied], "%7.2f"),                          //w
        nanify_vpts(profileBio[5 + iCopied], "%5.2f"),                          //ff
        nanify_vpts(profileBio[6 + iCopied], "%5.1f"),                          //dd
        nanify_vpts(profileBio[7 + iCopied],  "%6.2f"),                         //sd_vvp
        profileBio[8 + iCopied] == TRUE ? "TRUE" : "FALSE",                     // gap
        nanify_vpts(profileBio[11 + iCopied],  "%6.1f"),                        // eta
        nanify_vpts(profileBio[12 + iCopied],  "%6.2f"),                        // dens
        nanify_vpts(profileBio[9 + iCopied], "%6.2f"),                          // dbz
        nanify_vpts(profileAll[9 + iCopied], "%6.2f"),                          // DBZH
        nanify_vpts(profileBio[10 + iCopied],  "%5.f"),                         // n
        nanify_vpts(profileBio[13 + iCopied],  "%5.f"),                        // n_dbz
        nanify_vpts(profileAll[10 + iCopied],  "%5.f"),                        // n_all
        nanify_vpts(profileAll[13 + iCopied], "%5.f"),                         // n_dbz_all
        *rcs, *sd_vvp_thresh, *vcp, latitude, longitude, height, *wavelength, fileIn);
    }

    profileAll = NULL;
    profileBio = NULL;
    free((void*) profileAll);
    free((void*) profileBio);

    fclose(fp);
}

static void printCellProp(CELLPROP* cellProp, float elev, int nCells, int nCellsValid, vol2bird_t *alldata){
    
    // ---------------------------------------------------------- //
    // this function prints the cell properties struct to stderr  //
    // ---------------------------------------------------------- //
    
    vol2bird_err_printf("#Cell analysis for elevation %f:\n",elev*RAD2DEG);
    vol2bird_err_printf("#Minimum cell area in km^2     : %f\n",alldata->constants.areaCellMin);
    vol2bird_err_printf("#Threshold for mean dBZ cell   : %g dBZ\n",alldata->misc.cellDbzMin);
    vol2bird_err_printf("#Threshold for mean stdev cell : %g m/s\n",alldata->options.cellStdDevMax);
    vol2bird_err_printf("#Valid cells                   : %i/%i\n#\n",nCellsValid,nCells);
    vol2bird_err_printf("cellProp: .index .nGates .nGatesClutter   .Area .dbzAvg .texAvg .cv   .dbzMax .iRangOfMax .iAzimOfMax .drop\n");
    for (int iCell = 0; iCell < nCells; iCell++) {
        if (cellProp[iCell].drop == TRUE) {
            continue;
        }
        vol2bird_err_printf("cellProp: %6d %7d %14d %7.2f %7.2f %7.2f %5.2f %7.2f %11d %11d %5c\n",
            cellProp[iCell].index,
            cellProp[iCell].nGates,
            cellProp[iCell].nGatesClutter,
            cellProp[iCell].area,
            cellProp[iCell].dbzAvg,
            cellProp[iCell].texAvg,
            cellProp[iCell].cv,
            cellProp[iCell].dbzMax,
            cellProp[iCell].iRangOfMax,
            cellProp[iCell].iAzimOfMax,
            cellProp[iCell].drop == TRUE ? 'T' : 'F');
    }
}


static void printGateCode(char* flags, const unsigned int gateCode) {
    
    // --------------------------------------------------- //
    // this function prints the integer value gateCode as  //
    // the equivalent sequence of bits                     //
    // --------------------------------------------------- //

    int iFlag;
    int nFlagsNeeded;
    int nFlags;
    int nFlagsMax;
    
    
    if (gateCode <= 0) {
        nFlagsNeeded = 0;
    }
    else {
        nFlagsNeeded = (int) ceil(log(gateCode + 1)/log(2));
    }
    
    nFlagsMax = 9;
    if (nFlagsNeeded > nFlagsMax) {
        vol2bird_err_printf("There's only space for %d flags\n. Aborting",nFlagsMax);
        return;
    }
    
    nFlags = nFlagsMax;

    for (iFlag = nFlags-1; iFlag >= 0; iFlag--) {
    
        int iFlagIsEnabled = (gateCode & (1 << iFlag)) >> iFlag;
        
        if (iFlagIsEnabled == TRUE) {
            flags[nFlags-iFlag-1] = '1';
        }
        else {
            flags[nFlags-iFlag-1] = '0';
        };
    }
    
    flags[nFlags] = '\0';
    
    return;

} // printGateCode


void printImage(PolarScan_t* scan, const char* quantity) {
    int nAzim = PolarScan_getNrays(scan);
    int nRang = PolarScan_getNbins(scan);
    int iRang;
    int iAzim;
    int needsSignChar;
    int needsFloat;
    int maxValue;
    
    PolarScanParam_t* scanParam = PolarScan_getParameter(scan, quantity);
    
    if(scanParam == (PolarScanParam_t *) NULL){
        vol2bird_err_printf("warning::printImage: quantity %s not found in scan\n",quantity);
        return;
    }

    double thisValue;
    int nChars;
    char* formatStr;
    char* naStr;
    
    maxValue = 0;
    needsSignChar = FALSE;
    needsFloat = FALSE;
    
    RaveValueType type;
    
    // first, determine how many characters are needed to print array 'imageInt'
    for (iAzim = 0; iAzim < nAzim; iAzim++) { 
        
        for (iRang = 0; iRang < nRang; iRang++) {
            
            type = PolarScanParam_getValue(scanParam,iRang,iAzim,&thisValue);
            if (thisValue < 0) {
                needsSignChar = TRUE;
            }
            
            if (XABS(thisValue - (long) thisValue) >= 0.01) {
                needsFloat = TRUE;
            }

            thisValue = XABS(thisValue);

            if (thisValue > maxValue) {
                maxValue = thisValue;
            }
        }
    }


    nChars = (int) ceil(log(maxValue + 1)/log(10));

    if (needsSignChar) {
        nChars += 1;
    }
    
    if (!needsFloat){
        switch (nChars) {
            case 0 :
                formatStr = "  %1f";
                naStr=" NA";
                break; 
            case 1 :
                formatStr = "  %1f";
                naStr=" NA";
                break; 
            case 2 :
                formatStr = " %2f"; 
                naStr=" NA";
                break;
            case 3 :
                formatStr = " %3f"; 
                naStr=" NA ";
                break;
            case 4 :
                formatStr = " %4f"; 
                naStr=" NA  ";
                break;
            default :
                formatStr = " %8f"; 
                naStr=" NA      ";
                
        }        
    }
    else{
        switch (nChars) {
            case 0 :
                formatStr = " %1.2f";
                naStr=" NA  ";
                break; 
            case 1 :
                formatStr = " %1.2f";
                naStr=" NA  ";
                break; 
            case 2 :
                formatStr = " %2.2f"; 
                naStr=" NA   ";
                break;
            case 3 :
                formatStr = " %3.2f"; 
                naStr=" NA    ";
                break;
            case 4 :
                formatStr = " %4.2f"; 
                naStr=" NA     ";
                break;
            default :
                formatStr = " %8.2f"; 
                naStr=" NA         ";

        }
    }

    
    for (iAzim = 0; iAzim < nAzim; iAzim++) { 
        
        for (iRang = 0; iRang < nRang; iRang++) {
                
            type = PolarScanParam_getValue(scanParam,iRang,iAzim,&thisValue);
            
            if (type == RaveValueType_DATA){
                vol2bird_err_printf(formatStr,thisValue);
            }
            else{
                vol2bird_err_printf(naStr,thisValue);

            }
        }
        vol2bird_err_printf("\n");
    }
        
} // printImageInt



static int printMeta(PolarScan_t* scan, const char* quantity) {
    
    vol2bird_err_printf("scan->heig = %f\n",PolarScan_getHeight(scan));
    vol2bird_err_printf("scan->elev = %f\n",PolarScan_getElangle(scan));
    vol2bird_err_printf("scan->nRang = %ld\n",PolarScan_getNbins(scan));
    vol2bird_err_printf("scan->nAzim = %ld\n",PolarScan_getNrays(scan));
    vol2bird_err_printf("scan->rangeScale = %f\n",PolarScan_getRscale(scan));
    vol2bird_err_printf("scan->azimScale = %f\n",360.0f/PolarScan_getNrays(scan));
    
    PolarScanParam_t* scanParam = PolarScan_getParameter(scan,quantity);
    
    if(scanParam != (PolarScanParam_t*) NULL){
        vol2bird_err_printf("scan->%s->valueOffset = %f\n",quantity,PolarScanParam_getOffset(scanParam));
        vol2bird_err_printf("scan->%s->valueScale = %f\n",quantity,PolarScanParam_getGain(scanParam));
        vol2bird_err_printf("scan->%s->nodata = %f\n",quantity,PolarScanParam_getNodata(scanParam));
        vol2bird_err_printf("scan->%s->undetect = %f\n",quantity,PolarScanParam_getUndetect(scanParam));
    }
    
    return 0;

} // printMeta




static int selectCellsToDrop(CELLPROP *cellProp, int nCells, int dualpol, vol2bird_t* alldata){
    int output = 0;
    if(dualpol){
        output = selectCellsToDrop_dualPol(cellProp, nCells, alldata);
    }
    else{
        output = selectCellsToDrop_singlePol(cellProp, nCells, alldata);
    }
    return output;
}



static int selectCellsToDrop_singlePol(CELLPROP *cellProp, int nCells, vol2bird_t* alldata){
    // determine which blobs to drop from map based on low mean dBZ / high stdev /
    // small area / high percentage clutter
    for (int iCell = 0; iCell < nCells; iCell++) {
        int notEnoughGates = cellProp[iCell].area < alldata->constants.areaCellMin;
        int dbzTooLow = cellProp[iCell].dbzAvg < alldata->misc.cellDbzMin;
        int texTooHigh = cellProp[iCell].texAvg > alldata->options.cellStdDevMax;
        int tooMuchClutter = ((float) cellProp[iCell].nGatesClutter / cellProp[iCell].nGates) > alldata->constants.cellClutterFractionMax;
        
        if (notEnoughGates) {
            
            // this blob is too small too be a weather cell, more likely 
            // that these are birds. So, drop the blob from the record of 
            // weather cells 

            cellProp[iCell].drop = TRUE;
            
            continue;
        }

        if (dbzTooLow && texTooHigh) {

            // apparently, we are dealing a blob that is fairly large (it
            // passes the 'notEnoughGates' condition above, but it has both 
            // a low dbz and a high vrad texture. In contrast, weather cells
            // have high dbz and low vrad texture. It is therefore unlikely
            // that the blob is a weather cell.
            
            if (tooMuchClutter) {
                // pass
            }
            else {
                
                // So at this point we have established that we are 
                // dealing with a blob that is:
                //     1. fairly large
                //     2. has too low dbz to be precipitation
                //     3. has too high vrad texture to be precipitation
                //     4. has too little clutter to be attributed to clutter
                // Therefore we drop it from the record of known weather cells
                
                cellProp[iCell].drop = TRUE;
                
            }
        }
    }
    return 1;
}


static int selectCellsToDrop_dualPol(CELLPROP *cellProp, int nCells, vol2bird_t* alldata){
    // determine which blobs to drop from map based on small area
    for (int iCell = 0; iCell < nCells; iCell++) {
        int notEnoughGates = cellProp[iCell].area < alldata->constants.areaCellMin;

        if (notEnoughGates) {
            
            // this blob is too small too be a weather cell, more likely 
            // that these are birds. So, drop the blob from the record of 
            // weather cells 

            cellProp[iCell].drop = TRUE;
            continue;
        }
    }
    return 1;
}


static void sortCellsByArea(CELLPROP *cellProp, const int nCells) {

    // ---------------------------------------------------------------- //
    // Sorting of the cell properties based on cell area.               //
    // ---------------------------------------------------------------- // 

    int iCell;
    int iCellOther;
    CELLPROP tmp;

    //Sorting of data elements using straight insertion method.
    for (iCell = 1; iCell < nCells; iCell++) {

        tmp = cellProp[iCell];

        iCellOther = iCell - 1;

        while (iCellOther >= 0 && cellProp[iCellOther].nGates < tmp.nGates) {

            cellProp[iCellOther + 1] = cellProp[iCellOther];

            iCellOther--;
        }

        cellProp[iCellOther + 1] = tmp;

    } //for iCell

    return;
} // sortCellsByArea



radarDataFormat determineRadarFormat(char* filename){
    
#ifdef IRIS
    if (isIRIS(filename)==0){
        return radarDataFormat_ODIM;
    }
#endif

#ifdef RSL
    if(RSL_filetype(filename) != UNKNOWN){
        return radarDataFormat_RSL;
    }
#endif

    // try to load the file using Rave
    // unfortunately this loads the entire file into memory,
    // but no other file type check function available in Rave.
    RaveIO_t* raveio = RaveIO_open(filename, 0, NULL);

    // check that a valid RaveIO_t pointer was returned
    if (raveio != (RaveIO_t*) NULL){
        RAVE_OBJECT_RELEASE(raveio);
        return radarDataFormat_ODIM;
    }
    
    return radarDataFormat_UNKNOWN;
}


static int removeDroppedCells(CELLPROP *cellProp, const int nCells) {
    int iCell;
    int iCopy;
    int nCopied;
    CELLPROP* cellPropCopy;
    CELLPROP cellPropEmpty;
    
    cellPropEmpty.iRangOfMax = -1;
    cellPropEmpty.iAzimOfMax = -1;
    cellPropEmpty.nGates = -1;
    cellPropEmpty.nGatesClutter = -1;
    cellPropEmpty.dbzAvg = 0.0f;
    cellPropEmpty.texAvg = 0.0f;
    cellPropEmpty.dbzMax = 0.0f;
    cellPropEmpty.index = -1;
    cellPropEmpty.drop = TRUE;
    cellPropEmpty.cv = 0.0f;



    #ifdef FPRINTFON
    for (iCell = 0; iCell < nCells; iCell++) {
        vol2bird_err_printf("(%d/%d): index = %d, nGates = %d\n",iCell,nCells,cellProp[iCell].index,cellProp[iCell].nGates);
    }
    vol2bird_err_printf("end of list\n");
    #endif


    
    cellPropCopy = (CELLPROP*) malloc(sizeof(CELLPROP) * nCells);
    if (!cellPropCopy) {
        vol2bird_err_printf("Requested memory could not be allocated in removeDroppedCells!\n");
        return -1;
    }    

    iCopy = 0;
    
    for (iCell = 0; iCell < nCells; iCell++) {    
        
        if (cellProp[iCell].drop == TRUE) {
            // pass
        }  
        else {
            cellPropCopy[iCopy] = cellProp[iCell];
            iCopy += 1;
        }
    }
    
    nCopied = iCopy;
    
    for (iCopy = 0; iCopy < nCopied; iCopy++) {
        cellProp[iCopy] = cellPropCopy[iCopy];
    }

    for (iCell = nCopied; iCell < nCells; iCell++) {
        cellProp[iCell] = cellPropEmpty;
    }

    #ifdef FPRINTFON
    for (iCell = 0; iCell < nCells; iCell++) {
        vol2bird_err_printf("(%d/%d): copied = %c, index = %d, nGates = %d\n",iCell,nCells,iCell < nCopied ? 'T':'F',cellProp[iCell].index,cellProp[iCell].nGates);
    }
    #endif 
    
    free(cellPropCopy);
       
    return nCopied;
   
}




static void updateFlagFieldsInPointsArray(const float* yObs, const float* yFitted, const int* includedIndex, 
                                   const int nPointsIncluded, float* points_local, vol2bird_t* alldata) {
                                       
    // ----------------------------------------------------------------------------------- //
    // after the first svdfit to the selection of points, we want to identify gates that   //
    // deviate strongly from the fitted vrad value (we assume that they are outliers). So, //
    // this function marks the corresponding points in 'points' by enabling the            //
    // 'flagPositionVDifMax' bit in 'gateCode', such that outliers will be rejected during //
    // the second svdfit iteration                                                         //
    // ----------------------------------------------------------------------------------- //

    int iPointIncluded;
    int iPoint;
    unsigned int gateCode;

    for (iPointIncluded = 0; iPointIncluded < nPointsIncluded; iPointIncluded++) {

        float absVDif = fabs(yObs[iPointIncluded]-yFitted[iPointIncluded]);
        
        if (absVDif > alldata->constants.absVDifMax) {
            
            iPoint = includedIndex[iPointIncluded];
            gateCode = (unsigned int) points_local[iPoint * alldata->points.nColsPoints + alldata->points.gateCodeCol];
            points_local[iPoint * alldata->points.nColsPoints + alldata->points.gateCodeCol] = (float) (gateCode |= 1<<(alldata->flags.flagPositionVDifMax));

        }
    } 

} // updateFlagFieldsInPointsArray



static int updateMap(PolarScan_t* scan, CELLPROP *cellProp, const int nCells, vol2bird_t* alldata) {

    // ------------------------------------------------------------------------- //
    // This function updates cellImage by dropping cells and reindexing the map. //
    // Leaving index 0 unused, will be used for assigning cell fringes           //
    // ------------------------------------------------------------------------- //

    int iGlobal;
    int iCell;
    int iCellNew;
    int nCellsValid;
    int cellImageValue;

    PolarScanParam_t *cellParam = PolarScan_getParameter(scan, CELLNAME);
    int* cellImage = (int *) PolarScanParam_getData(cellParam);

    int nGlobal = PolarScanParam_getNbins(cellParam) * PolarScanParam_getNrays(cellParam);

    #ifdef FPRINTFON
    int minValue = cellImage[0];
    int maxValue = cellImage[0];
    for (iGlobal = 1;iGlobal < nGlobal;iGlobal++) {
        if (cellImage[iGlobal] < minValue) {
            minValue = cellImage[iGlobal];
        }
        if (cellImage[iGlobal] > maxValue) {
            maxValue = cellImage[iGlobal];
        }
    }
    vol2bird_err_printf("minimum value in cellImage array = %d.\n", minValue);
    vol2bird_err_printf("maximum value in cellImage array = %d.\n", maxValue);
    #endif

    nCellsValid = nCells;

    for (iGlobal = 0; iGlobal < nGlobal; iGlobal++) {

        if (cellImage[iGlobal] == -1) {
            continue;
        }

        cellImageValue = cellImage[iGlobal];

        if (cellImageValue > nCells - 1) {
            vol2bird_err_printf( "You just asked for the properties of cell %d, which does not exist.\n", cellImageValue);
            continue;
        }

        if (cellProp[cellImageValue].drop == TRUE) {
            cellImage[iGlobal] = -1;
        }
    }

    // label small cells so that 'removeDroppedCells()' will remove them (below)
    for (iCell = 0; iCell < nCells; iCell++) {
        if (cellProp[iCell].area < alldata->constants.areaCellMin) {
            cellProp[iCell].drop = TRUE;
        }
    }

    // remove all cell that have .drop == TRUE
    nCellsValid = removeDroppedCells(&cellProp[0],nCells);

    // sort the cells by area
    sortCellsByArea(&cellProp[0],nCells);


    #ifdef FPRINTFON
    vol2bird_err_printf("nCellsValid = %d\n",nCellsValid);
    vol2bird_err_printf("\n");
    #endif

    // replace the values in cellImage with newly calculated index values:
    for (iCell = 0; iCell < nCells; iCell++) {

        if (iCell < nCellsValid) {
            iCellNew = -1 * (iCell + 2 + 100);
        }
        else {
            iCellNew = -1;
        }

        #ifdef FPRINTFON
        vol2bird_err_printf("before: cellProp[%d].index = %d.\n",iCell,cellProp[iCell].index);
        vol2bird_err_printf("before: cellProp[%d].nGates = %d.\n",iCell,cellProp[iCell].nGates);
        vol2bird_err_printf("before: iCell = %d.\n",iCell);
        vol2bird_err_printf("before: iCellNew = %d.\n",iCellNew);
        vol2bird_err_printf("\n");
        #endif

        for (iGlobal = 0; iGlobal < nGlobal; iGlobal++) {
            if (cellImage[iGlobal] == cellProp[iCell].index) {
                cellImage[iGlobal] = iCellNew;
            }
        }
        // have the indices in cellProp match the re-numbering
        cellProp[iCell].index = iCellNew;

    } // (iCell = 0; iCell < nCells; iCell++)


    // once you've re-numbered everything, flip the sign back and
    // remove the offset of 100...
    for (iGlobal = 0; iGlobal < nGlobal; iGlobal++) {
        if (cellImage[iGlobal] == -1) {
            // do nothing
        }
        else {
            cellImage[iGlobal] = (-1 * cellImage[iGlobal]) - 100;
        }
    }
    // ...and make sure the indices in cellProp match that change
    for (iCell = 0; iCell < nCells; iCell++) {

        if (cellProp[iCell].index == -1) {
            // do nothing
        }
        else {
            cellProp[iCell].index = (-1 * cellProp[iCell].index) - 100;
        }

        #ifdef FPRINTFON
        vol2bird_err_printf("after: cellProp[%d].index = %d.\n",iCell,cellProp[iCell].index);
        vol2bird_err_printf("after: cellProp[%d].nGates = %d.\n",iCell,cellProp[iCell].nGates);
        vol2bird_err_printf("\n");
        #endif
    }
    RAVE_OBJECT_RELEASE(cellParam);
    return nCellsValid;
} // updateMap



void vol2birdCalcProfiles(vol2bird_t *alldata) {

  int nPasses;
  int iPoint;
  int iLayer;
  int iPass;
  int iProfileType;

  if (alldata->misc.initializationSuccessful == FALSE) {
    vol2bird_err_printf( "You need to initialize vol2bird before you can use it. Aborting.\n");
    return;
  }

  // calculate the profiles in reverse order, because you need the result
  // of iProfileType == 3 in order to check whether chi < stdDevMinBird
  // when calculating iProfileType == 1
  for (iProfileType = alldata->profiles.nProfileTypes; iProfileType > 0; iProfileType--) {
    // ------------------------------------------------------------- //
    //                        prepare the profile                    //
    //                                                               //
    // iProfileType == 1: birds                                      //
    // iProfileType == 2: non-birds                                  //
    // iProfileType == 3: birds+non-birds                            //
    // ------------------------------------------------------------- //

    // FIXME: we better get rid of ProfileType==2 altogether, instead of
    // skipping it here.
    if (iProfileType == 2)
      continue;

    alldata->profiles.iProfileTypeLast = iProfileType;

    // if the user does not require fitting a model to the observed
    // vrad values, we don't need a second pass to remove dealiasing outliers
    if (alldata->options.fitVrad == TRUE) {
      nPasses = 2;
    } else {
      nPasses = 1;
    }

    // set a flag that indicates if we want to keep earlier dealiased values
    int recycleDealias = FALSE;
    if (iProfileType < 3 && alldata->options.dealiasRecycle) {
      recycleDealias = TRUE;
    }

    // reset the flagPositionVDifMax bit before calculating each profile
    for (iPoint = 0; iPoint < alldata->points.nRowsPoints; iPoint++) {
      unsigned int gateCode = (unsigned int) alldata->points.points[iPoint * alldata->points.nColsPoints + alldata->points.gateCodeCol];
      alldata->points.points[iPoint * alldata->points.nColsPoints + alldata->points.gateCodeCol] = (float) (gateCode &= ~(1
          << (alldata->flags.flagPositionVDifMax)));
    }

    // reset the dealiased vrad value before calculating each profile
    if (!recycleDealias) {
      for (iPoint = 0; iPoint < alldata->points.nRowsPoints; iPoint++) {
        alldata->points.points[iPoint * alldata->points.nColsPoints + alldata->points.vraddValueCol] = alldata->points.points[iPoint
            * alldata->points.nColsPoints + alldata->points.vradValueCol];
      }
    }

    for (iLayer = 0; iLayer < alldata->options.nLayers; iLayer++) {

      // these variables are needed just outside of the iPass loop below
      float chi = NAN;
      int hasGap = TRUE;
      float birdDensity = NAN;

      for (iPass = 0; iPass < nPasses; iPass++) {

        const int iPointFrom = alldata->points.indexFrom[iLayer];
        const int nPointsLayer = alldata->points.nPointsWritten[iLayer];

        int iPointLayer;
        int iPointIncluded;
        int iPointIncludedZ;
        int nPointsIncluded;
        int nPointsIncludedZ;

        float parameterVector[] = { NAN, NAN, NAN };
        float avar[] = { NAN, NAN, NAN };

        float *pointsSelection = malloc(sizeof(float) * nPointsLayer * alldata->misc.nDims);
        float *yNyquist = malloc(sizeof(float) * nPointsLayer);
        float *yDealias = malloc(sizeof(float) * nPointsLayer);
        float *yObs = malloc(sizeof(float) * nPointsLayer);
        float *yFitted = malloc(sizeof(float) * nPointsLayer);
        int *includedIndex = malloc(sizeof(int) * nPointsLayer);

        float *yObsSvdFit = yObs;
        float dbzValue = NAN;
        float undbzValue = NAN;
        double undbzSum = 0.0;
        float undbzAvg = NAN;
        float dbzAvg = NAN;
        float reflectivity = NAN;
        float chisq = NAN;
        float hSpeed = NAN;
        float hDir = NAN;

        for (iPointLayer = 0; iPointLayer < nPointsLayer; iPointLayer++) {

          pointsSelection[iPointLayer * alldata->misc.nDims + 0] = 0.0f;
          pointsSelection[iPointLayer * alldata->misc.nDims + 1] = 0.0f;

          yNyquist[iPointLayer] = 0.0f;
          yDealias[iPointLayer] = 0.0f;
          yObs[iPointLayer] = 0.0f;
          yFitted[iPointLayer] = 0.0f;

          includedIndex[iPointLayer] = -1;

        };

        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 0] = (iLayer + 0.5) * alldata->options.layerThickness;
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 1] = alldata->options.layerThickness;
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 2] = NODATA;
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 3] = NODATA;
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 4] = NODATA;
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 5] = NODATA;
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 6] = NODATA;
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 7] = NODATA;
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 8] = NODATA;
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 9] = NODATA;
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 10] = NODATA;
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 11] = NODATA;
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 12] = NODATA;
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 13] = NODATA;

        //Calculate the average reflectivity Z of the layer
        iPointIncludedZ = 0;
        for (iPointLayer = iPointFrom; iPointLayer < iPointFrom + nPointsLayer; iPointLayer++) {

          unsigned int gateCode = (unsigned int) alldata->points.points[iPointLayer * alldata->points.nColsPoints + alldata->points.gateCodeCol];

          if (includeGate(iProfileType, 0, gateCode, alldata) == TRUE) {

            // get the dbz value at this [azimuth, elevation]
            dbzValue = alldata->points.points[iPointLayer * alldata->points.nColsPoints + alldata->points.dbzValueCol];
            // convert from dB scale to linear scale
            if (isnan(dbzValue) == TRUE) {
              undbzValue = 0;
            } else {
              undbzValue = (float) exp(0.1 * log(10) * dbzValue);
            }
            // sum the undbz in this layer
            undbzSum += undbzValue;
            // raise the counter
            iPointIncludedZ += 1;

          }
        } // endfor (iPointLayer = 0; iPointLayer < nPointsLayer; iPointLayer++) {
        nPointsIncludedZ = iPointIncludedZ;

        // calculate bird densities from undbzSum
        if (nPointsIncludedZ > alldata->constants.nPointsIncludedMin) {
          // when there are enough valid points, convert undbzAvg back to dB-scale
          undbzAvg = (float) (undbzSum / nPointsIncludedZ);
          dbzAvg = (10 * log(undbzAvg)) / log(10);
        } else {
          undbzAvg = UNDETECT;
          dbzAvg = UNDETECT;
        }

        // convert from Z (not dBZ) in units of mm^6/m^3 to
        // reflectivity eta in units of cm^2/km^3
        reflectivity = alldata->misc.dbzFactor * undbzAvg;

        if (iProfileType == 1) {
          // calculate bird density in number of birds/km^3 by
          // dividing the reflectivity by the (assumed) cross section
          // of one bird
          birdDensity = reflectivity / alldata->options.birdRadarCrossSection;
        } else {
          birdDensity = UNDETECT;
        }

        // birdDensity and reflectivity should also be UNDETECT when undbzAvg is
        if (undbzAvg == UNDETECT) {
          reflectivity = UNDETECT;
          birdDensity = UNDETECT;
        }

        //Prepare the arguments of svdfit
        iPointIncluded = 0;
        for (iPointLayer = iPointFrom; iPointLayer < iPointFrom + nPointsLayer; iPointLayer++) {

          unsigned int gateCode = (unsigned int) alldata->points.points[iPointLayer * alldata->points.nColsPoints + alldata->points.gateCodeCol];

          if (includeGate(iProfileType, 1, gateCode, alldata) == TRUE) {

            // copy azimuth angle from the 'points' array
            pointsSelection[iPointIncluded * alldata->misc.nDims + 0] = alldata->points.points[iPointLayer * alldata->points.nColsPoints
                + alldata->points.azimAngleCol];
            // copy elevation angle from the 'points' array
            pointsSelection[iPointIncluded * alldata->misc.nDims + 1] = alldata->points.points[iPointLayer * alldata->points.nColsPoints
                + alldata->points.elevAngleCol];
            // copy nyquist interval from the 'points' array
            yNyquist[iPointIncluded] = alldata->points.points[iPointLayer * alldata->points.nColsPoints + alldata->points.nyquistCol];
            // copy the observed vrad value at this [azimuth, elevation]
            yObs[iPointIncluded] = alldata->points.points[iPointLayer * alldata->points.nColsPoints + alldata->points.vradValueCol];
            // copy the dealiased vrad value at this [azimuth, elevation]
            yDealias[iPointIncluded] = alldata->points.points[iPointLayer * alldata->points.nColsPoints + alldata->points.vraddValueCol];
            // pre-allocate the fitted vrad value at this [azimuth,elevation]
            yFitted[iPointIncluded] = 0.0f;
            // keep a record of which index was just included
            includedIndex[iPointIncluded] = iPointLayer;
            // raise the counter
            iPointIncluded += 1;

          }
        } // endfor (iPointLayer = 0; iPointLayer < nPointsLayer; iPointLayer++) {
        nPointsIncluded = iPointIncluded;

        // check if there are directions that have almost no observations
        // (as this makes the svdfit result really uncertain)
        hasGap = hasAzimuthGap(&pointsSelection[0], nPointsIncluded, alldata);

        if (alldata->options.fitVrad == TRUE) {

          if (hasGap == FALSE) {

            // ------------------------------------------------------------- //
            //                  dealias radial velocities                    //
            // ------------------------------------------------------------- //

            // dealias velocities if requested by user
            // only dealias in first pass (later passes for removing dual-PRF dealiasing errors,
            // which show smaller offsets than (2*nyquist velocity) and therefore are not
            // removed by dealiasing routine)
            // The condition nyquistMinUsed<maxNyquistDealias enforces that if all scans
            // have a higher Nyquist velocity than maxNyquistDealias, dealiasing is suppressed
            if (alldata->options.dealiasVrad && iPass == 0 && !recycleDealias) {
#ifdef FPRINTFON
              vol2bird_err_printf("dealiasing %i points for profile %i, layer %i ...\n",nPointsIncluded,iProfileType,iLayer+1);
#endif
              int result = dealias_points(&pointsSelection[0], alldata->misc.nDims, &yNyquist[0], alldata->misc.nyquistMin, &yObs[0], &yDealias[0],
                  nPointsIncluded);
              // store dealiased velocities in points array (for re-use when iPass>0)
              for (int i = 0; i < nPointsIncluded; i++) {
                alldata->points.points[includedIndex[i] * alldata->points.nColsPoints + alldata->points.vraddValueCol] = yDealias[i];
              }

              if (result == 0) {
                vol2bird_err_printf( "Warning, failed to dealias radial velocities");
              }
            }

            //print the dealiased values to stderr
            if (alldata->options.printDealias == TRUE) {
              printDealias(&pointsSelection[0], alldata->misc.nDims, &yNyquist[0], &yObs[0], &yDealias[0], nPointsIncluded, iProfileType, iLayer + 1,
                  iPass + 1);
            }

            // yDealias is initialized to yObs, so we can always run svdfit
            // on yDealias, even when not running a dealiasing routine
            yObsSvdFit = yDealias;

            // ------------------------------------------------------------- //
            //                       do the svdfit                           //
            // ------------------------------------------------------------- //

            chisq = svdfit(&pointsSelection[0], alldata->misc.nDims, &yObsSvdFit[0], &yFitted[0], nPointsIncluded, &parameterVector[0], &avar[0],
                alldata->misc.nParsFitted);

            if (chisq < alldata->constants.chisqMin) {
              // the standard deviation of the fit is too low, as in the case of overfit
              // reset parameter vector array elements to NAN and continue with the next layer
              parameterVector[0] = NAN;
              parameterVector[1] = NAN;
              parameterVector[2] = NAN;
              // FIXME: if this happens, profile fields are not updated from UNDETECT to NODATA
            } else {

              chi = sqrt(chisq);
              hSpeed = sqrt(pow(parameterVector[0], 2) + pow(parameterVector[1], 2));
              hDir = (atan2(parameterVector[0], parameterVector[1]) * RAD2DEG);

              if (hDir < 0) {
                hDir += 360.0f;
              }

              // if the fitted vrad value is more than 'absVDifMax' away from the corresponding
              // observed vrad value, set the gate's flagPositionVDifMax bit flag to 1, excluding the
              // gate in the second svdfit iteration
              updateFlagFieldsInPointsArray(&yObsSvdFit[0], &yFitted[0], &includedIndex[0], nPointsIncluded, &(alldata->points.points[0]), alldata);

            }

          } // endif (hasGap == FALSE)

        }; // endif (fitVrad == TRUE)

        //---------------------------------------------//
        //         Fill the profile arrays             //
        //---------------------------------------------//

        // always fill below profile fields, these never have a NODATA or UNDETECT value.
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 0] = iLayer * alldata->options.layerThickness;
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 1] = (iLayer + 1) * alldata->options.layerThickness;
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 8] = (float) hasGap;
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 10] = (float) nPointsIncluded;
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 13] = (float) nPointsIncludedZ;

        // fill below profile fields when (1) VVP fit was not performed because of azimuthal data gap
        // and (2) layer contains range gates within the volume sampled by the radar.
        if (hasGap && nPointsIncludedZ > alldata->constants.nPointsIncludedMin) {
          alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 2] = UNDETECT;
          alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 3] = UNDETECT;
          alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 4] = UNDETECT;
          alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 5] = UNDETECT;
          alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 6] = UNDETECT;
          alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 7] = UNDETECT;
          alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 9] = dbzAvg;
          alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 11] = reflectivity;
          alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 12] = birdDensity;
        }
        // case of valid fit, fill profile fields with VVP fit parameters
        if (!hasGap) {
          alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 2] = parameterVector[0];
          alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 3] = parameterVector[1];
          alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 4] = parameterVector[2];
          alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 5] = hSpeed;
          alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 6] = hDir;
          alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 7] = chi;
          alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 9] = dbzAvg;
          alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 11] = reflectivity;
          alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 12] = birdDensity;
        }

        free((void*) yObs);
        free((void*) yFitted);
        free((void*) yNyquist);
        free((void*) yDealias);
        free((void*) pointsSelection);
        free((void*) includedIndex);

      } // endfor (iPass = 0; iPass < nPasses; iPass++)
      // You need some of the results of iProfileType == 3 in order
      // to calculate iProfileType == 1, therefore iProfileType == 3 is executed first
      if (iProfileType == 3) {
        // NOTE: chi can have NAN or numeric value at this point
        // when NAN, below condition evaluates to FALSE, i.e. scatterersAreNotBirds is set to FALSE
        if (chi < alldata->options.stdDevMinBird) {
          alldata->misc.scatterersAreNotBirds[iLayer] = TRUE;
        } else {
          alldata->misc.scatterersAreNotBirds[iLayer] = FALSE;
        }
      }
      if (iProfileType == 1) {
        // set the bird density to zero if radial velocity stdev below threshold:
        if (alldata->misc.scatterersAreNotBirds[iLayer] == TRUE) {
          alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 12] = 0.0;
        }
      }

    } // endfor (iLayer = 0; iLayer < nLayers; iLayer++)

    if (alldata->options.printProfileVar == TRUE) {
      printProfile(alldata);
    }
    if (iProfileType == 1 && alldata->options.exportBirdProfileAsJSONVar == TRUE) {
      exportBirdProfileAsJSON(alldata);
    }

    // ---------------------------------------------------------------- //
    //       this next section is a bit ugly but it does the job        //
    // ---------------------------------------------------------------- //

    int iCopied = 0;
    int iColProfile;
    int iRowProfile;
    for (iRowProfile = 0; iRowProfile < alldata->profiles.nRowsProfile; iRowProfile++) {
      for (iColProfile = 0; iColProfile < alldata->profiles.nColsProfile; iColProfile++) {
        switch (iProfileType) {
        case 1:
          alldata->profiles.profile1[iCopied] = alldata->profiles.profile[iCopied];
          break;
        case 2:
          alldata->profiles.profile2[iCopied] = alldata->profiles.profile[iCopied];
          break;
        case 3:
          alldata->profiles.profile3[iCopied] = alldata->profiles.profile[iCopied];
          break;
        default:
          vol2bird_err_printf( "Something is wrong this should not happen.\n");
        }
        iCopied += 1;
      }
    }

  } // endfor (iProfileType = nProfileTypes; iProfileType > 0; iProfileType--)

} // vol2birdCalcProfiles


int vol2birdGetNColsProfile(vol2bird_t *alldata) {

    if (alldata->misc.initializationSuccessful==FALSE) {
        vol2bird_err_printf("You need to initialize vol2bird before you can use it. Aborting.\n");
        return -1;
    }
    return alldata->profiles.nColsProfile;
} // vol2birdGetNColsProfile


int vol2birdGetNRowsProfile(vol2bird_t *alldata) {   

    if (alldata->misc.initializationSuccessful==FALSE) {
        vol2bird_err_printf("You need to initialize vol2bird before you can use it. Aborting.\n");
        return -1;
    }
    return alldata->profiles.nRowsProfile;
} // vol2birdGetNColsProfile


float* vol2birdGetProfile(int iProfileType, vol2bird_t *alldata) {
    if (alldata->misc.initializationSuccessful==FALSE) {
        vol2bird_err_printf("You need to initialize vol2bird before you can use it. Aborting.\n");
        return (float *) NULL;
    }

    switch (iProfileType) {
        case 1 : 
            return &(alldata->profiles.profile1[0]);
        case 2 : 
            return &(alldata->profiles.profile2[0]);
        case 3 : 
            return &(alldata->profiles.profile3[0]);
        default :
            vol2bird_err_printf( "Something went wrong; behavior not implemented for given iProfileType.\n");
    }

    return (float *) NULL;
}


void vol2birdPrintIndexArrays(vol2bird_t* alldata) {
    
    if (alldata->misc.initializationSuccessful==FALSE) {
        vol2bird_err_printf("You need to initialize vol2bird before you can use it. Aborting.\n");
        return;
    }
    
    int iLayer;

    vol2bird_err_printf( "iLayer  iFrom   iTo     iTo-iFrom nWritten\n");
    for (iLayer = 0; iLayer < alldata->options.nLayers; iLayer++) {
        vol2bird_err_printf( "%7d %7d %7d %10d %8d\n",
            iLayer, 
            alldata->points.indexFrom[iLayer], 
            alldata->points.indexTo[iLayer], 
            alldata->points.indexTo[iLayer] - alldata->points.indexFrom[iLayer], 
            alldata->points.nPointsWritten[iLayer]);
    }
} // vol2birdPrintIndexArrays



void vol2birdPrintOptions(vol2bird_t* alldata) {
    
    // ------------------------------------------------------- //
    // this function prints vol2bird's configuration to stderr //
    // ------------------------------------------------------- //
    
    if (alldata->misc.initializationSuccessful==FALSE) {
        vol2bird_err_printf("You need to initialize vol2bird before you can use it. Aborting.\n");
        return;
    }

    vol2bird_err_printf("\n\nvol2bird configuration:\n\n");

    vol2bird_err_printf("%-25s = %f\n","absVDifMax",alldata->constants.absVDifMax);
    vol2bird_err_printf("%-25s = %f\n","azimMax",alldata->options.azimMax);
    vol2bird_err_printf("%-25s = %f\n","azimMin",alldata->options.azimMin);
    vol2bird_err_printf("%-25s = %f\n","birdRadarCrossSection",alldata->options.birdRadarCrossSection);
    vol2bird_err_printf("%-25s = %f\n","cellClutterFractionMax",alldata->constants.cellClutterFractionMax);
    vol2bird_err_printf("%-25s = %f\n","cellEtaMin",alldata->options.cellEtaMin);
    vol2bird_err_printf("%-25s = %f\n","cellStdDevMax",alldata->options.cellStdDevMax);
    vol2bird_err_printf("%-25s = %f\n","chisqMin",alldata->constants.chisqMin);
    vol2bird_err_printf("%-25s = %f\n","clutterValueMin",alldata->options.clutterValueMin);
    vol2bird_err_printf("%-25s = %f\n","etaMax",alldata->options.etaMax);
    vol2bird_err_printf("%-25s = %f\n","dbzThresMin",alldata->options.dbzThresMin);
    vol2bird_err_printf("%-25s = %s\n","dbzType",alldata->options.dbzType);
    vol2bird_err_printf("%-25s = %f\n","elevMax",alldata->options.elevMax);
    vol2bird_err_printf("%-25s = %f\n","elevMin",alldata->options.elevMin);
    vol2bird_err_printf("%-25s = %d\n","fitVrad",alldata->options.fitVrad);
    vol2bird_err_printf("%-25s = %f\n","fringeDist",alldata->constants.fringeDist);
    vol2bird_err_printf("%-25s = %f\n","layerThickness",alldata->options.layerThickness);
    vol2bird_err_printf("%-25s = %f\n","minNyquist",alldata->options.minNyquist);
    vol2bird_err_printf("%-25s = %f\n","areaCellMin",alldata->constants.areaCellMin);
    vol2bird_err_printf("%-25s = %d\n","nAzimNeighborhood",alldata->constants.nAzimNeighborhood);
    vol2bird_err_printf("%-25s = %d\n","nBinsGap",alldata->constants.nBinsGap);
    vol2bird_err_printf("%-25s = %d\n","nCountMin",alldata->constants.nCountMin);
    vol2bird_err_printf("%-25s = %d\n","nLayers",alldata->options.nLayers);
    vol2bird_err_printf("%-25s = %d\n","nObsGapMin",alldata->constants.nObsGapMin);
    vol2bird_err_printf("%-25s = %d\n","nPointsIncludedMin",alldata->constants.nPointsIncludedMin);
    vol2bird_err_printf("%-25s = %d\n","nRangNeighborhood",alldata->constants.nRangNeighborhood);
    vol2bird_err_printf("%-25s = %f\n","radarWavelength",alldata->options.radarWavelength);
    vol2bird_err_printf("%-25s = %f\n","rangeMax",alldata->options.rangeMax);
    vol2bird_err_printf("%-25s = %f\n","rangeMin",alldata->options.rangeMin);
    vol2bird_err_printf("%-25s = %f\n","rCellMax",alldata->misc.rCellMax);
    vol2bird_err_printf("%-25s = %f\n","refracIndex",alldata->constants.refracIndex);
    vol2bird_err_printf("%-25s = %d\n","requireVrad",alldata->options.requireVrad);
    vol2bird_err_printf("%-25s = %f\n","stdDevMinBird",alldata->options.stdDevMinBird);
    vol2bird_err_printf("%-25s = %c\n","useClutterMap",alldata->options.useClutterMap == TRUE ? 'T' : 'F');
    vol2bird_err_printf("%-25s = %f\n","vradMin",alldata->constants.vradMin);
    
    vol2bird_err_printf("\n\n");

}  // vol2birdPrintOptions





void vol2birdPrintPointsArray(vol2bird_t* alldata) {
    
    // ------------------------------------------------- //
    // this function prints the 'points' array to stderr //
    // ------------------------------------------------- //
    
    if (alldata->misc.initializationSuccessful==FALSE) {
        vol2bird_err_printf("You need to initialize vol2bird before you can use it. Aborting.\n");
        return;
    }

    int iPoint;
    
    vol2bird_err_printf( "iPoint    range     azim    elev         dbz        vrad    cell    gateCode   flags           nyquist     vradd        clut\n");
    
    for (iPoint = 0; iPoint < alldata->points.nRowsPoints * alldata->points.nColsPoints; iPoint+=alldata->points.nColsPoints) {
        
            char gateCodeStr[10];  // 9 bits plus 1 position for the null character '\0'
            
            printGateCode(&gateCodeStr[0], (int) alldata->points.points[iPoint + alldata->points.gateCodeCol]);
        
            vol2bird_err_printf( "  %6d",    iPoint/alldata->points.nColsPoints);
            vol2bird_err_printf( "  %6.1f",  alldata->points.points[iPoint + alldata->points.rangeCol]);
            vol2bird_err_printf( "  %6.2f",  alldata->points.points[iPoint + alldata->points.azimAngleCol]);
            vol2bird_err_printf( "  %6.2f",  alldata->points.points[iPoint + alldata->points.elevAngleCol]);
            vol2bird_err_printf( "  %10.2f", alldata->points.points[iPoint + alldata->points.dbzValueCol]);
            vol2bird_err_printf( "  %10.2f", alldata->points.points[iPoint + alldata->points.vradValueCol]);
            vol2bird_err_printf( "  %6.0f",  alldata->points.points[iPoint + alldata->points.cellValueCol]);
            vol2bird_err_printf( "  %8.0f",  alldata->points.points[iPoint + alldata->points.gateCodeCol]);
            vol2bird_err_printf( "  %12s",   gateCodeStr);
            vol2bird_err_printf( "  %10.2f", alldata->points.points[iPoint + alldata->points.nyquistCol]);
            vol2bird_err_printf( "  %10.2f", alldata->points.points[iPoint + alldata->points.vraddValueCol]);
            vol2bird_err_printf( "  %10.2f", alldata->points.points[iPoint + alldata->points.clutValueCol]);
            vol2bird_err_printf( "\n");
    }    
} // vol2birdPrintPointsArray




void vol2birdPrintPointsArraySimple(vol2bird_t* alldata) {
    
    // ------------------------------------------------- //
    // this function prints the 'points' array to stderr //
    // ------------------------------------------------- //
    
    int iPoint;
    
    vol2bird_err_printf( "iPoint  azim    elev    dbz         vrad        cell     flags     nyquist vradd\n");
    
    for (iPoint = 0; iPoint < alldata->points.nRowsPoints * alldata->points.nColsPoints; iPoint+=alldata->points.nColsPoints) {
                
            vol2bird_err_printf( "  %6d",    iPoint/alldata->points.nColsPoints);
            vol2bird_err_printf( "  %6.2f",  alldata->points.points[iPoint + alldata->points.azimAngleCol]);
            vol2bird_err_printf( "  %6.2f",  alldata->points.points[iPoint + alldata->points.elevAngleCol]);
            vol2bird_err_printf( "  %10.2f", alldata->points.points[iPoint + alldata->points.dbzValueCol]);
            vol2bird_err_printf( "  %10.2f", alldata->points.points[iPoint + alldata->points.vradValueCol]);
            vol2bird_err_printf( "  %6.0f",  alldata->points.points[iPoint + alldata->points.cellValueCol]);
            vol2bird_err_printf( "  %8.0f",  alldata->points.points[iPoint + alldata->points.gateCodeCol]);
            vol2bird_err_printf( "  %10.2f", alldata->points.points[iPoint + alldata->points.nyquistCol]);
            vol2bird_err_printf( "  %10.2f", alldata->points.points[iPoint + alldata->points.vraddValueCol]);
            vol2bird_err_printf( "\n");
    }    
} // vol2birdPrintPointsArray






void printProfile(vol2bird_t* alldata) {
    
    if (alldata->misc.initializationSuccessful==FALSE) {
        vol2bird_err_printf("You need to initialize vol2bird before you can use it. Aborting.\n");
        return;
    }

    vol2bird_err_printf("\n\nProfile type: %d\n",alldata->profiles.iProfileTypeLast);

    vol2bird_err_printf("altmin-altmax: [u         ,v         ,w         ]; "
                   "hSpeed  , hDir    , chi     , hasGap  , dbzAvg  ,"
                   " nPoints, eta         , rhobird nPointsZ \n");

    int iLayer;
    
    for (iLayer = alldata->options.nLayers - 1; iLayer >= 0; iLayer--) {
        
        vol2bird_err_printf("%6.0f-%-6.0f: [%10.2f,%10.2f,%10.2f]; %8.2f, "
        "%8.1f, %8.1f, %8c, %8.2f, %7.0f, %12.2f, %8.2f %5.f\n",
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile +  0],
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile +  1],
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile +  2],
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile +  3],
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile +  4],
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile +  5],
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile +  6],
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile +  7],
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile +  8] == TRUE ? 'T' : 'F',
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile +  9],
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 10],
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 11],
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 12],
        alldata->profiles.profile[iLayer * alldata->profiles.nColsProfile + 13]);
    }

    
} // printProfile()

// reads a polar volume from file and returns it as a RAVE polar volume object
// remember to release the polar volume object when done with it
PolarVolume_t* vol2birdGetVolume(char* filenames[], int nInputFiles, float rangeMax, int small){
    
    PolarVolume_t* volume = NULL;
    
    #ifdef IRIS
    // test whether the file is in IRIS format
    if (isIRIS(filenames[0])==0){
        volume = vol2birdGetIRISVolume(filenames, nInputFiles);
        goto done;
    }
    #endif

    // not a rave complient file, attempt to read the file with the RSL library instead
    #ifdef RSL
    if(RSL_filetype(filenames[0]) != UNKNOWN){
        if (nInputFiles > 1){
            vol2bird_err_printf("Multiple input files detected in RSL format. Only single polar volume file import supported, using file %s only.\n", filenames[0]);
        }
        volume = vol2birdGetRSLVolume(filenames[0], rangeMax, small);
        goto done;
    }
    #endif
    
    volume = vol2birdGetODIMVolume(filenames, nInputFiles);

    if (volume != NULL) {
      PolarVolume_sortByElevations(volume,1);
    }
done:
    return volume;
}

#ifdef IRIS
PolarVolume_t* vol2birdGetIRISVolume(char* filenames[], int nInputFiles) {
    // initialize a polar volume to return
    PolarVolume_t* output = NULL;
    // initialize helper volume and scan to store intermediate file reads
    PolarVolume_t* volume = NULL;
    PolarScan_t* scan = NULL;
    
    int outputInitialised = FALSE;
    
    // initialize the rave object type of filename
    int rot = Rave_ObjectType_UNDEFINED;

    int ret = 0;
    
    file_element_s* file_element_p = NULL;

    for (int i=0; i<nInputFiles; i++){
        // read the iris file
        file_element_p = readIRIS(filenames[i]);

        if(file_element_p == NULL){
            vol2bird_err_printf( "Warning: failed to read file %s in IRIS format, ignoring.\n", filenames[i]);
            continue;
        }

        rot = objectTypeFromIRIS(file_element_p);
        
        if (rot == Rave_ObjectType_UNDEFINED){
            vol2bird_err_printf( "Warning: unknown object type while reading file %s in IRIS format, ignoring.\n", filenames[i]);
            free_IRIS(&file_element_p);
            continue;
        }
        
        // start a new output volume object if we do not have one yet
        if (output == NULL){
            output = RAVE_OBJECT_NEW(&PolarVolume_TYPE);
            if (output == NULL) {
                RAVE_CRITICAL0("Error: failed to create polarvolume instance");
                goto done;
            }
        }
        
        if (rot == Rave_ObjectType_PVOL) {
            volume = RAVE_OBJECT_NEW(&PolarVolume_TYPE);
            if (volume == NULL) {
                RAVE_CRITICAL0("Error: failed to create polarvolume instance");
                goto done;
            }
            
            // read iris data into rave polar volume object
            ret = populateObject((RaveCoreObject*) volume, file_element_p);

            if( ret != 0) {
                vol2bird_err_printf( "Error: could not populate IRIS data into a polar volume object\n");
                goto done;
            }
            
            if (!outputInitialised){
                RAVE_OBJECT_RELEASE(output); //may have been initialized earlier above
                output = RAVE_OBJECT_CLONE(volume);
                RAVE_OBJECT_RELEASE(volume);
                outputInitialised = TRUE;
                continue;
            }
            
            for (int j=0; j<PolarVolume_getNumberOfScans(volume); j++){
                scan = PolarVolume_getScan(volume, j);
                PolarVolume_addScan(output, scan);
                RAVE_OBJECT_RELEASE(scan);
            }
            
            free_IRIS(&file_element_p);
            RAVE_OBJECT_RELEASE(volume);
        }
    
        if (rot == Rave_ObjectType_SCAN) {
            scan = RAVE_OBJECT_NEW(&PolarScan_TYPE);
            if (scan == NULL) {
                RAVE_CRITICAL0("Error: failed to create polarscan instance");
                goto done;
            }
            
            // read iris data into rave polar volume object
            ret = populateObject((RaveCoreObject*) scan, file_element_p);

            if( ret != 0) {
                vol2bird_err_printf( "Error: could not populate IRIS data into a polar scan object\n");
                goto done;
            }
            
            if (!outputInitialised){
                // copy essential root metadata to volume
                PolarVolume_setDate(output, PolarScan_getDate(scan));
                PolarVolume_setTime(output, PolarScan_getTime(scan));
                PolarVolume_setLatitude(output, PolarScan_getLatitude(scan));
                PolarVolume_setLongitude(output, PolarScan_getLongitude(scan));
                PolarVolume_setHeight(output, PolarScan_getHeight(scan));
                PolarVolume_setSource(output, PolarScan_getSource(scan));
                outputInitialised = TRUE;
            }
            
            PolarVolume_addScan(output, scan);
            free_IRIS(&file_element_p);
            RAVE_OBJECT_RELEASE(scan);
        }
    
    }

    done:
    
        // clean up
        if(file_element_p != NULL) {
            free_IRIS(&file_element_p);
        }
        RAVE_OBJECT_RELEASE(volume);            
        RAVE_OBJECT_RELEASE(scan);
        
        return output;
}
#endif

PolarVolume_t* vol2birdGetODIMVolume(char* filenames[], int nInputFiles) {
    // initialize a polar volume to return
    PolarVolume_t* output = NULL;
    // initialize helper volume and scan to store intermediate file reads
    PolarVolume_t* volume = NULL;
    PolarScan_t* scan = NULL;
    
    int outputInitialised = FALSE;
        
    // initialize the rave object type of filename
    int rot = Rave_ObjectType_UNDEFINED;

    for (int i=0; i<nInputFiles; i++){
        // read the file
        RaveIO_t* raveio = RaveIO_open(filenames[i], 0, NULL);

        if(raveio == NULL){
            vol2bird_err_printf( "Warning: failed to read file %s in ODIM format, ignoring.\n", filenames[i]);
            continue;
        }
        
        rot = RaveIO_getObjectType(raveio);

        if (!(rot == Rave_ObjectType_PVOL || rot == Rave_ObjectType_SCAN)) {
            vol2bird_err_printf( "Warning: no scan or volume found when reading file %s in ODIM format, ignoring.\n", filenames[i]);
            RAVE_OBJECT_RELEASE(raveio);
            continue;
        }
        
        // start a new output volume object if we do not have one yet
        if (output == NULL){
            output = RAVE_OBJECT_NEW(&PolarVolume_TYPE);
            if (output == NULL) {
                RAVE_CRITICAL0("Error: failed to create polarvolume instance");
                goto done;
            }
        }
        
        if (rot == Rave_ObjectType_PVOL) {
            // REMOVED BY AHE. Will be overwritten when getting object from raveio
            // volume = RAVE_OBJECT_NEW(&PolarVolume_TYPE);
            //if (volume == NULL) {
            //    RAVE_CRITICAL0("Error: failed to create polarvolume instance");
            //    goto done;
            //}
            
            // read ODIM data into rave polar volume object
            volume = (PolarVolume_t*) RaveIO_getObject(raveio);

            if( volume == NULL) {
                RAVE_OBJECT_RELEASE(raveio)
                RAVE_CRITICAL0("Error: could not populate ODIM data into a polarvolume object");
                goto done;
            }
            
            if (!outputInitialised){
                RAVE_OBJECT_RELEASE(output); // Added by AHE. Otherwise will loose output
                output = RAVE_OBJECT_CLONE(volume);
                RAVE_OBJECT_RELEASE(volume);
                RAVE_OBJECT_RELEASE(raveio)
                outputInitialised = TRUE;
                continue;
            }
            
            for (int j=0; j<PolarVolume_getNumberOfScans(volume); j++){
                scan = PolarVolume_getScan(volume, j);
                PolarVolume_addScan(output, scan);
                RAVE_OBJECT_RELEASE(scan);
            }
            
            RAVE_OBJECT_RELEASE(raveio);
            RAVE_OBJECT_RELEASE(volume);
        }
    
        if (rot == Rave_ObjectType_SCAN) {
            // Removed by AHE. Overwritten when getting object from raveio
            //scan = RAVE_OBJECT_NEW(&PolarScan_TYPE);
            //if (scan == NULL) {
            //    RAVE_CRITICAL0("Error: failed to create polarscan instance");
            //    goto done;
            //}
            
            // read iris data into rave polar volume object
            scan = (PolarScan_t*) RaveIO_getObject(raveio);

            if (scan == NULL) {
                RAVE_CRITICAL0("Error: could not populate ODIM data into a polar scan object");
                RAVE_OBJECT_RELEASE(raveio)
                goto done;
            }
            
            if (!outputInitialised){
                // copy essential root metadata to volume
                PolarVolume_setDate(output, PolarScan_getDate(scan));
                PolarVolume_setTime(output, PolarScan_getTime(scan));
                PolarVolume_setLatitude(output, PolarScan_getLatitude(scan));
                PolarVolume_setLongitude(output, PolarScan_getLongitude(scan));
                PolarVolume_setHeight(output, PolarScan_getHeight(scan));
                PolarVolume_setSource(output, PolarScan_getSource(scan));
                outputInitialised = TRUE;
            }
            
            PolarVolume_addScan(output, scan);
            RAVE_OBJECT_RELEASE(raveio);
            RAVE_OBJECT_RELEASE(scan);
        }
        RAVE_OBJECT_RELEASE(raveio);
    }

    done:
        // clean up
        RAVE_OBJECT_RELEASE(volume);            
        RAVE_OBJECT_RELEASE(scan);

        return output;
}


RaveIO_t* vol2birdIO_open(const char* filename)
{
  RaveIO_t* result = NULL;

  if (filename == NULL) {
    goto done;
  }

  result = RAVE_OBJECT_NEW(&RaveIO_TYPE);
  if (result == NULL) {
    RAVE_CRITICAL0("Failed to create raveio instance");
    goto done;
  }

  if (!RaveIO_setFilename(result, filename)) {
    RAVE_CRITICAL0("Failed to set filename");
    RAVE_OBJECT_RELEASE(result);
    goto done;
  }

  if (!RaveIO_load(result, 0, NULL)) {
    RAVE_WARNING0("Failed to load file");
    RAVE_OBJECT_RELEASE(result);
    goto done;
  }

done:
  return result;
}


// loads configuration data in the alldata struct
#ifndef NOCONFUSE

int vol2birdLoadConfig(vol2bird_t* alldata, const char* optionsFile) {

    alldata->misc.loadConfigSuccessful = FALSE;

    const char * optsConfFilename = getenv(OPTIONS_CONF);
    if (optsConfFilename == NULL) {
         if(optionsFile == NULL){
            optsConfFilename = OPTIONS_FILE;
        }
        else{
            optsConfFilename = optionsFile;
        }
    }
    else{
        vol2bird_err_printf( "Searching user configuration file '%s' specified in environmental variable '%s'\n",optsConfFilename,OPTIONS_CONF);
    }

    if (readUserConfigOptions(&(alldata->cfg), optsConfFilename) != 0) {
        vol2bird_err_printf( "An error occurred while reading the user configuration file '%s'.\n", optsConfFilename);
        return -1; 
    }
   
    cfg_t** cfg = &(alldata->cfg);

    // ------------------------------------------------------------- //
    //              vol2bird options from options.conf               //
    // ------------------------------------------------------------- //

    alldata->options.azimMax = cfg_getfloat(*cfg, "AZIMMAX");
    alldata->options.azimMin = cfg_getfloat(*cfg, "AZIMMIN");
    alldata->options.layerThickness = cfg_getfloat(*cfg, "HLAYER");
    alldata->options.nLayers = cfg_getint(*cfg, "NLAYER");
    alldata->options.rangeMax = cfg_getfloat(*cfg, "RANGEMAX");
    alldata->options.rangeMin = cfg_getfloat(*cfg, "RANGEMIN");
    alldata->options.elevMax = cfg_getfloat(*cfg, "ELEVMAX");
    alldata->options.elevMin = cfg_getfloat(*cfg, "ELEVMIN");
    alldata->options.radarWavelength = cfg_getfloat(*cfg, "RADAR_WAVELENGTH_CM");
    alldata->options.useClutterMap = cfg_getbool(*cfg,"USE_CLUTTERMAP");
    alldata->options.clutterValueMin = cfg_getfloat(*cfg, "CLUTTERVALUEMIN");
    strcpy(alldata->options.clutterMap,cfg_getstr(*cfg,"CLUTTERMAP"));
    alldata->options.printDbz = cfg_getbool(*cfg,"PRINT_DBZ");
    alldata->options.printDealias = cfg_getbool(*cfg,"PRINT_DEALIAS");
    alldata->options.printVrad = cfg_getbool(*cfg,"PRINT_VRAD");
    alldata->options.printRhohv = cfg_getbool(*cfg,"PRINT_RHOHV");
    alldata->options.printTex = cfg_getbool(*cfg,"PRINT_TEXTURE");
    alldata->options.printCell = cfg_getbool(*cfg,"PRINT_CELL");
    alldata->options.printCellProp = cfg_getbool(*cfg,"PRINT_CELL_PROP");
    alldata->options.printClut = cfg_getbool(*cfg,"PRINT_CLUT");
    alldata->options.printOptions = cfg_getbool(*cfg,"PRINT_OPTIONS");
    alldata->options.printProfileVar = cfg_getbool(*cfg,"PRINT_PROFILE");
    alldata->options.printPointsArray = cfg_getbool(*cfg,"PRINT_POINTS_ARRAY");
    alldata->options.fitVrad = cfg_getbool(*cfg,"FIT_VRAD");
    alldata->options.exportBirdProfileAsJSONVar = cfg_getbool(*cfg,"EXPORT_BIRD_PROFILE_AS_JSON"); 
    alldata->options.minNyquist = cfg_getfloat(*cfg,"MIN_NYQUIST_VELOCITY");
    alldata->options.maxNyquistDealias = cfg_getfloat(*cfg,"MAX_NYQUIST_DEALIAS");
    alldata->options.birdRadarCrossSection = cfg_getfloat(*cfg,"SIGMA_BIRD");
    alldata->options.cellStdDevMax = cfg_getfloat(*cfg,"STDEV_CELL");
    alldata->options.stdDevMinBird = cfg_getfloat(*cfg,"STDEV_BIRD");
    alldata->options.etaMax = cfg_getfloat(*cfg,"ETAMAX");
    alldata->options.cellEtaMin = cfg_getfloat(*cfg,"ETACELL");
    strcpy(alldata->options.dbzType,cfg_getstr(*cfg,"DBZTYPE"));
    alldata->options.requireVrad = cfg_getbool(*cfg,"REQUIRE_VRAD");
    alldata->options.dealiasVrad = cfg_getbool(*cfg,"DEALIAS_VRAD");
    alldata->options.dealiasRecycle = cfg_getbool(*cfg,"DEALIAS_RECYCLE");
    alldata->options.dualPol = cfg_getbool(*cfg,"DUALPOL");
    alldata->options.singlePol = cfg_getbool(*cfg,"SINGLEPOL");
    alldata->options.dbzThresMin = cfg_getfloat(*cfg,"DBZMIN");
    alldata->options.rhohvThresMin = cfg_getfloat(*cfg,"RHOHVMIN");
    alldata->options.resample = cfg_getbool(*cfg,"RESAMPLE");
    alldata->options.resampleRscale = cfg_getfloat(*cfg,"RESAMPLE_RSCALE");
    alldata->options.resampleNbins = cfg_getint(*cfg,"RESAMPLE_NBINS");
    alldata->options.resampleNrays = cfg_getint(*cfg,"RESAMPLE_NRAYS");
    alldata->options.mistNetNElevs = cfg_size(*cfg, "MISTNET_ELEVS");
    for(int i=0; i<alldata->options.mistNetNElevs; i++){
        alldata->options.mistNetElevs[i] = cfg_getnfloat(*cfg, "MISTNET_ELEVS",i);
    }
    alldata->options.mistNetElevsOnly = cfg_getbool(*cfg, "MISTNET_ELEVS_ONLY");
    alldata->options.useMistNet = cfg_getbool(*cfg, "USE_MISTNET");
    strcpy(alldata->options.mistNetPath,cfg_getstr(*cfg,"MISTNET_PATH"));


    // ------------------------------------------------------------- //
    //              vol2bird options from constants.h                //
    // ------------------------------------------------------------- //

    alldata->constants.areaCellMin = AREACELL;
    alldata->constants.cellClutterFractionMax = CLUTPERCCELL;
    alldata->constants.chisqMin = CHISQMIN;
    alldata->constants.fringeDist = FRINGEDIST;
    alldata->constants.nBinsGap = NBINSGAP;
    alldata->constants.nPointsIncludedMin = NDBZMIN;
    alldata->constants.nNeighborsMin = NEIGHBORS;
    alldata->constants.nObsGapMin = NOBSGAPMIN;
    alldata->constants.nAzimNeighborhood = NTEXBINAZIM;
    alldata->constants.nRangNeighborhood = NTEXBINRANG;
    alldata->constants.nCountMin = NTEXMIN; 
    alldata->constants.refracIndex = REFRACTIVE_INDEX_OF_WATER;
    alldata->constants.absVDifMax = VDIFMAX;
    alldata->constants.vradMin = VRADMIN;

    // ------------------------------------------------------------- //
    //       some other variables, derived from user options         //
    // ------------------------------------------------------------- //

    alldata->misc.rCellMax = alldata->options.rangeMax + RCELLMAX_OFFSET;
    alldata->misc.nDims = 2;
    alldata->misc.nParsFitted = 3;

    // the following settings depend on wavelength, will be set in Vol2birdSetup
    alldata->misc.dbzFactor = NAN;
    alldata->misc.dbzMax = NAN;
    alldata->misc.cellDbzMin = NAN;

    alldata->misc.loadConfigSuccessful = TRUE;

    return 0;

}

#endif

char* get_radar_name(char* source) {
    char* radarName = NULL;
    char* p = strstr(source, "RAD:");
    if (p != NULL) {
        p += strlen("RAD:");
        radarName = p;
        char* end = strchr(p, ',');
        if (end != NULL) {
            *end = '\0';
        }
    }

    if (radarName == NULL) {
        radarName = "UNKNOWN";
    }

    return radarName;
}

//int vol2birdSetUp(PolarVolume_t* volume, cfg_t** cfg, vol2bird_t* alldata) {
int vol2birdSetUp(PolarVolume_t* volume, vol2bird_t* alldata) {
    
    alldata->misc.initializationSuccessful = FALSE;
    
    alldata->misc.vol2birdSuccessful = TRUE;

    vol2bird_printf("Running vol2birdSetUp\n");

    if (alldata->misc.loadConfigSuccessful == FALSE){
        vol2bird_err_printf("Vol2bird configuration not loaded. Run vol2birdLoadConfig prior to vol2birdSetup\n");
        return -1;
    }
 
    radarName = get_radar_name(PolarVolume_getSource(volume))
    alldata -> misc.radarName = radarName;

    // reading radar wavelength from polar volume attribute
    // if present, overwrite options.radarWavelength with the value found.
    double wavelength = PolarVolume_getWavelength(volume);
    if (wavelength > 0){
        alldata->options.radarWavelength = wavelength;
    }
    else{
        vol2bird_err_printf("Warning: radar wavelength not stored in polar volume. Using user-defined value of %f cm ...\n", alldata->options.radarWavelength);
    }
    
    // Now that we have the wavelength, calculate its derived quantities:
    // dbzFactor is the proportionality constant between reflectivity factor Z in mm^6/m^3
    // and reflectivity eta in cm^2/km^3, when radarWavelength in cm:
    alldata->misc.dbzFactor = (pow(alldata->constants.refracIndex,2) * 1000 * pow(PI,5))/pow(alldata->options.radarWavelength,4);
    alldata->misc.dbzMax = 10*log(alldata->options.etaMax / alldata->misc.dbzFactor)/log(10);
    alldata->misc.cellDbzMin = 10*log(alldata->options.cellEtaMin / alldata->misc.dbzFactor)/log(10);
    // if stdDevMinBird not set by STDEV_BIRD in options.conf, initialize it depending on wavelength:
    // was initialized to -FLT_MAX, i.e. negative
    if (alldata->options.stdDevMinBird < 0){
        if (alldata->options.radarWavelength < 7.5){
            //C-band default:
            alldata->options.stdDevMinBird = STDEV_BIRD;
        }
        else{
            //S-band default:
            alldata->options.stdDevMinBird = STDEV_BIRD_S;
        }
    }
    // Extract the vcp attribute if present (i.e. NEXRAD only)
    RaveAttribute_t *attr;
    long vcp;
    int result = 0;
    attr = PolarVolume_getAttribute(volume, "how/vcp");
    if (attr != (RaveAttribute_t *) NULL) result = RaveAttribute_getLong(attr, &vcp);
    if (result > 0){
        alldata->misc.vcp = (int) vcp;
    }
    else{
        alldata->misc.vcp = 0;
    }
    RAVE_OBJECT_RELEASE(attr);
    // ------------------------------------------------------------- //
    //                 determine which scans to use                  //
    // ------------------------------------------------------------- //
    
    vol2birdScanUse_t* scanUse=NULL;

    scanUse = determineScanUse(volume, alldata);

    if (!alldata->options.dealiasVrad && alldata->misc.nyquistMinUsed < alldata->options.maxNyquistDealias){   
       vol2bird_err_printf("Warning: Nyquist velocity below maxNyquistDealias threshold was found (%f<%f), consider dealiasing.\n",alldata->misc.nyquistMinUsed,alldata->options.maxNyquistDealias);
    }    
    if (alldata->options.dealiasVrad && alldata->misc.nyquistMinUsed > alldata->options.maxNyquistDealias){
       alldata->options.dealiasVrad=FALSE;
       //Maybe you want to give a warning here, but that generates a lot of warnings as this situation is common ...
    }
    if (alldata->options.dealiasVrad){
        vol2bird_err_printf("Warning: radial velocities will be dealiased...\n");
    }

    // ------------------------------------------------------------- //
    //     store all options and constants in task_args string       //
    // ------------------------------------------------------------- //

    // the radar wavelength setting is read from taken from the volume object
    // if a wavelength attribute is present. Therefore the task_args string is
    // set here and not in vol2birdLoadConfig(), which has no access to the volume    
    
    //FIXME: add mistNetNElevs (mistnet elevations) to the task_args string
    sprintf(alldata->misc.task_args,
        "azimMax=%f,azimMin=%f,layerThickness=%f,nLayers=%i,rangeMax=%f,"
        "rangeMin=%f,elevMax=%f,elevMin=%f,radarWavelength=%f,"
        "useClutterMap=%i,clutterMap=%s,fitVrad=%i,exportBirdProfileAsJSONVar=%i,"
        "minNyquist=%f,maxNyquistDealias=%f,birdRadarCrossSection=%f,stdDevMinBird=%f,"
        "cellEtaMin=%f,etaMax=%f,dbzType=%s,requireVrad=%i,"
        "dealiasVrad=%i,dealiasRecycle=%i,dualPol=%i,singlePol=%i,rhohvThresMin=%f,"
        "resample=%i,resampleRscale=%f,resampleNbins=%i,resampleNrays=%i,"
        "mistNetNElevs=%i,mistNetElevsOnly=%i,useMistNet=%i,mistNetPath=%s,"
    
        "areaCellMin=%f,cellClutterFractionMax=%f,"
        "chisqMin=%f,clutterValueMin=%f,dbzThresMin=%f,"
        "fringeDist=%f,nBinsGap=%i,nPointsIncludedMin=%i,nNeighborsMin=%i,"
        "nObsGapMin=%i,nAzimNeighborhood=%i,nRangNeighborhood=%i,nCountMin=%i,"
        "refracIndex=%f,cellStdDevMax=%f,absVDifMax=%f,vradMin=%f",

        alldata->options.azimMax,
        alldata->options.azimMin,
        alldata->options.layerThickness,
        alldata->options.nLayers,
        alldata->options.rangeMax,
        alldata->options.rangeMin,
        alldata->options.elevMax,
        alldata->options.elevMin,
        alldata->options.radarWavelength,
        alldata->options.useClutterMap,
        alldata->options.clutterMap,
        alldata->options.fitVrad,
        alldata->options.exportBirdProfileAsJSONVar,
        alldata->options.minNyquist,
        alldata->options.maxNyquistDealias,
        alldata->options.birdRadarCrossSection,
        alldata->options.stdDevMinBird,
        alldata->options.cellEtaMin,
        alldata->options.etaMax,
        alldata->options.dbzType,
        alldata->options.requireVrad,
        alldata->options.dealiasVrad,
        alldata->options.dealiasRecycle,
        alldata->options.dualPol,
	    alldata->options.singlePol,
        alldata->options.rhohvThresMin,
        alldata->options.resample,
        alldata->options.resampleRscale,
        alldata->options.resampleNbins,
        alldata->options.resampleNrays,
        alldata->options.mistNetNElevs,
        alldata->options.mistNetElevsOnly,
        alldata->options.useMistNet,
        alldata->options.mistNetPath,

        alldata->constants.areaCellMin,
        alldata->constants.cellClutterFractionMax,
        alldata->constants.chisqMin,
        alldata->options.clutterValueMin,
        alldata->options.dbzThresMin,
        alldata->constants.fringeDist,
        alldata->constants.nBinsGap,
        alldata->constants.nPointsIncludedMin,
        alldata->constants.nNeighborsMin,
        alldata->constants.nObsGapMin,
        alldata->constants.nAzimNeighborhood,
        alldata->constants.nRangNeighborhood,
        alldata->constants.nCountMin,
        alldata->constants.refracIndex,
        alldata->options.cellStdDevMax,
        alldata->constants.absVDifMax,
        alldata->constants.vradMin
    );
   
    if (scanUse == (vol2birdScanUse_t*) NULL){
        vol2bird_err_printf( "Error: no valid scans found in polar volume, aborting ...\n");
        return -1;
    }

    // Print warning missing rain specification
    if(!alldata->options.singlePol && !alldata->options.dualPol){
        vol2bird_err_printf("Warning: neither single- nor dual-polarization precipitation filter selected by user, continuing in SINGLE polarization mode\n");
		alldata->options.singlePol = TRUE;
    }

    // Disable single pol rain filtering for S-band data when dual-pol moments are available
    if(alldata->options.radarWavelength > 7.5 && alldata->options.singlePol && alldata->options.dualPol){
        vol2bird_err_printf("Warning: disabling single-polarization precipitation filter for S-band data, continuing in DUAL polarization mode\n");
		alldata->options.singlePol = FALSE;
    }

    // Print warning for S-band in single pol mode
    if(alldata->options.radarWavelength > 7.5 && !alldata->options.dualPol){
        vol2bird_err_printf("Warning: using experimental SINGLE polarization mode on S-band data, results may be unreliable!\n");
    }

    // Print warning for MistNet mode
    if(alldata->options.useMistNet && (alldata->options.dualPol || alldata->options.singlePol)){
        vol2bird_err_printf("Warning: using MistNet, disabling other segmentation methods\n");
        alldata->options.singlePol = FALSE;
        alldata->options.dualPol = FALSE;
    }

    // check that we are requesting the right number of elevation scans for MistNet segmentation model
    if(alldata->options.mistNetNElevs != MISTNET_N_ELEV){
        vol2bird_err_printf( "Error: MistNet segmentation model expects %i elevations, but %i are specified.\n", MISTNET_N_ELEV, alldata->options.mistNetNElevs);
        return -1;
    }
    
    // check that MistNet segmentation model can be found on disk
    if(alldata->options.useMistNet && !isRegularFile(alldata->options.mistNetPath)){
        vol2bird_err_printf( "Error: MistNet segmentation model '%s' not found.\n", alldata->options.mistNetPath);
        return -1;
    }

    // Print warning for mistnet mode
    if(alldata->options.useMistNet && alldata->options.radarWavelength < 7.5){
        vol2bird_err_printf("Warning: MistNet segmentation model has been trained on S-band data, results at other radar wavelengths may be unreliable!\n");
    }
	
    // ------------------------------------------------------------- //
    //             lists of indices into the 'points' array:         //
    //          where each altitude layer's data starts and ends     //
    // ------------------------------------------------------------- //

    int iLayer;
    
    // pre-allocate the list with start-from indexes for each 
    // altitude bin in the profile
    alldata->points.indexFrom = (int*) malloc(sizeof(int) * alldata->options.nLayers);
    if (alldata->points.indexFrom == NULL) {
        vol2bird_err_printf("Error pre-allocating array 'indexFrom'\n");
        return -1;
    }
    for (iLayer = 0; iLayer < alldata->options.nLayers; iLayer++) {
        alldata->points.indexFrom[iLayer] = 0;
    }

    // pre-allocate the list with end-before indexes for each 
    // altitude bin in the profile
    alldata->points.indexTo = (int*) malloc(sizeof(int) * alldata->options.nLayers);
    if (alldata->points.indexTo == NULL) {
        vol2bird_err_printf("Error pre-allocating array 'indexTo'\n");
        return -1;
    }
    for (iLayer = 0; iLayer < alldata->options.nLayers; iLayer++) {
        alldata->points.indexTo[iLayer] = 0;
    }

    // pre-allocate the list containing TRUE or FALSE depending on the 
    // results of calculating iProfileType == 3, which are needed when
    // calculating iProfileType == 1
    alldata->misc.scatterersAreNotBirds = (int*) malloc(sizeof(int) * alldata->options.nLayers);
    if (alldata->misc.scatterersAreNotBirds == NULL) {
        vol2bird_err_printf("Error pre-allocating array 'scatterersAreNotBirds'\n");
        return -1;
    }
    for (iLayer = 0; iLayer < alldata->options.nLayers; iLayer++) {
        alldata->misc.scatterersAreNotBirds[iLayer] = -1;
    }


    // for each altitude layer, you need to remember how many points 
    // were already written. This information is stored in the 
    // 'nPointsWritten' array
    alldata->points.nPointsWritten = (int*) malloc(sizeof(int) * alldata->options.nLayers);
    if (alldata->points.nPointsWritten == NULL) {
        vol2bird_err_printf("Error pre-allocating array 'nPointsWritten'\n");
        return -1;
    }
    for (iLayer = 0; iLayer < alldata->options.nLayers; iLayer++) {
        alldata->points.nPointsWritten[iLayer] = 0;
    }


    // ------------------------------------------------------------- //
    //               information about the 'points' array            //
    // ------------------------------------------------------------- //

    alldata->points.nColsPoints = 10;
    alldata->points.nRowsPoints = detSvdfitArraySize(volume, scanUse, alldata);

    alldata->points.rangeCol = 0;
    alldata->points.azimAngleCol = 1;
    alldata->points.elevAngleCol = 2;
    alldata->points.dbzValueCol = 3;
    alldata->points.vradValueCol = 4;
    alldata->points.cellValueCol = 5;
    alldata->points.gateCodeCol = 6;
    alldata->points.nyquistCol = 7;
    alldata->points.vraddValueCol = 8;
    alldata->points.clutValueCol = 9;

    // pre-allocate the 'points' array (note it has 'nColsPoints'
    // pseudo-columns)
    alldata->points.points = (float*) malloc(sizeof(float) * alldata->points.nRowsPoints * alldata->points.nColsPoints);
    if (alldata->points.points == NULL) {
        vol2bird_err_printf("Error pre-allocating array 'points'.\n");
        return -1;
    }

    int iRowPoints;
    int iColPoints;
        
    for (iRowPoints = 0; iRowPoints < alldata->points.nRowsPoints; iRowPoints++) {
        for (iColPoints = 0; iColPoints < alldata->points.nColsPoints; iColPoints++) {
            alldata->points.points[iRowPoints*alldata->points.nColsPoints + iColPoints] = NAN;
        }
    }

    // information about the flagfields of 'gateCode'
    
    alldata->flags.flagPositionStaticClutter = 0;
    alldata->flags.flagPositionDynamicClutter = 1;
    alldata->flags.flagPositionDynamicClutterFringe = 2;
    alldata->flags.flagPositionVradMissing = 3;
    alldata->flags.flagPositionDbzTooHighForBirds = 4;
    alldata->flags.flagPositionVradTooLow = 5;
    alldata->flags.flagPositionVDifMax = 6;
    alldata->flags.flagPositionAzimOutOfRange = 7;

    // segment precipitation using Mistnet deep convolutional neural net
#ifdef MISTNET
#ifdef VOL2BIRD_R
    if (check_mistnet_loaded_c()) {
#endif
      if(alldata->options.useMistNet){
        vol2bird_err_printf("Running segmentScansUsingMistnet.\n");
        int result = segmentScansUsingMistnet(volume, scanUse, alldata);
        if (result < 0) return -1;
      }
#ifdef VOL2BIRD_R      
    }
#endif    
#endif

    // construct the 'points' array
    constructPointsArray(volume, scanUse, alldata);

    // classify the gates based on the data in 'points'
    classifyGatesSimple(alldata);


    // ------------------------------------------------------------- //
    //              information about the 'profile' array            //
    // ------------------------------------------------------------- //

    alldata->profiles.nProfileTypes = 3;
    alldata->profiles.nRowsProfile = alldata->options.nLayers;
    alldata->profiles.nColsProfile = 14; 
    
    // pre-allocate the array holding any profiled data (note it has 
    // 'nColsProfile' pseudocolumns):
    alldata->profiles.profile = (float*) malloc(sizeof(float) * alldata->profiles.nRowsProfile * alldata->profiles.nColsProfile);
    if (alldata->profiles.profile == NULL) {
        vol2bird_err_printf("Error pre-allocating array 'profile'.\n");
        return -1;
    }

    // these next three variables are a quick fix
    alldata->profiles.profile1 = (float*) malloc(sizeof(float) * alldata->profiles.nRowsProfile * alldata->profiles.nColsProfile);
    if (alldata->profiles.profile1 == NULL) {
        vol2bird_err_printf("Error pre-allocating array 'profile1'.\n");
        return -1;
    }
    alldata->profiles.profile2 = (float*) malloc(sizeof(float) * alldata->profiles.nRowsProfile * alldata->profiles.nColsProfile);
    if (alldata->profiles.profile2 == NULL) {
        vol2bird_err_printf("Error pre-allocating array 'profile2'.\n");
        return -1;
    }
    alldata->profiles.profile3 = (float*) malloc(sizeof(float) * alldata->profiles.nRowsProfile * alldata->profiles.nColsProfile);
    if (alldata->profiles.profile3 == NULL) {
        vol2bird_err_printf("Error pre-allocating array 'profile3'.\n");
        return -1;
    }

    int iRowProfile;
    int iColProfile;
        
    for (iRowProfile = 0; iRowProfile < alldata->profiles.nRowsProfile; iRowProfile++) {
        for (iColProfile = 0; iColProfile < alldata->profiles.nColsProfile; iColProfile++) {
            alldata->profiles.profile[iRowProfile*alldata->profiles.nColsProfile + iColProfile] = NODATA;
            alldata->profiles.profile1[iRowProfile*alldata->profiles.nColsProfile + iColProfile] = NODATA;
            alldata->profiles.profile2[iRowProfile*alldata->profiles.nColsProfile + iColProfile] = NODATA;
            alldata->profiles.profile3[iRowProfile*alldata->profiles.nColsProfile + iColProfile] = NODATA;
        }
    }

    alldata->profiles.iProfileTypeLast = -1;


 
    // ------------------------------------------------------------- //
    //              initialising rave profile fields                 //
    // ------------------------------------------------------------- //

    alldata->vp = RAVE_OBJECT_NEW(&VerticalProfile_TYPE);
        
    alldata->misc.initializationSuccessful = TRUE;

    if (alldata->options.printOptions == TRUE) {
        vol2birdPrintOptions(alldata);
    }
    
    if (alldata->options.printPointsArray == TRUE) {

        vol2birdPrintIndexArrays(alldata);
        vol2birdPrintPointsArray(alldata);

    }

    if (scanUse != NULL) free(scanUse);

    return 0;

} // vol2birdSetUp


void vol2birdTearDown(vol2bird_t* alldata) {
    
    // ---------------------------------------------------------- //
    // free the memory that was previously allocated for vol2bird //
    // ---------------------------------------------------------- //

    if (alldata->misc.initializationSuccessful==FALSE) {
        vol2bird_err_printf("You need to initialize vol2bird before you can use it. Aborting.\n");
        return;
    }

    // free the points array, the indexes into it, the counters, as well
    // as the profile data array
    free((void*) alldata->points.points);
    free((void*) alldata->profiles.profile);
    free((void*) alldata->profiles.profile1);
    free((void*) alldata->profiles.profile2);
    free((void*) alldata->profiles.profile3);
    free((void*) alldata->points.indexFrom);
    free((void*) alldata->points.indexTo);
    free((void*) alldata->points.nPointsWritten);
    free((void*) alldata->misc.scatterersAreNotBirds);
   
    // free all rave fields
    RAVE_OBJECT_RELEASE(alldata->vp);
 
    // free the memory that holds the user configurable options
#ifndef NOCONFUSE    
    cfg_free(alldata->cfg);
#endif    
    // reset this variable to its initial value
    alldata->misc.initializationSuccessful = FALSE;
    alldata->misc.loadConfigSuccessful = FALSE;

} // vol2birdTearDown





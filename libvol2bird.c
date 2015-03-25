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


#include <stdio.h>
#include <confuse.h>
#include <stdlib.h>
#include <math.h>
#include <vertical_profile.h>

#include "polarvolume.h"
#include "polarscan.h"
#include "libvol2bird.h"
#include "libsvdfit.h"
#include "constants.h"


// non-public function prototypes (local to this file/translation unit)

static int analyzeCells(const unsigned char *dbzImage, const unsigned char *vradImage,
                        const unsigned char *texImage, const unsigned char *clutterImage, int *cellImage,
                        const SCANMETA *dbzMeta, const SCANMETA *vradMeta, const SCANMETA *texMeta, const SCANMETA *clutterMeta,
                        const int nCells, const int useStaticClutterData);

static float calcDist(const int range1, const int azim1, const int range2, const int azim2, const float rscale, const float ascale);

static void calcTexture(unsigned char *texImage, const unsigned char *vradImage, const unsigned char *dbzImage,
                        const SCANMETA *texMeta, const SCANMETA *vradMeta, const SCANMETA *dbzMeta);

static void classifyGatesSimple(void);

static int constructorInt(SCANMETA* meta, int* image, PolarScan_t* scan, const int nGlobal, const int initValue);

static int constructorUChar(SCANMETA* meta, unsigned char* image, PolarScan_t* scan, const int nGlobal, const unsigned char initValue);

static void constructPointsArray(PolarVolume_t* volume);

static int detNumberOfGates(const int iLayer, const float rangeScale, const float elevAngle,
                            const int nRang, const int nAzim, const float radarHeight);

static int detSvdfitArraySize(PolarVolume_t* volume);

static int findWeatherCells(const unsigned char *dbzImage, int *cellImage, const SCANMETA *dbzMeta);

static int findNearbyGateIndex(const int nAzimParent, const int nRangParent, const int iParent,
                               const int nAzimChild,  const int nRangChild,  const int iChild);

static void fringeCells(int *cellImage,int nRang, int nAzim, float aScale, float rScale);

static int getListOfSelectedGates(const SCANMETA* vradMeta, const unsigned char *vradImage,
                                  const SCANMETA* dbzMeta, const unsigned char *dbzImage,
                                  const int *cellImage,
                                  const float altitudeMin, const float altitudeMax,
                                  float* points, int iPoint);

static int hasAzimuthGap(const float *points, const int nPoints);

static int includeGate(const int iProfileType, const unsigned int gateCode);

static int mapDataFromRave(PolarScan_t* scan, SCANMETA *meta, 
                           unsigned char *values, char *paramStr);

static void mapDataToRave(void);

static void printGateCode(char* flags, const unsigned int gateCode);

static void printImageInt(const SCANMETA* meta, const int* imageInt);

static void printImageUChar(const SCANMETA* meta, const unsigned char* imageUChar);

static int printMeta(const SCANMETA* meta, const char* varName);

static void printProfile(void);

static void sortCells(CELLPROP *cellProp, const int nCells);

static void updateFlagFieldsInPointsArray(const float* yObs, const float* yFitted, const int* includedIndex, 
                                          const int nPointsIncluded, float* points);
static int updateMap(int *cellImage, const int nGlobal, CELLPROP *cellProp, const int nCells);





// ------------------------------------------------------------- //
//              vol2bird options from options.conf               //
// ------------------------------------------------------------- //

// the configuration options 
static cfg_t* cfg;

// the user can specify to exclude gates based on their azimuth;
// the maximum is set by azimMax
static float azimMax;

// the user can specify to exclude gates based on their azimuth;
// the minimum is set by azimMin
static float azimMin;

// the thickness in meters of a layer in the altitude profile
static float layerThickness;

// the number of layers in an altitude profile
static int nLayers;

// whether or not to print each scan's raw dbz values to stderr
static int printDbz;

// whether or not to print each scan's raw vrad values to stderr
static int printVrad;

// whether or not to print each scan's raw tex values to stderr
static int printTex;

// whether or not to print each scan's cell values to stderr
static int printCell;

// whether or not to print each scan's cell properties to stderr
static int printCellProp;

// whether or not to print each scan's clutter values to stderr
static int printClut;

// whether or not to print vol2bird's configuration options to stderr
static int printOptions;

// whether or not to print vol2bird's profiled data to stderr
static int printProfileVar;

// whether or not to print the 'points' array
static int printPointsArray;

// whether or not a model should be fitted to the observed vrad values
static int fitVrad;

// the range beyond which observations are excluded when constructing 
// the altitude profile
static float rangeMax;

// the range below which observations are excluded when constructing 
// the altitude profile
static float rangeMin;

// the wavelength of the radar in units of centimeter
static float radarWavelength;

// whether clutter data is used
static int useStaticClutterData;




// ------------------------------------------------------------- //
//              vol2bird options from constants.h                //
// ------------------------------------------------------------- //

// when analyzing cells, AREAMIN determines the minimum size of a 
// cell to be considered in the rest of the analysis
static int nGatesCellMin;

// Maximum percentage of clutter in a cell 
static float cellClutterFractionMax;

// when analyzing cells, only cells for which the average dbz is 
// more than DBZCELL are considered in the rest of the analysis
static float cellDbzMin;

// minimum quality of the fit
static float chisqMin;

// threshold dbz value for excluding gates as clutter (static clutter 
// only)
static float clutterValueMin;

// maximum dbz used in calculation of profile dbzAvg  
static float dbzMax;

// minimum dbz used in calculation of cell dbzAvg
static float dbzThresMin;

// each weather cell identified by findWeatherCells() is grown by a distance 
// equal to 'fringeDist' using a region-growing approach
static float fringeDist;

// when determining whether there are enough vrad observations in 
// each direction, use NBINSGAP sectors
static int nBinsGap;

// when calculating the altitude-layer averaged dbz, there should 
// be at least NDBZMIN valid data points
static int nPointsIncludedMin;

// the minimum number of direct neighbors with dbz value above 
// dbzThresMin as used in findWeatherCells()
static int nNeighborsMin;

// there should be at least NOBSGAPMIN vrad observations in each 
// sector
static int nObsGapMin;

// vrad's texture is calculated based on the local neighborhood. The
// neighborhood size in the azimuth direction is equal to NTEXBINAZIM
static int nAzimNeighborhood;

// vrad's texture is calculated based on the local neighborhood. The
// neighborhood size in the range direction is equal to NTEXBINRANG
static int nRangNeighborhood;

// the minimum number of neighbors for the texture value to be 
// considered valid, as used in calcTexture()
static int nCountMin; 

// the refractive index of water
static float refracIndex;

// the bird radar cross section
static float birdRadarCrossSection;

// when analyzing cells, only cells for which the stddev of vrad
// (aka the texture) is less than STDEVCELL are considered in the
// rest of the analysis
static float cellStdDevMax;

// Minimum radial velocity standard deviation to be labeled as birds
static float stdDevMinBird;

// after fitting the vrad data, throw out any vrad observations that 
// are more that VDIFMAX away from the fitted value, since these are
// likely outliers 
static float absVDifMax;

// When analyzing cells, radial velocities lower than VRADMIN 
// are treated as clutter
static float vradMin;


// ------------------------------------------------------------- //
//               information about the 'points' array            //
// ------------------------------------------------------------- //

// The data needed for calculating bird densities are collected
// in one big array, 'points'. Although this array is one 
// variable, it is partitioned into 'nLayers' parts. The parts 
// are not equal in size, therefore we need to keep track of  
// where the data pertaining to a certain altitude bin can be
// written. The valid range of indexes into 'points' are stored 
// in arrays 'indexFrom' and 'indexTo'.

// the 'points' array has this many pseudo-columns
static int nColsPoints;

// the 'points' array has this many rows
static int nRowsPoints;

// the psuedo-column in 'points' that holds the azimuth angle
static int azimAngleCol;

// the psuedo-column in 'points' that holds the elevation angle
static int elevAngleCol;

// the psuedo-column in 'points' that holds the dbz value
static int dbzValueCol;

// the psuedo-column in 'points' that holds the vrad value
static int vradValueCol;

// the psuedo-column in 'points' that holds the cell value
static int cellValueCol;

// the psuedo-column in 'points' that holds the gate classification code
static int gateCodeCol;

// the 'points' array itself
static float* points;

// ------------------------------------------------------------- //
//             lists of indices into the 'points' array:         //
//          where each altitude layer's data starts and ends     //
// ------------------------------------------------------------- //

// for a given altitude layer in the profile, only part of the 'points'
// array is relevant. The 'indexFrom' and 'indexTo' arrays keep track 
// which rows in 'points' pertains to a given layer
static int* indexFrom;
static int* indexTo;

// nPointsWritten stores the number of points that was copied from one 
// of the scan elevations to the 'points' array; it should therefore
// never exceed indexTo[i]-indexFrom[i]
static int* nPointsWritten;




// ------------------------------------------------------------- //
//          information about the flagfields of 'gateCode'       //
// ------------------------------------------------------------- //

// the bit in 'gateCode' that says whether this gate is true in the static
// clutter map (which we don't have yet TODO)
static int flagPositionStaticClutter;

// the bit in 'gateCode' that says whether this gate is part of the 
// calculated cluttermap (without fringe)
static int flagPositionDynamicClutter;

// the bit in 'gateCode' that says whether this gate is part of the 
// fringe of the calculated cluttermap
static int flagPositionDynamicClutterFringe;

// the bit in 'gateCode' that says whether this gate has reflectivity data 
// but no corresponding radial velocity data
static int flagPositionVradMissing;

// the bit in 'gateCode' the psuedo-columnsays whether this gate's dbz value is too
// high to be due to birds, it must be caused by something else
static int flagPositionDbzTooHighForBirds;

// the bit in 'gateCode' the psuedo-columnsays whether this gate's radial velocity is
// close to zero. These gates are all discarded to exclude ground 
// clutter, which often has a radial velocity near zero.
static int flagPositionVradTooLow;

// the bit in 'gateCode' that says whether this gate passed the VDIFMAX test
static int flagPositionVDifMax;

// the bit in 'gateCode' that says whether the gate's azimuth angle was too low
static int flagPositionAzimTooLow;

// the bit in 'gateCode' that says whether the gate's azimuth angle was too high
static int flagPositionAzimTooHigh;




// ------------------------------------------------------------- //
//              information about the 'profile' array            //
// ------------------------------------------------------------- //

// the number of different types of profile we're making
static int nProfileTypes;

// how many rows there are in a profile
static int nRowsProfile;

// columns in profile contain [altmin,altmax,u,v,w,hSpeed,hDir,chi,hasGap,dbzAvg,nPointsCopied,reflectivity,birdDensity]
static int nColsProfile; 

// the profile array itself
static float* profile;

// the type of profile that was last calculated
static int iProfileTypeLast;




// ------------------------------------------------------------- //
//                       some other variables                    //
// ------------------------------------------------------------- //

// rCellMax is defined as rangeMax + 5000.0f  in order to avoid 
// edge effects when calculating the fringe
static float rCellMax;

// the number of dimensions to describe the location of an echo 
// (azimuth and elevAngle) as used in the 'pointsSelection' array 
// that is passed to svdfit
static int nDims;

// how many parameters are fitted by the svdfit procedure
static int nParsFitted;

// the factor that is used when converting from Z to eta
static float dbzFactor;

// whether the vol2bird module has been initialized
static int initializationSuccessful = FALSE;

// during calculation of iProfileType == 3, use this array to store the
// result of (chi < stdDevMinBird), such that it can be used later during 
// calculation of iProfileType == 1
static int* scatterersAreNotBirds;





static int analyzeCells(const unsigned char *dbzImage, const unsigned char *vradImage,
        const unsigned char *texImage, const unsigned char *clutterImage, int *cellImage,
        const SCANMETA *dbzMeta, const SCANMETA *vradMeta, const SCANMETA *texMeta, const SCANMETA *clutterMeta,
        const int nCells, const int useStaticClutterData) {

    // ----------------------------------------------------------------------------------- // 
    //  This function analyzes the cellImage array found by the 'findWeatherCells'         //
    //  procedure. Small cells are rejected and the cells are re-numbered according        //
    //  to size. The final number of cells in cellImage is returned as an integer.         // 
    // ----------------------------------------------------------------------------------- //

    CELLPROP *cellProp;
    int iCell;
    int iGlobal;
    int nGlobal;
    int iRang;
    int iAzim;
    int nCellsValid;
    int nAzim;
    int nRang;
    float validArea;
    float dbzValue;
    float texValue;
    float vradValue;
    float clutterValue;

    nCellsValid = nCells;
    nRang = dbzMeta->nRang;
    nAzim = dbzMeta->nAzim;
    nGlobal = nAzim*nRang;
    nCellsValid = 0;

    if (nCells == 0) {
        iGlobal = 0;
        for (iAzim = 0; iAzim < nAzim; iAzim++) {
            for (iRang = 0; iRang < nRang; iRang++) {
                cellImage[iGlobal] = -1;
                iGlobal++;
            }
        }
        return nCellsValid;
    }

    // Allocating and initializing memory for cell properties.
    cellProp = (CELLPROP *)malloc(nCells*sizeof(CELLPROP));
    if (!cellProp) {
        fprintf(stderr,"Requested memory could not be allocated!\n");
        return -10;
    }
    for (iCell = 0; iCell < nCells; iCell++) {
        cellProp[iCell].iRangOfMax = -1;
        cellProp[iCell].iAzimOfMax = -1;
        cellProp[iCell].nGates = 0;
        cellProp[iCell].nGatesClutter = 0;
        cellProp[iCell].dbzAvg = 0;
        cellProp[iCell].texAvg = 0;
        cellProp[iCell].dbzMax = dbzMeta->valueOffset;
        cellProp[iCell].index = iCell;
        cellProp[iCell].drop = FALSE;
        cellProp[iCell].cv = 0;
    }

    // Calculation of cell properties.
    for (iAzim = 0; iAzim < nAzim; iAzim++) {
        for (iRang = 0; iRang < nRang; iRang++) {

            iGlobal = iRang + iAzim * nRang;

            dbzValue = dbzMeta->valueScale * (float) dbzImage[iGlobal] + dbzMeta->valueOffset;
            vradValue = vradMeta->valueScale * (float) vradImage[iGlobal] + vradMeta->valueOffset;
            clutterValue = clutterMeta->valueScale * (float) clutterImage[iGlobal] + clutterMeta->valueOffset;
            texValue = texMeta->valueScale * (float) texImage[iGlobal] + texMeta->valueOffset;

            iCell = cellImage[iGlobal];

            if (iCell<0) {
                continue;
            }

            #ifdef FPRINTFON
            fprintf(stderr,"dbzValue = %f; vradValue = %f; clutterValue = %f; texValue = %f\n",dbzValue,vradValue,clutterValue,texValue);
            fprintf(stderr,"iGlobal = %d, iCell = %d\n",iGlobal,iCell);
            #endif

            cellProp[iCell].nGates += 1;

            // low radial velocities are treated as clutter, not included in calculation cell properties
            if (vradValue < vradMin){

                cellProp[iCell].nGatesClutter += 1;

                #ifdef FPRINTFON
                fprintf(stderr,"iGlobal = %d: vrad too low...treating as clutter\n",iGlobal);
                #endif

                continue;
            }

            // pixels in clutter map not included in calculation cell properties
            if (useStaticClutterData == TRUE){
                if (clutterValue > clutterValueMin){
                    cellProp[iCell].nGatesClutter += 1;
                    continue;
                }
            }



            if (dbzValue > cellProp[iCell].dbzMax) {

                #ifdef FPRINTFON
                fprintf(stderr,"%d: new dbzMax value of %f found for this cell (%d).\n",iGlobal,dbzValue,iCell);
                #endif

                cellProp[iCell].dbzMax = dbzValue;
                cellProp[iCell].iRangOfMax = iGlobal%nRang;
                cellProp[iCell].iAzimOfMax = iGlobal/nRang;
            }
            cellProp[iCell].dbzAvg += dbzValue;
            cellProp[iCell].texAvg += texValue;
        } // for (iRang = 0; iRang < nRang; iRang++)
    } // for (iAzim = 0; iAzim < nAzim; iAzim++)


    for (iCell = 0; iCell < nCells; iCell++) {
        validArea = cellProp[iCell].nGates - cellProp[iCell].nGatesClutter;
        if (validArea > 0){
            cellProp[iCell].dbzAvg /= validArea;
            cellProp[iCell].texAvg /= validArea;
            cellProp[iCell].cv = cellProp[iCell].texAvg / cellProp[iCell].dbzAvg;
        }
    }

    // determine which cells to drop from map based on low mean dBZ / high stdev /
    // small area / high percentage clutter
    for (iCell = 0; iCell < nCells; iCell++) {
        if (cellProp[iCell].nGates < nGatesCellMin ||
            (cellProp[iCell].dbzAvg < cellDbzMin &&
             cellProp[iCell].texAvg > cellStdDevMax &&
             ((float) cellProp[iCell].nGatesClutter / cellProp[iCell].nGates) < cellClutterFractionMax )) { // FIXME this condition still looks wrong if you ask me
            // Terms 2,3 and 4 are combined with && to be conservative in labeling stuff as
            // bird migration --see discussion of issue #37 on GitHub.
            cellProp[iCell].drop = TRUE;
            cellProp[iCell].nGates = 0;
        }
    }

    // sorting cell properties according to cell area. Drop small cells from map
    nCellsValid = updateMap(cellImage,nGlobal,cellProp,nCells);

    // printing of cell properties to stderr
    if (printCellProp == TRUE) {
        fprintf(stderr,"#Cell analysis for elevation %f:\n",dbzMeta->elev);
        fprintf(stderr,"#Minimum cell area in pixels   : %i\n",nGatesCellMin);
        fprintf(stderr,"#Threshold for mean dBZ cell   : %g dBZ\n",cellDbzMin);
        fprintf(stderr,"#Threshold for mean stdev cell : %g dBZ\n",cellStdDevMax);
        fprintf(stderr,"#Valid cells                   : %i/%i\n#\n",nCellsValid,nCells);
        fprintf(stderr,"cellProp: .index .nGates .nGatesClutter .dbzAvg .texAvg .cv   .dbzMax .iRangOfMax .iAzimOfMax .drop\n");
        for (iCell = 0; iCell < nCells; iCell++) {
            if (cellProp[iCell].nGates == 0) {
                continue;
            }
            fprintf(stderr,"cellProp: %6d %7d %14d %7.2f %7.2f %5.2f %7.2f %11d %11d %5c\n",
                    cellProp[iCell].index,
                    cellProp[iCell].nGates,
                    cellProp[iCell].nGatesClutter,
                    cellProp[iCell].dbzAvg,
                    cellProp[iCell].texAvg,
                    cellProp[iCell].cv,
                    cellProp[iCell].dbzMax,
                    cellProp[iCell].iRangOfMax,
                    cellProp[iCell].iAzimOfMax,
                    cellProp[iCell].drop == TRUE ? 'T' : 'F');
        }
    } // endif (printCellProp == TRUE)

    free(cellProp);

    return nCellsValid;
} // analyzeCells





static float calcDist(const int iRang1, const int iAzim1, const int iRang2, const 
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





static void calcTexture(unsigned char* texImage, const unsigned char* vradImage,
        const unsigned char* dbzImage, const SCANMETA* texMeta, const SCANMETA* vradMeta,
        const SCANMETA* dbzMeta) {


    // --------------------------------------------------------------------------------------- //
    // This function computes a texture parameter based on a block of (nRangNeighborhood x     //
    // nAzimNeighborhood) pixels. The texture parameter equals the local standard deviation    //
    // in the radial velocity field                                                            //
    // --------------------------------------------------------------------------------------- //


    int iRang;
    int iAzim;
    int nRang;
    int nAzim;
    int iNeighborhood;
    int nNeighborhood;
    int count;
    unsigned char vradMissingValue;
    unsigned char dbzMissingValue;
    unsigned char texMissingValue;
    double vmoment1;
    double vmoment2;
    double dbz;
    double tex;
    int iGlobal;
    int iLocal;
    float texOffset;
    float texScale;
    float dbzOffset;
    float dbzScale;
    float vradOffset;
    float vradScale;
    float vRadDiff;

    nRang = vradMeta->nRang;
    nAzim = vradMeta->nAzim;

    dbzOffset = dbzMeta->valueOffset;
    dbzScale = dbzMeta->valueScale;
    dbzMissingValue = dbzMeta->missing;

    vradOffset = vradMeta->valueOffset;
    vradScale = vradMeta->valueScale;
    vradMissingValue = vradMeta->missing;

    texOffset = texMeta->valueOffset;
    texScale = texMeta->valueScale;
    texMissingValue = texMeta->missing;

    nNeighborhood = (int)(nRangNeighborhood * nAzimNeighborhood);

    for (iAzim = 0; iAzim < nAzim; iAzim++) {
        for (iRang = 0; iRang < nRang; iRang++) {

            iGlobal = iRang + iAzim * nRang;

            // count number of direct neighbors above threshold
            count = 0;
            vmoment1 = 0;
            vmoment2 = 0;

            dbz = 0;

            for (iNeighborhood = 0; iNeighborhood < nNeighborhood; iNeighborhood++) {

                iLocal = findNearbyGateIndex(nAzim,nRang,iGlobal,nAzimNeighborhood,nRangNeighborhood,iNeighborhood);

                #ifdef FPRINTFON
                fprintf(stderr, "iLocal = %d; ",iLocal);
                #endif

                if (iLocal < 0) {
                    // iLocal less than zero are error codes
                    continue;
                }

                if (vradImage[iLocal] == vradMissingValue || dbzImage[iLocal] == dbzMissingValue) {
                    continue;
                }

                vRadDiff = vradOffset + vradScale * (float) (vradImage[iGlobal] - vradImage[iLocal]);
                vmoment1 += vRadDiff;
                vmoment2 += SQUARE(vRadDiff);

                dbz += dbzOffset + dbzScale * (float) dbzImage[iLocal];

                count++;

            }

            vmoment1 /= count;
            vmoment2 /= count;
            dbz /= count;

            // when not enough neighbors, continue
            if (count < nCountMin) {
                texImage[iGlobal] = texMissingValue;
            }
            else {

                tex = sqrt(XABS(vmoment2-SQUARE(vmoment1)));

                double tmpTex = ROUND((tex - texOffset) / texScale);
                if (0 <= tmpTex && tmpTex <= 255) {
                    texImage[iGlobal] = (unsigned char) tmpTex;
                }
                else {
                    fprintf(stderr, "Error casting texture value of %f to unsigned char type at texImage[%d]. Aborting.\n",tmpTex,iGlobal);
                    return;
                }
                

                #ifdef FPRINTFON
                fprintf(stderr,
                        "\n(C) count = %d; nCountMin = %d; vmoment1 = %f; vmoment2 = %f; tex = %f; texBody[%d] = %d\n",
                        count, nCountMin, vmoment1, vmoment2, tex,
                        iGlobal, texImage[iGlobal]);
                #endif

            } //else
        } //for
    } //for
} // calcTexture





static void classifyGatesSimple(void) {
    
    int iPoint;
    
    for (iPoint = 0; iPoint < nRowsPoints; iPoint++) {
    
        const float azimValue = points[iPoint * nColsPoints + azimAngleCol];
        const float dbzValue = points[iPoint * nColsPoints + dbzValueCol];        
        const float vradValue = points[iPoint * nColsPoints + vradValueCol];
        const int cellValue = (int) points[iPoint * nColsPoints + cellValueCol];

        unsigned int gateCode = 0;
        
        if (FALSE) {
            // this gate is true in the static clutter map (which we don't have yet TODO)
            gateCode |= 1<<flagPositionStaticClutter;
        }

        if (cellValue == 1) {
            // this gate is part of the cluttermap (without fringe)
            gateCode |= 1<<flagPositionDynamicClutter;
        }

        if (cellValue == 2) {
            // this gate is part of the fringe of the cluttermap
            gateCode |= 1<<flagPositionDynamicClutterFringe;
        }

        if (FALSE) {
            // this gate has reflectivity data but no corresponding radial velocity data
            // TODO no condition for this yet
            gateCode |= 1<<flagPositionVradMissing;
        }

        if (dbzValue > dbzMax) {
            // this gate's dbz value is too high to be due to birds, it must be 
            // caused by something else
            gateCode |= 1<<flagPositionDbzTooHighForBirds;
        }

        if (vradValue < vradMin) {
            // this gate's radial velocity is too low to be due to actual scatterers; likely just noise
            gateCode |= 1<<flagPositionVradTooLow;
        }

        if (FALSE) {
            // this flag is set later on by updateFlagFieldsInPointsArray()
            // flagPositionVDifMax
        }

        if (azimValue < azimMin) {
            // the user can specify to exclude gates based on their azimuth;
            // this clause is for gates that have too low azimuth
            gateCode |= 1<<flagPositionAzimTooLow;
        }
        
        if (azimValue > azimMax) {
            // the user can specify to exclude gates based on their azimuth;
            // this clause is for gates that have too high azimuth
            gateCode |= 1<<flagPositionAzimTooHigh;
        }

        points[iPoint * nColsPoints + gateCodeCol] = (float) gateCode;
        
    }

    return;
    
    
};





static int constructorInt(SCANMETA* meta, int* image, PolarScan_t* scan, const int nGlobal, const int initValue) {

    int iGlobal;
    
    if (image == NULL) {
        fprintf(stderr,"Error allocating image");
        return -2;
    }
    for (iGlobal = 0; iGlobal < nGlobal; iGlobal++) {
        image[iGlobal] = initValue;
    }
    meta->heig = (float) PolarScan_getHeight(scan);
    meta->elev = (float) (360 * PolarScan_getElangle(scan) / 2 / PI );
    meta->nRang = (int) PolarScan_getNbins(scan);
    meta->nAzim = (int) PolarScan_getNrays(scan);
    meta->rangeScale = (float) PolarScan_getRscale(scan);
    meta->azimScale = 360.0f/meta->nAzim;   // FIXME not sure if this always works
    meta->valueOffset = 0.0f;
    meta->valueScale = 1.0f;
    meta->missing = (unsigned char) 255;  // FIXME this does not work as intended for type int
    
    return 0;

}




static int constructorUChar(SCANMETA* meta, unsigned char* image, PolarScan_t* scan, const int nGlobal, const unsigned char initValue) {

    int iGlobal;
    
    if (image == NULL) {
        fprintf(stderr,"Error allocating image");
        return -2;
    }
    for (iGlobal = 0; iGlobal < nGlobal; iGlobal++) {
        image[iGlobal] = initValue;
    }
    meta->heig = (float) PolarScan_getHeight(scan);
    meta->elev = (float) (360 * PolarScan_getElangle(scan) / 2 / PI );
    meta->nRang = (int) PolarScan_getNbins(scan);
    meta->nAzim = (int) PolarScan_getNrays(scan);
    meta->rangeScale = (float) PolarScan_getRscale(scan);
    meta->azimScale = 360.0f/meta->nAzim;   // FIXME not sure if this always works
    meta->valueOffset = 0.0f;
    meta->valueScale = 1.0f;
    meta->missing = (unsigned char) 255;
    
    return 0;

}



static void constructPointsArray(PolarVolume_t* volume) {
    
        // iterate over the scans in 'volume'
        int iScan;
        int nScans;
        
        // determine how many scan elevations the volume object contains
        nScans = PolarVolume_getNumberOfScans(volume);


        for (iScan = 0; iScan < nScans; iScan++) {
            
            // initialize the scan object
            PolarScan_t* scan = NULL;
        
            // extract the scan object from the volume object
            scan = PolarVolume_getScan(volume, iScan);
            
            // determine the number of array elements in the polar scan
            int nGlobal = (int) PolarScan_getNbins(scan) * PolarScan_getNrays(scan);
            
            // pre-allocate the dbz variables
            unsigned char* dbzImage = malloc(sizeof(unsigned char) * nGlobal);
            SCANMETA dbzMeta;
            constructorUChar(&dbzMeta, &dbzImage[0], scan, nGlobal, 0);
            
            // pre-allocate the vrad variables
            unsigned char* vradImage = malloc(sizeof(unsigned char) * nGlobal);
            SCANMETA vradMeta;
            constructorUChar(&vradMeta, &vradImage[0], scan, nGlobal, 0);
            
            // pre-allocate the tex variables
            unsigned char* texImage = malloc(sizeof(unsigned char) * nGlobal);
            SCANMETA texMeta;
            constructorUChar(&texMeta, &texImage[0], scan, nGlobal, 0);
            
            // pre-allocate the clutter variables
            unsigned char* clutterImage = malloc(sizeof(unsigned char) * nGlobal);
            SCANMETA clutterMeta;
            constructorUChar(&clutterMeta, &clutterImage[0], scan, nGlobal, 0);
            
            // pre-allocate the cell variables
            int* cellImage = malloc(sizeof(int) * nGlobal);
            SCANMETA cellMeta;
            constructorInt(&cellMeta, &cellImage[0], scan, nGlobal, -1);

            // populate the dbzMeta and dbzImage variables with data from 
            // the Rave scan object:
            int rcDbz = mapDataFromRave(scan, &dbzMeta, &dbzImage[0],"DBZH");
            if (rcDbz != 0) {
                fprintf(stderr, "Something went wrong while mapping DBZH data from RAVE to LIBVOL2BIRD.\n");
            }

            // populate the vradMeta and vradImage variables with data from  
            // the Rave scan object:
            int rcVrad = mapDataFromRave(scan, &vradMeta, &vradImage[0],"VRAD");
            if (rcVrad != 0) {
                fprintf(stderr, "Something went wrong while mapping VRAD data from RAVE to LIBVOL2BIRD.\n");
            }


            // ------------------------------------------------------------- //
            //                      calculate vrad texture                   //
            // ------------------------------------------------------------- //

            calcTexture(&texImage[0], &vradImage[0], &dbzImage[0], 
                        &texMeta, &vradMeta, &dbzMeta);


            // ------------------------------------------------------------- //
            //        find (weather) cells in the reflectivity image         //
            // ------------------------------------------------------------- //
            
            int nCells = findWeatherCells(&dbzImage[0], &cellImage[0], 
                                   &dbzMeta);
            
            if (printCellProp == TRUE) {
                fprintf(stderr,"(%d/%d): found %d cells.\n",iScan, nScans, nCells);
            }
            
            // ------------------------------------------------------------- //
            //                      analyze cells                            //
            // ------------------------------------------------------------- //

            int nCellsValid = analyzeCells(&dbzImage[0], &vradImage[0], &texImage[0], 
                &clutterImage[0], &cellImage[0], &dbzMeta, &vradMeta, &texMeta, 
                &clutterMeta, nCells, 
                useStaticClutterData);

            // ------------------------------------------------------------- //
            //                     calculate fringe                          //
            // ------------------------------------------------------------- //

            fringeCells(&cellImage[0], cellMeta.nRang, cellMeta.nAzim, 
                cellMeta.azimScale, cellMeta.rangeScale);
                

            // ------------------------------------------------------------- //
            //            print selected outputs to stderr                   //
            // ------------------------------------------------------------- //

            if (printDbz == TRUE) {
                printMeta(&dbzMeta,"dbzMeta");
                printImageUChar(&dbzMeta,&dbzImage[0]);
            }
            if (printVrad == TRUE) {
                printMeta(&vradMeta,"vradMeta");
                printImageUChar(&vradMeta,&vradImage[0]);
            }
            if (printTex == TRUE) {
                printMeta(&texMeta,"texMeta");
                printImageUChar(&texMeta,&texImage[0]);
            }
            if (printCell == TRUE) {
                printMeta(&cellMeta,"cellMeta");
                printImageInt(&cellMeta,&cellImage[0]);
            }
            if (printClut == TRUE) { 
                printMeta(&clutterMeta,"clutterMeta");
                printImageUChar(&clutterMeta,&clutterImage[0]);
            }
                        
            // ------------------------------------------------------------- //
            //    fill in the appropriate elements in the points array       //
            // ------------------------------------------------------------- //

            int iLayer;
            
            for (iLayer = 0; iLayer < nLayers; iLayer++) {
                
                float altitudeMin = iLayer * layerThickness;
                float altitudeMax = (iLayer + 1) * layerThickness;
                int iRowPoints = indexFrom[iLayer] + nPointsWritten[iLayer];

                int n = getListOfSelectedGates(&vradMeta, &vradImage[0],
                    &dbzMeta, &dbzImage[0], 
                    &cellImage[0], 
                    altitudeMin, altitudeMax, 
                    &points[0], iRowPoints);
                    
                nPointsWritten[iLayer] += n;

                if (indexFrom[iLayer] + nPointsWritten[iLayer] > indexTo[iLayer]) {
                    fprintf(stderr, "Problem occurred: writing over existing data\n");
                    return;
                }


            } // endfor (iLayer = 0; iLayer < nLayers; iLayer++)

            // ------------------------------------------------------------- //
            //                         clean up                              //
            // ------------------------------------------------------------- //


            // free previously malloc'ed arrays
            free((void*) dbzImage);
            free((void*) vradImage);
            free((void*) texImage);
            free((void*) cellImage);
            free((void*) clutterImage);
            
            RAVE_OBJECT_RELEASE(scan);
            
        } // endfor (iScan = 0; iScan < nScans; iScan++)

}





static int detNumberOfGates(const int iLayer, 
                     const float rangeScale, const float elevAngle,
                     const int nRang, const int nAzim,
                     const float radarHeight) {

    // Determine the number of gates that are within the limits set
    // by (rangeMin,rangeMax) as well as by (iLayer*layerThickness,
    // (iLayer+1)*layerThickness).

    int nGates;
    int iRang;

    float layerHeight;
    float range;
    float beamHeight;


    nGates = 0;
    layerHeight = (iLayer + 0.5) * layerThickness;
    for (iRang = 0; iRang < nRang; iRang++) {
        range = (iRang + 0.5) * rangeScale;
        if (range < rangeMin || range > rangeMax) {
            // the gate is too close to the radar, or too far away
            continue;
        }
        beamHeight = range * sin(elevAngle * DEG2RAD) + radarHeight;
        if (fabs(layerHeight - beamHeight) > 0.5*layerThickness) {
            // the gate is not close enough to the altitude layer of interest
            continue;
        }

        #ifdef FPRINTFON
        fprintf(stderr, "iRang = %d; range = %f; beamHeight = %f\n",iRang,range,beamHeight);
        #endif

        nGates += nAzim;

    } // for iRang

    return nGates;

} // detNumberOfGates()





static int detSvdfitArraySize(PolarVolume_t* volume) {
    
    int iScan;
    int nScans = PolarVolume_getNumberOfScans(volume);

    int iLayer;
    int nRowsPoints = 0;
    
    int* nGates = malloc(sizeof(int) * nLayers);
    int* nGatesAcc = malloc(sizeof(int) * nLayers);
    
    for (iLayer = 0; iLayer < nLayers; iLayer++) {
        nGates[iLayer] = 0;
        nGatesAcc[iLayer] = 0;
    }

    for (iScan = 0; iScan < nScans; iScan++) {
        for (iLayer = 0; iLayer < nLayers; iLayer++) {

            PolarScan_t* scan = PolarVolume_getScan(volume, iScan);
            
            int nRang = (int) PolarScan_getNbins(scan);
            int nAzim = (int) PolarScan_getNrays(scan);
            float elevAngle = (float) (360 * PolarScan_getElangle(scan) / 2 / PI );
            float rangeScale = (float) PolarScan_getRscale(scan);
            float radarHeight = (float) PolarScan_getHeight(scan);

            nGates[iLayer] += detNumberOfGates(iLayer, 
                rangeScale, elevAngle, nRang, 
                nAzim, radarHeight);
                
            RAVE_OBJECT_RELEASE(scan);

        }
    }

    for (iLayer = 0; iLayer < nLayers; iLayer++) {
        if (iLayer == 0) {
            nGatesAcc[iLayer] = nGates[iLayer];
        }
        else {
            nGatesAcc[iLayer] = nGatesAcc[iLayer-1] + nGates[iLayer];
        }
    }
    nRowsPoints = nGatesAcc[nLayers-1];

    for (iLayer = 0; iLayer < nLayers; iLayer++) {
        if (iLayer == 0) {
            indexFrom[iLayer] = 0;
        }
        else {
            indexFrom[iLayer] = nGatesAcc[iLayer-1];
        }
    }
    for (iLayer = 0; iLayer < nLayers; iLayer++) {
        indexTo[iLayer] = nGatesAcc[iLayer];
    }

    free((void*) nGates);
    free((void*) nGatesAcc);

    return nRowsPoints;
    
}  // detSvdfitArraySize()






static int findWeatherCells(const unsigned char *dbzImage, int *cellImage,
        const SCANMETA *dbzMeta) {

    //  ----------------------------------------------------------------------------- //
    //  This function detects the cells in 'dbzImage' using a threshold value of      //
    //  'dbzThresMin' and a non-recursive algorithm which looks for neighboring       //
    //  pixels above that threshold. On return the marked cells are contained by      //
    //  'cellImage'. The number of detected cells/highest index value is returned.    //
    //  ----------------------------------------------------------------------------- //


    int iCellIdentifier;
    int nCells;
    int iRang;
    int nRang;
    int iAzim;
    int nAzim;
    int iNeighborhood;
    int nNeighborhood;
    int count;
    int cellImageInitialValue;

    unsigned char dbzThres;

    int iGlobal;
    int iGlobalOther;
    int nGlobal;
    int iLocal;

    unsigned char dbzMissing;
    int dbznAzim;
    int dbznRang;
    float dbzValueOffset;
    float dbzValueScale;
    float dbzRangeScale;


    // These next two local variables 'nAzimNeighborhood' and 
    // 'nRangNeighborhood' shadow two static global variables of the 
    // same name, but that is in fact what you want. This is because the
    // value of 'nAzimNeighborHood' (function scope) is not equal to 
    // 'nAzimNeighborHood' (file scope) and because 'nRangNeighborHood' 
    // (function scope) is not equal to 'nRangNeighborHood' (file 
    // scope). The variable names are the same though, because they play
    // the same role in 'findNearbyGateIndex()'.
    int nAzimNeighborhood;
    int nRangNeighborhood;
    
    int nHalfNeighborhood;

    int cellIdentifierGlobal;
    int cellIdentifierGlobalOther;
    int iGlobalInner;


    #ifdef FPRINTFON
    int dbg = 0;
    #endif

    if (dbzImage==NULL) {
        fprintf(stderr,"Input argument dbzImage is NULL.\n");
        return -1;
    }

    dbzMissing = dbzMeta->missing;
    dbznAzim = dbzMeta->nAzim;
    dbznRang = dbzMeta->nRang;
    dbzValueOffset = dbzMeta->valueOffset;
    dbzValueScale = dbzMeta->valueScale;
    dbzRangeScale = dbzMeta->rangeScale;

    nAzim = dbznAzim;
    nRang = dbznRang;

    nGlobal = nAzim * nRang;
    
    // We use a neighborhood of 3x3, because in this function we are
    // interested in a cell's direct neighbors. See also comment at
    // variable definitions of 'nAzimNeighborhood' and 'nRangNeighborhood'
    // above. 
    nAzimNeighborhood = 3;
    nRangNeighborhood = 3;
    
    nNeighborhood = nAzimNeighborhood * nRangNeighborhood;
    nHalfNeighborhood = (nNeighborhood - 1)/2;


    if (dbzImage != NULL) {
        dbzThres = (unsigned char) ROUND((dbzThresMin - dbzValueOffset) / dbzValueScale);
    }

    cellImageInitialValue = -1;
    for (iGlobal = 0; iGlobal < nGlobal; iGlobal++) {
        cellImage[iGlobal] = cellImageInitialValue;
    }

    // If threshold value is equal to missing value, return.
    if (dbzThres == dbzMissing) {
        return -1;
    }

    // ----------------------------------------------------------------------- //
    // Labeling of groups of connected pixels using horizontal, vertical, and  //
    // diagonal connections. The algorithm is described in 'Digital Image      //
    // Processing' by Gonzales and Woods published by Addison-Wesley.          //
    // ----------------------------------------------------------------------- //

    iCellIdentifier = 0;

    for (iAzim = 0; iAzim < nAzim; iAzim++) {
        for (iRang = 0; iRang < nRang; iRang++) {

            iGlobal = iRang + iAzim * nRang;

            if ((float)(iRang + 1) * dbzRangeScale > rCellMax) {
                continue;
            }
            else {
                #ifdef FPRINTFON
                fprintf(stderr, "iGlobal = %d\niRang + 1 = %d\n"
                "dbzRangeScale = %f\n"
                "rCellMax = %f\n"
                "(iRang + 1) * dbzRangeScale = %f\n"
                "((iRang + 1) * dbzRangeScale > rCellMax) = %d\n"
                "dbg=%d\n",iGlobal,iRang + 1,dbzRangeScale,rCellMax,
                (iRang + 1) * dbzRangeScale,
                ((iRang + 1) * dbzRangeScale > rCellMax),dbg);
                
                dbg++;
                #endif
            }

            #ifdef FPRINTFON
            fprintf(stderr,"iGlobal = %d\n",iGlobal);
            #endif

            if (dbzImage[iGlobal] == (unsigned char) dbzMissing) {

                #ifdef FPRINTFON
                fprintf(stderr,"dbzImage[%d] == dbzMissing\n",iGlobal);
                #endif

                continue;
            }

            if (dbzImage[iGlobal] < (unsigned char) dbzThres) {
                continue;
            }

            // count number of direct neighbors above threshold
            count = 0;
            for (iNeighborhood = 0; iNeighborhood < nNeighborhood; iNeighborhood++) {

                iLocal = findNearbyGateIndex(nAzim,nRang,iGlobal,nAzimNeighborhood,nRangNeighborhood,iNeighborhood);

                if (iLocal < 0) {
                    // iLocal below zero are error codes
                    continue;
                }
                if (dbzImage[iLocal] > dbzThres) {
                    count++;
                }

            }
            // when not enough qualified neighbors, continue
            if (count - 1 < nNeighborsMin) {
                continue;
            }

            #ifdef FPRINTFON
            fprintf(stderr,"iGlobal = %d, count = %d\n",iGlobal,count);
            #endif


            // Looking for horizontal, vertical, forward diagonal, and backward diagonal connections.
            for (iNeighborhood = 0; iNeighborhood < nHalfNeighborhood; iNeighborhood++) {

                iLocal = findNearbyGateIndex(nAzim,nRang,iGlobal,nAzimNeighborhood,nRangNeighborhood,iNeighborhood);

                if (iLocal < 0) {
                    // iLocal less than zero are error codes
                    continue;
                }

                // no connection found, go to next pixel within neighborhood
                if (cellImage[iLocal] == cellImageInitialValue) {
                    continue;
                }

                // if pixel still unassigned, assign same iCellIdentifier as connection
                if (cellImage[iGlobal] == cellImageInitialValue) {
                    cellImage[iGlobal] = cellImage[iLocal];
                }
                else {
                    // if connection found but pixel is already assigned a different iCellIdentifier:
                    if (cellImage[iGlobal] != cellImage[iLocal]) {
                        // merging cells detected: replace all other occurences by value of connection:
                        for (iGlobalOther = 0; iGlobalOther < nGlobal; iGlobalOther++) {
                            if (cellImage[iGlobalOther] == cellImage[iGlobal]) {
                                cellImage[iGlobalOther] = cellImage[iLocal];
                            }
                            // note: not all iCellIdentifier need to be used eventually
                        }
                    }
                }
            }

            // When no connections are found, give a new number.
            if (cellImage[iGlobal] == cellImageInitialValue) {

                #ifdef FPRINTFON
                fprintf(stderr, "new cell found...assigning number %d\n",iCellIdentifier);
                #endif

                cellImage[iGlobal] = iCellIdentifier;
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
        iGlobalOther = findNearbyGateIndex(nAzim,nRang,iGlobal,3,3,1);

        #ifdef FPRINTFON
        fprintf(stderr,"iGlobal = %d, iGlobalOther = %d\n",iGlobal,iGlobalOther);
        #endif

        cellIdentifierGlobal = cellImage[iGlobal];
        cellIdentifierGlobalOther = cellImage[iGlobalOther];
        if (cellIdentifierGlobal != cellImageInitialValue && cellIdentifierGlobalOther != cellImageInitialValue ) {
            // adjacent gates, both part of a cell -> assign them the same identifier, i.e. assign
            // all elements of cellImage that are equal to cellImage[iGlobalOther] the value of
            // cellImage[iGlobal]

            for (iGlobalInner = 0; iGlobalInner < nGlobal; iGlobalInner++) {
                if (cellImage[iGlobalInner] == cellIdentifierGlobalOther) {
                    cellImage[iGlobalInner] = cellIdentifierGlobal;
                }
            }
        }
    }

    // Returning number of detected cells (including fringe/clutter)
    nCells = iCellIdentifier;

    return nCells;
} // findWeatherCells






static int findNearbyGateIndex(const int nAzimParent, const int nRangParent, const int iParent,
                        const int nAzimChild,  const int nRangChild,  const int iChild) {



    if (nRangChild%2 != 1) {

        #ifdef FPRINTFON
        fprintf(stderr, "nRangChild must be an odd integer number.\n");
        #endif

        return -1;
    }

    if (nAzimChild%2 != 1) {

        #ifdef FPRINTFON
        fprintf(stderr, "nAzimChild must be an odd integer number.\n");
        #endif

        return -2;
    }

    if (iChild > nAzimChild * nRangChild - 1) {

        #ifdef FPRINTFON
        fprintf(stderr, "iChild is outside the child window.\n");
        #endif

        return -3;
    }


    const int iAzimParent = iParent / nRangParent;
    const int iRangParent = iParent % nRangParent;

    const int iAzimChild = iChild / nRangChild;
    const int iRangChild = iChild % nRangChild;

    // the azimuth dimension is wrapped (polar plot); iAzim = 0 is adjacent to iAzim = nAzim-1:
    const int iAzimReturn = (iAzimParent - nAzimChild/2 + iAzimChild + nAzimParent) % nAzimParent;
    const int iRangReturn = iRangParent - nRangChild/2 + iRangChild;


    if (iRangReturn > nRangParent - 1) {

        #ifdef FPRINTFON
        fprintf(stderr, "iChild is outside the parent array on the right-hand side.\n");
        #endif

        return -4;
    }
    if (iRangReturn < 0) {

        #ifdef FPRINTFON
        fprintf(stderr, "iChild is outside the parent array on the left-hand side.\n");
        #endif

        return -5;
    }

    return iAzimReturn * nRangParent + iRangReturn;

} // findNearbyGateIndex








static void fringeCells(int *cellImage, int nRang, int nAzim, float aScale, float rScale) {

    // -------------------------------------------------------------------------- //
    // This function enlarges cells in cellImage by an additional fringe.         //
    // First a block around each pixel is searched for pixels within a distance   //
    // equal to 'fringeDist'.                                                     //
    // -------------------------------------------------------------------------- //

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


    rBlock = ROUND(fringeDist / rScale);
    for (iAzim = 0; iAzim < nAzim; iAzim++) {
        for (iRang = 0; iRang < nRang; iRang++) {

            iGlobal = iRang + iAzim * nRang;

            if (cellImage[iGlobal] <= 1) {
                continue; // with the next iGlobal; already fringe or not in cellImage
            }

            // determine whether current pixel is a pixel on the edge of a cell
            isEdge = FALSE;
            for (iNeighborhood = 0; iNeighborhood < 9 && isEdge == FALSE; iNeighborhood++) {

                iLocal = findNearbyGateIndex(nAzim,nRang,iGlobal,3,3,iNeighborhood);

                if (iLocal < 0) {
                    // iLocal less than zero are error codes
                    continue; // with the next iGlobal
                }

                if (cellImage[iLocal] <= 1) {
                    isEdge = TRUE;
                }

            }
            
            if (isEdge == FALSE) {
                continue; // with the next iGlobal
            }

            actualRange = (iRang+0.5) * rScale;
            circumferenceAtActualRange = 2 * PI * actualRange;
            aBlock = (fringeDist / circumferenceAtActualRange) * nAzim;


            #ifdef FPRINTFON
            fprintf(stderr, "actualRange = %f\n",actualRange);
            fprintf(stderr, "circumferenceAtActualRange = %f\n",circumferenceAtActualRange);
            fprintf(stderr, "fringeDist / circumferenceAtActualRange = %f\n",fringeDist / circumferenceAtActualRange);
            fprintf(stderr, "aBlock = %d\n", aBlock);
            fprintf(stderr, "rBlock = %d\n", rBlock);
            #endif

            nAzimChild = 2 * aBlock + 1;
            nRangChild = 2 * rBlock + 1;
            nNeighborhood = nAzimChild * nRangChild;

            for (iNeighborhood = 0; iNeighborhood < nNeighborhood; iNeighborhood++) {

                iLocal = findNearbyGateIndex(nAzim,nRang,iGlobal,nAzimChild,nRangChild,iNeighborhood);

                if (iLocal < 0) {
                    // iLocal less than zero are error codes
                    continue;
                }

                iAzimLocal = iLocal / nRang;
                iRangLocal = iLocal % nRang;

                // if not within range or already in cellImage or already a fringe, do nothing
                theDist = calcDist(iRang, iAzim, iRangLocal, iAzimLocal, rScale, aScale);
                if (theDist > fringeDist || cellImage[iLocal] >= 1) {
                    continue; // with the next iGlobal
                }
                // include pixel (iRangLocal,iAzimLocal) in fringe
                cellImage[iLocal] = 1;
                
            } // (iNeighborhood = 0; iNeighborhood < nNeighborhood; iNeighborhood++)
        } // (iRang = 0; iRang < nRang; iRang++)
    } // (iAzim = 0; iAzim < nAzim; iAzim++)

    return;

} // fringeCells




static int getListOfSelectedGates(const SCANMETA* vradMeta, const unsigned char *vradImage, 
                           const SCANMETA* dbzMeta, const unsigned char *dbzImage,
                           const int *cellImage,
                           const float altitudeMin, const float altitudeMax,
                           float* points, int iRowPoints) {

    // ------------------------------------------------------------------- //
    // Write combinations of an azimuth angle, an elevation angle, an      // 
    // observed vrad value, an observed dbz value, and a cell identifier   //
    // value into an external larger list.                                 //
    // ------------------------------------------------------------------- //

    int iAzim;
    int iRang;
    int iGlobal;
    int nRang;
    int nAzim;
    int nPointsWritten;
    const int nColsPoints = 6;

    float gateHeight;
    float gateRange;
    float gateAzim;
    float rangeScale;
    float azimuthScale;
    float elevAngle;
    float radarHeight;
    float vradValueOffset;
    float vradValueScale;
    float dbzValueOffset;
    float dbzValueScale;
    float vradValue;
    float dbzValue;
    
    nPointsWritten = 0;

    nRang = vradMeta->nRang;
    nAzim = vradMeta->nAzim;
    rangeScale = vradMeta->rangeScale;
    azimuthScale = vradMeta->azimScale;
    elevAngle = vradMeta->elev;
    radarHeight = vradMeta->heig;
    vradValueOffset = vradMeta->valueOffset;
    vradValueScale = vradMeta->valueScale;

    dbzValueOffset = dbzMeta->valueOffset;
    dbzValueScale = dbzMeta->valueScale;


    for (iRang = 0; iRang < nRang; iRang++) {

        // so gateRange represents a distance along the view direction (not necessarily horizontal)
        gateRange = ((float) iRang + 0.5f) * rangeScale;

        // note that "sin(elevAngle*DEG2RAD)" is equivalent to = "cos((90 - elevAngle)*DEG2RAD)":
        gateHeight = gateRange * (float) sin(elevAngle*DEG2RAD) + radarHeight;

        if (gateRange < rangeMin || gateRange > rangeMax) {
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

            iGlobal = iRang + iAzim * nRang;
            gateAzim = ((float) iAzim + 0.5f) * azimuthScale;
            vradValue = vradValueScale * (float) vradImage[iGlobal] + vradValueOffset;
            dbzValue = dbzValueScale * (float) dbzImage[iGlobal] + dbzValueOffset;

            // store the location as an azimuth angle, elevation angle combination
            points[iRowPoints * nColsPoints + 0] = gateAzim;
            points[iRowPoints * nColsPoints + 1] = elevAngle;

            // also store the dbz value --useful when estimating the bird density
            points[iRowPoints * nColsPoints + 2] = dbzValue;
            
            // store the corresponding observed vrad value
            points[iRowPoints * nColsPoints + 3] = vradValue;

            // store the corresponding cellImage value
            points[iRowPoints * nColsPoints + 4] = (float) cellImage[iGlobal];

            // set the gateCode to zero for now
            points[iRowPoints * nColsPoints + 5] = (float) 0;

            // raise the row counter by 1
            iRowPoints += 1;
            
            // raise number of points written by 1
            nPointsWritten += 1;


        }  //for iAzim
    } //for iRang

    return nPointsWritten;


} // getListOfSelectedGates




static int hasAzimuthGap(const float* points, const int nPoints) {

    int hasGap;
    int nObs[nBinsGap];
    int iPoint;
    int iBinGap;
    int iBinGapNext;
    float azimuth;

    hasGap = FALSE;

    // Initialize histogram
    for (iBinGap = 0; iBinGap < nBinsGap; iBinGap++) {
        nObs[iBinGap] = 0;
    }

    // Collect histogram data
    for (iPoint = 0; iPoint < nPoints; iPoint++) {
        azimuth = points[iPoint*nDims];
        iBinGap = ((int) floor((azimuth / 360.0) * nBinsGap)) % nBinsGap;
        nObs[iBinGap]++;
    }

    // Detect adjacent bins in which the number of azimuth observations 
    // is less than the minimum required number
    for (iBinGap = 0; iBinGap < nBinsGap; iBinGap++) {
        
        iBinGapNext = (iBinGap + 1) % nBinsGap;
        
        if (nObs[iBinGap] < nObsGapMin && nObs[iBinGapNext] < nObsGapMin) {
            hasGap = TRUE;
        }
    }

    return hasGap;
    
} // hasAzimuthGap






static int includeGate(const int iProfileType, const unsigned int gateCode) {
    
    int doInclude = TRUE;
    
    if (gateCode & 1<<flagPositionStaticClutter) {
        
        // i.e. flag 0 in gateCode is true
        // this gate is true in the static clutter map (which we don't have yet TODO)
        
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
                fprintf(stderr, "Something went wrong; behavior not implemented for given iProfileType.\n");
        }
    }
    
    if (gateCode & 1<<flagPositionDynamicClutter) {
        
        // i.e. flag 1 in gateCode is true
        // this gate is part of the cluttermap (without fringe)
        
        switch (iProfileType) {
            case 1 : 
                doInclude = FALSE;
                break;
            case 2 : 
                doInclude = TRUE;
                break;
            case 3 : 
                doInclude = TRUE;
                break;
            default :
                fprintf(stderr, "Something went wrong; behavior not implemented for given iProfileType.\n");
        }
    }
        
    if (gateCode & 1<<flagPositionDynamicClutterFringe) {
        
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
                doInclude = TRUE;
                break;
            default :
                fprintf(stderr, "Something went wrong; behavior not implemented for given iProfileType.\n");
        }
    }
    
    if (gateCode & 1<<flagPositionVradMissing) {
        
        // i.e. flag 3 in gateCode is true
        // this gate has reflectivity data but no corresponding radial velocity data
        
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
                fprintf(stderr, "Something went wrong; behavior not implemented for given iProfileType.\n");
        }
    }
    
    if (gateCode & 1<<flagPositionDbzTooHighForBirds) {
        
        // i.e. flag 4 in gateCode is true
        // this gate's dbz value is too high to be due to birds, it must be 
        // caused by something else

        
        switch (iProfileType) {
            case 1 : 
                doInclude = FALSE;
                break;
            case 2 : 
                doInclude = TRUE;
                break;
            case 3 : 
                doInclude = TRUE;
                break;
            default :
                fprintf(stderr, "Something went wrong; behavior not implemented for given iProfileType.\n");
        }
    }
    
    if (gateCode & 1<<flagPositionVradTooLow) {
        
        // i.e. flag 5 in gateCode is true
        // this gate's radial velocity is too low to be due to actual scatterers; likely just noise
        
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
                fprintf(stderr, "Something went wrong; behavior not implemented for given iProfileType.\n");
        }
    }
    
    if (gateCode & 1<<flagPositionVDifMax) {

        // i.e. flag 6 in gateCode is true
        // after the first svdfit, this gate's fitted vRad was more than 
        // VDIFMAX away from the observed vRad for that gate. It is therefore
        // considered an outlier
        
        switch (iProfileType) {
            case 1 : 
                doInclude = FALSE;
                break;
            case 2 : 
                // FIXME not sure if this should be FALSE for iProfileType==2
                doInclude = TRUE;
                break;
            case 3 : 
                // FIXME not sure if this should be FALSE for iProfileType==3
                doInclude = TRUE;
                break;
            default :
                fprintf(stderr, "Something went wrong; behavior not implemented for given iProfileType.\n");
        }
    }



    if (gateCode & 1<<flagPositionAzimTooLow) {

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
                fprintf(stderr, "Something went wrong; behavior not implemented for given iProfileType.\n");
        }
    }


    if (gateCode & 1<<flagPositionAzimTooHigh) {

        // i.e. flag 8 in gateCode is true
        // the user can specify to exclude gates based on their azimuth;
        // this clause is for gates that have too high azimuth
        
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
                fprintf(stderr, "Something went wrong; behavior not implemented for given iProfileType.\n");
        }
    }



    return doInclude;

} // includeGate





static int readUserConfigOptions(void) {


    cfg_opt_t opts[] = {
        CFG_FLOAT("HLAYER",200.0f, CFGF_NONE),
        CFG_INT("NLAYER",30, CFGF_NONE),
        CFG_FLOAT("RANGEMIN",5000.0f, CFGF_NONE),
        CFG_FLOAT("RANGEMAX",25000.0f, CFGF_NONE),
        CFG_FLOAT("AZIMMIN",0.0f, CFGF_NONE),
        CFG_FLOAT("AZIMMAX",360.0f, CFGF_NONE),
        CFG_FLOAT("RADAR_WAVELENGTH_CM",5.3f, CFGF_NONE),
        CFG_BOOL("USE_STATIC_CLUTTER_DATA",FALSE,CFGF_NONE),
        CFG_BOOL("VERBOSE_OUTPUT_REQUIRED",FALSE,CFGF_NONE),
        CFG_BOOL("PRINT_DBZ",TRUE,CFGF_NONE),
        CFG_BOOL("PRINT_VRAD",TRUE,CFGF_NONE),
        CFG_BOOL("PRINT_CELL",TRUE,CFGF_NONE),
        CFG_BOOL("PRINT_CELL_PROP",TRUE,CFGF_NONE),
        CFG_BOOL("PRINT_TEXTURE",TRUE,CFGF_NONE),
        CFG_BOOL("PRINT_CLUT",TRUE,CFGF_NONE),
        CFG_BOOL("PRINT_OPTIONS",TRUE,CFGF_NONE),
        CFG_BOOL("FIT_VRAD",TRUE,CFGF_NONE),
        CFG_BOOL("PRINT_PROFILE",TRUE,CFGF_NONE),
        CFG_BOOL("PRINT_POINTS_ARRAY",FALSE,CFGF_NONE),
        CFG_END()
    };
    
    cfg = cfg_init(opts, CFGF_NONE);
    
    if (cfg_parse(cfg, "options.conf") == CFG_PARSE_ERROR) {
        return 1;
    }
   
    return 0;

} // readUserConfigOptions
    





static int mapDataFromRave(PolarScan_t* scan, SCANMETA* meta, unsigned char* values, char* paramStr) {

    PolarScanParam_t* param = PolarScan_getParameter(scan,paramStr);

    int iRang;
    int iAzim;
    
    int nRang = PolarScan_getNbins(scan);
    int nAzim = PolarScan_getNrays(scan);
    int iGlobal;
    
    iGlobal = 0;    
    
    if (param != NULL) {

        // scan includes requested information, apparently
        
        meta->heig = (float) PolarScan_getHeight(scan);
        meta->elev = (float) (360 * PolarScan_getElangle(scan) / 2 / PI );
        meta->nRang = (int) PolarScan_getNbins(scan);
        meta->nAzim = (int) PolarScan_getNrays(scan);
        meta->rangeScale = (float) PolarScan_getRscale(scan);
        meta->azimScale = 360.0f/meta->nAzim;   // FIXME not sure if this always works
        meta->valueOffset = (float) PolarScanParam_getOffset(param);
        meta->valueScale = (float) PolarScanParam_getGain(param);
        meta->missing = (unsigned char) PolarScanParam_getNodata(param);

        double value;
        double* valuePtr = &value;
        unsigned char valueUChar;
        
        
        for (iAzim = 0; iAzim < nAzim; iAzim++) {
            for (iRang = 0; iRang < nRang; iRang++) {
            
                RaveValueType t = PolarScanParam_getValue(param, iRang, iAzim, valuePtr);
                
                if (0 <= value && value <=255) {
                    valueUChar = (unsigned char) value;
                }
                else {
                    fprintf(stderr,"Out of range value read at iRang = %d, iAzim = %d\n",iRang,iAzim);
                    return -1;
                }
                values[iGlobal] = valueUChar;
                iGlobal++;
            }
        }
    } 
    else {
        
        fprintf(stderr,"Requested parameter not available\n");
        return -1;

    }
    
    
    return 0;

} // mapDataFromRave




static void mapDataToRave(void) {
    
    // initialize the profile object
    VerticalProfile_t* profileRave = NULL;

    VerticalProfile_setLevels(profileRave,nLayers);
    VerticalProfile_setInterval(profileRave,layerThickness);

    RAVE_OBJECT_RELEASE(profileRave);

    
    return;
    
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
        fprintf(stderr,"There's only space for %d flags\n. Aborting",nFlagsMax);
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





static int printMeta(const SCANMETA* meta, const char* varName) {
    
    fprintf(stderr,"%s->heig = %f\n",varName,meta->heig);
    fprintf(stderr,"%s->elev = %f\n",varName,meta->elev);
    fprintf(stderr,"%s->nRang = %d\n",varName,meta->nRang);
    fprintf(stderr,"%s->nAzim = %d\n",varName,meta->nAzim);
    fprintf(stderr,"%s->rangeScale = %f\n",varName,meta->rangeScale);
    fprintf(stderr,"%s->azimScale = %f\n",varName,meta->azimScale);
    fprintf(stderr,"%s->valueOffset = %f\n",varName,meta->valueOffset);
    fprintf(stderr,"%s->valueScale = %f\n",varName,meta->valueScale);
    fprintf(stderr,"%s->missing = %d\n",varName,meta->missing);
    
    return 0;

} // printMeta





static void sortCells(CELLPROP *cellProp, const int nCells) {

    // ---------------------------------------------------------------- //
    // Sorting of the cell properties based on cell area.               //
    // Assume an area equal to zero for cells that are marked 'dropped' //
    // ---------------------------------------------------------------- // 

    int iCell;
    int iCellOther;
    CELLPROP tmp;

    /*Sorting of data elements using straight insertion method.*/
    for (iCell = 1; iCell < nCells; iCell++) {

        tmp = cellProp[iCell];

        iCellOther = iCell - 1;

        while (iCellOther >= 0 && cellProp[iCellOther].nGates * XABS(cellProp[iCellOther].drop - 1) < tmp.nGates * XABS(tmp.drop - 1)) {

            cellProp[iCellOther + 1] = cellProp[iCellOther];

            iCellOther--;
        }

        cellProp[iCellOther + 1] = tmp;

    } //for iCell

    return;
} // sortCells





static void updateFlagFieldsInPointsArray(const float* yObs, const float* yFitted, const int* includedIndex, 
                                   const int nPointsIncluded, float* points) {
                                       
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
        
        if (absVDif > absVDifMax) {
            
            iPoint = includedIndex[iPointIncluded];
            gateCode = (unsigned int) points[iPoint * nColsPoints + gateCodeCol];
            points[iPoint * nColsPoints + gateCodeCol] = (float) (gateCode |= 1<<flagPositionVDifMax);

        }
    } 

} // updateFlagFieldsInPointsArray




static int updateMap(int *cellImage, const int nGlobal, CELLPROP *cellProp, const int nCells) {

    // ------------------------------------------------------------------------- //
    // This function updates cellImage by dropping cells and reindexing the map. //
    // Leaving index 0 unused, will be used for assigning cell fringes           //
    // ------------------------------------------------------------------------- //

    int iGlobal;
    int iCell;
    int iCellNew;
    int nCellsValid;
    int cellImageValue;


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
    fprintf(stderr,"minimum value in cellImage array = %d.\n", minValue);
    fprintf(stderr,"maximum value in cellImage array = %d.\n", maxValue);
    #endif

    nCellsValid = nCells;

    for (iGlobal = 0; iGlobal < nGlobal; iGlobal++) {

        if (cellImage[iGlobal] == -1) {
            continue;
        }

        cellImageValue = cellImage[iGlobal];

        if (cellImageValue > nCells - 1) {
            fprintf(stderr, "You just asked for the properties of cell %d, which does not exist.\n", cellImageValue);
            continue;
        }

        if (cellProp[cellImageValue].drop == TRUE) {
            cellImage[iGlobal] = -1;
        }
    }


    // sort the cells by area and determine number of valid cells
    sortCells(cellProp, nCells);

    while (nCellsValid > 0 && cellProp[nCellsValid - 1].nGates < nGatesCellMin) {
        nCellsValid--;
    }

    #ifdef FPRINTFON
    fprintf(stderr,"nCellsValid = %d\n",nCellsValid);
    fprintf(stderr,"\n");
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
        fprintf(stderr,"before: cellProp[%d].index = %d.\n",iCell,cellProp[iCell].index);
        fprintf(stderr,"before: cellProp[%d].nGates = %d.\n",iCell,cellProp[iCell].nGates);
        fprintf(stderr,"before: iCell = %d.\n",iCell);
        fprintf(stderr,"before: iCellNew = %d.\n",iCellNew);
        fprintf(stderr,"\n");
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
        fprintf(stderr,"after: cellProp[%d].index = %d.\n",iCell,cellProp[iCell].index);
        fprintf(stderr,"after: cellProp[%d].nGates = %d.\n",iCell,cellProp[iCell].nGates);
        fprintf(stderr,"\n");
        #endif
    }

    return nCellsValid;
} // updateMap




void vol2birdCalcProfiles() {

    int nPasses;
    int iPoint;
    int iLayer;
    int iPass;
    int iProfileType;

    if (initializationSuccessful==FALSE) {
        fprintf(stderr,"You need to initialize vol2bird before you can use it. Aborting.\n");
        return;
    }
    
    // calculate the profiles in reverse order, because you need the result 
    // of iProfileType == 3 in order to check whether chi < stdDevMinBird  
    // when calculating iProfileType == 1
    for (iProfileType = nProfileTypes; iProfileType > 0; iProfileType--) {

        // ------------------------------------------------------------- //
        //                        prepare the profile                    //
        // ------------------------------------------------------------- //

        iProfileTypeLast = iProfileType;

        // if the user does not require fitting a model to the observed 
        // vrad values, one pass is enough 
        if (fitVrad == TRUE) {
            nPasses = 2;
        } 
        else {
            nPasses = 1;
        }

        // reset the flagPositionVDifMax bit before calculating each profile
        for (iPoint = 0; iPoint < nRowsPoints; iPoint++) {
            unsigned int gateCode = (unsigned int) points[iPoint * nColsPoints + gateCodeCol];
            points[iPoint * nColsPoints + gateCodeCol] = (float) (gateCode &= ~(1<<flagPositionVDifMax));
        }
        
        
        for (iLayer = 0; iLayer < nLayers; iLayer++) {

            // this variable is needed just outside of the iPass loop below 
            float chi = NAN;
        
            for (iPass = 0; iPass < nPasses; iPass++) {

                const int iPointFrom = indexFrom[iLayer];
                const int nPointsLayer = nPointsWritten[iLayer];

                int iPointLayer;
                int iPointIncluded;
                int nPointsIncluded;

                float parameterVector[] = {NAN,NAN,NAN};
                float avar[] = {NAN,NAN,NAN};
                
                float* pointsSelection = malloc(sizeof(float) * nPointsLayer * nDims);
                float* yObs = malloc(sizeof(float) * nPointsLayer);
                float* yFitted = malloc(sizeof(float) * nPointsLayer);
                int* includedIndex = malloc(sizeof(int) * nPointsLayer);
                
                float dbzValue = NAN;
                float undbzValue = NAN;
                double undbzSum = 0.0;
                float undbzAvg = NAN;
                float dbzAvg = NAN;
                float birdDensity = NAN;
                float reflectivity = NAN;
                int hasGap = TRUE;
                float chisq = NAN;
                float hSpeed = NAN;
                float hDir = NAN;


                for (iPointLayer = 0; iPointLayer < nPointsLayer; iPointLayer++) {

                    pointsSelection[iPointLayer * nDims + 0] = 0.0f;
                    pointsSelection[iPointLayer * nDims + 1] = 0.0f;
                    
                    yObs[iPointLayer] = 0.0f;
                    yFitted[iPointLayer] = 0.0f;
                    
                    includedIndex[iPointLayer] = -1;

                };

                profile[iLayer*nColsProfile +  0] = iLayer * layerThickness;
                profile[iLayer*nColsProfile +  1] = (iLayer + 1) * layerThickness;
                profile[iLayer*nColsProfile +  2] = NAN;
                profile[iLayer*nColsProfile +  3] = NAN;
                profile[iLayer*nColsProfile +  4] = NAN;
                profile[iLayer*nColsProfile +  5] = NAN;
                profile[iLayer*nColsProfile +  6] = NAN;
                profile[iLayer*nColsProfile +  7] = NAN;
                profile[iLayer*nColsProfile +  8] = NAN;
                profile[iLayer*nColsProfile +  9] = NAN;
                profile[iLayer*nColsProfile + 10] = NAN;
                profile[iLayer*nColsProfile + 11] = NAN;
                profile[iLayer*nColsProfile + 12] = NAN;

                iPointIncluded = 0;
                
                for (iPointLayer = iPointFrom; iPointLayer < iPointFrom + nPointsLayer; iPointLayer++) {
                    
                    unsigned int gateCode = (unsigned int) points[iPointLayer * nColsPoints + gateCodeCol];

                    if (includeGate(iProfileType,gateCode) == TRUE) {

                        // copy azimuth angle from the 'points' array
                        pointsSelection[iPointIncluded * nDims + 0] = points[iPointLayer * nColsPoints + azimAngleCol];
                        // copy elevation angle from the 'points' array
                        pointsSelection[iPointIncluded * nDims + 1] = points[iPointLayer * nColsPoints + elevAngleCol];
                        // get the dbz value at this [azimuth, elevation] 
                        dbzValue = points[iPointLayer * nColsPoints + dbzValueCol];
                        // convert from dB scale to linear scale 
                        undbzValue = (float) exp(0.1*log(10)*dbzValue);
                        // sum the undbz in this layer
                        undbzSum += undbzValue;
                        // copy the observed vrad value at this [azimuth, elevation] 
                        yObs[iPointIncluded] = points[iPointLayer * nColsPoints + vradValueCol];
                        // pre-allocate the fitted vrad value at this [azimuth,elevation]
                        yFitted[iPointIncluded] = 0.0f;
                        // keep a record of which index was just included
                        includedIndex[iPointIncluded] = iPointLayer;
                        
                        // raise the counter
                        iPointIncluded += 1;
                        
                    }
                } // endfor (iPointLayer = 0; iPointLayer < nPointsLayer; iPointLayer++) {
                nPointsIncluded = iPointIncluded;
                
                
                if (nPointsIncluded > nPointsIncludedMin) {   
                    // when there are enough valid points, convert undbzAvg back to dB-scale
                    undbzAvg = (float) (undbzSum/nPointsIncluded);
                    dbzAvg = (10*log(undbzAvg))/log(10);
                }
                else {
                    undbzAvg = NAN;
                    dbzAvg = NAN;
                }

                // convert from Z (not dBZ) in units of mm^6/m^3 to 
                // reflectivity eta in units of cm^2/km^3
                reflectivity = dbzFactor * undbzAvg;
                
                if (iProfileType == 1) {
                    // calculate bird density in number of birds/km^3 by
                    // dividing the reflectivity by the (assumed) cross section
                    // of one bird
                    birdDensity = reflectivity / birdRadarCrossSection;
                }
                else {
                    birdDensity = NAN;
                }

                // check if there are directions that have almost no observations
                // (as this makes the svdfit result really uncertain)  
                hasGap = hasAzimuthGap(&pointsSelection[0], nPointsIncluded);
                
                if (fitVrad == TRUE) {

                    if (hasGap==FALSE) {
                        
                        // ------------------------------------------------------------- //
                        //                       do the svdfit                           //
                        // ------------------------------------------------------------- //

                        chisq = svdfit(&pointsSelection[0], nDims, &yObs[0], &yFitted[0], 
                                nPointsIncluded, &parameterVector[0], &avar[0], nParsFitted);

                        if (chisq < chisqMin) {
                            // the fit was not good enough, reset parameter vector array 
                            // elements to NAN and continue with the next layer
                            parameterVector[0] = NAN;
                            parameterVector[1] = NAN;
                            parameterVector[2] = NAN;
                            continue;
                        } 
                        else {
                            
                            chi = sqrt(chisq);
                            hSpeed = sqrt(pow(parameterVector[0],2) + pow(parameterVector[1],2));
                            hDir = (atan2(parameterVector[0],parameterVector[1])*RAD2DEG);
                            
                            if (hDir < 0) {
                                hDir += 360.0f;
                            }
                            
                        }
                    } // endif (hasGap == FALSE)
                    
                    // if the fitted vrad value is more than 'absVDifMax' away from the corresponding
                    // observed vrad value, set the gate's flagPositionVDifMax bit flag to 1, excluding the 
                    // gate in the second svdfit iteration
                    updateFlagFieldsInPointsArray(&yObs[0], &yFitted[0], &includedIndex[0], nPointsIncluded,
                                                  &points[0]);

                }; // endif (fitVrad == TRUE)

                profile[iLayer*nColsProfile +  0] = iLayer * layerThickness;
                profile[iLayer*nColsProfile +  1] = (iLayer + 1) * layerThickness;
                profile[iLayer*nColsProfile +  2] = parameterVector[0];
                profile[iLayer*nColsProfile +  3] = parameterVector[1];
                profile[iLayer*nColsProfile +  4] = parameterVector[2];
                profile[iLayer*nColsProfile +  5] = hSpeed;
                profile[iLayer*nColsProfile +  6] = hDir;
                profile[iLayer*nColsProfile +  7] = chi;
                profile[iLayer*nColsProfile +  8] = (float) hasGap;
                profile[iLayer*nColsProfile +  9] = dbzAvg;
                profile[iLayer*nColsProfile + 10] = (float) nPointsIncluded;
                profile[iLayer*nColsProfile + 11] = reflectivity;
                profile[iLayer*nColsProfile + 12] = birdDensity;
                
                free((void*) yObs);
                free((void*) yFitted);
                free((void*) pointsSelection);
        
            } // endfor (iPass = 0; iPass < nPasses; iPass++)
            
            
            // You need some of the results of iProfileType == 3 in order 
            // to calculate iProfileType == 1
            if (iProfileType == 3) {
                if (chi < stdDevMinBird) {
                    scatterersAreNotBirds[iLayer] = TRUE;
                }
                else {
                    scatterersAreNotBirds[iLayer] = FALSE;
                }
            }
            
        } // endfor (iLayer = 0; iLayer < nLayers; iLayer++)

        if (printProfileVar ==  TRUE) {
            printProfile();
        }

    } // endfor (iProfileType = nProfileTypes; iProfileType > 0; iProfileType--)

} // vol2birdCalcProfiles




void printImageInt(const SCANMETA* meta, const int* imageInt) {


    int nRang = meta->nRang;
    int nAzim = meta->nAzim;
    int iRang;
    int iAzim;
    int iGlobal;
    int needsSignChar;
    int maxValue;
    

    int thisValue;
    int nChars;
    char* formatStr;
    
    iGlobal = 0;
    maxValue = 0;
    needsSignChar = FALSE;
    
    fprintf(stderr,"elevAngle = %f degrees\n",meta->elev);
    
    
    // first, determine how many characters are needed to print array 'imageInt'
    for (iAzim = 0; iAzim < nAzim; iAzim++) { 
        
        for (iRang = 0; iRang < nRang; iRang++) {
            
            thisValue = (int) XABS(imageInt[iGlobal]);
            if (imageInt[iGlobal] < 0) {
                needsSignChar = TRUE;
            }
            if (thisValue > maxValue) {
                maxValue = thisValue;
            } ;
            
            iGlobal += 1;
            
        }
    }


    nChars = (int) ceil(log(maxValue + 1)/log(10));

    if (needsSignChar) {
        nChars += 1;
    }
    
    switch (nChars) {
        case 0 :
            formatStr = " %1d";
            break; 
        case 1 :
            formatStr = " %1d";
            break; 
        case 2 :
            formatStr = " %2d"; 
            break;
        case 3 :
            formatStr = " %3d"; 
            break;
        case 4 :
            formatStr = " %4d"; 
            break;
        default :
            formatStr = " %8d"; 
    }
    
    
    iGlobal = 0;
    
    for (iAzim = 0; iAzim < nAzim; iAzim++) { 
        
        for (iRang = 0; iRang < nRang; iRang++) {
    
            iGlobal = iRang + iAzim * nRang;
            
            thisValue = (int) imageInt[iGlobal];
            
            fprintf(stderr,formatStr,thisValue);
            
            iGlobal += 1;
            
        }
        fprintf(stderr,"\n");
    }
        
} // printImageInt




void printImageUChar(const SCANMETA* meta, const unsigned char* imageUChar) {

    int nAzim;
    int iAzim;
    int nRang;
    int iRang;
    int iGlobal;

    nAzim = meta->nAzim;
    nRang = meta->nRang;
    
    int* imageInt = malloc(sizeof(int) * nRang * nAzim);
    iGlobal = 0;
    
    for (iAzim = 0; iAzim < nAzim; iAzim++) {
        for (iRang = 0; iRang < nRang; iRang++) {
            
            iGlobal = iRang + iAzim * nRang;
            imageInt[iGlobal] = (int) imageUChar[iGlobal];
            iGlobal += 1;
            
        }
    }     

    printImageInt(meta,imageInt);
    
    free(imageInt);
    
} // printImageUChar




void vol2birdPrintIndexArrays(void) {
    
    if (initializationSuccessful==FALSE) {
        fprintf(stderr,"You need to initialize vol2bird before you can use it. Aborting.\n");
        return;
    }
    
    int iLayer;

    fprintf(stderr, "iLayer  iFrom   iTo     iTo-iFrom nWritten\n");
    for (iLayer = 0; iLayer < nLayers; iLayer++) {
        fprintf(stderr, "%7d %7d %7d %10d %8d\n",
            iLayer, 
            indexFrom[iLayer], 
            indexTo[iLayer], 
            indexTo[iLayer] - indexFrom[iLayer], 
            nPointsWritten[iLayer]);
    }
} // vol2birdPrintIndexArrays



void vol2birdPrintOptions(void) {
    
    // ------------------------------------------------------- //
    // this function prints vol2bird's configuration to stderr //
    // ------------------------------------------------------- //
    
    if (initializationSuccessful==FALSE) {
        fprintf(stderr,"You need to initialize vol2bird before you can use it. Aborting.\n");
        return;
    }

    fprintf(stderr,"\n\nvol2bird configuration:\n\n");

    fprintf(stderr,"%-25s = %f\n","absVDifMax",absVDifMax);
    fprintf(stderr,"%-25s = %f\n","azimMax",azimMax);
    fprintf(stderr,"%-25s = %f\n","azimMin",azimMin);
    fprintf(stderr,"%-25s = %f\n","birdRadarCrossSection",birdRadarCrossSection);
    fprintf(stderr,"%-25s = %f\n","cellClutterFractionMax",cellClutterFractionMax);
    fprintf(stderr,"%-25s = %f\n","cellDbzMin",cellDbzMin);
    fprintf(stderr,"%-25s = %f\n","cellStdDevMax",cellStdDevMax);
    fprintf(stderr,"%-25s = %f\n","chisqMin",chisqMin);
    fprintf(stderr,"%-25s = %f\n","clutterValueMin",clutterValueMin);
    fprintf(stderr,"%-25s = %f\n","dbzFactor",dbzFactor);
    fprintf(stderr,"%-25s = %f\n","dbzMax",dbzMax);
    fprintf(stderr,"%-25s = %f\n","dbzThresMin",dbzThresMin);
    fprintf(stderr,"%-25s = %f\n","fringeDist",fringeDist);
    fprintf(stderr,"%-25s = %f\n","layerThickness",layerThickness);
    fprintf(stderr,"%-25s = %d\n","nGatesCellMin",nGatesCellMin);
    fprintf(stderr,"%-25s = %d\n","nAzimNeighborhood",nAzimNeighborhood);
    fprintf(stderr,"%-25s = %d\n","nBinsGap",nBinsGap);
    fprintf(stderr,"%-25s = %d\n","nCountMin",nCountMin);
    fprintf(stderr,"%-25s = %d\n","nLayers",nLayers);
    fprintf(stderr,"%-25s = %d\n","nObsGapMin",nObsGapMin);
    fprintf(stderr,"%-25s = %d\n","nPointsIncludedMin",nPointsIncludedMin);
    fprintf(stderr,"%-25s = %d\n","nRangNeighborhood",nRangNeighborhood);
    fprintf(stderr,"%-25s = %f\n","radarWavelength",radarWavelength);
    fprintf(stderr,"%-25s = %f\n","rangeMax",rangeMax);
    fprintf(stderr,"%-25s = %f\n","rangeMin",rangeMin);
    fprintf(stderr,"%-25s = %f\n","rCellMax",rCellMax);
    fprintf(stderr,"%-25s = %f\n","refracIndex",refracIndex);
    fprintf(stderr,"%-25s = %f\n","stdDevMinBird",stdDevMinBird);
    fprintf(stderr,"%-25s = %c\n","useStaticClutterData",useStaticClutterData == TRUE ? 'T' : 'F');
    fprintf(stderr,"%-25s = %f\n","vradMin",vradMin);
    
    fprintf(stderr,"\n\n");

}  // vol2birdPrintOptions





void vol2birdPrintPointsArray(void) {
    
    // ------------------------------------------------- //
    // this function prints the 'points' array to stderr //
    // ------------------------------------------------- //
    
    if (initializationSuccessful==FALSE) {
        fprintf(stderr,"You need to initialize vol2bird before you can use it. Aborting.\n");
        return;
    }

    int iPoint;
    
    fprintf(stderr, "iPoint  azim    elev    dbz         vrad        cell    gateCode  flags     \n");
    
    for (iPoint = 0; iPoint < nRowsPoints * nColsPoints; iPoint+=nColsPoints) {
        
            char gateCodeStr[10];  // 9 bits plus 1 position for the null character '\0'
            
            printGateCode(&gateCodeStr[0], (int) points[iPoint + gateCodeCol]);
        
            fprintf(stderr, "  %6d",    iPoint/nColsPoints);
            fprintf(stderr, "  %6.2f",  points[iPoint + azimAngleCol]);
            fprintf(stderr, "  %6.2f",  points[iPoint + elevAngleCol]);
            fprintf(stderr, "  %10.2f", points[iPoint + dbzValueCol]);
            fprintf(stderr, "  %10.2f", points[iPoint + vradValueCol]);
            fprintf(stderr, "  %6.0f",  points[iPoint + cellValueCol]);
            fprintf(stderr, "  %8.0f",  points[iPoint + gateCodeCol]);
            fprintf(stderr, "  %12s",   gateCodeStr);
            fprintf(stderr, "\n");
    }    
} // vol2birdPrintPointsArray




void printProfile(void) {
    
    if (initializationSuccessful==FALSE) {
        fprintf(stderr,"You need to initialize vol2bird before you can use it. Aborting.\n");
        return;
    }

    fprintf(stderr,"\n\nProfile type: %d\n",iProfileTypeLast);

    fprintf(stderr,"altmin-altmax: [u         ,v         ,w         ]; "
                   "hSpeed  , hDir    , chi     , hasGap  , dbzAvg  ,"
                   " nPoints, eta         , rhobird \n");

    int iLayer;
    
    for (iLayer = nLayers - 1; iLayer >= 0; iLayer--) {
        
        fprintf(stderr,"%6.0f-%-6.0f: [%10.2f,%10.2f,%10.2f]; %8.2f, "
        "%8.1f, %8.1f, %8c, %8.2f, %7.0f, %12.2f, %8.2f\n",
        profile[iLayer * nColsProfile +  0],
        profile[iLayer * nColsProfile +  1],
        profile[iLayer * nColsProfile +  2],
        profile[iLayer * nColsProfile +  3],
        profile[iLayer * nColsProfile +  4],
        profile[iLayer * nColsProfile +  5],
        profile[iLayer * nColsProfile +  6],
        profile[iLayer * nColsProfile +  7],
        profile[iLayer * nColsProfile +  8] == TRUE ? 'T' : 'F',
        profile[iLayer * nColsProfile +  9],
        profile[iLayer * nColsProfile + 10],
        profile[iLayer * nColsProfile + 11],
        profile[iLayer * nColsProfile + 12]);
    }

    
} // printProfile()




int vol2birdSetUp(PolarVolume_t* volume) {

    if (readUserConfigOptions() != 0) {
        fprintf(stderr, "An error occurred while reading the user configuration file 'options.conf'.\n");
        return -1; 
    }
   
    // ------------------------------------------------------------- //
    //              vol2bird options from options.conf               //
    // ------------------------------------------------------------- //

    azimMax = cfg_getfloat(cfg, "AZIMMAX");
    azimMin = cfg_getfloat(cfg, "AZIMMIN");
    layerThickness = cfg_getfloat(cfg, "HLAYER");
    nLayers = cfg_getint(cfg, "NLAYER");
    rangeMax = cfg_getfloat(cfg, "RANGEMAX");
    rangeMin = cfg_getfloat(cfg, "RANGEMIN");
    radarWavelength = cfg_getfloat(cfg, "RADAR_WAVELENGTH_CM");
    useStaticClutterData = cfg_getbool(cfg,"USE_STATIC_CLUTTER_DATA");
    printDbz = cfg_getbool(cfg,"PRINT_DBZ");
    printVrad = cfg_getbool(cfg,"PRINT_VRAD");
    printTex = cfg_getbool(cfg,"PRINT_TEXTURE");
    printCell = cfg_getbool(cfg,"PRINT_CELL");
    printCellProp = cfg_getbool(cfg,"PRINT_CELL_PROP");
    printClut = cfg_getbool(cfg,"PRINT_CLUT");
    printOptions = cfg_getbool(cfg,"PRINT_OPTIONS");
    printProfileVar = cfg_getbool(cfg,"PRINT_PROFILE");
    printPointsArray = cfg_getbool(cfg,"PRINT_POINTS_ARRAY");
    fitVrad = cfg_getbool(cfg,"FIT_VRAD");



    // ------------------------------------------------------------- //
    //              vol2bird options from constants.h                //
    // ------------------------------------------------------------- //

    nGatesCellMin = AREACELL;
    cellClutterFractionMax = CLUTPERCCELL;
    cellDbzMin = DBZCELL;
    chisqMin = CHISQMIN;
    clutterValueMin = DBZCLUTTER;
    dbzMax = DBZMAX;
    dbzThresMin = DBZMIN;
    fringeDist = FRINGEDIST;
    nBinsGap = NBINSGAP;
    nPointsIncludedMin = NDBZMIN;
    nNeighborsMin = NEIGHBORS;
    nObsGapMin = NOBSGAPMIN;
    nAzimNeighborhood = NTEXBINAZIM;
    nRangNeighborhood = NTEXBINRANG;
    nCountMin = NTEXMIN; 
    refracIndex = REFRACTIVE_INDEX_OF_WATER;
    birdRadarCrossSection = SIGMABIRD;
    cellStdDevMax = STDEVCELL;
    stdDevMinBird = STDEVBIRD;
    absVDifMax = VDIFMAX;
    vradMin = VRADMIN;




    // ------------------------------------------------------------- //
    //                       some other variables                    //
    // ------------------------------------------------------------- //

    rCellMax = rangeMax + 5000.0f;
    nDims = 2;
    nParsFitted = 3;
    dbzFactor = (pow(refracIndex,2) * 1000 * pow(PI,5))/pow(radarWavelength,4);




    // ------------------------------------------------------------- //
    //             lists of indices into the 'points' array:         //
    //          where each altitude layer's data starts and ends     //
    // ------------------------------------------------------------- //
    
    int iLayer;
    
    // pre-allocate the list with start-from indexes for each 
    // altitude bin in the profile
    indexFrom = (int*) malloc(sizeof(int) * nLayers);
    if (indexFrom == NULL) {
        fprintf(stderr,"Error pre-allocating array 'indexFrom'\n");
        return -1;
    }
    for (iLayer = 0; iLayer < nLayers; iLayer++) {
        indexFrom[iLayer] = 0;
    }

    // pre-allocate the list with end-before indexes for each 
    // altitude bin in the profile
    indexTo = (int*) malloc(sizeof(int) * nLayers);
    if (indexTo == NULL) {
        fprintf(stderr,"Error pre-allocating array 'indexTo'\n");
        return -1;
    }
    for (iLayer = 0; iLayer < nLayers; iLayer++) {
        indexTo[iLayer] = 0;
    }

    // pre-allocate the list containing TRUE or FALSE depending on the 
    // results of calculating iProfileType == 3, which are needed when
    // calculating iProfileType == 1
    scatterersAreNotBirds = (int*) malloc(sizeof(int) * nLayers);
    if (scatterersAreNotBirds == NULL) {
        fprintf(stderr,"Error pre-allocating array 'scatterersAreNotBirds'\n");
        return -1;
    }
    for (iLayer = 0; iLayer < nLayers; iLayer++) {
        scatterersAreNotBirds[iLayer] = -1;
    }


    // for each altitude layer, you need to remember how many points 
    // were already written. This information is stored in the 
    // 'nPointsWritten' array
    nPointsWritten = (int*) malloc(sizeof(int) * nLayers);
    if (nPointsWritten == NULL) {
        fprintf(stderr,"Error pre-allocating array 'nPointsWritten'\n");
        return -1;
    }
    for (iLayer = 0; iLayer < nLayers; iLayer++) {
        nPointsWritten[iLayer] = 0;
    }




    // ------------------------------------------------------------- //
    //               information about the 'points' array            //
    // ------------------------------------------------------------- //

    nColsPoints = 6;
    nRowsPoints = detSvdfitArraySize(volume);

    azimAngleCol = 0;
    elevAngleCol = 1;
    dbzValueCol = 2;
    vradValueCol = 3;
    cellValueCol = 4;
    gateCodeCol = 5;

    // pre-allocate the 'points' array (note it has 'nColsPoints'
    // pseudo-columns)
    points = (float*) malloc(sizeof(float) * nRowsPoints * nColsPoints);
    if (points == NULL) {
        fprintf(stderr,"Error pre-allocating array 'points'.\n"); 
        return -1;
    }

    int iRowPoints;
    int iColPoints;
        
    for (iRowPoints = 0; iRowPoints < nRowsPoints; iRowPoints++) {
        for (iColPoints = 0; iColPoints < nColsPoints; iColPoints++) {
            points[iRowPoints*nColsPoints + iColPoints] = NAN;
        }
    }

    // information about the flagfields of 'gateCode'
    
    flagPositionStaticClutter = 0;
    flagPositionDynamicClutter = 1;
    flagPositionDynamicClutterFringe = 2;
    flagPositionVradMissing = 3;
    flagPositionDbzTooHighForBirds = 4;
    flagPositionVradTooLow = 5;
    flagPositionVDifMax = 6;
    flagPositionAzimTooLow = 7;
    flagPositionAzimTooHigh = 8;

    // construct the 'points' array
    constructPointsArray(volume);

    // classify the gates based on the data in 'points'
    classifyGatesSimple();

    if (printPointsArray == TRUE) {

        vol2birdPrintIndexArrays();
        vol2birdPrintPointsArray();

    }

    // ------------------------------------------------------------- //
    //              information about the 'profile' array            //
    // ------------------------------------------------------------- //

    nProfileTypes = 3;
    nRowsProfile = nLayers;
    nColsProfile = 13; 
    
    // pre-allocate the array holding any profiled data (note it has 
    // 'nColsProfile' pseudocolumns):
    profile = (float*) malloc(sizeof(float) * nRowsProfile * nColsProfile);
    if (profile == NULL) {
        fprintf(stderr,"Error pre-allocating array 'profile'.\n"); 
        return -1;
    }

    int iRowProfile;
    int iColProfile;
        
    for (iRowProfile = 0; iRowProfile < nRowsProfile; iRowProfile++) {
        for (iColProfile = 0; iColProfile < nColsProfile; iColProfile++) {
            profile[iRowProfile*nColsProfile + iColProfile] = NAN;
        }
    }

    iProfileTypeLast = -1;

    initializationSuccessful = TRUE;

    if (printOptions == TRUE) {
        vol2birdPrintOptions();
    }

    return 0;

} // vol2birdSetUp




void vol2birdTearDown() {
    
    // ---------------------------------------------------------- //
    // free the memory that was previously allocated for vol2bird //
    // ---------------------------------------------------------- //

    if (initializationSuccessful==FALSE) {
        fprintf(stderr,"You need to initialize vol2bird before you can use it. Aborting.\n");
        return;
    }

    // free the points array, the indexes into it, the counters, as well
    // as the profile data array
    free((void*) points);
    free((void*) profile);
    free((void*) indexFrom);
    free((void*) indexTo);
    free((void*) nPointsWritten);
    
    
    // free the memory that holds the user configurable options
    cfg_free(cfg);
    
    // reset this variable to its initial value
    initializationSuccessful = FALSE;

} // vol2birdTearDown





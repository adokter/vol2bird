#include "cartesian.h"
#include "polarvolume.h"

Cartesian_t* polarVolumeToCartesian(PolarVolume_t* pvol, long dim, long res, double init);

double distance2height(double distance,double elev);

double distance2range(double distance,double elev);

double range2distance(double range,double elev);

double range2height(double range,double elev);

double*** init3DTensor(int dim1, int dim2, int dim3, double init);

float**** create4DTensor(float *array, int dim1, int dim2, int dim3, int dim4);

PolarVolume_t* PolarVolume_selectScansByElevation(PolarVolume_t* volume, float elevs[], int nElevs);

int polarVolumeTo3DTensor(PolarVolume_t* pvol, double ****tensor, int dim, long res, int nParam);

int fill3DTensor(double ***tensor, RaveObjectList_t* list, int dim1, int dim2, int dim3);

float* flatten3DTensor(double ***tensor, int dim1, int dim2, int dim3);

void free3DTensor(double ***tensor, int dim1, int dim2);

void free4DTensor(float ****tensor, int dim1, int dim2, int dim3);

#ifdef MISTNET 
PolarVolume_t* segmentScansUsingMistnet(PolarVolume_t* volume, vol2bird_t* alldata);
#endif


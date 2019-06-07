#include "cartesian.h"

Cartesian_t* polarVolumeToCartesian(PolarVolume_t* pvol, float elevs[], int nElevs, long dim, long res, double init);

double distance2height(double distance,double elev);

double distance2range(double distance,double elev);

double*** init3DTensor(int dim1, int dim2, int dim3, double init);

float**** create4DTensor(int dim1, int dim2, int dim3, int dim4);

int polarVolumeTo3DTensor(PolarVolume_t* pvol, double ****array, float elevs[], int nElevs, int dim, double res, int nParam);

int fill3DTensor(double ***array, RaveObjectList_t* list, int dim1, int dim2, int dim3);

float* flatten3DTensor(double ***array, int dim1, int dim2, int dim3);

void free1DTensor(float *array, int dim);

void free3DTensor(double ***array, int dim1, int dim2);

void free4DTensor(float ****array, int dim1, int dim2, int dim3);
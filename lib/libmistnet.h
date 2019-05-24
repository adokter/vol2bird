#include "cartesian.h"

Cartesian_t* polarVolumeToCartesian(PolarVolume_t* pvol, double elevs[], int nElevs, long dim, long res, double init);

double distance2height(double distance,double elev);

double distance2range(double distance,double elev);

void free3DArray(double ***array, int dim1, int dim2);

double*** init3DArray(int dim1, int dim2, int dim3, double init);

int mistnet3DArray(double ****array, PolarVolume_t* pvol, float elevs[], int nElevs, int dim, double res);

int fill3DArray(double ***array, Cartesian_t* cartesian, int dim1, int dim2, int dim3);
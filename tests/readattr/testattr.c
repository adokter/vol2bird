#include <rave_io.h>
#include <polarvolume.h>
#include <stdio.h>

int main(int argc, char** argv) {
  const char* filename = argv[1];
  double value = 0.0; 

  // read the input file and assign it to a generic rave object
  RaveIO_t* raveio = RaveIO_open(filename);
  if (raveio == (RaveIO_t*) NULL){
        fprintf(stderr,"cannot open file ...\n");
        return -1;
  }

  // extract polar volume object
  PolarVolume_t* pvol = (PolarVolume_t*) RaveIO_getObject(raveio);
  if (pvol == (PolarVolume_t*) NULL){
        fprintf(stderr,"cannot open volume ...\n");
        return -1;
  }  


  // extract wavelength attribute
  RaveAttribute_t* attr = PolarVolume_getAttribute(pvol, "how/RXloss");

    if (attr != (RaveAttribute_t *) NULL){
        RaveAttribute_getDouble(attr, &value);
        fprintf(stderr, "wavelength attribute found!\n");
    }
    else{
 	fprintf(stderr, "no attribute found ...\n");
    }

}

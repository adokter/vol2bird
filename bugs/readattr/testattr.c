#include <rave_io.h>
#include <polarvolume.h>
#include <stdio.h>

int main(int argc, char** argv) {
  const char* filename = argv[1];
  double value; 

  // read the input file and assign it to a generic rave object
  RaveIO_t* raveio = RaveIO_open(filename);

  // extract polar volume object
  PolarVolume_t* pvol = (PolarVolume_t*) RaveIO_getObject(raveio);

  // extract wavelength attribute
  RaveAttribute_t* attr = PolarVolume_getAttribute(pvol, "/how/wavelength");
    if (attr != (RaveAttribute_t *) NULL){
        RaveAttribute_getDouble(attr, &value);
        fprintf(stderr, "wavelength attribute found!\n");
    }
    else{
 	fprintf(stderr, "no attribute found ...\n");
    }

}

#include <rave_io.h>
#include <stdio.h>

int main(int argc, char** argv) {
  const char* filename = argv[1];

  // read the input file and assign it to a generic rave object
  RaveIO_t* raveio = RaveIO_open(filename);
  if (raveio == (RaveIO_t*) NULL){
        fprintf(stderr,"cannot open file ...\n");
        return -1;
  }
  else{
	fprintf(stderr,"file opened succesfully ...\n");
        return -1;
  }

}

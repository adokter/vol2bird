## vol2bird algorithm
Copyright 2015 Adriaan Dokter (University of Amsterdam) & Netherlands eScience Centre

If you are interested in collaborating and using this software, please contact me at a.m.dokter@uva.nl

This repository contains code to generate vertical profiles of birds (VPBs), following the methods described in:

*Bird migration flight altitudes studied by a network of operational weather radars*
Dokter A.M., Liechti F., Stark H., Delobbe L., Tabary P., Holleman I.
J. R. Soc. Interface, **8**, 30–43, 2011.
DOI: [10.1098/rsif.2010.0116](https://doi.org/10.1098/rsif.2010.0116)

* `lib` contains the main library 
* `pyvol2bird` is a python wrapper for the library

for compiling the code run `./configure --with-confuse=DIR`,
with `DIR` the root directory of your libconfuse installation (containing subfolders `DIR/include` and `DIR/lib` with the
libconfuse include and library files)

for running the code on a baltrad node add the paths listed in `lib/ldlib.txt` to `LD_LIBRARY_PATH`

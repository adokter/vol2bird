## vol2bird algorithm
Copyright 2015 Adriaan Dokter & Netherlands eScience Centre

If you want to collaborate and use this software, please contact me at a.m.dokter@uva.nl

This repository contains code to generate vertical profiles of birds (VPBs), following the methods described in:

*Bird migration flight altitudes studied by a network of operational weather radars*
Dokter A.M., Liechti F., Stark H., Delobbe L., Tabary P., Holleman I.
J. R. Soc. Interface, **8**, 30–43, 2011.
DOI: [10.1098/rsif.2010.0116](https://doi.org/10.1098/rsif.2010.0116)

1. `testbaltrad.c` is the program's main
2. `libvol2bird` is the library containing the routines for calculating VPBs
3. `libsvdfit` is a library used by libvol2bird to perform singular value decomposition, for fitting a VVP radial velocity model
4. `options.conf` for specifying configurable options
5. `constants.h` contains fixed constants

for compiling the code use the `Makefile`.

for running the code on a baltrad node add the paths listed in `ldlib.txt` to `LD_LIBRARY_PATH`

## vol2bird algorithm
Copyright 2017 Adriaan Dokter (University of Amsterdam, Cornell lab of ornithology) & Netherlands eScience Centre

If you are interested in using this software through collaboration, please contact me at a.m.dokter@uva.nl

This repository contains code to generate vertical profiles of birds (VPBs), following the methods described in:

*Bird migration flight altitudes studied by a network of operational weather radars*
Dokter A.M., Liechti F., Stark H., Delobbe L., Tabary P., Holleman I.
J. R. Soc. Interface, **8**, 30–43, 2011.
DOI: [10.1098/rsif.2010.0116](https://doi.org/10.1098/rsif.2010.0116)

* `data` contains two polar volume radar files for testing
* `doc` documentation on ODIM hdf5 data format
* `docker` scripts for generating and running a Docker container for vol2bird
* `etc` contains the configuration file `options.conf`. Put a copy in your working directory to load it on runtime
* `lib` contains the main library
* `pgfplugin` product generation framework (PGF) plugin for the BALTRAD system
* `pyvol2bird` is a python wrapper for the library
* `src` contains main executables for vol2bird (main program) and rsl2odim (converts NEXRAD to ODIM data format)
* `tests` unit tests for pgfplugin


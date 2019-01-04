## Vertical profiling of biological scatterers
**vol2bird** is an algorithm that calculates vertical profiles of birds and other biological scatterers from weather radar data. **vol2bird** is written in C.

See this publications for a description of the methodology:

[**Bird migration flight altitudes studied by a network of operational weather radars**](https://doi.org/10.1098/rsif.2010.0116)  
Dokter AM, Liechti F, Stark H, Delobbe L, Tabary P, Holleman I  
J. R. Soc. Interface, **8**, 30–43, 2011, DOI [10.1098/rsif.2010.0116](https://doi.org/10.1098/rsif.2010.0116)

Recent algorithm extensions and integration with the [bioRad R package](http://adokter.github.io/bioRad) are described here:

[**bioRad: biological analysis and visualization of weather radar data**](https://doi.org/10.1111/ecog.04028)  
Dokter AM, Desmet P, Spaaks JH, van Hoey S, Veen L, Verlinden L, Nilsson C, Haase G, Leijnse H, Farnsworth A, Bouten W, Shamoun-Baranes J.  
Ecography, 2018, DOI [10.1111/ecog.04028](https://doi.org/10.1111/ecog.04028)  

### repository contents

Browse source code at:
[https://github.com/adokter/vol2bird](https://github.com/adokter/vol2bird)

Report a bug at:
[https://github.com/adokter/vol2bird/issues](https://github.com/adokter/vol2bird/issues)

* `data` contains two polar volume radar files for testing
* `doc` documentation on ODIM hdf5 data format
* `docs` html pages of this documentation
* `docker` scripts for generating and running a Docker container for vol2bird
* `doxygen` doxygen configuration for generating this documentation
* `etc` contains the user configuration file `options.conf`. Put your own modified copy in your working directory to run vol2bird with non-default settings.
* `lib` contains the main library
* `pgfplugin` product generation framework (PGF) plugin for the BALTRAD system
* `pyvol2bird` is a python wrapper for the library
* `src` contains main executables for vol2bird (main program) and rsl2odim (converts NEXRAD to ODIM data format)
* `tests` unit tests for pgfplugin

Copyright 2010-2018 Adriaan M. Dokter (University of Amsterdam, Cornell lab of ornithology) & Netherlands eScience Centre

Contributors: Jurriaan H. Spaaks (NLeSC), Lourens Veen (NLeSC), Iwan Holleman (KNMI & Radboud University)

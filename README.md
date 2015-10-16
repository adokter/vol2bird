## Vol2Bird algorithm

This repository contains code to generate vertical profiles of birds (VPBs), following the methods described in:

*Bird migration flight altitudes studied by a network of operational weather radars*
Dokter A.M., Liechti F., Stark H., Delobbe L., Tabary P., Holleman I.
J. R. Soc. Interface, **8**, 30–43, 2011.
DOI: [10.1098/rsif.2010.0116](https://doi.org/10.1098/rsif.2010.0116)


Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.


1. `testbaltrad.c` is the program's main
2. `libvol2bird` is the library containing the routines for calculating VPBs
3. `libsvdfit` is a library used by libvol2bird to perform singular value decomposition, for fitting a VVP radial velocity model
4. `options.conf` for specifying configurable options
5. `constants.h` contains fixed constants

for compiling the code use the `Makefile`.

for running the code on a baltrad node add the paths listed in `ldlib.txt` to `LD_LIBRARY_PATH`

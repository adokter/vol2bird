# BALTRAD bird algorithm

This repository contains code to generate vertical profiles of birds (VPBs).

1. `testbaltrad.c` is the program's main
2. `libvol2bird` is the library containing the routines for calculating VPBs
3. `libsvdfit` is a library used by libvol2bird to perform singular value decomposition, for fitting a VVP radial velocity model
4. `options.conf` for specifying configurable options
5. `constants.h` contains fixed constants

for compiling the code use the `Makefile`.

for running the code on a baltrad node add the paths listed in `ldlib.txt` to `LD_LIBRARY_PATH`



To update on github:

```
git status
git add *
git commit -m 'some comment'
git push
```

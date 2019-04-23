# vol2bird 0.4.0

All issues included in this release can be found [here](https://github.com/adokter/vol2bird/milestone/2?closed=1).

* vol2bird now reads Vaisala Sigmet IRIS (Iris RAW) format, the native radar data format of Vaisala radar processors (e.g. used in Canada, Portugal, Finland) (#112)

* added functionality to vol2bird to read files containing single scans (sweeps) and merge them into polar volumes (#116)

* rsl2odim now converts Vaisala IRIS format and ODIM hdf5 format and can merge files containing scans / sweeps into a polar volume. 

* new input argument format that allows specifying multiple input files

* change default maximum range from 25 km to 35 km (#117)

* fixed a bug that in rare cases produced a negative reflectivity eta (#123)

* new determineRadarFormat() function to check whether a file is in ODIM, IRIS or RSL format

* added a documentation webpage gnerated by doxygen (available at http://adokter.github.io/vol2bird)

* added this NEWS.md file to document releases

# vol2bird 0.3.20 and older

All issues included in this release can be found [here](https://github.com/adokter/bioRad/milestone/3?closed=1).


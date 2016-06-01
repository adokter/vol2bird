Currently test is not passing on h5 file in the current repository folder:
```
./openfile 201503222304_IST_MON_YAZ_C.h5
```
leads to
```
HDF5-DIAG: Error detected in HDF5 (1.8.13) thread 0:
  #000: H5A.c line 1094 in H5Aread(): unable to read attribute
    major: Attribute
    minor: Read failed
  #001: H5A.c line 1154 in H5A_read(): unable to convert between src and dst datatypes
    major: Attribute
    minor: Feature is unsupported
  #002: H5T.c line 4535 in H5T_path_find(): no appropriate function for conversion path
    major: Datatype
    minor: Unable to initialize object
Segmentation fault: 11
```

## vol2bird Python wrapper
compile using Makefile

run in Python as
```
import _pyvol2bird, _raveio
polarvolume = _raveio.open("test.h5").object
v2b=_pyvol2bird.new(polarvolume)
# change a user option
v2b.constants_nGatesCellMin = 10
# calculate a vertical profile of birds
vpr=a.vol2bird(polarvolume, 0.0)
```
:

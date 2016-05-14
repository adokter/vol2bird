## vol2bird Python wrapper
compile using Makefile

run in Python as
```
import _vol2bird, _raveio
polarvolume = _raveio.open("test.h5").object
a=_vol2bird.new(polarvolume)
# change a user option
a.constants_nGatesCellMin = 10
# calculate a vertical profile of birds
vpr=a.vol2bird()
```
:

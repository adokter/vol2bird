## vol2bird Python wrapper
compile using Makefile

run in Python as
```
import _pyvol2bird, _raveio
polarvolume = _raveio.open("test.h5").object
v2b=_pyvol2bird.new(polarvolume)
# change a user option (proof of concept for one option only, not implemented fully)
v2b.options_cellEtaMin = 1000
# calculate a vertical profile of birds
vpr=v2b.vol2bird(polarvolume)
# write the profile to file
ios = _raveio.new()
ios.object = vpr
ios.filename = "testout.h5"
ios.save()
```


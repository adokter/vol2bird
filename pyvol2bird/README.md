## vol2bird Python wrapper
compile using Makefile

run in Python as
```
import _pyvol2bird, _raveio
polarvolume = _raveio.open("test.h5").object
v2b=_pyvol2bird.new(polarvolume)
# change a user option (proof of concept for one option only, not implemented fully)
v2b.constants_nGatesCellMin = 10
# calculate a vertical profile of birds
vpr=v2b.vol2bird(polarvolume0)
# write the profile to file
ios = _raveio.new()
ios.object = vpr
ios.filename = "testout.h5"
ios.save()
```

Python knows about location of shared libraries either via PYTHONPATH variable, or via site-packages directory.
```
>>> import site
>>> site.getsitepackages()
/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages
```
In the site-packages directory files pyvol2bird.pth, rave.pth are installed, listing the location of shared libraries.

'''
Copyright (C) 2010- Swedish Meteorological and Hydrological Institute (SMHI)
Copyright 2016 Netherlands eScience Center

This file is part of RAVE.

RAVE is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

RAVE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with RAVE.  If not, see <http://www.gnu.org/licenses/>.

'''
## Plugin for generating the vol2bird product generation that is initiated from the beast
## framework.
## Register in pgf with
## --name=eu.baltrad.beast.vol2bird -m rave_pgf_vol2bird_plugin -f generate
##
## @file
## @author Anders Henja, SMHI
## @author Jurriaan Spaaks, NLeSC
## @author Lourens Veen, NLeSC
## @date 2016-06-13

import _pyvol2bird
import _rave
import _raveio
import rave_tempfile
import odim_source
import logging
import rave_pgf_logger

logger = rave_pgf_logger.create_logger()

ravebdb = None
try:
  import rave_bdb
  ravebdb = rave_bdb.rave_bdb()
except:
  pass

import rave_dom_db

## Creates a dictionary from a rave argument list
#@param arglist the argument list
#@return a dictionary
def arglist2dict(arglist):
  result={}
  for i in range(0, len(arglist), 2):
    result[arglist[i]] = arglist[i+1]
  return result

## Creates a composite
#@param files the list of files to be used for generating the composite
#@param arguments the arguments defining the composite
#@return a temporary h5 file with the composite
def generate(files, arguments):
  args = arglist2dict(arguments)

  for fname in files:
    logger.info("Input file: %s" % fname)
    obj = None
    if ravebdb != None:
      obj = ravebdb.get_rave_object(fname)
    else:
      rio = _raveio.open(fname)
      obj = rio.object

  v2b = _pyvol2bird.new(obj)

  # calculate a vertical profile of birds
  vpr = v2b.vol2bird(obj)

  fileno, outfile = rave_tempfile.mktemp(suffix='.h5', close="True")
  logger.info("Output file: %s" % outfile)
  ios = _raveio.new()
  ios.object = vpr
  ios.filename = outfile
  ios.save()

  return outfile

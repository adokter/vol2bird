###########################################################################
# Copyright (C) 2011 Swedish Meteorological and Hydrological Institute, SMHI,
#
# This file is part of beamb.
#
# beamb is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# beamb is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with beamb.  If not, see <http://www.gnu.org/licenses/>.
# ------------------------------------------------------------------------
# 
# Main build file
# @file
# @author Anders Henja (Swedish Meteorological and Hydrological Institute, SMHI)
# @date 2011-11-10
###########################################################################

.PHONY:all
all: build

def.mk:
	+[ -f $@ ] || $(error You need to run ./configure)

.PHONY:build 
build: def.mk
	$(MAKE) -C lib
	$(MAKE) -C pyvol2bird

.PHONY:install
install: def.mk
	$(MAKE) -C lib install
	$(MAKE) -C pyvol2bird install
	@echo "################################################################"
	@echo "To run the binaries you will need to setup your library path to"
	@echo "LD_LIBRARY_PATH="`cat def.mk | grep LD_PRINTOUT | sed -e"s/LD_PRINTOUT=//"`
	@echo "################################################################"

.PHONY:doc
doc:
	$(MAKE) -C doxygen doc

.PHONY:test
test: def.mk
	@chmod +x ./tools/test_vol2bird.sh
	@./tools/test_vol2bird.sh

.PHONY:clean
clean:
	$(MAKE) -C lib clean
	$(MAKE) -C pyvol2bird clean
	@\rm -f *~

	

.PHONY:distclean
distclean:
	$(MAKE) -C lib distclean
	$(MAKE) -C pyvol2bird distclean
	@\rm -f *~ config.log config.status def.mk

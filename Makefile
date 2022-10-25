###########################################################################
# Copyright (C) 2016 Adriaan Dokter & Netherlands eScience Center
#
# This file is part of vol2bird.
#
# vol2bird is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# vol2bird is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with beamb.  If not, see <http://www.gnu.org/licenses/>.
# ------------------------------------------------------------------------
# 
# Main build file
###########################################################################

.PHONY:all
all: build

def.mk:
	+[ -f $@ ] || $(error You need to run ./configure)

.PHONY:build 
build: def.mk
	if [ -f "./libmistnet/Makefile" ]; then \
		$(MAKE) -C libmistnet; \
	fi
	$(MAKE) -C lib
	$(MAKE) -C src
	#$(MAKE) -C pyvol2bird
	#$(MAKE) -C pgfplugin

.PHONY:install
install: def.mk
	$(MAKE) -C lib install
	$(MAKE) -C src install
	#$(MAKE) -C pyvol2bird install
	#$(MAKE) -C pgfplugin install
	@if [ -f "./libmistnet/Makefile" ]; then \
		$(MAKE) -C libmistnet install; \
	fi
	@echo "################################################################"
	@echo "To run the binaries you will need to setup your library path to"
	@echo "LD_LIBRARY_PATH="`cat def.mk | grep LD_PRINTOUT | sed -e"s/LD_PRINTOUT=//"`
	@echo "################################################################"

.PHONY:doc
doc:
	$(MAKE) -C doxygen doc

.PHONY:test
test: def.mk
	$(MAKE) -C tests test

.PHONY:clean
clean:
	$(MAKE) -C lib clean
	$(MAKE) -C src clean
	#$(MAKE) -C pyvol2bird clean
	$(MAKE) -C tests clean
	@if [ -f "./libmistnet/Makefile" ]; then \
		$(MAKE) -C libmistnet clean; \
	fi
	@\rm -f *~

.PHONY:distclean
distclean:
	$(MAKE) -C lib distclean
	$(MAKE) -C src distclean
	#$(MAKE) -C pyvol2bird distclean
	$(MAKE) -C tests distclean
	@\rm -rf libmistnet/CMakeFiles libmistnet/Makefile libmistnet/configure libmistnet/cmake_install.cmake
	@\rm -f libmistnet/CMakeCache.txt libmistnet/install_manifest.txt libmistnet/libmistnet.so
	@\rm -f *~ config.log config.status def.mk vol2bird.sh rsl2odim.sh

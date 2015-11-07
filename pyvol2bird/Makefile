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
# beamb make file
# @file
# @author Anders Henja (Swedish Meteorological and Hydrological Institute, SMHI)
# @date 2011-11-14
###########################################################################
-include ../def.mk

# specific c flags
#
#CFLAGS= $(SHARED_FLAG) $(CCOPTS) $(DEFS) -I../lib -I. $(PYTHON_INCLUDE_FLAG)\
#			$(NUMPY_INCLUDE_FLAG) $(RAVE_INCLUDE_FLAG) $(RAVE_INCLUDE_FLAG)/python \
#			$(HLHDF_INCLUDE_FLAG) $(ZLIB_INCDIR) $(HDF5_INCDIR) $(PROJ_INCLUDE_FLAG)
CFLAGS= -I../lib -I. $(PYTHON_INCLUDE_FLAG) $(RAVE_MODULE_PYCFLAGS) -I/usr/include

# Linker flags
#
LDFLAGS= -L../lib -L. $(RAVE_MODULE_LDFLAGS) $(PROJ_LIBRARY_FLAG) -L/usr/lib

# We need to ensure that raves expat requirements are managed
#EXPAT_LIBRARY=
#ifeq ($(EXPAT_SUPPRESSED), no)
#CPPFLAGS+= $(EXPAT_INCLUDE_DIR) -DRAVE_XML_SUPPORTED
#CFLAGS+= $(EXPAT_INCLUDE_DIR) -DRAVE_XML_SUPPORTED
#LDFLAGS+= $(EXPAT_LIB_DIR)
#EXPAT_LIBRARY=-lexpat
#endif

#LIBRARIES=-lhlhdf -lproj  $(EXPAT_LIBRARY) -lravetransform -lravepyapi -lvol2bird -lrt

# --------------------------------------------------------------------
# Fixed definitions

VOL2BIRD_SOURCE= pyvol2bird.c
VOL2BIRD_OBJECTS= $(VOL2BIRD_SOURCE:.c=.o)
VOL2BIRD_TARGET= _vol2bird.so

MAKECDEPEND=$(CC) -MM $(CFLAGS) -MT '$(@D)/$(@F)' -o $(DF).d $<

DEPDIR=.dep
DF=$(DEPDIR)/$(*F)

# Ensures that the .dep directory exists
.PHONY=$(DEPDIR)
$(DEPDIR):
	+@[ -d $@ ] || mkdir -p $@

# And the rest of the make file targets
#
.PHONY=all
all:		$(VOL2BIRD_TARGET)

$(VOL2BIRD_TARGET): $(DEPDIR) $(VOL2BIRD_OBJECTS) ../lib/libvol2bird.so
	$(LDSHARED) -o $@ $(VOL2BIRD_OBJECTS) $(LDFLAGS) $(RAVE_MODULE_PYLIBRARIES)
	
.PHONY=install
install:
	@mkdir -p ${prefix}/share/vol2bird/pyvol2bird
	@cp -v -f $(VOL2BIRD_TARGET) ${prefix}/share/vol2bird/pyvol2bird
	@-echo "$(prefix)/share/vol2bird/pyvol2bird" > "$(SITEPACK_PYTHON)/pyvol2bird.pth"

.PHONY=clean
clean:
	@\rm -f *.o core *~
	@\rm -fr $(DEPDIR)

.PHONY=distclean		 
distclean:	clean
	@\rm -f $(VOL2BIRD_TARGET)

# --------------------------------------------------------------------
# Rules

# Contains dependency generation as well, so if you are not using
# gcc, comment out everything until the $(CC) statement.
%.o : %.c
	@$(MAKECDEPEND); \
	cp $(DF).d $(DF).P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		-e '/^$$/ d' -e 's/$$/ :/' < $(DF).d >> $(DF).P; \
	\rm -f $(DF).d
	$(CC) -c $(CFLAGS) $< -o $@

# NOTE! This ensures that the dependencies are setup at the right time so this should not be moved
-include $(VOL2BIRD_SOURCE:%.c=$(DEPDIR)/%.P)


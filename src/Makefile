#
# Copyright 2013 Netherlands eScience Center
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# 
-include ../def.mk 
 
CC      ?= gcc
#CFLAGS  += -fPIC -x c -DFPRINTFON
CFLAGS  += -I../lib -Wall -g
LDFLAGS += -L../lib

# define where the vol2bird stuff is
#SRC_VOL2BIRD_DIR         = /projects/baltrad/vol2bird/lib
#SRC_VOL2BIRD_DIR         = /Users/adriaan/git/vol2bird/lib
SRC_VOL2BIRD_DIR         = ./

all : vol2bird rsl2odim

../lib/libvol2bird.so :
	make -C ../lib libvol2bird.so

VOL2BIRD_DEPS = vol2bird.c ../lib/libvol2bird.h ../lib/constants.h
RSL2ODIM_DEPS = rsl2odim.c ../lib/libvol2bird.h

vol2bird : ../lib/libvol2bird.so $(VOL2BIRD_DEPS)
	#
	# ------------------------------------
	#       making vol2bird 
	# ------------------------------------
	#
	$(CC) $(CFLAGS) vol2bird.c \
	-I. \
	$(RAVE_MODULE_CFLAGS) \
	$(RAVE_MODULE_LDFLAGS) \
	$(PROJ_LIBRARY_FLAG) \
	$(RSL_LIBRARY_FLAG) \
	-lvol2bird $(RAVE_MODULE_LIBRARIES) -lm $(RSL_LIB) -o vol2bird $(LDFLAGS)
	
rsl2odim : ../lib/libvol2bird.so $(RSL2ODIM_DEPS)
	#
	# ------------------------------------
	#       making rsl2odim 
	# ------------------------------------
	#
	$(CC) $(CFLAGS) rsl2odim.c \
	-I. \
	$(RAVE_MODULE_CFLAGS) \
	$(RAVE_MODULE_LDFLAGS) \
	$(PROJ_LIBRARY_FLAG) \
	$(RSL_LIBRARY_FLAG) \
	-lvol2bird $(RAVE_MODULE_LIBRARIES) -lm $(RSL_LIB) -o rsl2odim $(LDFLAGS)
	
	#
	# (You may still have to change your LD_LIBRARY_PATH)
	#

.PHONY : install
install : 
	# ------------------------------------
	#  copy binaries to target directory
	# ------------------------------------
	install -d ${prefix}/bin
	install -m 755 vol2bird ${prefix}/bin/vol2bird
	install -m 755 rsl2odim ${prefix}/bin/rsl2odim

.PHONY : clean
clean : 
	# ------------------------------------
	#  cleaning up old library and binary
	# ------------------------------------
	@if [ -f "./vol2bird" ]; then \
		\rm vol2bird; \
	fi
	@if [ -f "./rsl2odim" ]; then \
		\rm rsl2odim; \
	fi
	@\rm -f *~

.PHONY : distclean
distclean:
	@\rm -f *~
	@\rm vol2bird
	@\rm rsl2odim



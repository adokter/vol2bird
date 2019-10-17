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
CXX     = /usr/bin/c++   #note: has to be same compiler as used for Pytorch
#CFLAGS  += -fPIC -x c -DFPRINTFON
CFLAGS  += -I../lib -Wall -g $(RSL_CFLAG) $(IRIS_CFLAG) $(MISTNET_CFLAG)
LDFLAGS += -L../lib

# define where the vol2bird stuff is
#SRC_VOL2BIRD_DIR         = /projects/baltrad/vol2bird/lib
#SRC_VOL2BIRD_DIR         = /Users/adriaan/git/vol2bird/lib
SRC_VOL2BIRD_DIR         = ./

all :

ifeq ($(RSL_CFLAG),-DRSL)
all : vol2bird.o vol2bird rsl2odim.o rsl2odim
else
all : vol2bird.o vol2bird
endif

../lib/libvol2bird.so :
	make -C ../lib libvol2bird.so

VOL2BIRD_DEPS = vol2bird.c ../lib/libvol2bird.h ../lib/constants.h
RSL2ODIM_DEPS = rsl2odim.c ../lib/libvol2bird.h

vol2bird.o : vol2bird.c
	#
	# ------------------------------------
	#       making vol2bird.o
	# ------------------------------------
	#
	$(CC) -c $(CFLAGS) vol2bird.c \
	-I. \
	$(RAVE_MODULE_CFLAGS) \

rsl2odim.o : rsl2odim.c
	#
	# ------------------------------------
	#       making rsl2odim.o
	# ------------------------------------
	#
	$(CC) -c $(CFLAGS) rsl2odim.c \
	-I. \
	$(RAVE_MODULE_CFLAGS) \


vol2bird : ../lib/libvol2bird.so $(VOL2BIRD_DEPS)
	#
	# ------------------------------------
	#       linking vol2bird
	# ------------------------------------
	#
	$(CXX) -o vol2bird vol2bird.o \
	$(RAVE_MODULE_LDFLAGS) \
	$(PROJ_LIBRARY_FLAG) \
	$(RSL_LIBRARY_FLAG) \
	$(MISTNET_INCLUDE_FLAG) \
	$(MISTNET_LIBRARY_FLAG) \
	-lvol2bird $(RAVE_MODULE_LIBRARIES) -lm $(RSL_LIB) $(IRIS_LIB) $(LDFLAGS) $(MISTNET_LIB)
	
rsl2odim : ../lib/libvol2bird.so $(RSL2ODIM_DEPS)
	#
	# ------------------------------------
	#       making rsl2odim 
	# ------------------------------------
	#
	$(CXX) -o rsl2odim rsl2odim.o \
	$(RAVE_MODULE_LDFLAGS) \
	$(PROJ_LIBRARY_FLAG) \
	$(RSL_LIBRARY_FLAG) \
	$(IRIS_LIBRARY_FLAG) \
	$(MISTNET_INCLUDE_FLAG) \
	$(MISTNET_LIBRARY_FLAG) \
	-lvol2bird $(RAVE_MODULE_LIBRARIES) -lm $(RSL_LIB) $(IRIS_LIB) $(LDFLAGS) $(MISTNET_LIB)
	
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
	@if [ -f "./rsl2odim" ]; then \
	        \install -m 755 rsl2odim ${prefix}/bin/rsl2odim; \
	fi

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
	@if [ -f "./vol2bird.o" ]; then \
		\rm vol2bird.o; \
	fi
	@\rm -f *~

.PHONY : distclean
distclean:
	@if [ -f "./vol2bird" ]; then \
		\rm -f vol2bird; \
	fi
	@if [ -f "./rsl2odim" ]; then \
		\rm -f rsl2odim; \
	fi
	@if [ -f "./vol2bird.o" ]; then \
		\rm -f volbird.o; \
	fi
	@\rm -f *~

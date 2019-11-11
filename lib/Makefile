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
 
# use gcc if not set already
CC      = /usr/bin/gcc
#CFLAGS  += $(RSL_CFLAG) -fPIC -x c -DFPRINTFON
CFLAGS  += $(RSL_CFLAG) $(IRIS_CFLAG) -fPIC -g -x c 
LDFLAGS += -shared -Wall

# set linker flag on OSX
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	LDFLAGS += -undefined dynamic_lookup
endif

# define where the vol2bird stuff is
SRC_VOL2BIRD_DIR = .

all : libvol2bird.so

LIBVOL2BIRD_DEPS = librender.h constants.h libsvdfit.h libdealias.h librsl.h libvol2bird.h librender.c libsvdfit.c libdealias.c librsl.c libvol2bird.c

libvol2bird.so : $(LIBVOL2BIRD_DEPS)
	# ------------------------------------
	#        making libvol2bird.so        
	# ------------------------------------
	$(CC) $(CFLAGS) $(MISTNET_CFLAG) \
	$(RAVE_MODULE_CFLAGS) \
	-I. \
	$(CONFUSE_INCLUDE_FLAG) \
	$(GSL_INCLUDE_FLAG) \
	$(RSL_INCLUDE_FLAG) \
	$(RAVE_MODULE_LDFLAGS) \
	$(PROJ_LIBRARY_FLAG) \
	$(CONFUSE_LIBRARY_FLAG) \
	$(GSL_LIBRARY_FLAG) \
	$(RSL_LIBRARY_FLAG) \
	$(SRC_VOL2BIRD_DIR)/libvol2bird.c \
	$(SRC_VOL2BIRD_DIR)/libsvdfit.c \
	$(SRC_VOL2BIRD_DIR)/libdealias.c \
	$(SRC_VOL2BIRD_DIR)/librsl.c \
	$(SRC_VOL2BIRD_DIR)/librender.c \
	$(LDFLAGS) \
	-Wall -o libvol2bird.so $(RAVE_MODULE_LIBRARIES) -lconfuse -lgsl -lgslcblas $(RSL_LIB) $(IRIS_LIB)

.PHONY : install
install : libvol2bird.so
	# ------------------------------------
	#  copy binaries to target directory
	# ------------------------------------
	install -d ${prefix}/lib
	install -m 644 libvol2bird.so ${prefix}/lib/libvol2bird.so

.PHONY : clean
clean : 
	# ------------------------------------
	#  cleaning up old library and binary
	# ------------------------------------
	@if [ -f "./libvol2bird.so" ]; then \
		\rm libvol2bird.so; \
	fi
	@\rm -f *~

.PHONY : distclean
distclean:
	@\rm -f libvol2bird.so
	@\rm -f *~

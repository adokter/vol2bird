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
 
CC      = gcc
#CFLAGS  = -fPIC -x c -DFPRINTFON
CFLAGS  = -fPIC -x c -I../lib
LDFLAGS = -shared -Wall -L../lib

# define where the vol2bird stuff is
#SRC_VOL2BIRD_DIR         = /projects/baltrad/vol2bird/lib
#SRC_VOL2BIRD_DIR         = /Users/adriaan/git/vol2bird/lib
SRC_VOL2BIRD_DIR         = ./

all : vol2bird

libvol2bird :
	make -C ../lib libvol2bird 

vol2bird : libvol2bird
	#
	#
	#
	#
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
	-lvol2bird $(RAVE_MODULE_LIBRARIES) -lm -o vol2bird $(LDFLAGS)
	
	#
	# (You may still have to change your LD_LIBRARY_PATH)
	#

install : vol2bird
	# ------------------------------------
	#  copy binaries to target directory
	# ------------------------------------
	install -d ${prefix}/bin
	install -m 755 vol2bird ${prefix}/bin/vol2bird

clean : 
	# ------------------------------------
	#  cleaning up old library and binary
	# ------------------------------------
	@if [ -f "./vol2bird" ]; then \
		\rm vol2bird; \
	fi
	@\rm -f *~

distclean:
	@\rm -f *~
	@\rm vol2bird



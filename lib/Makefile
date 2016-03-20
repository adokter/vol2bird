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
CFLAGS  = -fPIC -x c
LDFLAGS = -shared -Wall

# define where the vol2bird stuff is
#SRC_VOL2BIRD_DIR         = /projects/baltrad/vol2bird/lib
#SRC_VOL2BIRD_DIR         = /Users/adriaan/git/vol2bird/lib
SRC_VOL2BIRD_DIR         = ./

all : libvol2bird
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
	-L. \
	-L$(PWD) \
	-lvol2bird $(RAVE_MODULE_LIBRARIES) -lm -o vol2bird
	
	#
	# (You may still have to change your LD_LIBRARY_PATH)
	#

libvol2bird :
	# ------------------------------------
	#        making libvol2bird.so        
	# ------------------------------------
	$(CC) $(CFLAGS) \
	-I. \
	$(RAVE_MODULE_CFLAGS) \
	$(CONFUSE_INCLUDE_FLAG) \
	$(RAVE_MODULE_LDFLAGS) \
	$(PROJ_LIBRARY_FLAG) \
	$(CONFUSE_LIBRARY_FLAG) \
	$(SRC_VOL2BIRD_DIR)/libvol2bird.c \
	$(SRC_VOL2BIRD_DIR)/libsvdfit.c \
	$(LDFLAGS) \
	-Wall -o libvol2bird.so $(RAVE_MODULE_LIBRARIES) -lconfuse -lm


clean : 
	# ------------------------------------
	#  cleaning up old library and binary
	# ------------------------------------
	@if [ -f "./libvol2bird.so" ]; then \
		\rm libvol2bird.so; \
	fi
	@if [ -f "./vol2bird" ]; then \
		\rm vol2bird; \
	fi
	@\rm -f *~

distclean:
	@\rm -f libvol2bird.so
	@\rm -f *~
	@\rm vol2bird



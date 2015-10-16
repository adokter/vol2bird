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
 
 
# directory where the baltrad software was installed
BALTRAD_PREFIX          = /opt/baltrad

CC      = gcc
#CFLAGS  = -fPIC -x c -DFPRINTFON
CFLAGS  = -fPIC -x c
LDFLAGS = -shared -Wl,--no-undefined


# directories containing include files
INCLUDE_RAVE_DIR        = $(BALTRAD_PREFIX)/rave/include
INCLUDE_THIRD_PARTY_DIR = $(BALTRAD_PREFIX)/third_party/include
INCLUDE_CONFUSE         = /usr/include

# note: you have to add these LIB_*_DIRs to the LD_LIBRARY_PATH before 
# you can run the binary

LIB_RAVE_DIR            = $(BALTRAD_PREFIX)/rave/lib
LIB_HLHDF_DIR           = $(BALTRAD_PREFIX)/hlhdf/lib
LIB_PROJ_EXPAT_DIR      = $(BALTRAD_PREFIX)/third_party/lib
LIB_BBUFR_DIR           = $(BALTRAD_PREFIX)/bbufr/lib
LIB_CONFUSE             = /usr/lib/x86_64-linux-gnu

# define where the vol2bird stuff is
SRC_VOL2BIRD_DIR        = $(HOME)/github/enram/ncradar/lib/vol2bird




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
	-I$(INCLUDE_RAVE_DIR) \
	-I$(INCLUDE_THIRD_PARTY_DIR) \
	-I$(SRC_VOL2BIRD_DIR) \
	-L$(LIB_RAVE_DIR) \
	-L$(LIB_HLHDF_DIR) \
	-L$(LIB_PROJ_EXPAT_DIR) \
	-L$(LIB_BBUFR_DIR) \
	-L$(PWD) \
	-lravetoolbox -lhdf5 -lhlhdf -lproj -lexpat -lOperaBufr -lm -lvol2bird -o vol2bird
	
	#
	# (You may still have to change your LD_LIBRARY_PATH)
	#


libvol2bird :
	# ------------------------------------
	#        making libvol2bird.so        
	# ------------------------------------
	$(CC) $(CFLAGS) \
	-I$(INCLUDE_RAVE_DIR) \
	-I$(INCLUDE_THIRD_PARTY_DIR) \
	-I$(INCLUDE_CONFUSE) \
	-I$(SRC_VOL2BIRD_DIR) \
	-L$(LIB_RAVE_DIR) \
	-L$(LIB_PROJ_EXPAT_DIR) \
	-L$(LIB_CONFUSE) \
	$(SRC_VOL2BIRD_DIR)/libvol2bird.c \
	$(SRC_VOL2BIRD_DIR)/libsvdfit.c \
	-Wall -o libvol2bird.so -lravetoolbox -lconfuse -lm $(LDFLAGS)


clean : 
	# ------------------------------------
	#  cleaning up old library and binary
	# ------------------------------------
	if [ -f "./libvol2bird.so" ]; then \
		rm libvol2bird.so; \
	fi
	if [ -f "./vol2bird" ]; then \
		rm vol2bird; \
	fi





# Copyright 2019 Adriaan Dokter
#
-include ../def.mk 
 
# use gcc if not set already
CXX      = /usr/bin/c++
CXXFLAGS  +=  -std=gnu++11 
CFLAGS =  -fPIC -g
LDFLAGS += -shared -Wall

# FIXME these should be set in def.mk
MISTNET_INCLUDE_FLAG = -I/Users/amd427/Documents/radar/mistnet/libtorch/include
MISTNET_LIBRARY_FLAG = -L/Users/amd427/Documents/radar/mistnet/libtorch/lib/
MISTNET_LIB = -ltorch -lc10 -lcaffe2

all : librender.o libmistnet.so

LIBMISTNET_DEPS = libmistnet.h libmistnet.cpp librender.h librender.c

librender.o : librender.h librender.c
	# ------------------------------------
	#        making librender.o        
	# ------------------------------------
	$(CC) $(CFLAGS) librender.c \
	-I../lib \
	$(RAVE_MODULE_CFLAGS) \
	$(MISTNET_INCLUDE_FLAG) \
	 -o librender.o -c

libmistnet.so : $(LIBMISTNET_DEPS)
	# ------------------------------------
	#        making libmistnet.so        
	# ------------------------------------
	$(CXX) $(CXXFLAGS) libmistnet.cpp librender.o \
	-I../lib \
	-L../lib \
	$(MISTNET_INCLUDE_FLAG) \
	$(MISTNET_LIBRARY_FLAG) \
	$(RAVE_MODULE_LDFLAGS) \
	$(LDFLAGS) \
	-Wall -o libmistnet.so $(RAVE_MODULE_LIBRARIES)  $(MISTNET_LIB) -lm

.PHONY : install
install : libmistnet.so
	# ------------------------------------
	#  copy binaries to target directory
	# ------------------------------------
	install -d ${prefix}/lib
	install -m 644 libmistnet.so ${prefix}/lib/libmistnet.so

.PHONY : clean
clean : 
	# ------------------------------------
	#  cleaning up old library
	# ------------------------------------
	@if [ -f "./libmistnet.so" ]; then \
		\rm libmistnet.so; \
	fi
	@if [ -f "./librender.o" ]; then \
		\rm librender.o; \
	fi
	@\rm -f *~

.PHONY : distclean
distclean:
	@\rm -f libmistnet.so librender.o
	@\rm -f *~

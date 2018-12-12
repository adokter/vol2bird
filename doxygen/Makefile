############################################################
# Makefile for generating the documentation
# 
# Author(s):	Adriaan Dokter
#
# Date:		2018-12-11
############################################################
all:

.PHONY=doc
doc:
	@doxygen vol2bird_doxygen.cfg

.PHONY=clean
clean:
	@\rm -f *.o
	@\rm -f *~

.PHONY=distclean
distclean:
	@\rm -f *.o
	@\rm -f *~ core
	@\rm -fr doxygen

distribution:

install:


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
include ../def.mk 

print-% : ; @echo $* = $($*)

PREFIX = /opt/baltrad

test : HLHDFLIB = $(PREFIX)/hlhdf/lib
test : RAVELIB = $(PREFIX)/rave/lib
test : PYTHONLIB = $(PREFIX)/third_party/lib
test : RAVEPYTHON = $(PREFIX)/third_party/bin/python

test : VOL2BIRDLIB = ../lib
test : fixtures
	#
	#
	#
	# ------------------------------------
	#       running tests
	# ------------------------------------
	#
	LD_LIBRARY_PATH=$(PYTHONLIB):$(RAVELIB):$(HLHDFLIB):$(VOL2BIRDLIB) $(RAVEPYTHON) -m unittest discover


fixtures :
	git clone https://github.com/adokter/ODIM-hdf5-test fixtures

.PHONY: clean
clean:
	@\rm -rf fixtures
	@\rm -f *.h5
	@\rm -f *.pyc

.PHONY: distclean
distclean: clean


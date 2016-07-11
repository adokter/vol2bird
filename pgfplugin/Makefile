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
 
.PHONY : all
all :

.PHONY : install
install :
	# ------------------------------------
	#  copy binaries to target directory
	# ------------------------------------
	install -d ${prefix}/share/vol2bird
	install -m 644 rave_pgf_vol2bird_plugin.py ${prefix}/share/vol2bird
	install -d ${prefix}/etc
	@-echo "$(prefix)/share/vol2bird" > "${prefix}/etc/vol2bird-pgfplugin.pth"

.PHONY : clean
clean : 
	# ------------------------------------
	#  cleaning up old library and binary
	# ------------------------------------
	@\rm -f *.pyc
	@\rm -f *~

.PHONY : distclean
distclean:
	@\rm -f *.pyc
	@\rm -f *~



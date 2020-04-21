# Copyright (C) 2011-2012 Quantum ESPRESSO Foundation
#
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#

default : build_elpa 

all: build_elpa

build_elpa:
	if test -d src ; then \
	( cd src ; if test "$(MAKE)" = "" ; then make $(MFLAGS) ; \
	else $(MAKE) $(MFLAGS) ; fi ) ; fi

build_elpa_dev:
	if test -d devel ; then \
	( cd src ; if test "$(MAKE)" = "" ; then make $(MFLAGS) ; \
	else $(MAKE) $(MFLAGS) ; fi ) ; fi
	
clean:
	if test -d src ; then \
	( cd src ; if test "$(MAKE)" = "" ; then make $(MFLAGS) clean ; \
	else $(MAKE) $(MFLAGS) clean ; fi ) ; fi
	if test -d devel ; then \
	( cd devel ; if test "$(MAKE)" = "" ; then make $(MFLAGS) clean ; \
	else $(MAKE) $(MFLAGS) clean ; fi ) ; fi
	rm -f libelpa.a

veryclean: clean
	rm -f libelpa.a

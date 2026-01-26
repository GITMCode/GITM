#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission
#  For more information, see http://csem.engin.umich.edu/tools/swmf

SHELL =/bin/sh

include ../Makefile.def
include ../Makefile.conf
-include Makefile.DEPEND

MODULES = \
	ModCharSize.o\
	ModIE.o\
	ModIEDummy.o\
	ModIECoupling.o\

OBJECTS = \
	tools.o

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${MODULES}

#
#	Making executables and libraries
#
MY_LIB = libUAIE.a

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${LIBPREV} ${MF} ${MODULES} ${OBJECTS}
	cp -f ${LIBPREV} ${MY_LIB}
	${AR} ${MY_LIB} ${MODULES} ${OBJECTS}

clean: cleanfiles

distclean: clean
	rm -f Makefile.DEPEND

# keep this line
include ../srcMake/Makefile.conf
include ../Makefile.def
-include ./Makefile.DEPEND

#
#  Makefile for AB2D code:
#

MODULES = \
	AB_module.o\
	AB_XFER_module.o\
	AB_ERROR_module.o\
	AB_SPH_module.o\
	AB_COMM_module.o\
	AB_ARRAY_module.o

OBJECTS = \
	AB_XFER_array_util.o\
	AB_XFER_1blk_util.o

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl  ./ ${SEARCH} ${LIBDIR} ${OBJECTS} ${MODULES}


MY_LIB = ${LIBDIR}/libSphere.a

LIB:	DEPEND
	${MAKE} ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${MF} ${MODULES} ${OBJECTS}
	@echo 'Creating libSphere library'
	${AR} ${MY_LIB} ${MODULES} ${OBJECTS}

clean:
	rm -f *~ core *.o *.mod fort.* a.out *.exe *.a *.so *.protex

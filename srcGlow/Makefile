include ../Makefile.def
include ../srcMake/Makefile.conf
-include ./Makefile.DEPEND
#include ./Makefile.RULES

MODULES = Mod_GLOW.o\
	ModCIMPIT.o\
	ModCXPARS.o

OBJECTS = \
	egrid.o\
	ephoto.o\
	etrans.o\
	exsect.o\
	fieldm.o\
	gchem.o\
	geomag.o\
	glow.o\
	maxt.o\
	qback.o\
	rcolum.o\
	rout.o\
	solzen.o\
	ssflux.o\
	vquart.o\
	GL_library.o

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ./ ${SEARCH} ${LIBDIR} ${OBJECTS} ${MODULES}

EXE = ${BINDIR}/GLOW.exe

MYLIB = libGLOW.a

LIB:	DEPEND
	${MAKE} ${MYLIB}
	@echo
	@echo ${MYLIB} has been brought up to date.
	@echo

${MYLIB}: ${MODULES} ${OBJECTS}
	@echo 'Creating GLOW library'
	${AR} ${MYLIB} ${MODULES} ${OBJECTS}

MY_LIB=libUPTOGL.a

LIBADD: DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${MODULES} ${OBJECTS}
	${MAKE} ${MYLIB}
	# cp -f ${LIBPREV} ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS} ${MODULES}



GLOW: DEPEND
	 @make ${EXE}

glow: 	${OBJECTS} ${MODULES} 
	${COMPILE.f77} ${Cflag3} -o ${EXE} main.o ${OBJECTS} ${MODULES}

clean:
	rm -f *~ core *.o *.mod fort.* a.out *.exe *.a *.so *.protex

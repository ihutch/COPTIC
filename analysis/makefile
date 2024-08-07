COPTIC=coptic
########################################################################
# The root directory of the distribution (for now)
COPTIC_ROOT=..
SRC=src
#########################################################################
ACCISLIB=${COPTIC_ROOT}/accis/libaccisX.a
ifeq ("$(LIBPATH)","")
   LIBPATH:=-L/usr/X11R6/lib/ -L/usr/lib/mesa
endif
LIBRARIES = $(LIBPATH) -L${COPTIC_ROOT}/accis/ -laccisX -lX11 $(GLULIBS) -L. -lvtkwrite -lfftpack5.1 -lfftw3f -lanalcoptic
#Default accis driver choice
#Alternatives are vecx or vecglx, can be overriden by commandline option
VECX=vecglx
#    VECX=vecx
ifeq ("$(VECX)","vecglx")     
   GLULIBS= -lGL -lGLU
endif
######################################################################
# An option setting might override default compiler.
ifeq ("$(G77)","")
	override G77=$(shell cat ${COPTIC_ROOT}/${SRC}/compiler 2>/dev/null)
endif
# export this so it is inherited by sub-makes.
export G77
GFINAL=gcc-4.1 -v -pg -o $(COPTIC).prof $(COPTIC).o $(OBJECTS) -static-libgcc -lpthread_p -lm_p -lc -lg2c -lmpich -lrt -lfrtbegin  $(LIBRARIES)
#OPTIMIZE=-O3 -funroll-loops -finline-functions
OPTIMIZE=
COMPILE-SWITCHES=-I.
ifeq ("$(G77)","mpif77")	
  COMPILE-SWITCHES = -Wall $(OPTIMIZE) -I. -mcmodel=medium
  COMPILE-SWITCHES:=$(COMPILE-SWITCHES)
#-fallow-argument-mismatch this broke on engaging 4 Dec 2023.
  OPTIMIZE=-O3
endif	
COMPILE-SWITCHES := $(COMPILE-SWITCHES) -g -fbounds-check
##COMPILE-SWITCHES = -Wall -Wno-unused $(OPTIMIZE) -g -I.
NOBOUNDS= $(COMPILE-SWITCHES) -fno-bounds-check
NOGLOBALS= $(COMPILE-SWITCHES) -Wno-globals
HEADERS=${COPTIC_ROOT}/${SRC}/ptaccom.f
##########################################################################
TARGETS=partexamine denexamine phiexamine diagexamine fluxexamine philineread phioverlay phasespace phasereadplot

all: $(TARGETS)

# Work around the inability to have multiple definitions in modern gfortran
libcoptic.a : ${COPTIC_ROOT}/${SRC}/libcoptic.a
	echo Copying libcoptic.a from coptic root
	cp ${COPTIC_ROOT}/${SRC}/libcoptic.a libcoptic.a
	ar d libcoptic.a reduce.o

# We presume that all these routines are serial not parallel.
dummyreduce.o : ${COPTIC_ROOT}/${SRC}/dummyreduce.f
	$(G77) -c ${COPTIC_ROOT}/${SRC}/dummyreduce.f

% : %.f  makefile modesaving.f examdecl.f $(ACCISLIB) dummyreduce.o ${COPTIC_ROOT}/${SRC}/libcoptic.a $(HEADERS) libvtkwrite.a libfftpack5.1.a libanalcoptic.a libcoptic.a
	$(G77)  -o $* $(COMPILE-SWITCHES) $(PROFILING) $*.f dummyreduce.o libcoptic.a $(LIBRARIES)

% : %.f90  makefile modesaving.f examdecl.f $(ACCISLIB) dummyreduce.o libcoptic.a $(HEADERS) libvtkwrite.a libfftpack5.1.a libanalcoptic.a libcoptic.a
	$(G77)  -o $* $(COMPILE-SWITCHES) $(PROFILING) $*.f90 dummyreduce.o ${COPTIC_ROOT}/${SRC}/coptic.a $(LIBRARIES)

# Don't recompile accis every time the makefile is changed.
${COPTIC_ROOT}/accis/%.o : ${COPTIC_ROOT}/accis/%.f $(HEADERS)
	$(G77)  -c $(COMPILE-SWITCHES) $(PROFILING) $*.f

# Just putting the specials first ensures that the compile works.
%.o : %.f makefile $(HEADERS)
	$(G77)  -c $(COMPILE-SWITCHES) $(PROFILING) $*.f

%.o : %.f90 makefile 
	gfortran -c $(COMPILE-SWITCHES) $(PROFILING) $*.f90

%.o : %.c makefile
	cc -c $(PROFILING) $*.c -I.

${COPTIC_ROOT}/libcoptic.a : 
	make -C ${COPTIC_ROOT}/ libcoptic.a

libvtkwrite.a : vtkwritewrap.o visit_writer.o
	@ar -rs libvtkwrite.a vtkwritewrap.o visit_writer.o

libfftpack5.1.a : fftpack5.1.f
	$(G77) $(COMPILE-SWITCHES) -c -w -fallow-argument-mismatch fftpack5.1.f
	ar cr libfftpack5.1.a fftpack5.1.o
	rm -f fftpack5.1.o

libanalcoptic.a : modesaving.o diststat.o
	@ar -rs libanalcoptic.a modesaving.o diststat.o

modeplot : modeplot.f fitexponential.o modesaving.f examdecl.f $(ACCISLIB)
	$(G77)  -o modeplot $(COMPILE-SWITCHES) $(PROFILING) modeplot.f dummyreduce.o fitexponential.o ${COPTIC_ROOT}/${SRC}/libcoptic.a $(LIBRARIES)

modeshape : modeshape.f fitexponential.o modesaving.f examdecl.f $(ACCISLIB)
	$(G77)  -o modeshape $(COMPILE-SWITCHES) $(PROFILING) modeshape.f dummyreduce.o fitexponential.o ${COPTIC_ROOT}/${SRC}/libcoptic.a $(LIBRARIES)

########################################################
$(ACCISLIB) : ${COPTIC_ROOT}/accis/*.f ${COPTIC_ROOT}/accis/*.c ${COPTIC_ROOT}/accis/*.h
	make -C ${COPTIC_ROOT}/accis

######################################################
vecx :
	make clean
	make VECX=vecx -C ${COPTIC_ROOT}/accis
	make

#####################################################
clean :
	rm -f *.o $(TARGETS) *.html *.flx *.ph? *.den T*.0?? *.ps *.dat *.a


COPTIC=coptic
#########################################################################
ACCISLIB=../accis/libaccisX.a
LIBRARIES = -L/usr/X11R6/lib/ -L../accis/ -laccisX -lXt -lX11 $(GLULIBS)
#Default accis driver choice
#Alternatives are vecx or vecglx, can be overriden by commandline option
VECX=vecglx
#    VECX=vecx
ifeq ("$(VECX)","vecglx")     
   GLULIBS= -lGL -lGLU
endif
######################################################################
# An option setting might override default compiler.
G77=mpif77 -f77=g77 
#G77=mpif77
# For loki:
#G77=/opt/mvapich2.gbe/path/bin/mpif90
# export this so it is inherited by sub-makes.
export G77
GFINAL=gcc-4.1 -v -pg -o $(COPTIC).prof $(COPTIC).o $(OBJECTS) -static-libgcc -lpthread_p -lm_p -lc -lg2c -lmpich -lrt -lfrtbegin  $(LIBRARIES)
#OPTIMIZE=-O3 -funroll-loops -finline-functions
OPTIMIZE=-O3
COMPILE-SWITCHES = -Wall  $(OPTIMIZE)  -I. 
#COMPILE-SWITCHES = -Wall   $(OPTIMIZE) -I. -g -fbounds-check
##COMPILE-SWITCHES = -Wall -Wno-unused $(OPTIMIZE) -g -I.
NOBOUNDS= $(COMPILE-SWITCHES) -fno-bounds-check
NOGLOBALS= $(COMPILE-SWITCHES) -Wno-globals
HEADERS=../ptaccom.f
##########################################################################
TARGETS=partexamine denexamine phiexamine diagexamine fluxexamine philineread phioverlay


all: $(TARGETS)

# We presume that all these routines are serial not parallel.
dummyreduce.o : ../dummyreduce.f
	$(G77) -c ../dummyreduce.f

% : %.f  makefile $(ACCISLIB) dummyreduce.o ../libcoptic.a $(HEADERS)
	make -C ../ libcoptic.a
	$(G77)  -o $* $(COMPILE-SWITCHES) $(PROFILING) $*.f dummyreduce.o ../libcoptic.a $(LIBRARIES)

# Don't recompile accis every time the makefile is changed.
../accis/%.o : ../accis/%.f $(HEADERS)
	$(G77)  -c $(COMPILE-SWITCHES) $(PROFILING) $*.f

# Just putting the specials first ensures that the compile works.
%.o : %.f makefile $(HEADERS)
	$(G77)  -c $(COMPILE-SWITCHES) $(PROFILING) $*.f

%.o : %.c makefile
	cc -c $(PROFILING) $*.c

#../libcoptic.a : 
#	make -C ../ libcoptic.a

########################################################
$(ACCISLIB) : ../accis/*.f ../accis/*.c ../accis/*.h
	make -C ../accis

######################################################
vecx :
	make clean
	make VECX=vecx -C ../accis
	make

#####################################################
clean :
	rm -f *.o $(TARGETS) *.html *.flx *.ph? *.den T*.0?? *.ps *.dat


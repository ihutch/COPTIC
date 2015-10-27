COPTIC=coptic
########################################################################
# The root directory of the distribution (for now)
COPTIC_ROOT=..
SRC=src
#########################################################################
#########################################################################
ACCISLIB=${COPTIC_ROOT}/accis/libaccisX.a
LIBRARIES = -L/usr/X11R6/lib/ -L/usr/lib/mesa -L${COPTIC_ROOT}/accis/ -laccisX -lXt -lX11 $(GLULIBS) -L.
#Default accis driver choice
#Alternatives are vecx or vecglx, can be overriden by commandline option
VECX=vecglx
#    VECX=vecx
ifeq ("$(VECX)","vecglx")     
   GLULIBS= -lGL -lGLU
endif
######################################################################
# An option setting might override default compiler.
G77=mpif77
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
HEADERS=$(COPTIC_ROOT)/$(SRC)/ptaccom.f
##########################################################################
TARGETS=creintest mditeratetest mpibbdytest mditeratortest stresstest sortest 

all: $(TARGETS)

# We used to presume that all these routines are serial not parallel.
# But now we are using the full libcoptic without this replacement.
dummyreduce.o : $(COPTIC_ROOT)/$(SRC)/dummyreduce.f
	$(G77) -c $(COPTIC_ROOT)/$(SRC)/dummyreduce.f

% : %.f  makefile $(ACCISLIB) $(COPTIC_ROOT)/$(SRC)/libcoptic.a $(HEADERS)
# I don't remember why this explicit form was preferred to a dependency below.
#	make -C ../ libcoptic.a
	$(G77)  -o $* $(COMPILE-SWITCHES) $(PROFILING) $*.f $(COPTIC_ROOT)/$(SRC)/libcoptic.a $(LIBRARIES)

# Don't recompile accis every time the makefile is changed.
$(COPTIC_ROOT)/accis/%.o : $(COPTIC_ROOT)/accis/%.f $(HEADERS)
	$(G77)  -c $(COMPILE-SWITCHES) $(PROFILING) $*.f

# Just putting the specials first ensures that the compile works.
%.o : %.f makefile $(HEADERS)
	$(G77)  -c $(COMPILE-SWITCHES) $(PROFILING) $*.f

%.o : %.c makefile
	cc -c $(PROFILING) $*.c -I.

$(COPTIC_ROOT)/$(SRC)/libcoptic.a : 
	make -C $(COPTIC_ROOT)/$(SRC)/ libcoptic.a

########################################################
$(ACCISLIB) : $(COPTIC_ROOT)/accis/*.f $(COPTIC_ROOT)/accis/*.c $(COPTIC_ROOT)/accis/*.h
	make -C $(COPTIC_ROOT)/accis

######################################################
vecx :
	make clean
	make VECX=vecx -C $(COPTIC_ROOT)/accis
	make

#####################################################
clean :
	rm -f *.o $(TARGETS) *.html *.flx *.ph? *.den T*.0?? *.ps *.a 

testing :
#	echo $(TARGETS)
# Execute each target program with 1 second between each.
	$(patsubst %,./%;sleep 1;,$(TARGETS))

COPTIC=coptic
#########################################################################
ACCISLIB=./accis/libaccisX.a
LIBRARIES = -L/usr/X11R6/lib/ -L./accis/ -laccisX -lXt -lX11 $(GLULIBS)
#Default accis driver choice
#Alternatives are vecx or vecglx, can be overriden by commandline option
VECX=vecglx
#    VECX=vecx
ifeq ("$(VECX)","vecglx")     
   GLULIBS= -lGL -lGLU
endif
##########################################################################
# The reinjection choice:
#######################
#REINJECT=orbitinjnew.o extint.o
#GEOMFILE=geomsphere.dat
##################
#REINJECT=reinject.o
#GEOMFILE=geomsphere.dat
###################
REINJECT=cartreinject.o
GEOMFILE=geomcubic.dat
#GEOMFILE=geomz200x25.dat
##########################################################################
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
##########################################################################
FIXEDOBJECTS=sormpi.o sorrelaxgen.o mpibbdy.o  cijroutine.o cijplot.o 3dobjects.o mditerate.o svdsol.o padvnc.o chargetomesh.o slicesect.o randf.o reindiag.o pinit.o ccpicplot.o volint.o fluxdata.o stringsnames.o meshconstruct.o partwriteread.o checkcode.o stress.o average.o randc.o
#
SPECIALOBJECTS=bdyroutine.o reduce.o getfield.o interpolations.o 
# Things just needed for the test routine:
UTILITIES=udisplay.o
REGULAROBJECTS= $(FIXEDOBJECTS) ${REINJECT}
OBJECTS=$(SPECIALOBJECTS) $(REGULAROBJECTS)
HEADERS=bbdydecl.f meshcom.f objcom.f 3dcom.f partcom.f rancom.f ran1com.f creincom.f ptaccom.f colncom.f examdecl.f griddecl.f ptchcom.f mditcom.f
TARGETS=mpibbdytest mditeratetest sormpitest fieldtest
##########################################################################
# If this rule does not seem to recognize the file you are trying to make,
# then run 'make' to completion first. It is something to do with the
# match-anything rules and prerequisites. I think that the rule is being
# interpreted as "terminal" which means it does not apply unless its
# prerequisites exist. By my interpretation of the info, this ought not
# to be happening (requires a double colon), but on horace, it is. 
#% : %.f  makefile $(OBJECTS) $(UTILITIES) $(ACCISLIB)
#	$(G77)  -o $* $(COMPILE-SWITCHES) $(PROFILING) $*.f $(OBJECTS) $(UTILITIES) $(LIBRARIES)

% : %.f  makefile libcoptic.a $(ACCISLIB)
	$(G77)  -o $* $(COMPILE-SWITCHES) $(PROFILING) $*.f libcoptic.a $(LIBRARIES)

# Don't recompile accis every time the makefile is changed.
./accis/%.o : ./accis/%.f $(HEADERS)
	$(G77)  -c $(COMPILE-SWITCHES) $(PROFILING) $*.f

# Just putting the specials first ensures that the compile works.
%.o : %.f makefile $(HEADERS)
	$(G77)  -c $(COMPILE-SWITCHES) $(PROFILING) $*.f

%.o : %.c makefile
	cc -c $(PROFILING) $*.c

#default target
# Problem when using geometry that can't do smt check. 
smt.out : $(COPTIC) ccpicgeom.dat
	if [ -f "$(GEOMFILE)" ] ; then ln -s -f $(GEOMFILE) ccpicgeom.dat ; fi
	@if [ -f smt.out ] ; then mv smt.out smt.prev ; fi
	./$(COPTIC)
	@if [ -f smt.prev ] ;then if [ -f smt.out ] ;then diff smt.prev smt.out ;else touch smt.out ;fi ;fi

libcoptic.a : makefile $(OBJECTS) $(UTILITIES)
	ar -rs libcoptic.a $(OBJECTS) $(UTILITIES)

#mpi checking target
mpicheck : $(COPTIC)
	mpiexec -l -n 2 ./$(COPTIC) >mpicheck.out
	if [ -f mpicheck.prev ] ; then diff mpicheck.prev mpicheck.out ; else mv mpicheck.out mpicheck.prev  ; fi
	mv mpicheck.out mpicheck.prev

#####################################################
# Things to compile with non-standard switches
# We make one of these the first thing in objects to force the header
# dependence to be reported, not just ignored by make on pattern rule.
interpolations.o : interpolations.f makefile $(HEADERS)
	$(G77)  -c $(NOBOUNDS) $(PROFILING) $*.f

getfield.o : getfield.f makefile $(HEADERS)
	$(G77)  -c $(NOBOUNDS) $(PROFILING) $*.f

reduce.o : reduce.f makefile $(HEADERS)
	$(G77) -c $(NOGLOBALS) $(PROFILING) $*.f

bdyroutine.o : bdyroutine.f makefile $(HEADERS)
	$(G77) -c $(NOGLOBALS) $(PROFILING) $*.f

#####################################################
# Main program explicit to avoid make bugs:
$(COPTIC) : $(COPTIC).f makefile $(ACCISLIB) $(OBJECTS) $(UTILITIES) libcoptic.a
	@echo "      rjscheme="\'$(REINJECT)\'" " > REINJECT.f
	$(G77)  -o $(COPTIC) $(COMPILE-SWITCHES) $(PROFILING) $(COPTIC).f  libcoptic.a $(LIBRARIES)

$(ACCISLIB) : ./accis/*.f ./accis/*.c ./accis/*.h
	make -C accis

######################################################
testing : testing/mpibbdytest testing/fieldtest testing/stresstest
	@echo Made tests in directory testing. Run them to test.

#mpibbdytest : mpibbdytest.o udisplay.o mpibbdy.o mditerate.o reduce.o
#	$(G77) -o mpibbdytest  mpibbdytest.f mpibbdy.o udisplay.o  mditerate.o reduce.o

testing/mpibbdytest : testing/mpibbdytest.o udisplay.o mpibbdy.o mditerate.o reduce.o makefile
	$(G77) -o testing/mpibbdytest  testing/mpibbdytest.f mpibbdy.o udisplay.o  mditerate.o reduce.o

testing/fieldtest : testing/fieldtest.f makefile libcoptic.a
	$(G77)  -o testing/fieldtest $(COMPILE-SWITCHES) $(PROFILING) testing/fieldtest.f libcoptic.a $(LIBRARIES)

testing/stresstest : testing/stresstest.f stress.o $(ACCISLIB)
	$(G77) -o testing/stresstest testing/stresstest.f stress.o $(LIBRARIES)

vecx :
	make clean
	make VECX=vecx -C accis
	make

#####################################################
clean :
	rm -f *.o $(TARGETS) *.html *.flx *.ph? *.den T*.0?? *.ps libcoptic.a
	make -C accis mproper

ftnchek :
	./ftnchekrun "$(COPTIC).f $(OBJECTS)"
	@echo To view do: konqueror CallTree.html &

coptic.prof : makefile $(OBJECTS) 
	make clean
	make PROFILING=-pg coptic
	make PROFILING=-pg coptic.o
	$(GFINAL)
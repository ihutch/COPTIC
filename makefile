GLULIBS= -lGL -lGLU
ACCISLIB=./accis/libaccisX.a
LIBRARIES = -L/usr/X11R6/lib/ -L./accis/ -laccisX -lXt -lX11 $(GLULIBS)
COPTIC=coptic
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
FIXEDOBJECTS=sormpi.o sorrelaxgen.o mpibbdy.o  cijroutine.o cijplot.o 3dobjects.o mditerate.o svdsol.o padvnc.o chargetomesh.o slicesect.o randf.o reindiag.o pinit.o ccpicplot.o volint.o fluxdata.o stringsnames.o meshconstruct.o partwriteread.o checkcode.o stress.o average.o bdyroutine.o reduce.o getfield.o interpolations.o
# Things just needed for the test routine:
UTILITIES=udisplay.o
OBJECTS=$(FIXEDOBJECTS) ${REINJECT}
HEADERS=bbdydecl.f meshcom.f objcom.f 3dcom.f partcom.f rancom.f ran1com.f creincom.f ptaccom.f
TARGETS=mpibbdytest mditeratetest sormpitest fieldtest
##########################################################################
G77=mpif77 -f77=g77 
#G77=mpif77
# export this so it is inherited by sub-makes.
export G77
GFINAL=gcc-4.1 -v -pg -o $(COPTIC).prof $(COPTIC).o $(OBJECTS) -static-libgcc -lpthread_p -lm_p -lc -lg2c -lmpich -lrt -lfrtbegin  $(LIBRARIES)
#OPTIMIZE=-O3 -funroll-loops -finline-functions
OPTIMIZE=-O3
COMPILE-SWITCHES = -Wall  $(OPTIMIZE)  -I. 
#COMPILE-SWITCHES = -Wall   $(OPTIMIZE) -I. -g -fbounds-check
#COMPILE-SWITCHES = -Wall   $(OPTIMIZE) -I. -g 
##COMPILE-SWITCHES = -Wall -Wno-unused $(OPTIMIZE) -g -I.
NOBOUNDS= $(COMPILE-SWITCHES) -fno-bounds-check
NOGLOBALS= $(COMPILE-SWITCHES) -Wno-globals
#PROFILING=-pg
#PROFILING= -pg -static-libgcc -lpthread_p -lm_p
##########################################################################
# If this rule does not seem to recognize the file you are trying to make,
# then run 'make' to completion first. It is something to do with the
# match-anything rules and prerequisites. I think that the rule is being
# interpreted as "terminal" which means it does not apply unless its
# prerequisites exist. By my interpretation of the info, this ought not
# to be happening (requires a double colon), but on horace, it is. 
% : %.f  makefile $(OBJECTS) $(UTILITIES) $(ACCISLIB)
	$(G77)  -o $* $(COMPILE-SWITCHES) $(PROFILING) $*.f $(OBJECTS) $(UTILITIES) $(LIBRARIES)

# Don't recompile accis every time the makefile is changed.
./accis/%.o : ./accis/%.f $(HEADERS)
	$(G77)  -c $(COMPILE-SWITCHES) $(PROFILING) $*.f

%.o : %.f makefile $(HEADERS)
	$(G77)  -c $(COMPILE-SWITCHES) $(PROFILING) $*.f

%.o : %.c makefile
	cc -c $(PROFILING) $*.c

#default target
# Problem when using geometry that can't do smt check. 
smt.out : $(COPTIC) ccpicgeom.dat
	@if [ -f smt.out ] ; then mv smt.out smt.prev ; fi
	./$(COPTIC)
	@if [ -f smt.prev ] ;then if [ -f smt.out ] ;then diff smt.prev smt.out ;else touch smt.out ;fi ;fi

#mpi checking target
mpicheck : $(COPTIC)
	mpiexec -l -n 2 ./$(COPTIC) >mpicheck.out
	if [ -f mpicheck.prev ] ; then diff mpicheck.prev mpicheck.out ; else mv mpicheck.out mpicheck.prev  ; fi
	mv mpicheck.out mpicheck.prev

#####################################################
# Things to compile with non-standard switches
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
$(COPTIC) : $(COPTIC).f makefile $(ACCISLIB) $(OBJECTS) $(UTILITIES)
	if [ -f "$(GEOMFILE)" ] ; then rm -f ccpicgeom.dat ; ln -s $(GEOMFILE) ccpicgeom.dat ; fi
	echo "      rjscheme="\'$(REINJECT)\'" " > REINJECT.f
	$(G77)  -o $(COPTIC) $(COMPILE-SWITCHES) $(PROFILING) $(COPTIC).f  $(OBJECTS) $(UTILITIES) $(LIBRARIES)

$(ACCISLIB) : ./accis/*.f ./accis/*.c ./accis/*.h
	make -C accis

######################################################
testing : testing/mpibbdytest testing/fieldtest testing/stresstest
	@echo Made tests in directory testing. Run them to test.

#mpibbdytest : mpibbdytest.o udisplay.o mpibbdy.o mditerate.o reduce.o
#	$(G77) -o mpibbdytest  mpibbdytest.f mpibbdy.o udisplay.o  mditerate.o reduce.o

testing/mpibbdytest : testing/mpibbdytest.o udisplay.o mpibbdy.o mditerate.o reduce.o makefile
	$(G77) -o testing/mpibbdytest  testing/mpibbdytest.f mpibbdy.o udisplay.o  mditerate.o reduce.o

testing/fieldtest : testing/fieldtest.f makefile $(OBJECTS) $(ACCISLIB)
	$(G77)  -o testing/fieldtest $(COMPILE-SWITCHES) $(PROFILING) testing/fieldtest.f $(OBJECTS) $(UTILITIES) $(LIBRARIES)

testing/stresstest : testing/stresstest.f stress.o $(ACCISLIB)
	$(G77) -o testing/stresstest testing/stresstest.f stress.o $(LIBRARIES)

#####################################################
clean :
	rm -f *.o $(TARGETS) *.html *.flx *.phi *.phiave *.den T*.0?? *.ps
	make -C accis mproper

ftnchek :
	./ftnchekrun "$(COPTIC).f $(OBJECTS)"
	@echo To view do: konqueror CallTree.html &

final : 
	$(GFINAL)

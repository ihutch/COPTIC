COPTIC=coptic
#########################################################################
ACCISLIB=./accis/libaccisX.a
LIBRARIES = -L/usr/X11R6/lib/ -L./accis/ -laccisX -lXt -lX11 $(GLULIBS)
#Default accis driver choice
#Alternatives are vecx or vecglx, can be overriden by commandline option
ifeq ("$(VECX)","")
	VECX=vecglx
endif	
# Alternative is    VECX=vecx
ifeq ("$(VECX)","vecglx")     
   GLULIBS= -lGL -lGLU
endif
##########################################################################
# The reinjection choice:
#######################
# This does not work with vaccheck because of outer boundary alteration.
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
# In g77 -Wno-globals silences spurious type messages on reduce.f
# This is unrecognized by gfortan. For which no-unused is better.
NGW=-Wno-unused
ifeq ("$(G77)","")
	G77=$(shell cat compiler 2>/dev/null)
	ifeq ("$(G77)","")	
		G77=mpif77 -f77=g77
	endif
	ifeq ("$(G77)","mpif77 -f77=g77")	
		NGW=-Wno-globals
	endif
endif
#
# export this so it is inherited by sub-makes.
export G77
##########################################################################
GFINAL=gcc-4.1 -v -pg -o $(COPTIC).prof $(COPTIC).o $(OBJECTS) -static-libgcc -lpthread_p -lm_p -lc -lg2c -lmpich -lrt -lfrtbegin  $(LIBRARIES)
#OPTIMIZE=-O3 -funroll-loops -finline-functions
OPTIMIZE=-O3
COMPILE-SWITCHES = -Wall  $(OPTIMIZE)  -I. 
#COMPILE-SWITCHES = -Wall   $(OPTIMIZE) -I. -g -fbounds-check
##COMPILE-SWITCHES = -Wall -Wno-unused $(OPTIMIZE) -g -I.
# Noboundscheck switches are not compatible with e.g. pathscale compiler:
ifeq ($(findstring g77,$(G77)),g77)
	NOBOUNDS= $(COMPILE-SWITCHES) -fno-bounds-check
else
	NOBOUNDS= $(COMPILE-SWITCHES)
endif
NOGLOBALS= $(COMPILE-SWITCHES) $(NGW)
##########################################################################
##########################################################################
FIXEDOBJECTS=sormpi.o sorrelaxgen.o mpibbdy.o cijroutine.o cijplot.o 3dobjects.o mditerate.o  padvnc.o chargetomesh.o slicesect.o randf.o reindiag.o pinit.o phisoluplot.o orbit3plot.o volint.o fluxdata.o stringsnames.o meshconstruct.o partwriteread.o partaccum.o checkcode.o stress.o average.o objplot.o randc.o cmdline.o
#
SPECIALOBJECTS=bdyroutine.o faddu.o reduce.o getfield.o interpolations.o 
# Things just needed for the test routine:
UTILITIES=udisplay.o
SOLOBJECTS= cijroutine.o mditerate.o mpibbdy.o sormpi.o sorrelaxgen.o meshconstruct.o getfield.o interpolations.o cijplot.o phisoluplot.o slicesect.o 3dobjects.o bdysetsol.o faddu.o
REGULAROBJECTS= $(FIXEDOBJECTS) ${REINJECT}
OBJECTS=$(SPECIALOBJECTS) $(REGULAROBJECTS)
HEADERS=bbdydecl.f meshcom.f objcom.f 3dcom.f partcom.f rancom.f ran1com.f creincom.f ptaccom.f colncom.f examdecl.f griddecl.f ptchcom.f mditcom.f sectcom.f plascom.f slpcom.f myidcom.f
SOLHEADERS= bbdydecl.f meshcom.f objcom.f 3dcom.f accis/world3.h 
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

##########################################
#default target
# Problem when using geometry that can't do smt check. 
smt.out : $(COPTIC) copticgeom.dat
	@if [ -f smt.out ] ; then mv smt.out smt.prev ; fi
	@if [ -f T1e0v000P020L1e0z005x05.phi ] ; then mv T1e0v000P020L1e0z005x05.phi prior.phi ; echo "Created prior.phi" ; fi
	./$(COPTIC)
	@if [ -f smt.prev ] ;then if [ -f smt.out ] ;then diff smt.prev smt.out ;else touch smt.out ;fi ;fi
	@if [ -f prior.phi ] && [ -f T1e0v000P020L1e0z005x05.phi ] ; then if ! diff T1e0v000P020L1e0z005x05.phi prior.phi ; then echo "**** RESULT CHANGED" ; fi; rm prior.phi; else echo "File[s] lacking to compare result."; fi 

# For now we are using a big hammer to ensure libcoptic is clean.
libcoptic.a : makefile $(OBJECTS) $(UTILITIES)
	rm -f libcoptic.a
	ar -rs libcoptic.a $(OBJECTS) $(UTILITIES)

libcopsol.a : makefile $(SOLOBJECTS) $(UTILITIES) $(SOLHEADERS)
	rm -f libcopsol.a
	ar -rs libcopsol.a $(SOLOBJECTS) $(UTILITIES)

copticgeom.dat : $(GEOMFILE)
	if [ -f "$(GEOMFILE)" ] ; then ln -s -f $(GEOMFILE) copticgeom.dat ; fi

#mpi checking target
mpicheck : $(COPTIC)
	mpiexec -l -n 2 ./$(COPTIC) >mpicheck.out
	if [ -f mpicheck.prev ] ; then diff mpicheck.prev mpicheck.out ; else mv mpicheck.out mpicheck.prev  ; fi
	mv mpicheck.out mpicheck.prev

# Attempt to implement makefile configure for compiler.
compiler :
	@if which mpif77 >/dev/null;\
 then echo -n "MPI system. ";\
  if which g77 >/dev/null ;then  echo -n "Force g77. ";GHERE="mpif77 -f77=g77";\
  else GHERE=mpif77 ; fi\
 else echo -n "Not MPI System. ";\
  if which g77 >/dev/null ; then GHERE="g77";fi\
 fi;\
 echo "G77="$${GHERE}; G77=$${GHERE}; echo $${G77} > compiler
# The question then is how to use this compiler setting.
# A recursive way is:
#	@make 
#	@rm compiler
# But this does not seem to inherit G77. So you have to make compiler;make.
## @if mpif77 | grep gfortran >/dev/null; then echo gfortran also available; fi

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
	echo G77=$G77
	@echo "      rjscheme="\'$(REINJECT)\'" " > REINJECT.f
	$(G77)  -o $(COPTIC) $(COMPILE-SWITCHES) $(PROFILING) $(COPTIC).f  libcoptic.a $(LIBRARIES)

# Sorserial links nonmpibbdy.o explicitly, so that none of those routines
# are linked from the main libcopsol that need mpi.
sorserial : sortest.f makefile $(ACCISLIB) $(SOLOBJECTS) nonmpibbdy.o
	$(G77) -o sorserial  $(COMPILE-SWITCHES) $(PROFILING) sortest.f nonmpibbdy.o libcopsol.a $(LIBRARIES)

$(ACCISLIB) : ./accis/*.f ./accis/*.c ./accis/*.h
	@echo "******************* Making accis with VECX=${VECX} **********"
	make -C accis

######################################################
testing : testing/mpibbdytest testing/fieldtest testing/stresstest
	@echo Made tests in directory testing. Run them to test.

testing/mpibbdytest : testing/mpibbdytest.o udisplay.o mpibbdy.o mditerate.o reduce.o makefile
	$(G77) -o testing/mpibbdytest  testing/mpibbdytest.f mpibbdy.o udisplay.o  mditerate.o reduce.o

testing/fieldtest : testing/fieldtest.f makefile libcoptic.a
	$(G77)  -o testing/fieldtest $(COMPILE-SWITCHES) $(PROFILING) testing/fieldtest.f libcoptic.a $(LIBRARIES)

testing/stresstest : testing/stresstest.f stress.o $(ACCISLIB)
	$(G77) -o testing/stresstest testing/stresstest.f stress.o $(LIBRARIES)

vecx :
	make clean
	make VECX=vecx

#####################################################
clean :
	rm -f *.o $(TARGETS) *.html *.flx *.ph? *.den T*.* *.ps *.aux *.log *.out *.toc libcoptic.a
	make -C accis mproper

ftnchek :
	./ftnchekrun "$(COPTIC).f $(OBJECTS)"
	@echo To view do: google-chrome CallTree.html

tree :
	./ftnchekrun "-nocheck $(COPTIC).f $(OBJECTS)"
	@echo To view do: 
	google-chrome CallTree.html

vcg :
	./ftnchekrun "-nocheck -vcg $(COPTIC).f $(OBJECTS)"

fordocu :
	testing/fordocu.sh "$(COPTIC).f $(OBJECTS)"
	google-chrome html/index.html

coptic.prof : makefile $(OBJECTS) 
	make clean
	make PROFILING=-pg coptic
	make PROFILING=-pg coptic.o
	$(GFINAL)

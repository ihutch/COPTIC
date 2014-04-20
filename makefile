
# This makefile assumes the shell is Bash.
########################################################################
# The eventual target:
COPTIC=coptic
#########################################################################
# To get to compile with X, you might need to supplement this path.
LIBPATH= -L./accis/ -L/usr/lib/mesa
#########################################################################
# Test whether X libraries are found. Null => yes.
 TESTGL:=$(shell ld  $(LIBPATH) -lGLU -lGL -o /dev/null none.o 2>&1 | grep GL)
 TESTX11:=$(shell ld $(LIBPATH) -lXt -lX11 -o /dev/null none.o 2>&1 | grep X)
##########################################################################
# The reinjection choice:
#######################
# This does not work with vaccheck because of outer boundary alteration.
#REINJECT=orbitinjnew.o extint.o
#GEOMFILE=geometry/geomsphere.dat
##################
#REINJECT=reinject.o
#GEOMFILE=geometry/geomsphere.dat
###################
REINJECT=cartreinject.o
GEOMFILE=geometry/geomcubic.dat
########################################################################
# Decide accis driver choice. Alternatives are vec4014 vecx or vecglx. 
# Automatic choice can be overriden by commandline option e.g. 
# make VECX=vec4014
# or in accis if libraries are unfound. 
# But we need to be able to tell which ACCISDRV to use based upon 
# accis configuration. So we decide here. 
ifeq ("$(VECX)","vec4014")
 ACCISDRV=accis
 LIBRARIES = $(LIBPATH) -l$(ACCISDRV)
else
 ifeq ("$(VECX)","vecx")
   ifeq ("$(TESTX11)","")
     ACCISDRV=accisX
     LIBRARIES = $(LIBPATH) -l$(ACCISDRV) -lXt -lX11 $(GLULIBS)
   else
# Wanted vecx but could not have it:
     ACCISDRV=accis
     LIBRARIES = $(LIBPATH) -l$(ACCISDRV)
     VECX:=vec4014
   endif
 else
# VECX not vec4014 or vecx
   ifeq ("$(TESTGL)","")
     ACCISDRV=accisX
     GLULIBS= -lGL -lGLU
     LIBRARIES = $(LIBPATH) -l$(ACCISDRV) -lXt -lX11 $(GLULIBS)
     VECX=vecglx
   else
    ifeq ("$(TESTX11)","")
     ACCISDRV=accisX
     GLULIBS=
     LIBRARIES = $(LIBPATH) -l$(ACCISDRV) -lXt -lX11 $(GLULIBS)
     VECX=vecx
    else
     ACCISDRV=accis
     LIBRARIES = $(LIBPATH) -l$(ACCISDRV)
     VECX:=vec4014
    endif
   endif
 endif
endif
ACCISLIB=./accis/lib$(ACCISDRV).a
# For submakes:
export VECX
##########################################################################
# Decide which compiler to use.
ifeq ("$(G77)","")
# I don't know why this has to be overridden. 
# But within this section of code G77 is not set without an override.
	override G77=$(shell cat compiler 2>/dev/null)
#	ifeq ("$(G77)","")
# Default compiler to be used if a strange make target is used 
# on the very first call. Not now set. Depend on compiler tests.
# After that, compiler ought to be set on disk and used.
#		override G77=mpif77
#	endif
endif
# In g77 -Wno-globals silences spurious type messages on reduce.f
# This is unrecognized by gfortan. For which no-unused is better.
ifeq ("$(G77)","mpif77 -f77=g77")	
  NGW=-Wno-globals
endif
ifeq ("$(findstring gfortran,$(G77))","")
else
  NGW=-Wno-unused	
endif
# export this so it is inherited by sub-makes.
export G77
##########################################################################
# A couple of special compile/link cases, not usually used.
GFINAL=gcc-4.1 -v -pg -o $(COPTIC).prof $(COPTIC).o $(OBJECTS) -static-libgcc -lpthread_p -lm_p -lc -lg2c -lmpich -lrt -lfrtbegin  $(LIBRARIES)
GCURR=gcc -v -pg -o $(COPTIC).prof $(COPTIC).o $(OBJECTS) -static-libgcc -lpthread_p -lm_p -lc -lg2c -lmpich -lrt -lfrtbegin  $(LIBRARIES)
#OPTIMIZE=-O3 -funroll-loops -finline-functions
OPTIMIZE=
COMPILE-SWITCHES=-I.
ifeq ("$(G77)","mpif77 -f77=g77")	
  COMPILE-SWITCHES = -Wall $(OPTIMIZE) -I.
  OPTIMIZE=-O3
endif	
#COMPILE-SWITCHES = -Wall   $(OPTIMIZE) -I. -g -fbounds-check
##COMPILE-SWITCHES = -Wall -Wno-unused $(OPTIMIZE) -g -I.
# Noboundscheck switches are not compatible with e.g. pathscale compiler:
ifeq ($(findstring g77,"$(G77)"),g77)
	NOBOUNDS= $(COMPILE-SWITCHES) -fno-bounds-check
else
	NOBOUNDS= $(COMPILE-SWITCHES)
endif
NOGLOBALS= $(COMPILE-SWITCHES) $(NGW)
##########################################################################
##########################################################################
FIXEDOBJECTS=sormpi.o sorrelaxgen.o cijroutine.o cijplot.o 3dobjects.o mditerate.o padvnc.o chargetomesh.o slicesect.o randf.o reindiag.o pinit.o phisoluplot.o orbit3plot.o volint.o fluxdata.o stringsnames.o meshconstruct.o partwriteread.o partaccum.o checkcode.o stress.o average.o objplot.o cmdline.o fsects.o bdyshare.o randc.o toms659.o
ifeq ("$(findstring mpi,"$(G77)")","")
# non MPI compiler (e.g. gfortran) is used
       MPIOBJECTS:=dummyreduce.o nonmpibbdy.o
else
       MPIOBJECTS:=reduce.o mpibbdy.o
endif
SPECIALOBJECTS=bdyroutine.o faddu.o getfield.o interpolations.o 
# Things just needed for the test routine:
UTILITIES=udisplay.o
SOLOBJECTS= cijroutine.o mditerate.o mpibbdy.o sormpi.o sorrelaxgen.o meshconstruct.o getfield.o interpolations.o cijplot.o phisoluplot.o slicesect.o 3dobjects.o bdysetsol.o faddu.o fsects.o cmdline.o bdyshare.o
REGULAROBJECTS= $(MPIOBJECTS) $(FIXEDOBJECTS) ${REINJECT} 
OBJECTS=$(SPECIALOBJECTS) $(REGULAROBJECTS)
HEADERS=bbdydecl.f meshcom.f objcom.f 3dcom.f partcom.f rancom.f ran1com.f creincom.f ptaccom.f colncom.f griddecl.f ptchcom.f mditcom.f sectcom.f plascom.f slpcom.f myidcom.f facebcom.f cdistcom.f ndimsdecl.f
# Nothing in root directory now depends on  examdecl.f
SOLHEADERS= bbdydecl.f meshcom.f objcom.f 3dcom.f accis/world3.h 
TARGETS=mpibbdytest mditeratetest sormpitest fieldtest
##########################################################################
# If this rule does not seem to recognize the file you are trying to make,
# then run 'make' to completion first. It is something to do with the
# match-anything rules and prerequisites. I think that the rule is being
# interpreted as "terminal" which means it does not apply unless its
# prerequisites exist.

% : %.f  makefile libcoptic.a $(ACCISLIB)
	$(G77)  -o $* $(NOGLOBALS) $(PROFILING) $*.f libcoptic.a $(LIBRARIES)

%.o : %.f makefile $(HEADERS)
	$(G77)  -c $(NOGLOBALS) $(PROFILING) $*.f

%.o : %.c makefile
	cc -c $(PROFILING) $*.c

# The testing target pattern used for files in geometry/
# Ensure we are using a standard allocated array size.
# Run coptic on the dat file. If it works and gives phi output, sum it.
# Compare with the old cks file and see if it's the same.
%.cks : %.dat $(COPTIC)
	./setdimens 64 64 128
	rm -f T*.phi
	./${COPTIC} $*.dat
	@if [ -f T*.phi ] ; then sum T*.phi >checksum ; else echo NO .phi FILE GENERATED; ls T*.phi 2>/dev/null; exit 1; fi
	@if diff checksum $*.cks >diffout 2>&1; then echo; echo "        Case $*.cks: OK. No differences"; touch $*.cks; else cat diffout; echo '******** Failed geometry test on $*.cks *********'; cat diffout >> GeometryTests; fi; rm -f diffout checksum
	@if [ -f $*.cks ] ; then echo ; else echo "******** $*.cks not present. Creating it."; sum *.phi >$*.cks; fi
	@echo -----------------------------------------------------------------

##########################################
# Default target compiler must always be the first dependency.
# Problem when using geometry that can't do smt check. 
smt.out : compiler $(COPTIC) copticgeom.dat
	@if [ -f smt.out ] ; then mv smt.out smt.prev ; fi
	@if [ -f T1e0v000P200L1e0z005x05.phi ] ; then mv T1e0v000P200L1e0z005x05.phi prior.phi ; echo "Created prior.phi" ; fi
	./$(COPTIC)
	@if [ -f smt.prev ] ;then if [ -f smt.out ] ;then diff smt.prev smt.out ;else touch smt.out ;fi ;fi
	@if [ -f prior.phi ] && [ -f T1e0v000P200L1e0z005x05.phi ] ; then if ! diff T1e0v000P200L1e0z005x05.phi prior.phi ; then echo "**** RESULT CHANGED" ; rm prior.phi; fi; else echo "File[s] lacking to compare result.\nProbably you've just made coptic for the first time."; fi 
	@if [ "$(G77)" = "gfortran" ] ; then echo "Compiled serial coptic. make clean; make for MPI version if MPI available." ; fi

# For now we are using a big hammer to ensure libcoptic is clean.
libcoptic.a : compiler makefile $(OBJECTS) $(UTILITIES)
#	@echo $(OBJECTS)
	rm -f libcoptic.a
	@ar -rs libcoptic.a $(OBJECTS) $(UTILITIES)

libcopsol.a : compiler makefile $(SOLOBJECTS) $(UTILITIES) $(SOLHEADERS)
	rm -f libcopsol.a
	ar -rs libcopsol.a $(SOLOBJECTS) $(UTILITIES)

copticgeom.dat : $(GEOMFILE)
	@if [ -f "$(GEOMFILE)" ] ; then ln -s -f $(GEOMFILE) copticgeom.dat ; fi

#mpi checking target
mpicheck : $(COPTIC)
	mpiexec -l -n 2 ./$(COPTIC) >mpicheck.out
	if [ -f mpicheck.prev ] ; then diff mpicheck.prev mpicheck.out ; else mv mpicheck.out mpicheck.prev  ; fi
	mv mpicheck.out mpicheck.prev

# Configure compiler. Mostly one long continued bash script.
compiler : makefile
	@echo -n Compiler tests. $${G77}
	@\
 if which $${G77} 2>/dev/null >/dev/null ; then\
  GHERE="$${G77}";\
 else\
 if which mpif77 >/dev/null;\
 then echo -n " MPI system. "; GHERE=mpif77;\
  if which g77 >/dev/null ; then\
    if [ -z "`mpif77 --version | grep Portland`" ] ; then\
	 echo -n " Force g77. "; GHERE="mpif77 -f77=g77";\
    fi\
  fi\
 else echo -n " Not MPI System. ";\
  if which g77 >/dev/null ; then\
     GHERE="g77";\
  else if which f77 >/dev/null ; then GHERE="f77";else\
          echo "$${G77} NO COMPILER found! Specify G77= ..." ; exit 1;\
       fi\
  fi;\
 fi;\
 fi;\
 echo "Chosen G77="$${GHERE}; G77=$${GHERE}; echo $${G77} > compiler;
# To obtain this information, one has to make a second time.
	@echo "*********** Remaking COPTIC with chosen G77 ****************"
	@export MAKEFLAGS=; export G77=$${GHERE}; $(MAKE) coptic
#	@echo "*********** Terminating recursive make. Not an error *******"
#	exit 1

#####################################################
# Things to compile with non-standard switches
# We make one of these the first thing in objects to force the header
# dependence to be reported, not just ignored by make on pattern rule.
interpolations.o : interpolations.f compiler makefile $(HEADERS)
	$(G77)  -c $(NOBOUNDS) $(PROFILING) $*.f

getfield.o : getfield.f compiler makefile $(HEADERS)
	$(G77)  -c $(NOBOUNDS) $(PROFILING) $*.f

reduce.o : reduce.f compiler makefile $(HEADERS)
	$(G77) -c $(NOGLOBALS) $(PROFILING) $*.f

bdyroutine.o : bdyroutine.f compiler makefile $(HEADERS)
	$(G77) -c $(NOGLOBALS) $(PROFILING) $*.f

#####################################################
# Main program explicit to avoid make bugs:
$(COPTIC) : compiler $(COPTIC).f makefile $(ACCISLIB) $(OBJECTS) $(UTILITIES) libcoptic.a
	@echo "      rjscheme="\'$(REINJECT)\'" " > REINJECT.f
	$(G77) -o $(COPTIC) $(COMPILE-SWITCHES) $(PROFILING) $(COPTIC).f libcoptic.a $(LIBRARIES)
#	@echo "TestGL:$(TESTGL), TestX11:$(TESTX11), vecx=$(VECX), G77=$(G77)"

# Sorserial links nonmpibbdy.o explicitly, so that none of those routines
# are linked from the main libcopsol that need mpi.
sorserial : libcopsol.a sortest.f compiler makefile $(ACCISLIB) $(SOLOBJECTS) nonmpibbdy.o
	$(G77) -o sorserial $(COMPILE-SWITCHES) $(PROFILING) sortest.f nonmpibbdy.o libcopsol.a $(LIBRARIES)

$(ACCISLIB) : ./accis/*.f ./accis/*.c ./accis/*.h
	@echo "******************* Making accis with VECX=${VECX} **********"
	make -C accis
	@if [ -f ./accis/$(VECX) ] ; then echo ; else echo "Failed making accis with $(VECX). Might need to specify a different driver."; fi

######################################################
testing : compiler $(COPTIC).f makefile $(ACCISLIB) $(OBJECTS) $(UTILITIES) libcoptic.a 
	make -C testing
	@echo Made tests in directory testing. Now running them to test.
	make -C testing testing

vecx :
	make clean
	make VECX=vecx

#####################################################
geometry : geometry/*.cks
	rm -f T1*
	date >>GeometryTests
	cat GeometryTests
#####################################################
clean :
	rm -f *.o $(TARGETS) *.html *.flx *.ph? *.den T*.* *.ps *.aux *.log *.out *.toc *.prev *.tlg *.synctex.gz ftnchek.output libcoptic.a storedgeom.dat
	make -C testing clean
	make -C accis mproper

mproper :
	rm -f compiler REINJECT.f coptic copticgeom.dat storedgeom.dat GeometryTests
	make clean
	make -C analysis clean

ftnchek :
	./ftnchekrun "$(COPTIC).f $(OBJECTS)" >ftnchek.output
	@echo To view do: google-chrome CallTree.html
	less ftnchek.output

tree :
	./ftnchekrun "-nocheck $(COPTIC).f $(OBJECTS)"
	@echo To view do: 
	firefox CallTree.html

vcg :
	./ftnchekrun "-nocheck -vcg $(COPTIC).f $(OBJECTS)"

fordocu :
	testing/fordocu.sh "$(COPTIC).f $(OBJECTS)"
	firefox html/index.html

coptic.prof : compiler makefile $(OBJECTS) 
	make clean
	make PROFILING=-pg coptic
	make PROFILING=-pg coptic.o
	$(GCURR)

help :
	@echo Targets: clean mproper ftnchek tree coptic.prof vecx
	@echo Tests:   geometry testing testanal

testanal : 
	make clean
	make -C analysis clean
	make
	make -C analysis
	analysis/partexamine -vtk T1e0v000P200L1e0z005x05
	analysis/phiexamine T1e0v000P200L1e0z005x05.pha -w
	analysis/fluxexamine -q T1e0v000P200L1e0z005x05.flx
	@echo "******* Completed tests with no obvious analysis errors."
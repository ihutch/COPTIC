# The shell assumed in this makefile
SHELL=/bin/bash
########################################################################
# The eventual target:
COPTIC=coptic
########################################################################
# The root directory of the distribution (for now)
COPTIC_ROOT=..
#########################################################################
# To get to compile with X, you might need to supplement this path.
LIBPATH= -L$(COPTIC_ROOT)/accis/ -L/usr/lib/mesa
#########################################################################
# Test whether X libraries are found. Null => yes.
 TESTGL:=$(shell ld  $(LIBPATH) -lGLU -lGL -o /dev/null none.o 2>&1 | grep GL)
 TESTX11:=$(shell ld $(LIBPATH) -lX11 -o /dev/null none.o 2>&1 | grep X)
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
GEOMFILE=$(COPTIC_ROOT)/geometry/geomcubic.dat
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
     LIBRARIES = $(LIBPATH) -l$(ACCISDRV) -lX11 $(GLULIBS)
   else
# Wanted vecx but could not have it:
     ACCISDRV=accis
     LIBRARIES=$(LIBPATH) -l$(ACCISDRV)
     VECX:=vec4014
   endif
 else
# VECX not vec4014 or vecx
   ifeq ("$(TESTGL)","")
     ACCISDRV=accisX
     GLULIBS=-lGL -lGLU
     LIBRARIES=$(LIBPATH) -l$(ACCISDRV) -lX11 $(GLULIBS)
     VECX=vecglx
   else
    ifeq ("$(TESTX11)","")
     ACCISDRV=accisX
     GLULIBS=
     LIBRARIES=$(LIBPATH) -l$(ACCISDRV) -lX11 $(GLULIBS)
     VECX=vecx
    else
     ACCISDRV=accis
     LIBRARIES=$(LIBPATH) -l$(ACCISDRV)
     VECX:=vec4014
    endif
   endif
 endif
endif
ACCISLIB=${COPTIC_ROOT}/accis/lib$(ACCISDRV).a
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
# export this so it is inherited by sub-makes.
export G77
##########################################################################
#OPTIMIZE=-O3 -funroll-loops -finline-functions
OPTIMIZE=-O3
COMPILE-SWITCHES =-Wall $(OPTIMIZE) -I.
#COMPILE-SWITCHES =-Wall   $(OPTIMIZE) -I. -g -fbounds-check
##########################################################################
# Decide on various compiler options based on compiler name.
# -fno-bounds-check is not recognized reliably except by g77 and gfortran
# In g77 -Wno-globals silences spurious type messages on reduce.f
# This is unrecognized by gfortran. no-unused-dummy... is better.
ifeq ("$(findstring g77,$(G77))","")
 NOBOUNDS=$(COMPILE-SWITCHES)
 ifeq ("$(findstring mpif77,$(G77))","")
  ifeq ("$(findstring gfortran,$(G77))","")
  else
# Allow unused dummy arguments to be warned with plain gfortran
    NGW=
  endif
 else
  ifeq ("$(findstring gfortran,$(shell mpif77 -show))","")
  else
    NGW=-Wno-unused-dummy-argument
    NOBOUNDS=$(COMPILE-SWITCHES) -fno-bounds-check
# Not now needed -Wno-conversion 
  endif
 endif
else
  NGW=-Wno-globals
  NOBOUNDS=$(COMPILE-SWITCHES) -fno-bounds-check
endif
NOGLOBALS:=$(COMPILE-SWITCHES) $(NGW)
##########################################################################
# A couple of special compile/link cases, not usually used.
# Probably good only for g77
#GCURR=gcc -v -pg -o $(COPTIC).prof $(COPTIC).o $(OBJECTS) -static-libgcc -lpthread_p -lm_p -lc -lg2c -lmpich -lrt -lfrtbegin  $(LIBRARIES)
GCURR=$(G77) -v -pg -o $(COPTIC).prof $(COPTIC).o $(OBJECTS) $(LIBRARIES)
##########################################################################
FIXEDOBJECTS=sormpi.o sorrelaxgen.o cijroutine.o cijplot.o 3dobjects.o mditerate.o padvnc.o chargetomesh.o slicesect.o reindiag.o pinit.o phisoluplot.o orbit3plot.o volint.o fluxdata.o stringsnames.o meshconstruct.o partwriteread.o partaccum.o checkcode.o stress.o average.o objplot.o cmdline.o fsects.o bdyshare.o toms659.o ranlux.o intersects.o trapinit.o
ifeq ("$(G77)","")
# First time through with undefined compiler don't define MPIOBJECTS.
else
 ifeq ("$(findstring mpi,"$(G77)")","")
# non MPI compiler (e.g. gfortran) is used
       MPIOBJECTS=dummyreduce.o nonmpibbdy.o
 else
       MPIOBJECTS=reduce.o mpibbdy.o
 endif
endif
SPECIALOBJECTS=bdyroutine.o faddu.o getfield.o interpolations.o 
# Things just needed for the test routine:
UTILITIES=udisplay.o
SOLOBJECTS=cijroutine.o mditerate.o mpibbdy.o sormpi.o sorrelaxgen.o meshconstruct.o getfield.o interpolations.o cijplot.o phisoluplot.o slicesect.o 3dobjects.o bdysetsol.o faddu.o fsects.o cmdline.o bdyshare.o
REGULAROBJECTS=$(MPIOBJECTS) $(FIXEDOBJECTS) ${REINJECT} 
OBJECTS=$(SPECIALOBJECTS) $(REGULAROBJECTS)
HEADERS=bbdydecl.f meshcom.f objcom.f 3dcom.f partcom.f rancom.f creincom.f ptaccom.f colncom.f griddecl.f ptchcom.f mditcom.f sectcom.f plascom.f slpcom.f myidcom.f facebcom.f cdistcom.f ndimsdecl.f reincom.f
# Nothing in root directory now depends on  examdecl.f
SOLHEADERS=bbdydecl.f meshcom.f objcom.f 3dcom.f $(COPTIC_ROOT)/accis/world3.h 
TARGETS=mpibbdytest mditeratetest sormpitest fieldtest
##########################################################################
# If this rule does not seem to recognize the file you are trying to make,
# then run 'make' to completion first. It is something to do with the
# match-anything rules and prerequisites. I think that the rule is being
# interpreted as "terminal" which means it does not apply unless its
# prerequisites exist.

% : %.f  makefile libcoptic.a $(ACCISLIB)
	$(G77) -o $* $(NOGLOBALS) $(PROFILING) $*.f libcoptic.a $(LIBRARIES)

%.o : %.f makefile $(HEADERS)
	$(G77) -c $(NOGLOBALS) $(PROFILING) $*.f

%.o : %.c makefile
	cc -c $(PROFILING) $*.c

# The testing target pattern used for files in geometry/
# Ensure we are using a standard allocated array size.
# Run coptic on the dat file. If it works and gives phi output, sum it.
# Compare with the old cks file and see if it's the same.
%.cks : %.dat $(COPTIC)	
	@./setdimens 64 64 128 >makegeom.log
	@rm -f T*.phi
	@./${COPTIC} $*.dat >>makegeom.log
	@if [ -f T*.phi ] ; then sum T*.phi >checksum ;\
	   else echo NO .phi FILE GENERATED; ls T*.phi 2>/dev/null; exit 1; fi
	@if diff checksum $*.cks >diffout 2>&1; then\
	  echo "        Case $*.cks: OK. No differences" >&2; touch $*.cks;\
	 else cat diffout;\
	  echo '******** Failed geometry test on $*.cks *********' >&2;\
	  echo "$*.cks:" >> GeometryTests; cat diffout >> GeometryTests;\
         fi; rm -f diffout checksum
	@if [ -f $*.cks ] ; then echo >/dev/null ;\
	 else echo "******** $*.cks not present. Creating it.";\
	  sum T*.phi >$*.cks;\
	 fi
#	@echo -----------------------------------------------------------------

##########################################
# Default target compiler must always be the first dependency.
tested : compiler $(COPTIC) copticgeom.dat
	@if [ -f T1e0v000P200L1e0z005x05.phi ] ; then mv T1e0v000P200L1e0z005x05.phi prior.phi ; echo "Created prior.phi" ; fi
	./$(COPTIC)
	@date > tested
	@if [ -f prior.phi ] && [ -f T1e0v000P200L1e0z005x05.phi ] ; then \
if [ `ls -s prior.phi | sed -e "s/[ ].*//"` -ne \
    `ls -s T1e0v000P200L1e0z005x05.phi | sed -e "s/[ ].*//"` ] ;then \
 echo "**** SIZE CHANGED. You probably ran a non-standard case just prior." ; \
else if ! diff T1e0v000P200L1e0z005x05.phi prior.phi ; then \
 echo "**** RESULT CHANGED" ; rm prior.phi; fi; fi;\
 else echo "File[s] lacking to compare result.";\
 echo "Probably you've just made coptic for the first time."; fi 
	@if [ "$(G77)" = "gfortran" ] ; then echo "Compiled serial coptic. make mproper; make for MPI version if MPI available." ; fi
	@make GeometryTests >/dev/null

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

#mpi checking target doesn't work well
#mpicheck : $(COPTIC)
#	mpiexec -n 2 ./$(COPTIC) >mpicheck.out
#	@if [ -f mpicheck.prev ] ; then diff mpicheck.prev mpicheck.out; else echo "No prior file to compare"  ; fi
#	mv mpicheck.out mpicheck.prev

# Configure compiler. Mostly one long continued bash script.
# Currently not forcing f77=g77 by default. 
# Preference order mpif77, g77, f77, gfortran
compiler : makefile
	@echo Compiler tests. $${G77}
	@\
 if which `echo $${G77} | sed -e "s/ .*$$//"` 2>/dev/null >/dev/null ; then\
  GHERE="$${G77}";\
 else\
 if which mpif77 >/dev/null;\
 then echo -n " MPI system. "; GHERE=mpif77;\
 else echo -n " Not MPI System. ";\
  if which g77 >/dev/null ; then\
     GHERE="g77";\
  else if which f77 >/dev/null ; then GHERE="f77";else\
          if which gfortran >/dev/null; then GHERE="gfortran"; else\
             echo "$${G77} NO COMPILER found! Specify G77= ..." ; exit 1;\
          fi\
       fi\
  fi;\
 fi;\
 fi;\
 echo "Chosen G77="$${GHERE}; G77=$${GHERE}; echo $${G77} > compiler;
# To obtain this information, one has to make a second time.
	@echo "*********** Remaking COPTIC with chosen G77 ****************"
	@export MAKEFLAGS=; export G77=$${GHERE}; $(MAKE) coptic

# Make the binary writes access='stream'. This has greater portability
# between compilers but is g77 incompatible.
streamset : 
	for file in *.f; do sed -i -e "/form=[']unformatted['][,)]/s/[']unformatted[']/\'unformatted\'\n     \$$  ,access=\'stream\'/" $${file}; done
	make mproper
	@echo
	@echo "**** Remake with FORTRAN 2003 compatible compiler e.g. G77=gfortran"

# Undo the access='stream' setting. Note tricky $$ escaping.
streamunset :
	for file in *.f ; do sed -i -e "/^.*form=[']unformatted[']$$/N;s/\n     [$$]  ,access=[']stream[']//" $${file}; done
	make

#####################################################
# Things to compile with non-standard switches
# We make one of these the first thing in objects to force the header
# dependence to be reported, not just ignored by make on pattern rule.
interpolations.o : interpolations.f compiler makefile $(HEADERS)
	$(G77) -c $(NOBOUNDS) $(NGW) $(PROFILING) $*.f

getfield.o : getfield.f compiler makefile $(HEADERS)
	$(G77) -c $(NOBOUNDS) $(NGW) $(PROFILING) $*.f

reduce.o : reduce.f compiler makefile $(HEADERS)
	$(G77) -c $(NOGLOBALS) $(PROFILING) $*.f

bdyroutine.o : bdyroutine.f compiler makefile $(HEADERS)
	$(G77) -c $(NOGLOBALS) $(PROFILING) $*.f

#####################################################
# Main program explicit to avoid make bugs:
$(COPTIC) : compiler $(COPTIC).f makefile $(ACCISLIB) $(OBJECTS) $(UTILITIES) libcoptic.a
	@echo "      rjscheme="\'$(REINJECT)\'" " > REINJECT.f
	$(G77) -o $(COPTIC) $(COMPILE-SWITCHES) $(PROFILING) $(COPTIC).f libcoptic.a $(LIBRARIES)

# Sorserial links nonmpibbdy.o explicitly, so that none of those routines
# are linked from the main libcopsol that need mpi.
sorserial : libcopsol.a sortest.f compiler makefile $(ACCISLIB) $(SOLOBJECTS) nonmpibbdy.o
	$(G77) -o sorserial $(COMPILE-SWITCHES) $(PROFILING) sortest.f nonmpibbdy.o libcopsol.a $(LIBRARIES)

$(ACCISLIB) : $(COPTIC_ROOT)/accis/*.f $(COPTIC_ROOT)/accis/*.c $(COPTIC_ROOT)/accis/*.h
	@echo "******************* Making accis with VECX=${VECX} **********"
	make -C $(COPTIC_ROOT)/accis
	@if [ -f $(COPTIC_ROOT)/accis/$(VECX) ] ; then echo ; else echo "Failed making accis with $(VECX). Might need to specify a different driver."; fi

vecx :
	make clean
	make VECX=vecx

######################################################
testing : compiler $(COPTIC).f makefile $(ACCISLIB) $(OBJECTS) $(UTILITIES) libcoptic.a 
	make -C $(COPTIC_ROOT)/testing
	@echo Made tests in directory testing. Now running them to test.
	make -C $(COPTIC_ROOT)/testing testing
	@echo If all programs finished without crashing. Thats good enough.

testanal : 
	make clean
	make -C $(COPTIC_ROOT)/analysis clean
	make
	make -C $(COPTIC_ROOT)/analysis
	$(COPTIC_ROOT)/analysis/partexamine -vtk T1e0v000P200L1e0z005x05
	$(COPTIC_ROOT)/analysis/phiexamine T1e0v000P200L1e0z005x05.pha -w
	$(COPTIC_ROOT)/analysis/fluxexamine -q T1e0v000P200L1e0z005x05.flx
	@echo "******* Completed tests with no obvious analysis errors."

geometry : $(COPTIC_ROOT)/geometry/*.cks
	@echo 'Geometry Tests completed or up to date' >&2
	rm -f T1*
	date >>GeometryTests
	cat GeometryTests

regeom :
	@echo 'Regenerating geometry/*.cks files' >&2
	for file in $(COPTIC_ROOT)/geometry/*.dat ; do rm $${file%.dat}.cks; make $${file%.dat}.cks; done

# Use compiler (which depends on makefile) as a test of major updates
GeometryTests : compiler
	@if [ "`sed -n '/g77/p' compiler`" == "" ]; then\
  echo ' ' >&2;\
  echo 'Starting background tests. Recompile only when completed.' >&2;\
  make geometry;\
else echo "Skipping background GeometryTests" >&2 ; fi &
	@sleep 1
	@echo

#####################################################
clean :
	rm -f *.o $(TARGETS) *.html *.flx *.ph? *.den T*.* *.ps *.aux *.log *.out *.toc *.prev *.tlg *.synctex.gz ftnchek.output libcoptic.a storedgeom.dat
	make -C $(COPTIC_ROOT)/testing clean
	make -C $(COPTIC_ROOT)/accis mproper

mproper :
	rm -f compiler REINJECT.f coptic copticgeom.dat storedgeom.dat GeometryTests tested
	make clean
	make -C $(COPTIC_ROOT)/analysis clean

ftnchek :
	ftnchek  -nopure\
 -truncation=no-significant-figures\
 -usage=no-com-block-unused,no-var-set-unused,no-arg-array-alias,no-ext-unused\
 -arguments=no-arrayness\
 -noarray -noextern\
 -calltree=tree,no-prune\
 -nonovice\
 -mkhtml=vars-and-consts-lowercase,keywords-lowercase\
 -quiet\
 `echo "$(COPTIC).f $(OBJECTS)" | sed -e "s/[.]o/[.]f/g"` >ftnchek.output
	@echo To view do: firefox CallTree.html
	@less ftnchek.output

tree :
	ftnchek  -nopure\
 -truncation=no-significant-figures\
 -usage=no-com-block-unused,no-var-set-unused,no-arg-array-alias,no-ext-unused\
 -arguments=no-arrayness\
 -noarray -noextern\
 -calltree=tree,no-prune\
 -nonovice\
 -mkhtml=vars-and-consts-lowercase,keywords-lowercase\
 -quiet -nocheck\
 `echo "$(COPTIC).f $(OBJECTS)" | sed -e "s/[.]o/[.]f/g"` >ftnchek.output
	@echo To view do: 
	firefox CallTree.html

vcg :
	ftnchek  -nopure\
 -truncation=no-significant-figures\
 -usage=no-com-block-unused,no-var-set-unused,no-arg-array-alias,no-ext-unused\
 -arguments=no-arrayness\
 -noarray -noextern\
 -calltree=tree,no-prune\
 -nonovice\
 -mkhtml=vars-and-consts-lowercase,keywords-lowercase\
 -quiet -nocheck -vcg\
 `echo "$(COPTIC).f $(OBJECTS)" | sed -e "s/[.]o/[.]f/g"`

fordocu :
	$(COPTIC_ROOT)/testing/fordocu.sh "$(COPTIC).f $(OBJECTS)"
	firefox html/index.html

coptic.prof : compiler makefile $(OBJECTS) 
	make clean
	make PROFILING=-pg coptic
	make PROFILING=-pg coptic.o
	$(GCURR)

help :
	@echo Targets: clean mproper ftnchek tree coptic.prof vecx
	@echo Tests:   geometry testing testanal regeom
	@echo Setup:   streamset streamunset

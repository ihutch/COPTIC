# The shell assumed in this makefile
SHELL=/bin/bash
########################################################################
# The eventual target:
COPTIC=coptic
########################################################################
# The root directory of the distribution (for now)
COPTIC_ROOT=..
# Defining ACCISPARENT prevents construction of accis elsewhere.
ACCISPARENT:=$(realpath $(COPTIC_ROOT))
export ACCISPARENT
include ACCIS.mk
# To get to compile with X, you might need to supplement this path.
LIBPATH:= -L$(COPTIC_ROOT)/accis/ -L/usr/lib/mesa
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
#########################################################################
ifeq ("$(G77)","")
# Configure compiler. Mostly one long continued bash script.
# Preference order mpif77, ftn, g77, f77, gfortran
 COMPILER:=\
$(shell \
 if which mpif77 >/dev/null 2>&1; then echo -n "mpif77";else\
  if which ftn >/dev/null 2>&1 ; then echo -n "ftn";else\
    if which g77 >/dev/null 2>&1 ; then echo -n "g77";else\
     if which f77 >/dev/null 2>&1 ; then echo -n "f77";else\
      if which gfortran >/dev/null 2>&1; then echo -n "gfortran";else\
	echo "Unable to decide compiler. Specify via G77=..." >&2; exit 1;\
      fi;\
     fi;\
    fi;\
  fi;\
 fi;\
)
G77:=$(COMPILER)
else
COMPILER:=$(G77)
endif
export G77
#########################################################################
# In the Cray environment we need to tell explicitly to link -dynamic.
ifeq ("$(G77)","ftn")
LIBPATH:=-dynamic $(LIBPATH)
endif
export LIBPATH	
##########################################################################
# Test whether X libraries are found. Null => yes.
 TESTGL:=$(shell $(COMPILER)  $(LIBPATH) -lGLU -lGL -o /dev/null 2>&1 | grep GL)
 TESTX11:=$(shell $(COMPILER) $(LIBPATH) -lX11 -o /dev/null 2>&1 | grep X)
# Test if fast fourier transform libraries are available. 
 TESTFFTW:=$(shell $(COMPILER) $(LIBPATH) -lfftw3 -o /dev/null 2>&1 | grep fftw)
########################################################################
# Decide accis driver choice. Alternatives are vec4014 vecx or vecglx. 
# Automatic choice can be overriden by commandline option e.g. 
# make VECX=vec4014. Default order now: vecx, vec4014.
# This avoids dependency headaches with libGL and remote displays.
ifeq ("$(VECX)","vec4014")
 ACCISDRV=accis
 LIBRARIES = $(LIBPATH) -l$(ACCISDRV)
else
 ifeq ("$(VECX)","vecglx")
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
 else
# VECX not vec4014 or vecglx. Default trying vecx.
   ifeq ("$(TESTX11)","")
     ACCISDRV=accisX
     LIBRARIES = $(LIBPATH) -l$(ACCISDRV) -lX11
     VECX:=vecx
   else
# Wanted vecx but could not have it:
     ACCISDRV=accis
     LIBRARIES=$(LIBPATH) -l$(ACCISDRV)
     VECX:=vec4014
   endif
 endif
endif
FFTWSOLVE=fftdummy.o
ifeq ("$(TESTFFTW)","")
	LIBRARIES:=$(LIBRARIES) -lfftw3
#	FFTWSOLVE:=fftsolve.o
	FFTWSOLVE:=ffttrid.o
endif
ACCISLIB=${COPTIC_ROOT}/accis/lib$(ACCISDRV).a
# For submakes:
export VECX
##########################################################################
# Cray compilers don't recognize -Wall and optimize by default
ifneq ("$(G77)","ftn")
	OPTIMIZE=-O3
        COMPILE-SWITCHES:=$(COMPILE-SWITCHES) -Wall $(OPTIMIZE) -I.
else
# -fPIC Not needed. All that's needed is the  -dynamic switch for link.
endif
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
  ifneq ("$(findstring gfortran,$(shell mpif77 -show))","")
    NGW=-Wno-unused-dummy-argument
    NOBOUNDS=$(COMPILE-SWITCHES) -fno-bounds-check
# Not now needed -Wno-conversion 
  endif
 endif
else
  NGW=-Wno-globals
  NOBOUNDS=$(COMPILE-SWITCHES) -fno-bounds-check
endif
PROFILE= -mcmodel=medium $(PROFILING)
##########################################################################
FIXEDOBJECTS=sormpi.o sorrelaxgen.o cijroutine.o cijplot.o 3dobjects.o mditerate.o padvnc.o chargetomesh.o slicesect.o reindiag.o pinit.o phisoluplot.o orbit3plot.o volint.o fluxdata.o stringsnames.o meshconstruct.o partwriteread.o partaccum.o checkcode.o stress.o average.o objplot.o cmdline.o fsects.o bdyshare.o toms659.o ranlux.o intersects.o trapinit.o mainroutines.o phasespace.o qinit.o $(FFTWSOLVE) xofcum.o
MPIOBJECTS=reduce.o mpibbdy.o
DUMMYMPIOBJECTS=nonmpibbdy.o dummyreduce.o
export DUMMYMPIOBJECTS
ifeq ("$(findstring mpi,"$(G77)")","")
   ifeq ("$(findstring ftn,"$(G77)")","")
# Assume Non-MPI compiler (e.g. gfortran) is used
       MPIOBJECTS=dummyreduce.o nonmpibbdy.o
       NGW=-Wno-unused-dummy-argument
   endif
endif
SPECIALOBJECTS=bdyroutine.o faddu.o getfield.o interpolations.o 
# Things just needed for the test routine:
UTILITIES=udisplay.o
SOLOBJECTS=cijroutine.o mditerate.o mpibbdy.o sormpi.o sorrelaxgen.o meshconstruct.o getfield.o interpolations.o cijplot.o phisoluplot.o slicesect.o 3dobjects.o bdysetsol.o faddu.o fsects.o cmdline.o bdyshare.o
REGULAROBJECTS=$(MPIOBJECTS) $(FIXEDOBJECTS) ${REINJECT} 
OBJECTS=$(SPECIALOBJECTS) $(REGULAROBJECTS)
HEADERS=bbdydecl.f meshcom.f objcom.f 3dcom.f partcom.f rancom.f creincom.f ptaccom.f colncom.f griddecl.f ptchcom.f mditcom.f sectcom.f plascom.f slpcom.f myidcom.f facebcom.f cdistcom.f ndimsdecl.f reincom.f fvcom.f
# Nothing in root directory now depends on  examdecl.f
SOLHEADERS=bbdydecl.f meshcom.f objcom.f 3dcom.f $(COPTIC_ROOT)/accis/world3.h 
TARGETS=mpibbdytest mditeratetest sormpitest fieldtest
NOGLOBALS:=$(COMPILE-SWITCHES) $(NGW)
##########################################################################
# If this rule does not seem to recognize the file you are trying to make,
# then run 'make' to completion first. It is something to do with the
# match-anything rules and prerequisites. I think that the rule is being
# interpreted as "terminal" which means it does not apply unless its
# prerequisites exist.

% : %.f  makefile libcoptic.a $(ACCISLIB)
	$(G77) -o $* $(NOGLOBALS) $(PROFILE) $*.f libcoptic.a $(LIBRARIES)

%.o : %.f makefile $(HEADERS)
	$(G77) -c $(NOGLOBALS) $(PROFILE) $*.f

%.o : %.f90  makefile $(HEADERS)
	$(G77) -c $(NOGLOBALS) $(PROFILE) $*.f90

%.o : %.c makefile
	cc -c $(PROFILE) $*.c

# The testing target pattern used for files in geometry/
# Ensure we are using a standard allocated array size.
# Run coptic on the dat file. If it works and gives phi output, sum it.
# Compare with the old cks file and see if it's the same.
%.cks : %.dat $(COPTIC)
	@rm -f T*.phi
	@./${COPTIC} $*.dat >>makegeom.log
	@if [ -f T*.phi ] ; then sum T*.phi >checksum ;\
	   else echo NO .phi FILE GENERATED; ls T*.phi 2>/dev/null; exit 1; fi
	@if diff checksum $*.cks >diffout 2>&1; then\
	  touch $*.cks;\
	 else cat diffout;\
	  echo '******** Failed geometry test on $*.cks *********' >&2;\
	  echo "$*.cks:" >> GeometryTests; cat diffout >> GeometryTests;\
         fi; rm -f diffout checksum
	@if [ -f $*.cks ] ; then echo >/dev/null ;\
	 else echo "******** $*.cks not present. Creating it.";\
	  sum T*.phi >$*.cks;\
	 fi

##########################################
# Default Target
tested : $(COPTIC) copticgeom.dat
	@if [ -f T1P-200L1x5y5z5.phi ] ; then mv T1P-200L1x5y5z5.phi prior.phi ; echo "Created prior.phi" ; fi
	./$(COPTIC)
	@date > tested
	@if [ -f prior.phi ] && [ -f T1P-200L1x5y5z5.phi ] ; then \
if [ `ls -s prior.phi | sed -e "s/[ ].*//"` -ne \
    `ls -s T1P-200L1x5y5z5.phi | sed -e "s/[ ].*//"` ] ;then \
 echo "**** SIZE CHANGED. You probably ran a non-standard case just prior." ; \
else if ! diff T1P-200L1x5y5z5.phi prior.phi ; then \
 echo "**** RESULT CHANGED" ; rm prior.phi; fi; fi;\
 else echo "File[s] lacking to compare result.";\
 echo "Probably you've just made coptic for the first time."; fi 
	@if [ "$(G77)" = "gfortran" ] ; then echo "Compiled serial coptic. make mproper; make for MPI version if MPI available." ; fi
	@make --no-print-directory NOACCISCHECK=true GeometryTests

# For now we are using a big hammer to ensure libcoptic is clean.
libcoptic.a : makefile $(OBJECTS) $(UTILITIES)
#	@echo $(OBJECTS)
	rm -f libcoptic.a
	@ar -rs libcoptic.a $(OBJECTS) $(UTILITIES)

libcopsol.a : makefile $(SOLOBJECTS) $(UTILITIES) $(SOLHEADERS)
	rm -f libcopsol.a
	ar -rs libcopsol.a $(SOLOBJECTS) $(UTILITIES)

# Create griddecl and copticgeom.dat only if they don't exist.
griddecl.f :
	./setdimens ../geometry/geomcubic.dat

copticgeom.dat : $(GEOMFILE)
	@if [ -f "$(GEOMFILE)" ] ; then ln -s -f $(GEOMFILE) copticgeom.dat ; fi

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
interpolations.o : interpolations.f makefile $(HEADERS)
	$(G77) -c $(NGW) $(PROFILE) $(NOBOUNDS) $*.f

getfield.o : getfield.f makefile $(HEADERS)
	$(G77) -c $(NGW) $(PROFILE) $(NOBOUNDS) $*.f

reduce.o : reduce.f makefile $(HEADERS)
	$(G77) -c $(NOGLOBALS) $(PROFILE) $*.f

bdyroutine.o : bdyroutine.f makefile $(HEADERS)
	$(G77) -c $(NOGLOBALS) $(PROFILE) $*.f

#####################################################
# Main program explicit to avoid make bugs:
$(COPTIC) : $(COPTIC).f makefile $(ACCISLIB) $(OBJECTS) $(UTILITIES) libcoptic.a
	@echo $(G77) > compiler
	@echo "      rjscheme="\'$(REINJECT)\'" " > REINJECT.f
	$(G77) -o $(COPTIC) $(COMPILE-SWITCHES) $(PROFILE) $(COPTIC).f libcoptic.a $(LIBRARIES)

# Sorserial links nonmpibbdy.o explicitly, so that none of those routines
# are linked from the main libcopsol that need mpi.
sorserial : libcopsol.a sortest.f makefile $(ACCISLIB) $(SOLOBJECTS) nonmpibbdy.o
	$(G77) -o sorserial $(COMPILE-SWITCHES) $(PROFILE) sortest.f nonmpibbdy.o libcopsol.a $(LIBRARIES)

$(ACCISLIB) : $(COPTIC_ROOT)/accis/*.f $(COPTIC_ROOT)/accis/*.c $(COPTIC_ROOT)/accis/*.h
#	@echo "TESTGL=$(TESTGL). TESTX11=$(TESTX11)."
	@echo "G77=$(G77). COMPILER=$(COMPILER)."
	@echo "******************* Making accis with VECX=${VECX} **********"
	make -C $(COPTIC_ROOT)/accis
	@if [ -f $(COPTIC_ROOT)/accis/$(VECX) ] ; then echo ; else echo "Failed making accis with $(VECX). Might need to specify a different driver."; fi

phasespace.o : phasespace.f makefile $(HEADERS) phasecom.f
	$(G77) -c $(NOGLOBALS) $(PROFILE) phasespace.f

vecx :
	make clean
	make VECX=vecx

analysis : $(DUMMYMPIOBJECTS)
	make -C $(COPTIC_ROOT)/analysis

serial : $(COPTIC) $(DUMMYMPIOBJECTS)
	$(G77) -o $(COPTIC).serial $(DUMMYMPIOBJECTS) $(COMPILE-SWITCHES) $(PROFILE) $(COPTIC).f libcoptic.a $(LIBRARIES)
######################################################
testing : $(COPTIC).f makefile $(ACCISLIB) $(OBJECTS) $(UTILITIES) libcoptic.a 
	make -C $(COPTIC_ROOT)/testing
	@echo Made tests in directory testing. Now running them to test.
	make -C $(COPTIC_ROOT)/testing testing
	@echo If all programs finished without crashing. Thats good enough.

testanal : 
	make clean
	make -C $(COPTIC_ROOT)/analysis clean
	make $(COPTIC)
	make -C $(COPTIC_ROOT)/analysis
	./$(COPTIC)
	$(COPTIC_ROOT)/analysis/partexamine -vtk T1e0v000P200L1e0z005x05
	$(COPTIC_ROOT)/analysis/phiexamine T1e0v000P200L1e0z005x05.pha -w
	$(COPTIC_ROOT)/analysis/fluxexamine -q T1e0v000P200L1e0z005x05.flx
	@echo "******* Completed tests with no obvious analysis errors."

geometry : $(COPTIC_ROOT)/geometry/*.cks
	@echo 'Geometry Tests completed or up to date' >&2
	@rm -f T1*
	@if [ -z "`tail -n1 GeometryTests | grep EST`" ] ; then echo "****** Differences were found. Check file GeometryTests ********"; else echo "All OK"; fi
	@date >>GeometryTests
#	tail GeometryTests
#	cd $(COPTIC_ROOT)/src/; ./coptic; cd -

regeom :
	@echo 'Regenerating geometry/*.cks files' >&2
	for file in $(COPTIC_ROOT)/geometry/*.dat ; do rm $${file%.dat}.cks; make $${file%.dat}.cks; done

GeometryTests :  $(COPTIC)
	@if ./setdimens 64 64 128 NO ; then\
  if [ "$(G77)" != "mpif77" ] ; then echo "Not doing geometry tests" >&2;\
  else\
   if $(G77) --version 2>&1 | grep "GNU Fortran" ; then\
    echo '****Starting background tests. Recompile only when completed.' >&2;\
    sleep 1; make geometry &\
   fi; \
  fi; \
else\
  echo Grid changes nececessary for GeometryTests. So not done.;\
  echo If you want to do tests, first do make newgrid.;\
fi

newgrid : griddeclnew.f
	mv griddeclnew.f griddecl.f
	date > newgrid

#####################################################
clean :
	rm -f *.o $(TARGETS) *.html *.flx *.ph? *.den T*.* *.ps *.aux *.log *.out *.toc *.prev *.tlg *.synctex.gz ftnchek.output libcoptic.a storedgeom.dat
	make -C $(COPTIC_ROOT)/testing clean
	make -C $(COPTIC_ROOT)/accis mproper

mproper :
	rm -f REINJECT.f coptic copticgeom.dat storedgeom.dat GeometryTests tested
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

coptic.prof : makefile $(OBJECTS) 
	make clean
	make PROFILING=-pg coptic
	make PROFILING=-pg coptic.o
	$(G77) -v -pg -o $(COPTIC).prof $(COPTIC).o $(OBJECTS) $(LIBRARIES)

help :
	@echo
	@echo Targets: clean mproper ftnchek tree coptic.prof vecx analysis serial
	@echo Tests:   GeometryTests testing testanal regeom
	@echo Setup:   streamset streamunset
	@echo Profiling etc do:    make PROFILING='<switches>'
	@echo

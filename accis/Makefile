########################################################################
# Makefile for accis routines.
# The shell assumed in this Makefile
SHELL=/bin/bash
#########################################################################
# To get to compile with X, you might need to supplement this path.
LIBPATH:=-L.
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
	echo "Unable to decide compiler. Specify via G77=..."; exit 1;\
      fi;\
     fi;\
    fi;\
  fi;\
 fi;\
)
G77:=$(COMPILER)
else
COMPILER=$(G77)
endif
ifeq ("$(G77)","ftn")
LIBPATH:=-dynamic $(LIBPATH)
endif	
G90=gfortran
#######################################################################
ifeq ("$(NOBACKSLASH)","")
NOBACKSLASH:=\
$(shell \
if [ -n "`$(COMPILER) --version 2>&1 | grep GNU`" ] ; then\
 echo "-w -fno-backslash";else \
if [ -n "`$(COMPILER) --version 2>&1 | grep PathScale`" ] ; then\
 echo "-backslash";else \
if [ -n "`$(COMPILER) --version 2>&1 | grep Portland`" ] ; then\
 echo "-Mbackslash"; \
fi fi fi\
)
endif
##########################################################################
# Test whether X libraries are found. Null => yes.
 TESTGL:=$(shell $(COMPILER)  $(LIBPATH) -lGLU -lGL -o /dev/null none.o 2>&1 | grep GL)
 TESTX11:=$(shell $(COMPILER) $(LIBPATH) -lX11 -o /dev/null 2>&1 none.o | grep X)
##########################################################################
ifneq ("$(VECX)","")
# VECX explicitly set. Use the tests to convey the choice.
 ifneq ("$(VECX)","vecglx")
  TESTGL=NO
 endif
 ifneq ("$(VECX)","vecx")
  TESTX11=NO
 endif
endif
##########################################################################
# VECX choice. Preference order vecglx, vecx, vec4014.
ifeq ("$(TESTGL)","")
    VECX=vecglx
    ACCISDRV=accisX
    libraries=-lX11 -lGL -lGLU $(LIBPATH)
else
    ifeq ("$(TESTX11)","")
     VECX=vecx
     ACCISDRV=accisX
     libraries=-lX11 $(LIBPATH)
    else
     VECX=vec4014
     ACCISDRV=accis
     libraries=-L.
    endif
endif
#########################################################################
ifeq ("$(G77)","gfortran")
	WSWITCHES=-Wall -Wno-conversion
endif
#########################################################################
standard_object_files= \
autopl.o \
contour.o \
fitrange.o \
axis.o \
axlabels.o \
axregion.o \
charsize.o \
lautopl.o \
pltinit.o \
vecw.o \
polyline.o \
labeline.o \
minmax.o \
folwcont.o \
scale.o \
vec4014.o \
ticset.o \
hid.o \
vec3.o \
examine.o \
axon.o \
arrowplot.o\
initiald.o\
sysfont.o\
slicing.o\
fillgrad.o\
boxcarave.o\
surfaces.o

noback_object_files = drwstr.o fontdata.o vecnp.o 

object_files=$(standard_object_files) $(noback_object_files)
root_files=$(patsubst %.o,%,$(object_files))

#header include files
headers= hidcom.h plotcom.h world3.h

#The real Makefile
MAKEFILE=Makefile

#######################################################################
# Start of targets
#######################################################################

all : $(MAKEFILE) lib$(ACCISDRV).a

RefManual.html : RefManual.tex RefManual.pdf
	if which tth ; then tth -e2 RefManual ; fi
	rm -f *.log *.out *.dvi *.toc *.ilg *.idx *.synctex.gz

RefManual.pdf : RefManual.tex
	pdflatex RefManual
	rm -f *.log *.out *.dvi *.toc *.ilg *.idx *.synctex.gz

#update the libraries.
libaccis.a : $(object_files)
	@echo "Updating libaccis. For $(G77), $(VECX), $(ACCISDRV)."
	ar -rs libaccis.a $(object_files)

libaccisX.a : libaccis.a $(VECX)

# Headers must be updated, but this section seems to override the
# implicit rule unless the header dependency is :: double.
# Objects depend on headers
$(standard_object_files) :: $(headers) $(MAKEFILE)

#vecx is the C linkage to X-window. Threading is not (yet) completed.
vecx.o : vecx.c $(MAKEFILE)
	$(CC)  $(THREADING) -c vecx.c

vecglx.o : vecglx.c
	$(CC) -c $(THREADING) vecglx.c
#	$(CC) -c $(WSWITCHES) $(THREADING) vecglx.c

# The file drwstr.f must be compiled with this switch which disables
# for gnu compilers the interpretation of backslashes as escapes.
# drwstr needs the NOBACKSLASH switch (and can't be in :: rule above). 
$(noback_object_files) : drwstr.f fontdata.f vecnp.f $(headers) $(MAKEFILE)
	$(G77) -c $(NOBACKSLASH) $(WSWITCHES) $*.f

# Specific test programs of one sort and another should be put here 
# only if they need special switches.
testing/plottest90 : testing/plottest90.f90 interface.f90
	$(G90) $(WSWITCHES) -o testing/plottest90 testing/plottest90.f90 -l$(ACCISDRV) $(libraries)

#pattern rule, compile using the external definitions of commons
%.o : %.f ;
	$(G77) -c $(WSWITCHES) $*.f

# For fortran 90 executables.
% : %.f90 lib$(ACCISDRV).a $(VECX)
	$(G90) $(WSWITCHES) -o $* $*.f90 -l$(ACCISDRV)  $(libraries)

# The main f77 executable pattern.
% : %.f lib$(ACCISDRV).a $(VECX)
	$(G77) $(WSWITCHES) -o $* $*.f -l$(ACCISDRV)  $(libraries)

%.4014 : %.f libaccis.a
	$(G77) $(WSWITCHES) -o $*.4014 $*.f -L. -laccis

vecx : libaccis.a vecx.o
	cp libaccis.a libaccisX.a
	ar -d libaccisX.a vec4014.o
	ar -q libaccisX.a vecx.o
	date >vecx
	rm -f vecglx

vecglx : libaccis.a vecglx.o
	cp libaccis.a libaccisX.a
	ar -d libaccisX.a vec4014.o
	ar -q libaccisX.a vecglx.o
	date >vecglx
	rm -f vecx

vec4014 : libaccis.a
	date >vec4014

# Fortran 90 interface
convert : convert.f95
# Really convert.f is a f90 program. But named so not accidentally removed.
	$(G90) -o convert convert.f95
	@rm *.mod

# There are problems with Metcalf's automatic converter and parameters that
# define dimensions. Parameters are not included. So the interfaces that it
# generates are broken. I've hacked the code to include them.

interface.f90 : convert
	@echo "! Interfaces to accis routines generated from fortran code">interface.f90
	@echo making interface.f90, using convert code.
	@for file in $(root_files); do echo $${file} 0 0 t t | ./convert>/dev/null; done
	@for file in $(root_files); do sed -i $${file}.f90 -e"/BLOCKDATA/,/BLOCKDATA/d" ; done
	@for file in $(root_files); do echo "! $${file}.f">>interface.f90; cat $${file}.f90>>interface.f90; echo " ">>interface.f90; done
# If you wish to keep individual interfaces. Comment the following out.
	@for file in $(root_files); do rm $${file}.f90; done 

# Synchronization of versions.
sync : syncsource synccoptic syncsceptic

syncsilas : lib$(ACCISDRV).a RefManual.html
	rsync -u -e ssh  --copy-links -v *.h *.f *.c RefManual.* Makefile silas:~/accis/
	date > syncsilas

syncsource : lib$(ACCISDRV).a RefManual.html
	cd ~/src/accis/ ; git pull origin
	date > syncsource

synccoptic : lib$(ACCISDRV).a RefManual.html syncsource
	cd ~/src/accis/ ; make mproper
	rsync -av --exclude '.git' ~/src/accis/ ~/src/coptic/accis/ 
	date > synccoptic

syncsceptic : lib$(ACCISDRV).a RefManual.html
	cd ~/src/accis/ ; make mproper
	rsync -av --exclude '.git' ~/src/accis/ ~/src/sceptic/accis/ 
	date > synsceptic

tests : 
	make
	make `ls testing/*.f | sed -e "s/[.]f//g"`
	@echo Run executables in testing/

cleantest :
	rm -f `ls testing/*.f | sed -e "s/[.]f//g"`

mproper : clean
	rm -f *.a *~ eye.dat sync*
	rm -f vecx vecglx vec4014
	rm -f convert interface.f90

clean : cleantest
	rm -f *.o
	rm -f *.ida
	rm -f plot*.ps plot*.pg*

help :
	@echo Targets: mproper clean cleantest tests interface.f90 sync

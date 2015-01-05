# Makefile for accis routines.
G90=gfortran

# These defaults generally set in ./configure, but we allow null configure.
#Default compiler
ifeq ("$(G77)","")
    G77=gfortran
else
    COMPILER=$(G77)
endif
#Default xlib path
ifeq ("$(XLIB)","")
    XLIB=/usr/lib/mesa/
endif
#Default driver choice
ifeq ("$(VECX)","")
    VECX=vecx
endif
ACCISDRV=accisX
ifeq ("$(VECX)","vec4014")
    ACCISDRV=accis
endif

WSWITCHES=-Wall
ifeq ("$(G77)","gfortran")
	WSWITCHES=-Wall -Wno-conversion
endif

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

ifeq ("$(VECX)","vec4014")
	libraries = -L.
else
	libraries = -L$(XLIB) -L. -lX11 $(GLLIBS)
endif

ifeq ("$(VECX)","vecglx")
	GLLIBS= -lGL -lGLU
else
	GLLIBS=
endif

#header include files
headers= hidcom.h plotcom.h world3.h

#The real makefile
MAKEFILE=makefile

#######################################################################
# Start of targets
#######################################################################

all : $(MAKEFILE) lib$(ACCISDRV).a

# A way of making a general makefile. By default the file 
# makefile is used if present, but if not, Makefile is used. 
# So we start as Makefile. The configure writes a new one: makefile.
# However, if it doesn't we have an infinite loop, which is bad.
# So we make the second call to an explicit file.
makefile : Makefile configure
	@echo Configuring the Makefile for this platform. 
	@export COMPILER="$(COMPILER)"; ./configure
	@echo Now running make again using the new Makefile.
# The problem with this recursion is that it does not override any command
# argument such as G77= Thus G77 might be forced blank inside the submake
# because by default command arguments override internal settings.
# It is not enough to give MAKEFILE= in the make arguments. You must do this:
	export MAKEFLAGS=; $(MAKE) -f $(MAKEFILE) $(MAKECMDGOALS)
	@echo =================================================================
	@echo Completed make with new makefile. [Ignore any following error.]
	@touch $(VECX)
# Touch vecx to prevent a libaccisX remake with the outermost VECX definition,
# For which vecx is standard.

RefManual.html : RefManual.tex RefManual.pdf
	if which tth ; then tth -e2 RefManual ; fi
	rm -f *.log *.out *.dvi *.toc *.ilg *.idx *.synctex.gz

RefManual.pdf : RefManual.tex
	pdflatex RefManual
	rm -f *.log *.out *.dvi *.toc *.ilg *.idx *.synctex.gz

#update the libraries.
libaccis.a : $(object_files)
	echo "Updating libaccis."
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
	$(CC) -c $(WSWITCHES) vecglx.c

# The file drwstr.f must be compiled with this switch which disables
# for gnu compilers the interpretation of backslashes as escapes.
# drwstr needs the NOBACKSLASH switch (and can't be in :: rule above). 
$(noback_object_files) : drwstr.f fontdata.f vecnp.f $(headers) $(MAKEFILE)
	$(G77) -c $(NOBACKSLASH) $(WSWITCHES) $*.f

# Specific test programs of one sort and another should be put here 
# only if they need special switches.
testing/plottest90 : testing/plottest90.f90 interface.f90
	$(G90) $(WSWITCHES) -o testing/plottest90 testing/plottest90.f90  -l$(ACCISDRV) $(libraries)

#pattern rule, compile using the external definitions of commons
%.o : %.f ;
	$(G77) -c $(WSWITCHES) $*.f

# For fortran 90 executables.
% : %.f90 lib$(ACCISDRV).a $(VECX)
	$(G90) $(WSWITCHES) -o $* $*.f90  -l$(ACCISDRV) $(libraries)

# The main f77 executable pattern.
% : %.f lib$(ACCISDRV).a $(VECX)
	$(G77) $(WSWITCHES) -o $* $*.f  -l$(ACCISDRV) $(libraries)

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

# There are major with Metcalf's automatic converter and parameters that
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
	rsync -u -e ssh  --copy-links -v *.h *.f *.c RefManual.* configure Makefile silas:~/accis/
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
	rm -f makefile
	rm -f vecx vecglx vec4014
	rm -f convert interface.f90

clean : cleantest
	rm -f *.o
	rm -f *.ida
	rm -f plot*.ps plot*.pg*


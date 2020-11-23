#######################################################################
# Makefile for accis routines.
# The shell assumed in this Makefile
SHELL=/bin/bash
#########################################################################
# Include the make commands to define compiler and locations.
include ACCIS.mk
#########################################################################
# Cross-compiling definitions obviously need these executables.
#FORTRAN=i686-w64-mingw32-gfortran -H -mwindows -mconsole --static
#--verbose --static
#GCC=i686-w64-mingw32-gcc -H -mwindows -mconsole
#AR=i686-w64-mingw32-gcc-ar
#########################################################################
# To get to compile with X, you might need to supplement this path.
LIBPATH:=-L. $(LIBPATH)
ifeq ("$(FORTRAN)","ftn")
LIBPATH:=-dynamic $(LIBPATH)
endif
G90=gfortran
#######################################################################
ifeq ("$(NOBACKSLASH)","")
NOBACKSLASH:=\
$(shell \
if [ -n "`$(FORTRAN) --version 2>&1 | grep PathScale`" ] ; then\
 echo "-backslash";else \
if [ -n "`$(FORTRAN) --version 2>&1 | grep Portland`" ] ; then\
 echo "-Mbackslash";else \
if [ -n "`$(FORTRAN) --version 2>&1 | grep GNU`" ] ; then\
 echo "-w -fno-backslash"; \
fi fi fi\
)
endif
##########################################################################
ifneq ("$(VECX)","")
# VECX explicitly set. Use the tests to convey the choice.
 ifeq ("$(TESTX11)","")
# If vecx is available, then if vecglx not the explicit set rule it out.
   ifneq ("$(VECX)","vecglx")
    TESTGL=NO
   endif
 endif
# And if the explicit setting is not vecx, rule it out.
 ifneq ("$(VECX)","vecx")
   TESTX11=NO
 endif
endif
##########################################################################
# VECX choice. Preference order vecx, vecglx, vec4014.
# Normally if vecx fails then vecglx will also fail. 
ifeq ("$(TESTX11)","")
  VECX=vecx
  ACCISDRV=accisX
  libraries= $(LIBPATH) -lX11
else
  ifeq ("$(TESTGL)","")
    VECX=vecglx
    ACCISDRV=accisX
    libraries=$(LIBPATH) -lX11 -lGL -lGLU 
  else
    VECX=vec4014
    ACCISDRV=accis
    libraries=-L.
  endif
endif
#########################################################################
ifeq ("$(FORTRAN)","gfortran")
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
sysfont.o\
slicing.o\
fillgrad.o\
boxcarave.o\
surfaces.o\
acpath.o\
shapes.o

noback_object_files = drwstr.o fontdata.o vecnp.o 
object_files=$(standard_object_files) $(noback_object_files)
root_files=$(patsubst %.o,%,$(object_files))
#header include files
headers= hidcom.h plotcom.h world3.h
#The real Makefile
MAKEFILE=Makefile ACCIS.mk
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
	@echo "Updating libaccis. For $(FORTRAN), $(VECX), $(ACCISDRV)."
	$(AR) -rs libaccis.a $(object_files)

# Headers must be updated, but this section seems to override the
# implicit rule unless the header dependency is :: double.
# Objects depend on headers
$(standard_object_files) :: $(headers) $(MAKEFILE)

#vecx is the C linkage to X-window. Threading is not (yet) completed.
vecx.o : vecx.c $(MAKEFILE)
	$(CC)  $(THREADING) -c vecx.c

vecglx.o : vecglx.c
	$(CC) $(THREADING) -c vecglx.c
#	$(CC) -c $(WSWITCHES) $(THREADING) vecglx.c

#vecwin is the C MS windows driver only made by cross-compiling
vecwin.o : vecwin.c
	$(GCC) $(libraries) -Wall -c vecwin.c

# The file drwstr.f must be compiled with this switch which disables
# for gnu compilers the interpretation of backslashes as escapes.
# drwstr needs the NOBACKSLASH switch (and can't be in :: rule above). 
$(noback_object_files) : drwstr.f fontdata.f vecnp.f $(headers) $(MAKEFILE)
	$(FORTRAN) -c $(NOBACKSLASH) $(COMPILE-SWITCHES) $(WSWITCHES) $*.f

# Specific test programs of one sort and another should be put here 
# only if they need special switches.
testing/plottest90 : testing/plottest90.f90 interface.f90
	$(G90) $(WSWITCHES) -o testing/plottest90 testing/plottest90.f90 -l$(ACCISDRV) $(libraries)

#pattern rule, compile using the external definitions of commons
%.o : %.f ;
	$(FORTRAN) -c $(COMPILE-SWITCHES) $(WSWITCHES) $*.f

# For fortran 90 executables.
% : %.f90 lib$(ACCISDRV).a $(VECX)
	$(G90) $(COMPILE-SWITCHES) $(WSWITCHES) -o $* $*.f90 -l$(ACCISDRV)  $(libraries)

# The main f77 executable pattern.
% : %.f lib$(ACCISDRV).a $(VECX)
	$(FORTRAN) $(COMPILE-SWITCHES) $(WSWITCHES) -o $* $*.f -l$(ACCISDRV)  $(libraries)

%.4014 : %.f libaccis.a
	$(FORTRAN) $(COMPILE-SWITCHES) $(WSWITCHES) -o $*.4014 $*.f -L. -laccis

vecx libaccisX.a : libaccis.a vecx.o
	cp libaccis.a libaccisX.a
	$(AR) -d libaccisX.a vec4014.o
	$(AR) -q libaccisX.a vecx.o
	date >vecx
	rm -f vecglx

vecglx : libaccis.a vecglx.o
	cp libaccis.a libaccisX.a
	$(AR) -d libaccisX.a vec4014.o
	$(AR) -q libaccisX.a vecglx.o
	date >vecglx
	rm -f vecx

vec4014 : libaccis.a
	$(AR) -d libaccis.a noscreen.o
	$(AR) -q libaccis.a vec4014.o
	date >vec4014
	rm -f noscreen

# Make the windows driver and a test executable.
vecwin libaccisWin.a : 
	make clean
	make GCC="i686-w64-mingw32-gcc -H -mwindows -mconsole" vecwin.o
	make FORTRAN="i686-w64-mingw32-gfortran -H -mwindows -mconsole --static" GCC="i686-w64-mingw32-gcc -H -mwindows -mconsole" libaccis.a
	cp libaccis.a libaccisWin.a
	i686-w64-mingw32-gcc-ar -d libaccisWin.a vec4014.o
	i686-w64-mingw32-gcc-ar -q libaccisWin.a vecwin.o
# Use this as the template for windows executables
	i686-w64-mingw32-gfortran -H -mwindows -mconsole --static -o testing/plottest.exe testing/plottest.f -L. -laccisWin
	rm -f *.o libaccis.a

noscreen : libaccis.a
	sed -e "s/data i14no\/0\//data i14no\/1\//" vec4014.f >noscreen.f
	$(FORTRAN) -c $(WSWITCHES) noscreen.f
	$(AR) -d libaccis.a vec4014.o
	$(AR) -q libaccis.a noscreen.o
	date >noscreen
	rm -f vec4014

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
sync : gitcommit syncsource synccoptic syncsceptic

gitcommit :
	@if [ -n "`git status | grep nothing`" ] ; then \
 git push; else echo "GIT COMMIT CHANGES FIRST. NOTHING DONE. git status:";\
 git status ; exit 1 ; fi

syncsilas : lib$(ACCISDRV).a RefManual.html
	rsync -u -e ssh  --copy-links -v *.h *.f *.c RefManual.* *.png Makefile silas:~/accis/
	date > syncsilas

syncsource : lib$(ACCISDRV).a RefManual.html
	cd ~/src/accis/ ; git fetch origin; git reset --hard origin/master
	date > syncsource

synccoptic : lib$(ACCISDRV).a RefManual.html syncsource
	cd ~/src/coptic/accis/ ; make mproper
	rsync -av --exclude '.git' ~/src/accis/ ~/src/coptic/accis/ 
	date > synccoptic

syncsceptic : lib$(ACCISDRV).a RefManual.html syncsource
	cd ~/src/sceptic/accis/ ; make mproper
	rsync -av --exclude '.git' ~/src/accis/ ~/src/sceptic/accis/ 
	date > syncsceptic

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
	rm -f noscreen*

help :
	@echo
	@echo Targets: mproper clean cleantest tests interface.f90 sync
	@echo
	@echo 'Additional makefile definitions. Include as $$ make <def> ...'
	@echo Compiler: FORTRAN=gfortran
	@echo 'Drivers, Chosen: VECX=vecx'
	@echo 'Possible choices:'
	@echo '         VECX=vecx       plain X11 driver'
	@echo '         VECX=vecglx     OpenGL driver'
	@echo '         VECX=vec4014    Tektronix 4014 driver'
	@echo '         VECX=noscreen   No Screen driver. ps output only.'
	@echo To specify this local version use:  ACCISPARENT=../
	@tail -n 13 ACCIS.mk

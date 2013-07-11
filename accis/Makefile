# Makefile for accis routines.

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
winset.o \
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
boxcarave.o

noback_object_files = drwstr.o fontdata.o vecnp.o 

object_files=$(standard_object_files) $(noback_object_files)

ifeq ("$(VECX)","vec4014")
	libraries = -L.
else
	libraries = -L$(XLIB) -L. -lXt -lX11 $(GLLIBS)
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

all : $(MAKEFILE) lib$(ACCISDRV).a accissum.html

# A way of making a general makefile. By default the file 
# makefile is used if present, but if not, Makefile is used. 
# So we start as Makefile. The configure writes a new one: makefile.
# However, if it doesn't we have an infinite loop, which is bad.
# So we make the second call to an explicit file.
makefile : Makefile configure
	@echo Configuring the Makefile for this platform. 
	export COMPILER="$(COMPILER)"; export VECX=$(VECX); ./configure
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

accissum.html : accissum.tex
	tth -e2 accissum
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

#vecx is the C linkage to X-window. If told -DPTHREADS it uses threads.
vecx.o : vecx.c $(MAKEFILE)
	$(CC)  $(THREADING) -c vecx.c

vecglx.o : vecglx.c
	$(CC) -c vecglx.c

# The file drwstr.f must be compiled with this switch which disables
# for gnu compilers the interpretation of backslashes as escapes.
# drwstr needs the NOBACKSLASH switch (and can't be in :: rule above). 
$(noback_object_files) : drwstr.f fontdata.f vecnp.f $(headers) $(MAKEFILE)
	$(G77) -c $(NOBACKSLASH) $*.f

#pattern rule, compile using the external definitions of commons
%.o : %.f ;
	$(G77) -c $*.f

# Specific test programs of one sort and another should be put here 
# only if they need special switches.

# The main executable pattern.
% : %.f lib$(ACCISDRV).a $(VECX)
	$(G77) -o $* $*.f  -l$(ACCISDRV) $(libraries)

%.4014 : %.f libaccis.a
	$(G77) -o $*.4014 $*.f -L. -laccis

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

# Synchronization of versions on other machines.
sync : syncsource synccoptic syncsceptic

syncsceptic : lib$(ACCISDRV).a accissum.html
	rsync -u -e ssh --copy-links -v *.h *.f *.c accissum.* configure Makefile ~/src/sceptic/accis/
	date > synchorace

syncsilas : lib$(ACCISDRV).a accissum.html
	rsync -u -e ssh  --copy-links -v *.h *.f *.c accissum.* configure Makefile silas:~/accis/
	date > syncsilas

syncunity : lib$(ACCISDRV).a accissum.html
	rsync -u -e ssh  --copy-links -v *.h *.f *.c accissum.* configure Makefile unity.psfc.mit.edu:~/accis/
	date > syncunity

syncsource : lib$(ACCISDRV).a accissum.html
	rsync -u -e ssh  --copy-links -v *.h *.f *.c accissum.* configure Makefile ~/src/accis/
	date > syncsource

synccoptic : lib$(ACCISDRV).a accissum.html
	rsync -u -e ssh  --copy-links -v *.h *.f *.c accissum.* configure Makefile ~/src/coptic/accis/
	date > synccoptic

mproper : clean
	rm -f *.a
	rm -f makefile
	rm -f vecx vecglx vec4014

clean :
	rm -f *.o
	rm -f *.ida
	rm -f plot*.ps
	rm -f *test fontshow



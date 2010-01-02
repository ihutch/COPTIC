# Makefile for accis routines. CCPIC version.

#Default xlib path (32 bit).
ifeq ("$(XLIB)","")
    XLIB=/usr/X11R6/lib/
endif
#Default driver choice
#Alternatives are vecx or vecglx
ifeq ("$(VECX)","")
    VECX=vecglx
endif

F77=$(G77)
# Use correct interpretation and suppress g77 warnings for font etc.
NOBACKSLASH=-fno-backslash -w

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
fillgrad.o

noback_object_files = drwstr.o fontdata.o vecnp.o 

object_files=$(standard_object_files) $(noback_object_files)

# For 64-bit machines this ought to be lib64/
#libraries = -L$(XLIB) -L. -lXt -lX11 $(GLLIBS)

libraries = -L$(XLIB) -L. -lXt -lX11 $(GLLIBS)

GLLIBS= -lGL -lGLU


#header include files
headers= hidcom.h plotcom.h world3.h

#The real makefile
MAKEFILE=Makefile

all : $(MAKEFILE) libaccisX.a accissum.html

accissum.html : accissum.tex
	tth -e2 accissum
	purge.tex

#update the libraries.
libaccis.a : $(object_files) $(headers)
	echo "Updating libaccis."
	ar -rs libaccis.a $(object_files)

libaccisX.a : libaccis.a $(VECX)

#	cp libaccis.a libaccisX.a
#	ar -d libaccisX.a vec4014.o
#	ar -q libaccisX.a vecx.o

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
	$(F77) -c $(NOBACKSLASH) $*.f

#pattern rule, compile using the external definitions of commons
%.o : %.f ;
	$(F77) -c $*.f

# Specific test programs of one sort and another should be put here 
# only if they need special switches.

# The main executable pattern.
% : %.f libaccisX.a $(VECX)
	$(F77) -o $* $*.f  -laccisX $(libraries)

%.4014 : %.f libaccis.a
	$(F77) -o $*.4014 $*.f -L. -laccis

vecx : libaccis.a vecx.o
	cp libaccis.a libaccisX.a
	ar -d libaccisX.a vec4014.o
	ar -q libaccisX.a vecx.o
	date >vecx

vecglx : libaccis.a vecglx.o
	cp libaccis.a libaccisX.a
	ar -d libaccisX.a vec4014.o
	ar -q libaccisX.a vecglx.o
	date >vecglx

# Synchronization of versions on other machines.
sync : synchorace syncunity syncsilas

synchorace : libaccisX.a accissum.html
	rsync -u -e ssh --copy-links -v *.h *.f *.c accissum.* configure Makefile tex2html.net:~/accis/
	date > synchorace

syncsilas : libaccisX.a accissum.html
	rsync -u -e ssh  --copy-links -v *.h *.f *.c accissum.* configure Makefile silas:~/accis/
	date > syncsilas

syncunity : libaccisX.a accissum.html
	rsync -u -e ssh  --copy-links -v *.h *.f *.c accissum.* configure Makefile unity.psfc.mit.edu:~/accis/
	date > syncunity

mproper : clean
	rm -f *.a
	rm -f makefile

clean :
	rm -f *.o
	rm -f *.ida
	rm -f plot*.ps
	rm -f *test fontshow


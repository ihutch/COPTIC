LIBRARIES = -L/usr/X11R6/lib/ -L/home/hutch/accis/ -laccisX -lXt -lX11 
# Things just needed for the test routine:
UTILITIES=udisplay.o
# The sormpi system.
OBJECTS=sormpi.o sorrelaxgen.o mpibbdy.o  cijroutine.o cijplot.o 3dobjects.o mditerate.o interpolations.o svdsol.o getfield.o padvnc.o chargetomesh.o slicesect.o randf.o randc.o reinject.o pinit.o bbdyroutine.o ccpicplot.o
HEADERS=bbdydecl.f meshcom.f objcom.f 3dcom.f partcom.f
TARGETS=mpibbdytest mditeratetest sormpitest fieldtest
G77=mpif77
#COMPILE-SWITCHES = -Wall -O2  -I. 
#COMPILE-SWITCHES = -Wall  -O2 -I. -g -ffortran-bounds-check
COMPILE-SWITCHES = -Wall  -O2 -I.
NOBOUNDS= -Wall -O2 -I.
PROFILING= -pg
#PROFILING=

# If this rule does not seem to recognize the file you are trying to make,
# then run 'make' to completion first. It is something to do with the
# match-anything rules and prerequisites. I think that the rule is being
# interpreted as "terminal" which means it does not apply unless its
# prerequisites exist. By my interpretation of the info, this ought not
# to be happening (requires a double colon), but on horace, it is. 
% : %.f  makefile $(OBJECTS) $(UTILITIES)
	$(G77)  -o $* $(COMPILE-SWITCHES) $(PROFILING) $*.f $(OBJECTS) $(UTILITIES) $(LIBRARIES)

%.o : %.f makefile $(HEADERS)
	$(G77)  -c $(COMPILE-SWITCHES) $(PROFILING) $*.f

#default target
smt.out : ccpic
	if [ -f smt.out ] ; then mv smt.out smt.prev ; fi
	./ccpic -p
	diff smt.prev smt.out

# Things to compile without the standard switches

interpolations.o : interpolations.f makefile $(HEADERS)
	$(G77)  -c $(NOBOUNDS) $(PROFILING) $*.f

getfield.o : getfield.f makefile $(HEADERS)
	$(G77)  -c $(NOBOUNDS) $(PROFILING) $*.f

# Main program explicit to avoid make bugs:
ccpic : ccpic.f makefile $(OBJECTS) $(UTILITIES)
	$(G77)  -o ccpic $(COMPILE-SWITCHES) ccpic.f  $(OBJECTS) $(UTILITIES) $(LIBRARIES)

clean :
	rm -f *.o $(TARGETS)


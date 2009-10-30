#GLULIBS= -lGL -lGLU
LIBRARIES = -L/usr/X11R6/lib/ -L/home/hutch/accis/ -laccisX -lXt -lX11 $(GLULIBS)
# Things just needed for the test routine:
UTILITIES=udisplay.o
# The reinjection choice:
#REINJECT=reinject.o
REINJECT=orbitinjnew.o extint.o
# The sormpi system.
FIXEDOBJECTS=sormpi.o sorrelaxgen.o mpibbdy.o  cijroutine.o cijplot.o 3dobjects.o mditerate.o interpolations.o svdsol.o getfield.o padvnc.o chargetomesh.o slicesect.o randf.o randc.o reindiag.o pinit.o bbdyroutine.o ccpicplot.o volint.o fluxdata.o stringsnames.o meshconstruct.o partwriteread.o checkcode.o reduce.o stress.o average.o

OBJECTS=$(FIXEDOBJECTS) ${REINJECT}
HEADERS=bbdydecl.f meshcom.f objcom.f 3dcom.f partcom.f rancom.f ran1com.f
TARGETS=mpibbdytest mditeratetest sormpitest fieldtest
G77=mpif77
#COMPILE-SWITCHES = -Wall -O2  -I. 
#COMPILE-SWITCHES = -Wall  -O2 -I. -g -fbounds-check
#COMPILE-SWITCHES = -Wall  -O2 -I. -g 
COMPILE-SWITCHES = -Wall -Wno-unused -O2 -I.
#COMPILE-SWITCHES = -Wall -Wno-unused -I.
NOBOUNDS= -Wall -Wno-unused -O2 -I.
#PROFILING= -pg
PROFILING=

# If this rule does not seem to recognize the file you are trying to make,
# then run 'make' to completion first. It is something to do with the
# match-anything rules and prerequisites. I think that the rule is being
# interpreted as "terminal" which means it does not apply unless its
# prerequisites exist. By my interpretation of the info, this ought not
# to be happening (requires a double colon), but on horace, it is. 
% : %.f  makefile $(OBJECTS) $(UTILITIES) /home/hutch/accis/libaccisX.a
	$(G77)  -o $* $(COMPILE-SWITCHES) $(PROFILING) $*.f $(OBJECTS) $(UTILITIES) $(LIBRARIES)

%.o : %.f makefile $(HEADERS)
	$(G77)  -c $(COMPILE-SWITCHES) $(PROFILING) $*.f

%.o : %.c makefile
	cc -c $(PROFILING) $*.c

#default target
smt.out : ccpic ccpicgeom.dat
	if [ -f smt.out ] ; then mv smt.out smt.prev ; fi
	./ccpic
	if [ -f smt.prev ] ; then diff smt.prev smt.out ; fi

fieldtest : fieldtest.f makefile $(OBJECTS) /home/hutch/accis/libaccisX.a
	$(G77)  -o fieldtest $(COMPILE-SWITCHES) $(PROFILING) fieldtest.f $(OBJECTS) $(UTILITIES) $(LIBRARIES)
#	$(G77)  -o fieldtest $(COMPILE-SWITCHES) $(PROFILING) fieldtest.f $(FIXEDOBJECTS) reinject.o $(UTILITIES) $(LIBRARIES)

#mpi checking target
mpicheck : ccpic
	mpiexec -l -n 2 ./ccpic >mpicheck.out
	if [ -f mpicheck.prev ] ; then diff mpicheck.prev mpicheck.out ; else mv mpicheck.out mpicheck.prev  ; fi
	mv mpicheck.out mpicheck.prev

# Things to compile without the standard switches

interpolations.o : interpolations.f makefile $(HEADERS)
	$(G77)  -c $(NOBOUNDS) $(PROFILING) $*.f

getfield.o : getfield.f makefile $(HEADERS)
	$(G77)  -c $(NOBOUNDS) $(PROFILING) $*.f

# Main program explicit to avoid make bugs:
ccpic : ccpic.f makefile $(OBJECTS) $(UTILITIES) /home/hutch/accis/libaccisX.a
	$(G77)  -o ccpic $(COMPILE-SWITCHES) $(PROFILING) ccpic.f  $(OBJECTS) $(UTILITIES) $(LIBRARIES)

mpibbdytest : mpibbdytest.o udisplay.o mpibbdy.o mditerate.o reduce.o
	$(G77) -o mpibbdytest  mpibbdytest.f mpibbdy.o udisplay.o  mditerate.o reduce.o

stresstest : stresstest.f stress.o /home/hutch/accis/libaccisX.a
	$(G77) -o stresstest stresstest.f stress.o $(LIBRARIES)

clean :
	rm -f *.o $(TARGETS) *.html *.flx *.phi T*.0?? *.ps

ftnchek :
	./ftnchekrun "ccpic.f $(OBJECTS)"
	@echo To view do: konqueror CallTree.html &

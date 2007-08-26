LIBRARIES = -L/usr/X11R6/lib/ -L/home/hutch/accis/ -laccisX -lXt -lX11 
# Things just needed for the test routine:
UTILITIES=udisplay.o
# The sormpi system.
OBJECTS=sormpi.o sorrelaxgen.o mpibbdy.o  cijroutine.o cijplot.o 3dobjects.o mditerate.o interpolations.o svdsol.o getfield.o padvnc.o chargetomesh.o slicesect.o
HEADERS=bbdydecl.f meshcom.f objcom.f 3dcom.f partcom.f
TARGETS=mpibbdytest mditeratetest sormpitest
G77=mpif77
#COMPILE-SWITCHES = -Wall -O2  -I. 
COMPILE-SWITCHES = -Wall  -I. -g  -ffortran-bounds-check



% : %.f makefile $(OBJECTS) $(UTILITIES)
	$(G77)  -o $* $(COMPILE-SWITCHES) $*.f $(OBJECTS) $(UTILITIES) $(LIBRARIES)

%.o : %.f makefile $(HEADERS)
	$(G77)  -c $(COMPILE-SWITCHES) $*.f

#default target
smt.out : sormpitest
	if [ -f smt.out ] ; then mv smt.out smt.prev ; fi
#	mpiexec -n 8 ./sormpitest -p
#	mpiexec -n 1 ./sormpitest -p
	./sormpitest -p
#	diff testing/smt.fixed smt.out
	diff smt.prev smt.out

# Things to compile without the standard switches

interpolations.o : interpolations.f makefile $(HEADERS)
	$(G77)  -c -Wall $*.f

getfield.o : getfield.f makefile $(HEADERS)
	$(G77)  -c -Wall $*.f

sormpitest : sormpitest.f makefile $(OBJECTS) $(UTILITIES)
	$(G77)  -o sormpitest $(COMPILE-SWITCHES) sormpitest.f  $(OBJECTS) $(UTILITIES) $(LIBRARIES)

clean :
	rm -f *.o $(TARGETS)


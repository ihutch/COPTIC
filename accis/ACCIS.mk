SHELL=/bin/bash
AR=ar
#########################################################################
ifeq ("$(FORTRAN)","")
# Configure compiler. Mostly one long continued bash script.
# Preference order mpif90, ftn, gfortran, f77
 FORTRAN:=\
$(shell \
 if which mpif90 >/dev/null 2>&1; then echo -n "mpif90";else\
  if which ftn >/dev/null 2>&1 ; then echo -n "ftn";else\
    if which gfortran >/dev/null 2>&1; then echo -n "gfortran";else\
     if which f77 >/dev/null 2>&1 ; then echo -n "f77";else\
	echo "Unable to decide compiler. Specify via FORTRAN=..." >&2; exit 1;\
     fi;\
    fi;\
  fi;\
 fi;\
)
endif
export FORTRAN
#########################################################################
# Define the directories, variables, defaults that ACCIS uses
ACCISPARENT= $(HOME)/src/
ACCISHOME=${ACCISPARENT}accis/
ACCISX=$(ACCISHOME)libaccisX.a
LIBPATH= -L$(ACCISHOME) -L.
LIBRARIES = -laccisX -lX11
LIBDEPS = $(ACCISHOME)libaccisX.a
COMPILE-SWITCHES = -Wall -O2
# -fbounds-check
#########################################################################
# Always check that the accis library is available and make it,
# unless we are in the ACCISHOME directory doing things explicitly.
ACCISCHECK:=\
$(shell if [ "${PWD}/" != "${ACCISHOME}" ];\
 then echo -n >&2 "Checking accis library ... ";\
 if [ -f "${ACCISX}" ] ; then echo>&2 "Library ${ACCISX} exists."; else\
   if [ -d "${ACCISPARENT}" ] ; then echo>&2 -n "src directory exists. ";\
     else mkdir ${ACCISPARENT} ; fi;\
   if [ -d "${ACCISHOME}" ] ; then echo>&2 -n "accis directory exists. ";\
     else if cd ${ACCISPARENT}; then\
	git clone https://github.com/ihutch/accis.git; cd - >/dev/null; fi;fi;\
     if cd ${ACCISHOME}; then make >&2; cd - >/dev/null; fi;\
   if [ -f "${ACCISX}" ] ; then echo>&2 "Made ${ACCISX}";\
     else echo>&2 "Error making ${ACCISX}"; fi;\
 fi;fi;\
)
#########################################################################
# This dependency should be included in the application makefile.
#$(ACCISX) : $(ACCISHOME)Makefile
#	@echo "$(ACCISCHECK)"
#	cd $(ACCISHOME); make; cd -
#########################################################################

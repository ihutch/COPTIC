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
ifeq ("$(ACCISPARENT)","") 
  ACCISPARENT:= $(HOME)/src
endif
ACCISHOME:=$(realpath $(ACCISPARENT))/accis
ACCISX=$(ACCISHOME)/libaccisX.a
LIBPATH= -L $(ACCISHOME) -L.
LIBRARIES = -laccisX -lX11
LIBDEPS = $(ACCISHOME)/libaccisX.a
#Allow higher makefile to overrule compile switches:
COMPILE-SWITCHES := -Wall -O2 $(COMPILE-SWITCHES)
#COMPILE-SWITCHES = -g -fbounds-check
#########################################################################
# Always check that the accis library is available and make it,
# unless we are in the ACCISHOME directory doing things explicitly.
ACCISCHECK:=\
$(shell if [ "${CURDIR}" != "$(ACCISHOME)" ];\
 then	echo -n >&2 "Checking accis library ... ";\
 if [ -f "${ACCISX}" ] ; then echo>&2 "${ACCISX} exists."; else\
   if [ -d "${ACCISPARENT}" ] ; then echo>&2 -n "parent directory exists. ";\
     else mkdir ${ACCISPARENT} ; fi;\
   if [ -d "${ACCISHOME}" ] ; then echo>&2 "accis directory exists. ";\
     else if cd ${ACCISPARENT}; then\
	git clone https://github.com/ihutch/accis.git; cd - >/dev/null; fi;fi;\
   if cd ${ACCISHOME}; then pwd >&2; make >&2; cd - >/dev/null; fi;\
   if [ -f "${ACCISX}" ] ; then echo>&2 "Made ${ACCISX}";\
     else echo>&2 "Error making ${ACCISX}"; fi;\
 fi;\
else echo >&2 "${CURDIR}" is the ACCISHOME: "$(ACCISHOME)". No tests. ; fi;\
)
##########################################################################
ifneq ("$(FORTRAN)","")
# Detect what system we are on and whether X11 is installed.
TESTMACOS:=$(shell ls /Users 2>&1 | grep "No such")
ifeq ("$(TESTMACOS)","") 
#MacOS with xquartz requires:
 CC=gcc -I /opt/X11/include
 LIBPATH:=-L /opt/X11/lib $(LIBPATH)
# Test whether X libraries are found. Null => yes.
 TESTX11:=$(shell $(FORTRAN) $(LIBPATH) -lX11 -o /dev/null 2>&1 | grep X)
   ifneq ("$(TESTX11)","")
     XNOTFOUND:=$(shell \
	echo "No X11 libraries found! On MacOS: brew cask install xquartz" >&2;\
	echo "or                           sudo port install xorg-libX11" >&2;)
   endif
else
# Test whether X libraries are found. Null => yes.
  TESTX11:=$(shell $(FORTRAN) $(LIBPATH) -lX11 -o /dev/null 2>&1 | grep X)
   ifneq ("$(TESTX11)","")
     XNOTFOUND:=$(shell echo\
	 "No X11 libraries found! Install package  libx11-dev" >&2;)
   endif
endif
TESTGL:=$(shell $(FORTRAN)  $(LIBPATH) -lGLU -lGL -o /dev/null 2>&1 | grep GL)
endif
#########################################################################
# To satisfy dependencies by building accis in the standard place, copy
# the file ACCIS.mk to your make directory; insert (before the first target)
#       include ACCIS.mk
# in the definitions of your makefile; and add the dependencies 
#       $(LIBDEPS), and  $(LIBPATH) $(LIBRARIES) to executables.
# These will cause the accis library to be cloned from github and compiled.
# If accis should be made in some non-standard $(ACCISHOME) location then
# the makefile can specify the parent directory by containing
#       ACCISPARENT:=$(realpath <ParentDirectory>)
#	export ACCISPARENT
# just before the   include ACCIS.mk   line. 
#########################################################################

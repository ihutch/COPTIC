#!/bin/bash
FFILES=`echo $* | sed -e "s/[.]o/[.]f/g"`
export FORDOCUROOT=/home/hutch/Download/Fordocu/Fordocu
perl ${FORDOCUROOT}/fordocu.pl -fixed -1 $FFILES

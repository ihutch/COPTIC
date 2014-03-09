#!/bin/bash
# Test printing of postscript and displaying with gv
make printtest
./printtest
gv plot0002.ps

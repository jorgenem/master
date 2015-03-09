#!/bin/bash
rm hermass_rev.exe
gfortran -o hermass_rev.exe hermass_rev.f herwig6521.o -L/home/jorgenem/tools/cernlib/2006b/x86_64-slc5-gcc43-opt/lib/ -lkernlib -lpacklib
./hermass_rev.exe
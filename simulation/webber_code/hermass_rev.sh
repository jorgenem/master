#!/bin/bash
rm hermass_rev.exe hermass_rev_orig.o
gfortran -c herwig6510.f
gfortran -c hermass_rev_orig.f
gfortran -o hermass_rev.exe hermass_rev_orig.o herwig6510.o /home/jorgenem/tools/minuit/libminuit.a
./hermass_rev.exe
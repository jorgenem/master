#!/bin/bash
rm hermass_rev_orig_write_events.exe hermass_rev_orig_write_events.o
gfortran -c herwig6510.f
gfortran -c hermass_rev_orig_write_events.f
gfortran -o hermass_rev_orig_write_events.exe hermass_rev_orig_write_events.o herwig6510.o /home/jorgenem/tools/minuit/libminuit.a
./hermass_rev_orig_write_events.exe
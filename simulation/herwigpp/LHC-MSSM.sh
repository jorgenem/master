#!/bin/bash
make
# Herwig++ read LHC-MSSM.in
Herwig++ run LHC-MSSM.run -N$1 -d1
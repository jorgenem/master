#!/bin/bash
rm main9.3
g++ -std=c++11 minimization9.3_fit_mass_differences_with_dilepton_edge_multiple_combos_summed.cpp -o main9.3 -O3 -I/usr/local/include  -I/usr/local/include -larmadillo
./main9.3